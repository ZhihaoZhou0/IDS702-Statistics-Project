# 01_modeling.R — Data preparation & logistic regression
library(data.table)
library(dplyr)
library(tidyr)
library(knitr)
library(broom)
library(emmeans)
options(timeout = 600)

# 1. Load data 
data_dir <- file.path("data", "bts_2026")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

fields_needed <- c(
  "Year", "Month", "DayofMonth", "DayOfWeek", "FlightDate",
  "Reporting_Airline", "Flight_Number_Reporting_Airline",
  "Origin", "Dest", "CRSDepTime", "DepTime", "DepDelay",
  "DepDel15", "Cancelled", "Diverted", "Distance",
  "CRSElapsedTime", "ArrDel15", "ArrDelay"
)

read_rdu_month <- function(month_number) {
  file_name <- paste0(
    "On_Time_Reporting_Carrier_On_Time_Performance_1987_present_2026_",
    month_number, ".zip"
  )
  zip_path <- file.path(data_dir, file_name)
  url <- paste0("https://transtats.bts.gov/PREZIP/", file_name)
  
  if (!file.exists(zip_path)) {
    download.file(url, destfile = zip_path, mode = "wb", quiet = TRUE)
  }
  
  zip_contents <- unzip(zip_path, list = TRUE)$Name
  csv_files <- zip_contents[grepl("\\.csv$", zip_contents, ignore.case = TRUE)]
  if (length(csv_files) < 1) {
    stop("No CSV found in BTS ZIP for month ", month_number)
  }
  
  # Work with one month's data at a time to keep memory use lower.
  extract_dir <- file.path(tempdir(), paste0("bts_month_", month_number))
  dir.create(extract_dir, showWarnings = FALSE, recursive = TRUE)
  csv_path <- unzip(zip_path, files = csv_files[1], exdir = extract_dir)[1]
  
  # Check that the required columns exist
  actual_fields <- names(
    data.table::fread(csv_path, nrows = 0)
  )
  
  missing_fields <- setdiff(fields_needed, actual_fields)
  
  if (length(missing_fields) > 0) {
    stop(
      "Missing BTS columns: ",
      paste(missing_fields, collapse = ", ")
    )
  }
  
  flight_data <- data.table::fread(
    csv_path, select = fields_needed,
    na.strings = c("", "NA"), showProgress = FALSE
  )
  unlink(extract_dir, recursive = TRUE)
  
  as_tibble(flight_data) %>%
    filter(Origin == "RDU")
}

rdu_all <- bind_rows(lapply(c(7, 8), read_rdu_month))

if (nrow(rdu_all) < 300) {
  stop("The RDU subset has fewer than 300 records; please review the source files.")
}

rdu_eda <- rdu_all %>%
  mutate(
    departure_hour = floor(as.numeric(CRSDepTime) / 100),
    departure_period = case_when(
      between(departure_hour, 5, 11) ~ "Morning",
      between(departure_hour, 12, 17) ~ "Afternoon",
      between(departure_hour, 18, 23) ~ "Evening",
      TRUE ~ NA_character_
    ),
    departure_period = factor(
      departure_period, levels = c("Morning", "Afternoon", "Evening")
    )
  )

rdu_delay <- rdu_eda %>%
  filter(Cancelled == 0,
         !is.na(DepDel15),
         !is.na(departure_period)) %>%
  mutate(delayed_15 = as.integer(DepDel15))

if (nrow(rdu_delay) < 300) {
  stop("The eligible RDU delay sample has fewer than 300 flights.")
}

# 2. Sample flow (for Technical Appendix)
sample_flow <- tibble(
  Step = c(
    "All RDU departures, Jul-Aug 2026",
    "Excluded: canceled",
    "Excluded: operated but missing delay indicator",
    "Excluded: outside study periods",
    "Final analytic sample"
  ),
  Flights = c(
    nrow(rdu_eda),
    sum(rdu_eda$Cancelled == 1, na.rm = TRUE),
    sum(rdu_eda$Cancelled == 0 & is.na(rdu_eda$DepDel15)),
    sum(rdu_eda$Cancelled == 0 & !is.na(rdu_eda$DepDel15) &
          is.na(rdu_eda$departure_period)),
    nrow(rdu_delay)
  )
)
print(sample_flow)

# 3. Check airline group sizes
# Total Number of Flights by Each Airline
rdu_delay %>%
  count(Reporting_Airline, sort = TRUE) %>%
  print(n = Inf)

# Airline × Time Period: Number of Flights and Number of Delays
rdu_delay %>%
  group_by(Reporting_Airline, departure_period) %>%
  summarise(flights = n(),
            delayed = sum(delayed_15),
            .groups = "drop") %>%
  pivot_wider(names_from = departure_period,
              values_from = c(flights, delayed),
              values_fill = 0) %>%
  print(n = Inf, width = Inf)

# 4. Recode variables for modeling
main_airlines <- c("DL", "AA", "WN", "YX", "UA", "F9", "B6")

model_data <- rdu_delay %>%
  mutate(
    airline = if_else(Reporting_Airline %in% main_airlines,
                      Reporting_Airline, "Other"),
    airline = factor(
      airline,
      levels = c(main_airlines, "Other"),   # First = Reference Group (DL)
      labels = c("Delta", "American", "Southwest", "Republic",
                 "United", "Frontier", "JetBlue", "Other")
    ),
    distance100 = Distance / 100,
    month = factor(Month, levels = c(7, 8),
                   labels = c("July", "August")),
    day_of_week = factor(DayOfWeek, levels = 1:7,
                         labels = c("Monday", "Tuesday", "Wednesday",
                                    "Thursday", "Friday", "Saturday",
                                    "Sunday"))
  )

table(model_data$airline, model_data$departure_period)

# 5. Fit models
# m0: main effects only
m0 <- glm(delayed_15 ~ departure_period + airline + distance100 +
            month + day_of_week,
          family = binomial, data = model_data)

# m1: with departure period x airline interaction (primary model)
m1 <- glm(delayed_15 ~ departure_period * airline + distance100 +
            month + day_of_week,
          family = binomial, data = model_data)

summary(m1)

# Likelihood ratio test for the interaction
anova(m0, m1, test = "Chisq")
AIC(m0, m1)

# 6. Odds ratio table (exponentiated) 
or_table <- tidy(m1, exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(
    term = term %>%
      gsub("departure_period", "", .) %>%
      gsub("airline", "", .) %>%
      gsub("day_of_week", "", .) %>%
      gsub("month", "", .) %>%
      gsub(":", " × ", .) %>%
      sub("distance100", "Distance (per 100 miles)", .) %>%
      sub("\\(Intercept\\)", "Intercept (baseline odds)", .),
    p.value = format.pval(p.value, digits = 2, eps = 0.001)
  ) %>%
  select(Term = term,
         `Odds ratio` = estimate,
         `SE (log-odds scale)` = std.error,
         `95% CI lower` = conf.low,
         `95% CI upper` = conf.high,
         `p-value` = p.value)

kable(or_table, digits = c(0, 2, 3, 2, 2, 0))

# 7. Airline-specific ORs: each period vs morning
emm <- emmeans(m1, ~ departure_period | airline)

period_by_airline <- contrast(emm, method = "trt.vs.ctrl",
                              ref = 1, adjust = "none") %>%
  summary(infer = TRUE, type = "response") %>%
  as_tibble() %>%
  transmute(Airline = airline,
            Comparison = gsub(" / ", " vs ", contrast),
            `Odds ratio` = odds.ratio,
            `95% CI lower` = asymp.LCL,
            `95% CI upper` = asymp.UCL,
            `p-value` = p.value) %>%
  mutate(`p-value` = format.pval(`p-value`, digits = 2, eps = 0.001))

kable(period_by_airline, digits = 2)