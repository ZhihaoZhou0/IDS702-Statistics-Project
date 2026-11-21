# Run scripts/01_modeling.R first in the same R session.
library(dplyr)

# Check the airline-by-period groups used in the interaction model.
cell_counts <- model.frame(m1) |>
  group_by(airline, departure_period, .drop = FALSE) |>
  summarise(
    flights = n(),
    delayed = sum(delayed_15 == 1),
    under_15_min = sum(delayed_15 == 0),
    .groups = "drop"
  ) |>
  arrange(flights)

print(cell_counts, n = Inf)

# Flag groups with fewer than five flights in either outcome category.
small_cells <- cell_counts |>
  filter(delayed < 5 | under_15_min < 5)

print(small_cells, n = Inf)

# Check whether a nonlinear distance term improves model fit.
m_distance_spline <- update(
  m1,
  . ~ . - distance100 + splines::ns(distance100, df = 3),
  data = model_data
)

print(c(linear = nobs(m1), spline = nobs(m_distance_spline)))
print(anova(m1, m_distance_spline, test = "Chisq"))
print(AIC(m1, m_distance_spline))

# Compare predicted and observed delay rates in ten risk groups.
calibration_table <- tibble::tibble(
  predicted = fitted(m1),
  delayed = model.response(model.frame(m1))
) |>
  mutate(risk_group = ntile(predicted, 10)) |>
  group_by(risk_group) |>
  summarise(
    flights = n(),
    predicted_pct = 100 * mean(predicted),
    observed_pct = 100 * mean(delayed),
    .groups = "drop"
  ) |>
  mutate(gap_pp = observed_pct - predicted_pct)

print(calibration_table, n = Inf)

# Sensitivity of evening-to-morning odds ratios to influential flights.
fit_data <- model.frame(m1)
cook_d <- cooks.distance(m1)
top_five <- head(order(cook_d, decreasing = TRUE), 5)

influence_table <- fit_data[top_five, c("airline", "departure_period", "delayed_15")]
influence_table$cooks_d <- round(unname(cook_d[top_five]), 5)
print(influence_table)

m_without_five <- update(m1, data = fit_data[-top_five, , drop = FALSE])

# Factor order in the fitted model is Morning, Afternoon, Evening.
stopifnot(identical(levels(fit_data$departure_period),
                    c("Morning", "Afternoon", "Evening")))

evening_or <- function(model) {
  means <- emmeans::emmeans(model, ~ departure_period | airline)
  comparison <- emmeans::contrast(
    means, method = list("Evening vs Morning" = c(-1, 0, 1))
  )
  as.data.frame(summary(comparison, type = "response")) |>
    select(airline, odds.ratio)
}

influence_comparison <- evening_or(m1) |>
  rename(full_OR = odds.ratio) |>
  left_join(
    rename(evening_or(m_without_five), reduced_OR = odds.ratio),
    by = "airline"
  ) |>
  mutate(change_pct = 100 * (reduced_OR / full_OR - 1)) |>
  arrange(desc(abs(change_pct)))

print(c(full = nobs(m1), without_five = nobs(m_without_five)))
print(as.data.frame(influence_comparison), digits = 4, row.names = FALSE)


# Allow flights on the same date to have correlated outcomes.
stopifnot(
  nobs(m1) == nrow(model_data),
  !anyNA(model_data$FlightDate)
)

date_vcov <- sandwich::vcovCL(
  m1, cluster = model_data$FlightDate, type = "HC1"
)

# Compare evening vs morning within each airline using both sets of SEs.
evening_ci <- function(covariance) {
  means <- emmeans::emmeans(
    m1, ~ departure_period | airline, vcov. = covariance
  )
  comparison <- emmeans::contrast(
    means, method = list("Evening vs Morning" = c(-1, 0, 1))
  )
  as.data.frame(confint(comparison, level = 0.95, type = "response")) |>
    select(airline, odds.ratio, asymp.LCL, asymp.UCL)
}

ordinary_ci <- evening_ci(vcov(m1)) |>
  rename(OR = odds.ratio, ordinary_low = asymp.LCL, ordinary_high = asymp.UCL)

clustered_ci <- evening_ci(date_vcov) |>
  select(airline, clustered_low = asymp.LCL, clustered_high = asymp.UCL)

ci_comparison <- left_join(ordinary_ci, clustered_ci, by = "airline")
print(as.data.frame(ci_comparison), digits = 4, row.names = FALSE)

# Joint Wald test of all 14 airline-by-period interaction coefficients.
interaction_rows <- grep(":", names(coef(m1)), fixed = TRUE)
restrictions <- diag(length(coef(m1)))[interaction_rows, , drop = FALSE]
colnames(restrictions) <- names(coef(m1))

print(car::linearHypothesis(
  m1, restrictions, vcov. = date_vcov, test = "Chisq"
))


# Check convergence and whether all coefficients were estimated.
fit_checks <- data.frame(
  converged = m1$converged,
  model_rank = m1$rank,
  total_coefficients = length(coef(m1)),
  missing_coefficients = sum(is.na(coef(m1)))
)
print(fit_checks)

# Check collinearity among the predictors before adding interactions.
vif_table <- car::vif(m0)
print(round(vif_table, 3))