library(ggplot2)
source("simulation.R")
source("analysis.R")


# Scenario name
scenario_name <- "baseline"


# Parameters
p <- make_params(
  herd_size = 100,
  horizon_years = 40,
  
  # 1. Corrective management
  corrective_management = NULL,
  
  # 2. Test-and-remove
  testing_enabled = FALSE,
  testing_interval_years = 0.25,
  testing_target_group = "all",
  testing_remove_policy = "high_only",
  testing_sensitivity = 0.95,
  
  # 3. Vaccination
  vaccination_enabled = TRUE,
  vaccination_groups = c("heifer", "cow"),
  vaccination_coverage = 0.90,
  vaccine_efficacy = 0.70,
)


# Run simulations

n_runs <- 10
master_seed <- 12345L

results <- run_strategy(
  p,
  n_runs = n_runs,
  master_seed = master_seed
)


# outcome
final_summary <- summarise_final_results(results)

print(final_summary)


# Main intervention outcomes: one row per scenario
intervention_summary <- summarise_intervention_outcomes(results)

intervention_summary$scenario <- scenario_name

intervention_summary <- intervention_summary[
  c(
    "scenario",
    "n_runs",
    "eradication_probability",
    "eradication_ci_low",
    "eradication_ci_high",
    "mean_final_infected",
    "mean_final_population",
    "mean_final_prevalence",
    "median_final_prevalence",
    "median_time_to_eradication",
    "mean_tests",
    "mean_test_positive_removed",
    "mean_vaccinations",
    "mean_external_purchases"
  )
]

print(intervention_summary)
View(intervention_summary)

write.csv(
  intervention_summary,
  paste0("intervention_summary_", scenario_name, ".csv"),
  row.names = FALSE
)



# Simple mean final outcomes across all replicates
final_mean <- data.frame(
  scenario = scenario_name,
  n_runs = n_runs,
  
  mean_final_calf = mean(final_summary$final_calf),
  mean_final_heifer = mean(final_summary$final_heifer),
  mean_final_cow = mean(final_summary$final_cow),
  
  mean_final_population = mean(final_summary$final_population),
  mean_final_infected = mean(final_summary$final_infected),
  mean_final_prevalence = mean(final_summary$final_prevalence),
  
  eradication_probability = mean(final_summary$eradicated)
)

print(final_mean)
View(final_mean)
