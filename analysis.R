# Analysis is separated from the model. It has no automatic runs or plots.
source("simulation.R")

summarise_runs <- function(x) data.frame(
  replicate = seq_along(x$runs),
  seed = vapply(x$runs, function(r) r$seed, integer(1)),
  eradicated = vapply(x$runs, function(r) r$eradicated, logical(1)),
  final_population = vapply(x$runs, function(r) sum(unlist(r$final_population)), integer(1))
)


annual_records_table <- function(run_result) {
    
    do.call(
      rbind,
      lapply(run_result$year_records, function(x) {
        
        data.frame(
          year = x$year,
          time = x$time,
          
          calf = x$population_by_age$calf,
          heifer = x$population_by_age$heifer,
          cow = x$population_by_age$cow,
          
          total_population = sum(
            x$population_by_age$calf,
            x$population_by_age$heifer,
            x$population_by_age$cow
          ),
          
          infected_calf =
            x$population_by_age_state$calf$I_L +
            x$population_by_age_state$calf$I_H,
          
          infected_heifer =
            x$population_by_age_state$heifer$I_L +
            x$population_by_age_state$heifer$I_H,
          
          infected_cow =
            x$population_by_age_state$cow$I_L +
            x$population_by_age_state$cow$I_H,

          total_infected =
            x$population_by_age_state$calf$I_L +
            x$population_by_age_state$calf$I_H +
            x$population_by_age_state$heifer$I_L +
            x$population_by_age_state$heifer$I_H +
            x$population_by_age_state$cow$I_L +
            x$population_by_age_state$cow$I_H,

          prevalence = if (
            sum(
              x$population_by_age$calf,
              x$population_by_age$heifer,
              x$population_by_age$cow
            ) > 0
          ) {
            (
              x$population_by_age_state$calf$I_L +
              x$population_by_age_state$calf$I_H +
              x$population_by_age_state$heifer$I_L +
              x$population_by_age_state$heifer$I_H +
              x$population_by_age_state$cow$I_L +
              x$population_by_age_state$cow$I_H
            ) / sum(
              x$population_by_age$calf,
              x$population_by_age$heifer,
              x$population_by_age$cow
            )
          } else {
            NA_real_
          },
          
          births_female = x$births_female,
          births_male = x$births_male,
          
          natural_deaths = sum(
            unlist(x$deaths_natural)
          ),
          
          deaths_blv = sum(
            unlist(x$dead_blv)
          ),
          
          culled_cow = x$culled_cow,
          external_purchases = x$external_purchases,
          heifers_sold_surplus =
            x$heifers_sold_surplus,

          tested_calf = x$tested_by_age$calf,
          tested_heifer = x$tested_by_age$heifer,
          tested_cow = x$tested_by_age$cow,
          tested_total = sum(unlist(x$tested_by_age)),

          removed_test_positive_calf = x$test_positive_removed_by_age$calf,
          removed_test_positive_heifer = x$test_positive_removed_by_age$heifer,
          removed_test_positive_cow = x$test_positive_removed_by_age$cow,
          removed_test_positive_total = sum(unlist(x$test_positive_removed_by_age)),

          vaccinated_calf = x$vaccinated_by_age$calf,
          vaccinated_heifer = x$vaccinated_by_age$heifer,
          vaccinated_cow = x$vaccinated_by_age$cow,
          vaccinated_total = x$vaccinated_total
        )
      })
    )
  }

# One row per replicate, matching the final-summary table used for the
# original BLV_control1308.R implementation.
summarise_final_results <- function(results) {
  rows <- lapply(seq_along(results$runs), function(i) {
    one_run <- results$runs[[i]]
    last_index <- length(one_run$year_records)
    last_record <- one_run$year_records[[last_index]]

    infected_by_age <- vapply(AGE_CLASSES, function(age_class) {
      states <- last_record$population_by_age_state[[age_class]]
      states$I_L + states$I_H
    }, numeric(1))

    data.frame(
      replicate = i,
      seed = one_run$seed,
      eradicated = one_run$eradicated,
      final_calf = last_record$population_by_age$calf,
      final_heifer = last_record$population_by_age$heifer,
      final_cow = last_record$population_by_age$cow,
      final_population = sum(unlist(last_record$population_by_age)),
      final_infected = sum(infected_by_age),
      stringsAsFactors = FALSE
    )
  })

  output <- do.call(rbind, rows)
  output$final_prevalence <- ifelse(
    output$final_population > 0,
    output$final_infected / output$final_population,
    NA_real_
  )
  output
}

# One row per strategy.  Eradication is the model's operational surveillance
# outcome (three clean surveillance rounds), not guaranteed infection absence.
summarise_intervention_outcomes <- function(results) {
  final <- summarise_final_results(results)
  n_runs <- nrow(final)
  n_eradicated <- sum(final$eradicated)
  eradication_ci <- stats::binom.test(n_eradicated, n_runs)$conf.int

  eradication_times <- vapply(results$runs, function(one_run) {
    if (isTRUE(one_run$eradicated)) one_run$eradication_time else NA_real_
  }, numeric(1))

  burden <- do.call(rbind, lapply(results$runs, function(one_run) {
    annual <- annual_records_table(one_run)
    data.frame(
      tests = sum(annual$tested_total),
      test_positive_removed = sum(annual$removed_test_positive_total),
      vaccinations = sum(annual$vaccinated_total),
      external_purchases = sum(annual$external_purchases)
    )
  }))

  data.frame(
    n_runs = n_runs,
    eradication_probability = n_eradicated / n_runs,
    eradication_ci_low = unname(eradication_ci[1]),
    eradication_ci_high = unname(eradication_ci[2]),
    mean_final_prevalence = mean(final$final_prevalence, na.rm = TRUE),
    median_final_prevalence = stats::median(final$final_prevalence, na.rm = TRUE),
    mean_final_infected = mean(final$final_infected),
    mean_final_population = mean(final$final_population),
    median_time_to_eradication = if (all(is.na(eradication_times))) NA_real_ else
      stats::median(eradication_times, na.rm = TRUE),
    mean_tests = mean(burden$tests),
    mean_test_positive_removed = mean(burden$test_positive_removed),
    mean_vaccinations = mean(burden$vaccinations),
    mean_external_purchases = mean(burden$external_purchases),
    stringsAsFactors = FALSE
  )
}
  
