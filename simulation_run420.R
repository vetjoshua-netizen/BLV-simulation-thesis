# Run the full BLV intervention factorial batch.
#
# Design: 3 herd sizes x 5 corrective-management settings x 7 testing
# settings x 4 vaccination coverages = 420 scenarios.
# Each scenario is simulated 100 times with matched replicate seeds.
#
# Output folder (created automatically):
#   results/intervention_factorial_420_100_runs/
#     scenario_grid.csv                 complete design matrix
#     initial_state_used.csv            starting-state settings actually used
#     scenario_001.rds ...              one complete saved result per scenario
#     scenario_summary_partial.csv      refreshed after each completed scenario
#     replicate_final_results_partial.csv
#     progress.csv                      status, elapsed time, and any errors
#     scenario_summary.csv              final scenario summary when complete
#     replicate_final_results.csv       final 420 x 100 run-level table
#
# To resume after closing R or an interruption, re-run this file with
# overwrite_existing <- FALSE. Completed compatible scenario .rds files are
# reused; incomplete scenarios are run again.

script_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_file))
  setwd(dirname(normalizePath(sub("^--file=", "", script_file[1]))))

source("simulation.R")
source("analysis.R")

# -------------------------------------------------------------------------
# Batch settings
# -------------------------------------------------------------------------
n_runs <- 100L
master_seed <- 12345L
horizon_years <- 40L
overwrite_existing <- FALSE

# Keep the original 12.5% calf, 37.5% heifer, 50% cow composition generated
# by make_params().  The burn-in is used for infection prevalence only, not
# for the starting age structure.
initial_state_mode <- "original_parameter_structure"

# Use the median endemic prevalence among burn-in runs that were still
# infected at year 100, calculated separately for each herd size. The same
# prevalence is assigned within each starting age group so that it represents
# whole-herd prevalence, rather than cow-only prevalence.
initial_prevalence_scope <- "all_age_groups"
initial_LH_ratio <- 0.70

# Four coverage levels give the requested 420-scenario factorial design.
# The 70% level is a provisional intermediate coverage. Replace it with the
# final C50-derived coverage BEFORE running if the vaccine dose-response work
# identifies a different value.
vaccination_levels <- c(0, 0.50, 0.70, 1.00)
herd_sizes <- c(100L, 200L, 800L)
management_levels <- c("CM_0", "CM_25", "CM_50", "CM_90", "CM_99")

testing_levels <- data.frame(
  testing_setting = c(
    "none",
    "high_only_3_months", "high_only_6_months", "high_only_12_months",
    "all_positive_3_months", "all_positive_6_months", "all_positive_12_months"
  ),
  testing_enabled = c(FALSE, rep(TRUE, 6L)),
  testing_remove_policy = c("all_positive", rep("high_only", 3L), rep("all_positive", 3L)),
  testing_interval_years = c(0.25, 0.25, 0.50, 1.00, 0.25, 0.50, 1.00),
  stringsAsFactors = FALSE
)

results_dir <- file.path("results", "intervention_factorial_420_100_runs")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# Select and record the burn-in endemic prevalence and starting structure
# -------------------------------------------------------------------------
burnin_file <- file.path(
  "results", "introduction_burnin_replacement_100", "introduction_summary.csv"
)

if (!file.exists(burnin_file))
  stop("Cannot find replacement-rule burn-in summary: ", burnin_file)

burnin <- read.csv(burnin_file, stringsAsFactors = FALSE)
stable_prevalence <- aggregate(
  mean_prevalence_active_runs ~ herd_size,
  data = burnin,
  FUN = median
)
names(stable_prevalence)[2] <- "initial_prevalence"

initial_state_used <- data.frame(
  herd_size = herd_sizes,
  initial_state_mode = initial_state_mode,
  initial_prevalence_scope = initial_prevalence_scope,
  init_calves = round(0.125 * herd_sizes),
  init_heifers = round(0.375 * herd_sizes),
  init_cows = herd_sizes - round(0.125 * herd_sizes) - round(0.375 * herd_sizes),
  initial_prevalence = stable_prevalence$initial_prevalence[match(herd_sizes, stable_prevalence$herd_size)],
  initial_LH_ratio = initial_LH_ratio,
  stringsAsFactors = FALSE
)
write.csv(initial_state_used, file.path(results_dir, "initial_state_used.csv"), row.names = FALSE)

# -------------------------------------------------------------------------
# Build the 420-scenario design matrix
# -------------------------------------------------------------------------
scenario_grid <- expand.grid(
  herd_size = herd_sizes,
  corrective_management = management_levels,
  testing_setting = testing_levels$testing_setting,
  vaccination_coverage = vaccination_levels,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
scenario_grid <- merge(scenario_grid, testing_levels, by = "testing_setting", all.x = TRUE, sort = FALSE)
scenario_grid$corrective_management <- factor(scenario_grid$corrective_management, levels = management_levels)
scenario_grid$testing_setting <- factor(scenario_grid$testing_setting, levels = testing_levels$testing_setting)
scenario_grid <- scenario_grid[order(
  scenario_grid$herd_size,
  scenario_grid$corrective_management,
  scenario_grid$testing_setting,
  scenario_grid$vaccination_coverage
), ]
scenario_grid$corrective_management <- as.character(scenario_grid$corrective_management)
scenario_grid$testing_setting <- as.character(scenario_grid$testing_setting)
scenario_grid$scenario_id <- sprintf("scenario_%03d", seq_len(nrow(scenario_grid)))
scenario_grid$vaccination_enabled <- scenario_grid$vaccination_coverage > 0
scenario_grid$initial_state_mode <- initial_state_mode
scenario_grid$initial_prevalence_scope <- initial_prevalence_scope
scenario_grid <- scenario_grid[, c(
  "scenario_id", "herd_size", "initial_state_mode", "initial_prevalence_scope", "corrective_management",
  "testing_setting", "testing_enabled", "testing_remove_policy",
  "testing_interval_years", "vaccination_enabled", "vaccination_coverage"
)]

stopifnot(nrow(scenario_grid) == 420L)
write.csv(scenario_grid, file.path(results_dir, "scenario_grid.csv"), row.names = FALSE)

model_version <- "v4_global_foi_replacement_22_24m_endemic_prevalence_all_ages_intervention_420"
progress_file <- file.path(results_dir, "progress.csv")
progress <- data.frame(
  scenario_id = scenario_grid$scenario_id,
  status = "pending",
  started_at = NA_character_,
  completed_at = NA_character_,
  elapsed_minutes = NA_real_,
  message = NA_character_,
  stringsAsFactors = FALSE
)
if (file.exists(progress_file)) {
  previous_progress <- read.csv(progress_file, stringsAsFactors = FALSE)
  keep <- intersect(names(progress), names(previous_progress))
  progress[, keep] <- previous_progress[match(progress$scenario_id, previous_progress$scenario_id), keep]
}
write_progress <- function() write.csv(progress, progress_file, row.names = FALSE)
write_progress()

completed_summaries <- list()
completed_finals <- list()

save_partial_tables <- function() {
  if (length(completed_summaries)) {
    write.csv(do.call(rbind, completed_summaries),
              file.path(results_dir, "scenario_summary_partial.csv"), row.names = FALSE)
    write.csv(do.call(rbind, completed_finals),
              file.path(results_dir, "replicate_final_results_partial.csv"), row.names = FALSE)
  }
}

for (i in seq_len(nrow(scenario_grid))) {
  scenario <- scenario_grid[i, , drop = FALSE]
  result_file <- file.path(results_dir, paste0(scenario$scenario_id, ".rds"))
  progress_index <- match(scenario$scenario_id, progress$scenario_id)

  use_existing <- FALSE
  if (file.exists(result_file) && !overwrite_existing) {
    existing <- tryCatch(readRDS(result_file), error = function(e) NULL)
    use_existing <- !is.null(existing) &&
      identical(existing$model_version, model_version) &&
      identical(existing$n_runs, n_runs) &&
      identical(existing$master_seed, master_seed) &&
      identical(existing$scenario, scenario)
  }

  if (use_existing) {
    message("Using existing result: ", scenario$scenario_id, " (", i, "/", nrow(scenario_grid), ")")
    completed_summaries[[length(completed_summaries) + 1L]] <- existing$summary
    completed_finals[[length(completed_finals) + 1L]] <- existing$final_results
    progress$status[progress_index] <- "complete_existing"
    progress$completed_at[progress_index] <- as.character(Sys.time())
    progress$message[progress_index] <- "Compatible saved result reused."
    write_progress()
    save_partial_tables()
    next
  }

  message("Running ", scenario$scenario_id, " (", i, "/", nrow(scenario_grid), ")")
  scenario_started <- Sys.time()
  progress$status[progress_index] <- "running"
  progress$started_at[progress_index] <- as.character(scenario_started)
  progress$completed_at[progress_index] <- NA_character_
  progress$elapsed_minutes[progress_index] <- NA_real_
  progress$message[progress_index] <- NA_character_
  write_progress()

  state_row <- initial_state_used[match(scenario$herd_size, initial_state_used$herd_size), ]
  params_args <- list(
    herd_size = scenario$herd_size,
    horizon_years = horizon_years,
    init_prevalence = state_row$initial_prevalence,
    init_LH_ratio = initial_LH_ratio,
    initial_prevalence_scope = initial_prevalence_scope,
    corrective_management = if (scenario$corrective_management == "CM_0") NULL else scenario$corrective_management,
    testing_enabled = scenario$testing_enabled,
    testing_interval_years = scenario$testing_interval_years,
    testing_target_group = "all",
    testing_remove_policy = scenario$testing_remove_policy,
    vaccination_enabled = scenario$vaccination_enabled,
    vaccination_groups = c("heifer", "cow"),
    vaccination_coverage = scenario$vaccination_coverage
  )
  one_result <- tryCatch({
    params <- do.call(make_params, params_args)
    simulation_results <- run_strategy(params, n_runs = n_runs, master_seed = master_seed)
    scenario_summary <- cbind(scenario, master_seed = master_seed,
                              horizon_years = horizon_years,
                              summarise_intervention_outcomes(simulation_results))
    final_results <- cbind(scenario_id = scenario$scenario_id,
                           summarise_final_results(simulation_results))
    list(params = params, simulation_results = simulation_results,
         summary = scenario_summary, final_results = final_results)
  }, error = function(e) e)

  elapsed_minutes <- as.numeric(difftime(Sys.time(), scenario_started, units = "mins"))
  if (inherits(one_result, "error")) {
    progress$status[progress_index] <- "error"
    progress$completed_at[progress_index] <- as.character(Sys.time())
    progress$elapsed_minutes[progress_index] <- elapsed_minutes
    progress$message[progress_index] <- conditionMessage(one_result)
    write_progress()
    message("ERROR in ", scenario$scenario_id, ": ", conditionMessage(one_result))
    next
  }

  scenario_result <- list(
    model_version = model_version,
    scenario = scenario,
    n_runs = n_runs,
    master_seed = master_seed,
    horizon_years = horizon_years,
    initial_state_used = state_row,
    params = one_result$params,
    summary = one_result$summary,
    final_results = one_result$final_results,
    simulation_results = one_result$simulation_results
  )
  saveRDS(scenario_result, result_file, compress = "gzip")
  completed_summaries[[length(completed_summaries) + 1L]] <- one_result$summary
  completed_finals[[length(completed_finals) + 1L]] <- one_result$final_results
  progress$status[progress_index] <- "complete"
  progress$completed_at[progress_index] <- as.character(Sys.time())
  progress$elapsed_minutes[progress_index] <- elapsed_minutes
  progress$message[progress_index] <- "Saved successfully."
  write_progress()
  save_partial_tables()
  message("Completed ", scenario$scenario_id, " in ", round(elapsed_minutes, 1), " minutes.")
}

save_partial_tables()
if (length(completed_summaries)) {
  scenario_summary <- do.call(rbind, completed_summaries)
  replicate_final_results <- do.call(rbind, completed_finals)
  write.csv(scenario_summary, file.path(results_dir, "scenario_summary.csv"), row.names = FALSE)
  write.csv(replicate_final_results, file.path(results_dir, "replicate_final_results.csv"), row.names = FALSE)
}

write_progress()
message("Finished batch check. Completed scenarios: ", sum(progress$status %in% c("complete", "complete_existing")), " / ", nrow(scenario_grid))
message("Results saved in: ", normalizePath(results_dir))

