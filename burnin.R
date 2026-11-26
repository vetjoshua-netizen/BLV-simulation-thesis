
source("simulation.R")

herd_sizes <- c(100L, 200L, 800L)
index_ages <- c("calf", "heifer", "cow")
index_states <- c("I_L", "I_H")
n_runs <- getOption("blv_burnin_n_runs", 100L)
horizon_years <- 100L  
master_seed <- 12345L
selected_scenario_ids <- getOption("blv_burnin_selected_scenario_ids", NULL)
output_dir <- getOption("blv_burnin_output_dir",
                        file.path("results", "introduction_burnin_replacement"))
cache_dir <- getOption("blv_burnin_cache_dir", NULL)
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

scenario_grid <- expand.grid(
  herd_size=herd_sizes,
  index_age=index_ages,
  index_state=index_states,
  KEEP.OUT.ATTRS=FALSE,
  stringsAsFactors=FALSE
)
scenario_grid <- scenario_grid[order(scenario_grid$herd_size,
                                     match(scenario_grid$index_age, index_ages),
                                     match(scenario_grid$index_state, index_states)), ]
scenario_grid$scenario_id <- sprintf("scenario_%03d", seq_len(nrow(scenario_grid)))
scenario_grid <- scenario_grid[, c("scenario_id", "herd_size", "index_age", "index_state")]
stopifnot(nrow(scenario_grid) == 18L)
write.csv(scenario_grid, file.path(output_dir, "scenario_grid.csv"), row.names=FALSE)
selected_indices <- if (is.null(selected_scenario_ids))
  seq_len(nrow(scenario_grid)) else match(selected_scenario_ids, scenario_grid$scenario_id)
if (anyNA(selected_indices))
  stop("Unknown scenario ID in selected_scenario_ids.")


set.seed(master_seed)
replicate_seeds <- sample.int(.Machine$integer.max, n_runs)
model_md5 <- unname(tools::md5sum("simulation.R"))

extract_year <- function(record) {
  counts <- record$population_by_age_state
  age_values <- unlist(lapply(AGE_CLASSES, function(ac) {
    setNames(
      c(record$population_by_age[[ac]],
        vapply(STATES, function(s) counts[[ac]][[s]], integer(1))),
      c(paste0(ac, "_population"), paste0(ac, "_", STATES))
    )
  }))

  population <- sum(unlist(record$population_by_age))
  low <- sum(vapply(AGE_CLASSES, function(ac) counts[[ac]][["I_L"]], integer(1)))
  high <- sum(vapply(AGE_CLASSES, function(ac) counts[[ac]][["I_H"]], integer(1)))
  event_values <- c(
    births_female=record$births_female,
    deaths_natural_calf=record$deaths_natural$calf,
    deaths_natural_heifer=record$deaths_natural$heifer,
    deaths_natural_cow=record$deaths_natural$cow,
    deaths_blv_calf=record$dead_blv$calf,
    deaths_blv_heifer=record$dead_blv$heifer,
    deaths_blv_cow=record$dead_blv$cow,
    new_infected_calf=record$new_infected$calf,
    new_infected_heifer=record$new_infected$heifer,
    new_infected_cow=record$new_infected$cow,
    cows_culled=record$culled_cow,
    external_cow_purchases=record$external_purchases,
    surplus_heifers_sold=record$heifers_sold_surplus
  )
  as.data.frame(as.list(c(
    year=as.integer(record$year),
    population=population,
    infected_low=low,
    infected_high=high,
    prevalence_total=(low+high)/population,
    prevalence_low=low/population,
    prevalence_high=high/population,
    age_values,
    event_values
  )), check.names=FALSE)
}

draw_scenario_plot <- function(trajectories, scenario, years) {
  par(mfrow=c(2,2), mar=c(4,4,3,1), oma=c(0,0,3,0))
  plot_prevalence <- function(column, label, color) {
    values <- sapply(trajectories, function(x) x[[column]])
    plot(range(years), c(0,100), type="n", xlab="Year", ylab="Prevalence (%)",
         main=label)
    matlines(years, 100*values, lty=1,
             col=adjustcolor(color, alpha.f=.25))
    lines(years, 100*apply(values, 1, median), col="black", lwd=2)
  }
  plot_prevalence("prevalence_total", "All infected (I_L + I_H)", "#246A73")
  plot_prevalence("prevalence_low", "Low infectiousness (I_L)", "#2C8C72")
  plot_prevalence("prevalence_high", "High infectiousness (I_H)", "#BD5B50")
  active <- sapply(trajectories, function(x)
    x$infected_low + x$infected_high > 0)
  plot(years, rowSums(active), type="l", lwd=2, col="#334E68",
       ylim=c(0,length(trajectories)), xlab="Year",
       ylab="Runs with infection", main="Infection persistence")
  mtext(sprintf("%s | %s cattle | one %s %s | %d runs",
                scenario$scenario_id, scenario$herd_size,
                scenario$index_age, scenario$index_state,
                length(trajectories)), outer=TRUE, cex=1.1)
}

summary_file <- file.path(output_dir, "introduction_summary.csv")
summary_table <- if (file.exists(summary_file))
  read.csv(summary_file, stringsAsFactors=FALSE) else NULL
years <- 0:horizon_years

for (i in selected_indices) {
  scenario <- scenario_grid[i, ]
  size <- scenario$herd_size
  composition <- switch(as.character(size),
                        "100"=c(calf=12L, heifer=38L, cow=50L),
                        "200"=c(calf=25L, heifer=75L, cow=100L),
                        "800"=c(calf=100L, heifer=300L, cow=400L))
  params <- make_params(
    herd_size=size,
    horizon_years=horizon_years,
    init_calves=unname(composition["calf"]),
    init_heifers=unname(composition["heifer"]),
    init_cows=unname(composition["cow"]),
    init_prevalence=0,
    initial_infected_age=scenario$index_age,
    initial_infected_state=scenario$index_state,
    corrective_management=NULL,
    testing_enabled=FALSE,
    vaccination_enabled=FALSE,
    surveillance_interval_years=NULL,
    external_replacements_tested=TRUE
  )

  scenario_dir <- file.path(output_dir, scenario$scenario_id)
  dir.create(scenario_dir, showWarnings=FALSE)
  trajectories <- vector("list", n_runs)

  for (j in seq_len(n_runs)) {
    run_file <- file.path(scenario_dir, sprintf("run_%03d.rds", j))
    existing_run_file <- run_file
    if (!file.exists(existing_run_file) && !is.null(cache_dir)) {
      cached_file <- file.path(cache_dir, scenario$scenario_id,
                               sprintf("run_%03d.rds", j))
      if (file.exists(cached_file)) existing_run_file <- cached_file
    }
    if (file.exists(existing_run_file)) {
      saved <- readRDS(existing_run_file)
      if (!identical(saved$params, params) ||
          !identical(saved$model_md5, model_md5) ||
          !identical(saved$seed, replicate_seeds[j]))
        stop("Saved run does not match the current model or parameters: ", existing_run_file)
      trajectory <- saved$trajectory
      if (!identical(existing_run_file, run_file))
        saveRDS(saved, run_file)
      message(scenario$scenario_id, " (", i, "/18), run ", j, "/", n_runs,
              " loaded")
    } else {
      elapsed <- system.time({
        run <- run_one(params, replicate_seeds[j])  #run the simulation
        trajectory <- do.call(rbind, lapply(run$year_records, extract_year))
      })[["elapsed"]]
      # The final flush may repeat the horizon snapshot; keep its final state.
      trajectory <- trajectory[!duplicated(trajectory$year, fromLast=TRUE), ]
      if (length(trajectory$year) != length(years) ||
          any(trajectory$year != years))
        stop("Unexpected yearly records for ", scenario$scenario_id, " run ", j)
      if (trajectory$infected_low[1] + trajectory$infected_high[1] != 1L)
        stop("The run must start with exactly one infected animal.")
      saveRDS(list(params=params, model_md5=model_md5,
                   seed=replicate_seeds[j], trajectory=trajectory), run_file)
      rm(run)
      message(scenario$scenario_id, " (", i, "/18), run ", j, "/", n_runs,
              " completed in ", round(elapsed, 1), " s")
    }
    trajectory$replicate <- j
    trajectory$seed <- replicate_seeds[j]
    trajectories[[j]] <- trajectory
  }

  yearly <- do.call(rbind, trajectories)
  final <- yearly[yearly$year == horizon_years, ]
  active <- final$infected_low + final$infected_high > 0
  summary_row <- data.frame(
    scenario_id=scenario$scenario_id,
    herd_size=size,
    index_age=scenario$index_age,
    index_state=scenario$index_state,
    horizon_years=horizon_years,
    n_runs=n_runs,
    active_runs_at_horizon=sum(active),
    extinct_runs_at_horizon=sum(!active),
    mean_prevalence_all_runs=mean(final$prevalence_total),
    median_prevalence_all_runs=median(final$prevalence_total),
    mean_prevalence_active_runs=if (any(active))
      mean(final$prevalence_total[active]) else NA_real_,
    median_final_calves=median(final$calf_population),
    median_final_heifers=median(final$heifer_population),
    median_final_cows=median(final$cow_population),
    mean_final_calves=mean(final$calf_population),
    mean_final_heifers=mean(final$heifer_population),
    mean_final_cows=mean(final$cow_population),
    mean_external_cow_purchases_per_year=
      mean(yearly$external_cow_purchases),
    mean_surplus_heifers_sold_per_year=
      mean(yearly$surplus_heifers_sold)
  )
  if (!is.null(summary_table))
    summary_table <- summary_table[summary_table$scenario_id != scenario$scenario_id, ]
  summary_table <- rbind(summary_table, summary_row)
  summary_table <- summary_table[order(match(summary_table$scenario_id,
                                              scenario_grid$scenario_id)), ]
  saveRDS(list(scenario=scenario, params=params,
               replicate_seeds=replicate_seeds,
               yearly=yearly, summary=summary_row),
          file.path(output_dir, paste0(scenario$scenario_id, ".rds")))
  write.csv(yearly, file.path(output_dir, paste0(scenario$scenario_id, "_yearly.csv")),
            row.names=FALSE)
  write.csv(summary_table, summary_file, row.names=FALSE)

  png(file.path(output_dir, paste0(scenario$scenario_id, ".png")),
      width=1600, height=1200, res=150)
  draw_scenario_plot(trajectories, scenario, years)
  dev.off()
  pdf(file.path(output_dir, paste0(scenario$scenario_id, ".pdf")),
      width=11, height=8.5)
  draw_scenario_plot(trajectories, scenario, years)
  dev.off()
  print(summary_row)
}

options(
  blv_burnin_n_runs=100L,
  blv_burnin_output_dir=file.path("results", "introduction_burnin_replacement_100"),
  blv_burnin_cache_dir=file.path("results", "introduction_burnin_replacement_pilot_20")
)
source("burnin_initialization.R")
