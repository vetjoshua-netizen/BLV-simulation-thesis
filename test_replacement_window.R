# Small checks for the 22–24 month replacement-heifer policy.
source("simulation.R")

make_empty_sim <- function() {
  p <- make_params(
    herd_size=1L,
    init_calves=0L,
    init_heifers=0L,
    init_cows=0L,
    target_cows=1L,
    horizon_years=1,
    female_calf_entry_rate=0,
    cow_culling_rate=0,
    calf_mortality_rate=0,
    heifer_mortality_rate=0,
    blv_mortality_rate=c(calf=0, heifer=0, cow=0),
    surveillance_interval_years=NULL
  )
  Simulation$new(p, seed=1L)
}

# Only the 23-month heifer is eligible; the 21-month heifer must remain a heifer.
sim <- make_empty_sim()
young <- sim$herd$new_animal(-21/12, "heifer", "S")
eligible <- sim$herd$new_animal(-23/12, "heifer", "S")
stopifnot(identical(eligible_replacement_heifer_ids(sim), eligible$id))
invisible(replace_cow(sim, "test"))
stopifnot(sim$herd$n_age("cow") == 1L,
          sim$herd$get_animal(young$id)$age_class == "heifer",
          sim$herd$get_animal(eligible$id)$age_class == "cow",
          sim$acc$external_purchases == 0L)

# With no eligible heifer, the missing cow is bought externally.
sim <- make_empty_sim()
invisible(sim$herd$new_animal(-21/12, "heifer", "S"))
invisible(replace_cow(sim, "test"))
stopifnot(sim$herd$n_age("cow") == 1L,
          sim$herd$n_age("heifer") == 1L,
          sim$acc$external_purchases == 1L)

cat("Replacement-window checks passed.\n")
