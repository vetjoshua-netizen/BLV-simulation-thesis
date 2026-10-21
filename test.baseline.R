source("simulation.R")


# Animal creation; addition/removal; age and disease indexes.

h <- Herd$new()

a <- h$new_animal(
  0,
  "calf",
  "S"
)

stopifnot(
  a$alive,
  h$n_age("calf") == 1L,
  herd_indexes_consistent(h)
)

h$set_state(
  a,
  "I_L"
)

h$set_age_class(
  a,
  "heifer"
)

stopifnot(
  h$n_state("heifer", "I_L") == 1L,
  herd_indexes_consistent(h)
)

h$remove_animal(a)

stopifnot(
  !a$alive,
  h$n_age("heifer") == 0L
)


# Scheduled events are chronological and FIFO on ties.

q <- EventHeap$new()

q$push(
  2,
  "later"
)

q$push(
  1,
  "first"
)

q$push(
  1,
  "second"
)

got <- vapply(
  replicate(
    3,
    q$pop(),
    simplify = FALSE
  ),
  function(x) x$kind,
  character(1)
)

stopifnot(
  identical(
    got,
    c("first", "second", "later")
  )
)


# Birth / vertical transmission copy the dam's state under probability one.

p <- make_params(
  init_calves = 0,
  init_heifers = 0,
  init_cows = 1,
  init_prevalence = 1,
  init_LH_ratio = 1,
  prop_female = 1,
  vertical_transmission_rate = 1,
  birth_rate = 0,
  cow_culling_rate = 0,
  seed = 1
)

s <- Simulation$new(p)

do_birth_from_dam(
  s,
  random_member(s, "cow")
)

stopifnot(
  s$herd$n_state("calf", "I_L") == 1L,
  herd_indexes_consistent(s$herd)
)


# Natural death and one-ID background culling (the old sample(ids, k) hazard).

s <- Simulation$new(
  make_params(
    init_calves = 1,
    init_heifers = 0,
    init_cows = 0,
    birth_rate = 0,
    seed = 2
  )
)

execute_stochastic(
  s,
  list(
    list(
      rate = 1,
      kind = "death_natural",
      ac = "calf"
    )
  ),
  1
)

stopifnot(
  s$herd$n_age("calf") == 0L
)

s <- Simulation$new(
  make_params(
    init_calves = 0,
    init_heifers = 0,
    init_cows = 1,
    birth_rate = 0,
    seed = 3
  )
)

execute_stochastic(
  s,
  list(
    list(
      rate = 1,
      kind = "cull_background",
      ac = "cow"
    )
  ),
  1
)

stopifnot(
  s$herd$n_age("cow") == 1L,
  s$acc$external_purchases == 1L
)


# Horizontal infection preserves the current population-mix rule.

s <- Simulation$new(
  make_params(
    init_calves = 0,
    init_heifers = 0,
    init_cows = 0,
    seed = 4
  )
)

s$herd$new_animal(
  0,
  "heifer",
  "S"
)

s$herd$new_animal(
  0,
  "heifer",
  "I_L"
)

do_infection(
  s,
  "heifer"
)

stopifnot(
  s$herd$n_state("heifer", "I_L") == 2L
)


# Ageing and current oldest-heifer replacement rule.

s <- Simulation$new(
  make_params(
    init_calves = 0,
    init_heifers = 0,
    init_cows = 0,
    seed = 5
  )
)

x <- s$herd$new_animal(
  -.5,
  "calf",
  "S"
)

handle_scheduled(
  s,
  "age_to_heifer",
  list(id = x$id)
)

stopifnot(
  x$age_class == "heifer"
)

young <- s$herd$new_animal(
  -1.2,
  "heifer",
  "S"
)

old <- s$herd$new_animal(
  -1.8,
  "heifer",
  "S"
)

replace_cow(
  s,
  "normal_cull"
)

stopifnot(
  old$age_class == "cow",
  young$age_class == "heifer"
)


# Horizon, reproducibility, actual seed recording, and zero-population edge case.

p <- make_params(
  herd_size = 5,
  horizon_years = .2,
  seed = 10
)

r1 <- run_one(
  p,
  99
)

r2 <- run_one(
  p,
  99
)

stopifnot(
  r1$seed == 99L,
  identical(r1, r2),
  tail(r1$year_records, 1)[[1]]$time == .2
)

e <- run_one(
  make_params(
    init_calves = 0,
    init_heifers = 0,
    init_cows = 0,
    horizon_years = .1,
    seed = 11
  )
)

stopifnot(
  sum(unlist(e$final_population)) == 0L
)

cat("All baseline tests passed.\n")
















#Load model
source("simulation.R")



#1.Small test helpers
pass <- function(name) {
  cat("[PASS]", name, "\n")
}

fail <- function(name, e) {
  cat("[FAIL]", name, "\n")
  cat("       ", conditionMessage(e), "\n")
}

run_test <- function(name, expr) {
  tryCatch(
    {
      force(expr)
      pass(name)
      TRUE
    },
    error = function(e) {
      fail(name, e)
      FALSE
    }
  )
}



# 2. Common test parameters
make_test_params <- function(...) {
  make_params(
    herd_size = 100,
    horizon_years = 1,
    init_calves = 20,
    init_heifers = 30,
    init_cows = 50,
    seed = 1,
    ...
  )
}



# TEST 1
# Initial population and state indexes

test_initialisation <- run_test(
  "Initialisation and index consistency",
  {
    p <- make_test_params()
    
    sim <- Simulation$new(p, seed = 1)
    
    # Total herd size
    total <- sum(
      vapply(
        AGE_CLASSES,
        sim$herd$n_age,
        integer(1)
      )
    )
    
    stopifnot(total == 100)
    
    # Expected age structure
    stopifnot(sim$herd$n_age("calf") == 20)
    stopifnot(sim$herd$n_age("heifer") == 30)
    stopifnot(sim$herd$n_age("cow") == 50)
    
    # Initial BLV prevalence
    stopifnot(
      sim$herd$n_state("cow", "I_L") == 5
    )
    
    stopifnot(
      sim$herd$n_state("cow", "I_H") == 5
    )
    
    # Index consistency
    stopifnot(
      herd_indexes_consistent(sim$herd)
    )
  }
)



# TEST 2
# Singleton sampling

test_singleton <- run_test(
  "Singleton sampling",
  {
    s <- IndexedSet$new()
    
    s$add(100L)
    
    x <- s$random_one()
    
    stopifnot(identical(x, 100L))
  }
)



# TEST 3
# Stochastic rates are actually generated (This catches the bug where stochastic_rates() returned list()


test_stochastic_rates <- run_test(
  "Stochastic event rates are generated",
  {
    p <- make_test_params(
      birth_rate = 1,
      calf_mortality_rate = 0.05,
      heifer_mortality_rate = 0.03,
      cow_culling_rate = 0.25
    )
    
    sim <- Simulation$new(p, seed = 1)
    
    rates <- stochastic_rates(sim)
    
    # Must not be empty
    stopifnot(length(rates) > 0)
    
    # Extract event names
    kinds <- vapply(
      rates,
      function(x) x$kind,
      character(1)
    )
    
    # Birth must exist
    stopifnot("birth" %in% kinds)
    
    # Birth rate = 50 cows * 1
    birth_index <- which(kinds == "birth")[1]
    
    stopifnot(
      rates[[birth_index]]$rate == 50
    )
  }
)



# TEST 4
# Direct birth
test_birth <- run_test(
  "Birth creates calves",
  {
    p <- make_test_params(
      prop_female = 1,
      calf_mortality_rate = 0,
      heifer_mortality_rate = 0,
      cow_culling_rate = 0
    )
    
    sim <- Simulation$new(p, seed = 1)
    
    before <- sim$herd$n_age("calf")
    
    # Force 10 births
    for (i in seq_len(10)) {
      do_birth(sim)
    }
    
    after <- sim$herd$n_age("calf")
    
    stopifnot(
      after == before + 10
    )
  }
)



# TEST 5
# Vertical transmission
test_vertical_transmission <- run_test(
  "Vertical transmission",
  {
    p <- make_params(
      herd_size = 1,
      horizon_years = 1,
      init_calves = 0,
      init_heifers = 0,
      init_cows = 1,
      init_prevalence = 1,
      init_LH_ratio = 1,
      prop_female = 1,
      vertical_transmission_rate = 1,
      calf_mortality_rate = 0,
      heifer_mortality_rate = 0,
      cow_culling_rate = 0,
      seed = 1
    )
    
    sim <- Simulation$new(p, seed = 1)
    
    do_birth(sim)
    
    # Infected dam + transmission probability 1
    # -> infected calf
    infected_calf <-
      sim$herd$n_state("calf", "I_L") +
      sim$herd$n_state("calf", "I_H")
    
    stopifnot(
      infected_calf == 1
    )
  }
)



# TEST 6
# Horizontal infection
test_infection <- run_test(
  "Horizontal infection",
  {
    p <- make_test_params(
      calf_mortality_rate = 0,
      heifer_mortality_rate = 0,
      cow_culling_rate = 0
    )
    
    sim <- Simulation$new(p, seed = 1)
    
    before <-
      sim$herd$n_state("cow", "I_L") +
      sim$herd$n_state("cow", "I_H")
    
    # Force one infection event
    do_infection(sim, "cow")
    
    after <-
      sim$herd$n_state("cow", "I_L") +
      sim$herd$n_state("cow", "I_H")
    
    stopifnot(
      after == before + 1
    )
    
    stopifnot(
      herd_indexes_consistent(sim$herd)
    )
  }
)



# TEST 7
# Natural death
test_natural_death <- run_test(
  "Natural death removes animal",
  {
    p <- make_params(
      herd_size = 1,
      horizon_years = 1,
      init_calves = 1,
      init_heifers = 0,
      init_cows = 0,
      calf_mortality_rate = 1,
      heifer_mortality_rate = 0,
      cow_culling_rate = 0,
      birth_rate = 0,
      seed = 1
    )
    
    sim <- Simulation$new(p, seed = 1)
    
    before <- sim$herd$n_age("calf")
    
    # Force a natural death event
    rates <- list(
      list(
        rate = 1,
        kind = "death_natural",
        ac = "calf"
      )
    )
    
    execute_stochastic(
      sim,
      rates,
      total = 1
    )
    
    after <- sim$herd$n_age("calf")
    
    stopifnot(
      before == 1
    )
    
    stopifnot(
      after == 0
    )
    
    stopifnot(
      herd_indexes_consistent(sim$herd)
    )
  }
)



# TEST 8
# Cow culling + replacement
test_culling_replacement <- run_test(
  "Cow culling and replacement",
  {
    p <- make_params(
      herd_size = 2,
      horizon_years = 1,
      init_calves = 0,
      init_heifers = 1,
      init_cows = 1,
      cow_culling_rate = 1,
      birth_rate = 0,
      calf_mortality_rate = 0,
      heifer_mortality_rate = 0,
      cap_herd_at_target = FALSE,
      seed = 1
    )
    
    sim <- Simulation$new(p, seed = 1)
    
    before_cows <- sim$herd$n_age("cow")
    before_heifers <- sim$herd$n_age("heifer")
    
    rates <- list(
      list(
        rate = 1,
        kind = "cull_background",
        ac = "cow"
      )
    )
    
    execute_stochastic(
      sim,
      rates,
      total = 1
    )
    
    after_cows <- sim$herd$n_age("cow")
    after_heifers <- sim$herd$n_age("heifer")
    
    # One cow removed, one heifer becomes cow
    stopifnot(before_cows == 1)
    stopifnot(before_heifers == 1)
    
    stopifnot(after_cows == 1)
    stopifnot(after_heifers == 0)
    
    stopifnot(
      sim$acc$culled_cow == 1
    )
    
    stopifnot(
      herd_indexes_consistent(sim$herd)
    )
  }
)



# TEST 9
# BLV death + replacement
test_blv_death <- run_test(
  "BLV death and replacement",
  {
    p <- make_params(
      herd_size = 2,
      horizon_years = 1,
      init_calves = 0,
      init_heifers = 1,
      init_cows = 1,
      init_prevalence = 1,
      init_LH_ratio = 1,
      cow_culling_rate = 0,
      birth_rate = 0,
      calf_mortality_rate = 0,
      heifer_mortality_rate = 0,
      cap_herd_at_target = FALSE,
      seed = 1
    )
    
    sim <- Simulation$new(p, seed = 1)
    
    before_infected <-
      sim$herd$n_state("cow", "I_L") +
      sim$herd$n_state("cow", "I_H")
    
    rates <- list(
      list(
        rate = 1,
        kind = "death_blv",
        ac = "cow"
      )
    )
    
    execute_stochastic(
      sim,
      rates,
      total = 1
    )
    
    after_infected <-
      sim$herd$n_state("cow", "I_L") +
      sim$herd$n_state("cow", "I_H")
    
    stopifnot(
      before_infected == 1
    )
    
    stopifnot(
      after_infected == 0
    )
    
    stopifnot(
      sim$acc$dead_blv$cow == 1
    )
    
    stopifnot(
      herd_indexes_consistent(sim$herd)
    )
  }
)



# TEST 10
# Wane: R calf -> S
test_waning <- run_test(
  "Maternal waning R -> S",
  {
    p <- make_params(
      herd_size = 1,
      horizon_years = 1,
      init_calves = 0,
      init_heifers = 0,
      init_cows = 0,
      maternal_waning_rate = 1,
      seed = 1
    )
    
    sim <- Simulation$new(p, seed = 1)
    
    # Create an R calf manually
    a <- sim$herd$new_animal(
      birth_time = 0,
      age_class = "calf",
      state = "R"
    )
    
    stopifnot(
      sim$herd$n_state("calf", "R") == 1
    )
    
    rates <- list(
      list(
        rate = 1,
        kind = "wane",
        ac = "calf"
      )
    )
    
    execute_stochastic(
      sim,
      rates,
      total = 1
    )
    
    stopifnot(
      sim$herd$n_state("calf", "R") == 0
    )
    
    stopifnot(
      sim$herd$n_state("calf", "S") == 1
    )
    
    stopifnot(
      herd_indexes_consistent(sim$herd)
    )
  }
)



# TEST 11
# Aging: calf -> heifer -> cow
test_aging <- run_test(
  "Age transitions",
  {
    p <- make_params(
      herd_size = 1,
      horizon_years = 5,
      init_calves = 0,
      init_heifers = 0,
      init_cows = 0,
      cap_herd_at_target = FALSE,
      seed = 1
    )
    
    sim <- Simulation$new(p, seed = 1)
    
    a <- sim$herd$new_animal(
      birth_time = 0,
      age_class = "calf",
      state = "S"
    )
    
    # Calf -> heifer
    handle_scheduled(
      sim,
      "age_to_heifer",
      list(id = a$id)
    )
    
    stopifnot(
      sim$herd$n_age("calf") == 0
    )
    
    stopifnot(
      sim$herd$n_age("heifer") == 1
    )
    
    # Heifer -> cow
    handle_scheduled(
      sim,
      "age_to_cow",
      list(id = a$id)
    )
    
    stopifnot(
      sim$herd$n_age("heifer") == 0
    )
    
    stopifnot(
      sim$herd$n_age("cow") == 1
    )
    
    stopifnot(
      herd_indexes_consistent(sim$herd)
    )
  }
)



# TEST 12
# Reproducibility
test_reproducibility <- run_test(
  "Reproducibility with identical seed",
  {
    p <- make_params(
      herd_size = 50,
      horizon_years = 2,
      seed = 123
    )
    
    r1 <- run_one(
      p,
      seed = 99
    )
    
    r2 <- run_one(
      p,
      seed = 99
    )
    
    # Same seed should give exactly same result
    stopifnot(
      identical(r1, r2)
    )
    
    stopifnot(
      r1$seed == 99L
    )
    
    stopifnot(
      r2$seed == 99L
    )
  }
)



# TEST 13
# Reproducible replicate seeds generated from one master seed
test_multiple_runs <- run_test(
  "Master seed generates unique and reproducible replicate seeds",
  {
    p <- make_params(
      herd_size = 50,
      horizon_years = 2,
      seed = 1
    )
    
    out <- run_strategy(
      p,
      n_runs = 5,
      master_seed = 12345
    )
    out_again <- run_strategy(
      p,
      n_runs = 5,
      master_seed = 12345
    )
    
    seeds <- vapply(
      out$runs,
      function(x) x$seed,
      integer(1)
    )
    
    stopifnot(
      length(unique(seeds)) == 5L,
      identical(seeds, out$replicate_seeds),
      identical(seeds, out_again$replicate_seeds),
      identical(out$master_seed, 12345L)
    )
  }
)



# TEST 14
# Full simulation smoke test
test_full_simulation <- run_test(
  "Full simulation smoke test",
  {
    p <- make_params(
      herd_size = 100,
      horizon_years = 5,
      seed = 1
    )
    
    out <- run_one(
      p,
      seed = 1
    )
    
    # Result fields exist
    stopifnot(
      all(
        c(
          "seed",
          "params_seed",
          "year_records",
          "eradicated",
          "eradication_time",
          "ever_infected",
          "final_population"
        ) %in%
          names(out)
      )
    )
    
    # Expected number of yearly records:
    # year 0, 1, 2, 3, 4, final year 5
    stopifnot(
      length(out$year_records) == 6
    )
    
    # All years should be non-decreasing
    years <- vapply(
      out$year_records,
      function(x) x$year,
      numeric(1)
    )
    
    stopifnot(
      all(diff(years) >= 0)
    )
    
    # Check every record
    for (rec in out$year_records) {
      
      # Age-class counts
      age_counts <- vapply(
        AGE_CLASSES,
        function(ac) rec$population_by_age[[ac]],
        integer(1)
      )
      
      stopifnot(
        all(age_counts >= 0)
      )
      
      # State counts
      for (ac in AGE_CLASSES) {
        
        state_counts <- vapply(
          STATES,
          function(s)
            rec$population_by_age_state[[ac]][[s]],
          integer(1)
        )
        
        # State counts must sum to age-class count
        stopifnot(
          sum(state_counts) == age_counts[[ac]]
        )
        
        # No negative state counts
        stopifnot(
          all(state_counts >= 0)
        )
      }
      
      # Total population must be non-negative
      stopifnot(
        sum(age_counts) >= 0
      )
    }
  }
)


# TEST 15
# Larger simulation: 1000 animals, 40 years, 10 runs
test_large_simulation <- run_test(
  "Large simulation: 1000 animals, 40 years, 10 runs",
  {
    p <- make_params(
      herd_size = 1000,
      horizon_years = 40,
      seed = 1
    )
    
    out <- run_strategy(
      p,
      n_runs = 10,
      master_seed = 12345
    )
    
    # 10 replicates
    stopifnot(
      length(out$runs) == 10
    )
    
    # Seeds are unique and reproducible from the master seed.
    seeds <- vapply(
      out$runs,
      function(x) x$seed,
      integer(1)
    )
    
    stopifnot(
      length(unique(seeds)) == 10L,
      identical(seeds, out$replicate_seeds),
      identical(out$master_seed, 12345L)
    )
    
    # Check all replicates
    for (run in out$runs) {
      
      stopifnot(
        length(run$year_records) == 41
      )
      
      for (rec in run$year_records) {
        
        age_counts <- vapply(
          AGE_CLASSES,
          function(ac)
            rec$population_by_age[[ac]],
          integer(1)
        )
        
        stopifnot(
          all(age_counts >= 0)
        )
        
        for (ac in AGE_CLASSES) {
          
          state_counts <- vapply(
            STATES,
            function(s)
              rec$population_by_age_state[[ac]][[s]],
            integer(1)
          )
          
          stopifnot(
            sum(state_counts) == age_counts[[ac]]
          )
          
          stopifnot(
            all(state_counts >= 0)
          )
        }
      }
    }
  }
)


