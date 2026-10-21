# BLV individual-based simulation. No plotting dependency.
# Hybrid model: exact SSA stochastic events + chronological scheduled events.
# TODOs preserved: confirm receiving-age-class I_L/I_H assignment, dam-state
# vertical assignment, immediate oldest-heifer replacement, and non-strict cap.
# "eradicated" means undetected under surveillance, not infection-free.
# All epidemiological defaults below are PLACEHOLDERS.  The intervention
# defaults implement the agreed hypothetical vaccination scenario: 70% vaccine
# take, all-or-nothing protection, one lifetime dose, no waning and no effect
# on animals already infected with BLV.

if (!requireNamespace("R6", quietly=TRUE)) 
  stop("Package 'R6' is required.")

#1 constants
STATES <- c("S","R","I_L","I_H");
AGE_CLASSES <- c("calf","heifer","cow")
default_initial_composition <- function(n) { 
  x <- round(c(calf=.2,heifer=.3)*n); 
  list(calves=x[1],heifers=x[2],cows=n-sum(x)) 
}

#2 parameters
make_params <- function(
    herd_size=100L,
    horizon_years=40,
    init_calves=NULL,
    init_heifers=NULL,
    init_cows=NULL,
    init_prevalence=.2,
    init_LH_ratio=.5,
    birth_rate=1,
    prop_female=.5,
    calf_mortality_rate=.05,
    heifer_mortality_rate=.03,
    cow_culling_rate=.25,
    calf_to_heifer_age=.5,
    heifer_to_cow_age=2,
    vertical_transmission_rate=.1,
    maternal_waning_rate=NULL,
    beta_L=c(calf=0,heifer=.3,cow=.5),
    beta_H=c(calf=0,heifer=.6,cow=1),
    blv_mortality_rate=c(calf=0,heifer=.01,cow=.05),
    external_replacements_tested=TRUE,
    external_prevalence=.15,
    external_LH_ratio=.5,
    surveillance_interval_years=.25,
    surveillance_sensitivity=.95,
    n_clean_samples_required=3L,
    cap_herd_at_target=TRUE,

    # Corrective management reduces the horizontal force of infection by the
    # stated efficacy.  NULL means no corrective management.
    corrective_management=NULL,

    corrective_management_efficacy=c(low=.15, medium=.35, high=.60),

    # Test-and-remove intervention.
    testing_enabled=FALSE,

    testing_interval_years=.5,
    testing_target_group="all",
    testing_remove_policy="all_positive",
    testing_sensitivity=.95,

    # Preventive vaccination intervention.  A vaccinated S animal has a
    # vaccine_efficacy probability of becoming R; I_L and I_H are unchanged.
    vaccination_enabled=FALSE,

    vaccination_groups=c("heifer", "cow"),
    vaccination_coverage=1,
    vaccine_efficacy=.70,
    vaccination_initial_campaign=TRUE,
    seed=12345L
) {d<-default_initial_composition(herd_size); 
if(is.null(init_calves))
  init_calves<-d$calves; 
if(is.null(init_heifers))
  init_heifers<-d$heifers; 
if(is.null(init_cows))
  init_cows<-d$cows; 
as.list(environment())
}

#3 parameter validation
validate_params <- function(p) {
  if(
    any(
      c(
        p$horizon_years,p$init_calves,p$init_heifers,p$init_cows
      )<0
    )
  )
    stop("Invalid negative size/horizon.")

  if (!is.null(p$corrective_management) &&
      !(p$corrective_management %in% names(p$corrective_management_efficacy)))
    stop("corrective_management must be NULL, 'low', 'medium', or 'high'.")

  if (!is.numeric(p$vaccination_coverage) || length(p$vaccination_coverage) != 1L ||
      is.na(p$vaccination_coverage) || p$vaccination_coverage < 0 || p$vaccination_coverage > 1)
    stop("vaccination_coverage must be a single value from 0 to 1.")

  if (!is.numeric(p$vaccine_efficacy) || length(p$vaccine_efficacy) != 1L ||
      is.na(p$vaccine_efficacy) || p$vaccine_efficacy < 0 || p$vaccine_efficacy > 1)
    stop("vaccine_efficacy must be a single value from 0 to 1.")

  if (any(!(p$vaccination_groups %in% AGE_CLASSES)))
    stop("vaccination_groups must contain only calf, heifer, and/or cow.")

  if (!is.numeric(p$testing_interval_years) || length(p$testing_interval_years) != 1L ||
      p$testing_interval_years <= 0)
    stop("testing_interval_years must be positive.")

  if (!(p$testing_target_group %in% c("all", AGE_CLASSES)))
    stop("testing_target_group must be 'all', 'calf', 'heifer', or 'cow'.")

  if (!(p$testing_remove_policy %in% c("all_positive", "high_only")))
    stop("testing_remove_policy must be 'all_positive' or 'high_only'.")

  if (!is.numeric(p$testing_sensitivity) || length(p$testing_sensitivity) != 1L ||
      is.na(p$testing_sensitivity) || p$testing_sensitivity < 0 || p$testing_sensitivity > 1)
    stop("testing_sensitivity must be a single value from 0 to 1.")
}

# Backwards-compatible name used by older baseline scripts.
validate_baseline_params <- validate_params

#4 animal class
Animal<-R6::R6Class(
  "Animal",
  public=list(
    id=NULL,
    birth_time=NULL,
    age_class=NULL,
    state=NULL,
    r_source=NULL,
    vaccinated=FALSE,
    alive=TRUE,
    initialize=function(
    id,
    birth_time,
    age_class,
    state,
    r_source=NULL
    )
    {
      self$id<-id;
      self$birth_time<-birth_time;
      self$age_class<-age_class;
      self$state<-state;self$r_source<-r_source
      }
    )
  )

#5 index function
IndexedSet<-R6::R6Class(
  "IndexedSet",
  public=list(
    ids=NULL,
    initialize=function()
      self$ids<-integer(),
    add=function(id){
      if(!(id%in%self$ids))
        self$ids<-c(self$ids,id)
      },
    discard=function(id)
      self$ids<-self$ids[self$ids!=id],
    size=function()length(self$ids),
    as_vector=function()
      self$ids,
    random_one=function(){
      if(!length(self$ids))
        stop("empty set");
      self$ids[sample.int(length(self$ids),1L)]
      }
    )
  )

#6 herd class
Herd <- R6::R6Class(
  "Herd",
  public=list(
    
    animals=NULL,
    members=NULL,
    by_age=NULL,
    id_counter=0L,
    
    initialize=function() {
      self$animals <- new.env(hash=TRUE, parent=emptyenv())
      
      self$members <- setNames(
        lapply(
          AGE_CLASSES,
          function(a)
            setNames(
              lapply(STATES, function(s) IndexedSet$new()),
              STATES
            )
        ),
        AGE_CLASSES
      )
      
      self$by_age <- setNames(
        lapply(AGE_CLASSES, function(a) IndexedSet$new()),
        AGE_CLASSES
      )
    },
    
    add_animal=function(a) {
      assign(
        as.character(a$id),
        a,
        self$animals
      )
      
      self$members[[a$age_class]][[a$state]]$add(a$id)
      self$by_age[[a$age_class]]$add(a$id)
    },
    
    new_animal=function(
    birth_time,
    age_class,
    state,
    r_source=NULL
    ) {
      self$id_counter <- self$id_counter + 1L
      
      a <- Animal$new(
        self$id_counter,
        birth_time,
        age_class,
        state,
        r_source
      )
      
      self$add_animal(a)
      a
    },
    
    get_animal=function(id) {
      k <- as.character(id)
      
      if (!exists(
        k,
        self$animals,
        inherits=FALSE
      ))
        return(NULL)
      
      get(
        k,
        self$animals,
        inherits=FALSE
      )
    },
    
    remove_animal=function(a) {
      if (a$alive) {
        self$members[[a$age_class]][[a$state]]$discard(a$id)
        self$by_age[[a$age_class]]$discard(a$id)
        a$alive <- FALSE
      }
    },
    
    set_state=function(
    a,
    state,
    r_source=NULL
    ) {
      self$members[[a$age_class]][[a$state]]$discard(a$id)
      a$state <- state
      a$r_source <- r_source
      self$members[[a$age_class]][[state]]$add(a$id)
    },
    
    set_age_class=function(a, age_class) {
      self$members[[a$age_class]][[a$state]]$discard(a$id)
      self$by_age[[a$age_class]]$discard(a$id)
      
      a$age_class <- age_class
      
      self$members[[age_class]][[a$state]]$add(a$id)
      self$by_age[[age_class]]$add(a$id)
    },
    
    n_age=function(ac)
      self$by_age[[ac]]$size(),
    
    n_state=function(ac, s)
      self$members[[ac]][[s]]$size()
  )
)

#7 eventheap
EventHeap <- R6::R6Class(
  "EventHeap",
  public=list(
    events=NULL,
    sequence=0L,
    
    initialize=function()
      self$events <- list(),
    
    size=function()
      length(self$events),
    
    peek_time=function()
      if (length(self$events))
        self$events[[1]]$time
    else
      Inf,
    
    push=function(time, kind, meta=list()) {
      self$sequence <- self$sequence + 1L
      
      self$events[[length(self$events)+1L]] <- list(
        time=time,
        kind=kind,
        meta=meta,
        sequence=self$sequence
      )
      
      self$events <- self$events[
        order(
          vapply(
            self$events,
            function(x) x$time,
            numeric(1)
          ),
          vapply(
            self$events,
            function(x) x$sequence,
            integer(1)
          )
        )
      ]
    },
    
    pop=function() {
      if (!length(self$events))
        return(NULL)
      
      x <- self$events[[1]]
      self$events <- self$events[-1]
      x
    }
  )
)

#8 event accumulators
blank_acc <- function() {
  list(
    births_female = 0L,
    births_male = 0L,
    
    deaths_natural = setNames(
      as.list(rep(0L, 3)),
      AGE_CLASSES
    ),
    
    dead_blv = setNames(
      as.list(rep(0L, 3)),
      AGE_CLASSES
    ),
    
    new_infected = setNames(
      as.list(rep(0L, 3)),
      AGE_CLASSES
    ),
    
    external_purchases = 0L,
    external_purchases_for_test_positive = 0L,
    culled_cow = 0L,
    heifers_sold_surplus = 0L,
    
    tested_by_age = setNames(
      as.list(rep(0L, 3)),
      AGE_CLASSES
    ),

    test_positive_removed_by_age = setNames(
      as.list(rep(0L, 3)),
      AGE_CLASSES
    ),
    
    vaccinated_total = 0L,
    
    vaccinated_by_age = setNames(
      as.list(rep(0L, 3)),
      AGE_CLASSES
    )
  )
}

#9 simulation class
Simulation <- R6::R6Class(
  "Simulation",
  public = list(
    
    p = NULL,
    seed = NULL,
    herd = NULL,
    t = 0,
    heap = NULL,
    acc = NULL,
    year_records = NULL,
    processed_events = NULL,
    ever_infected = FALSE,
    consecutive_clean_samples = 0L,
    eradicated = FALSE,
    eradication_time = NULL,
    
    initialize = function(params, seed = NULL) {
      validate_params(params)
      
      self$p <- params
      self$seed <- as.integer(
        if (is.null(seed)) params$seed else seed
      )
      
      set.seed(self$seed)
      
      self$herd <- Herd$new()
      self$heap <- EventHeap$new()
      self$acc <- blank_acc()
      self$year_records <- list()
      self$processed_events <- list()
      
      init_population(self)
      schedule_recurring_events(self)
    },
    
    run = function() {
      run_ssa(self)
    }
  )
)

#10 aging
schedule_age_up <- function(sim, a) {
  if (a$age_class == "calf")
    sim$heap$push(
      a$birth_time + sim$p$calf_to_heifer_age,
      "age_to_heifer",
      list(id = a$id)
    )
  
  if (a$age_class == "heifer")
    sim$heap$push(
      a$birth_time + sim$p$heifer_to_cow_age,
      "age_to_cow",
      list(id = a$id)
    )
}

#11 initial population
init_population <- function(sim) {
  p <- sim$p
  h <- sim$herd
  
  for (i in seq_len(p$init_calves)) {
    a <- h$new_animal(
      -runif(1, 0, p$calf_to_heifer_age),
      "calf",
      "S"
    )
    
    schedule_age_up(sim, a)
  }
  
  for (i in seq_len(p$init_heifers)) {
    a <- h$new_animal(
      -runif(
        1,
        p$calf_to_heifer_age,
        p$heifer_to_cow_age
      ),
      "heifer",
      "S"
    )
    
    schedule_age_up(sim, a)
  }
  
  ni <- round(
    p$init_cows * p$init_prevalence
  )
  
  nl <- round(
    ni * p$init_LH_ratio
  )
  
  z <- c(
    rep("S", p$init_cows - ni),
    rep("I_L", nl),
    rep("I_H", ni - nl)
  )
  
  if (length(z) > 1)
    z <- z[sample.int(length(z))]
  
  for (s in z) {
    h$new_animal(
      -runif(
        1,
        p$heifer_to_cow_age,
        p$heifer_to_cow_age + 6
      ),
      "cow",
      s
    )
    
    if (s %in% c("I_L", "I_H"))
      sim$ever_infected <- TRUE
  }
}

#12 recur a scheduled event
schedule_recurring_events <- function(sim) {
  # At programme start, offer one dose to the initially eligible animals.
  # Later cohorts are offered vaccination when they enter an eligible class.
  if (isTRUE(sim$p$vaccination_enabled) && isTRUE(sim$p$vaccination_initial_campaign))
    sim$heap$push(0, "vaccination_campaign")

  sim$heap$push(
    0,
    "snapshot"
  )
  
  if (
    !is.null(sim$p$surveillance_interval_years) &&
    sim$p$surveillance_interval_years > 0
  )
    sim$heap$push(
      sim$p$surveillance_interval_years,
      "surveillance"
    )

  if (isTRUE(sim$p$testing_enabled))
    sim$heap$push(
      sim$p$testing_interval_years,
      "control_test"
    )
  
  sim$heap$push(
    sim$p$horizon_years,
    "end"
  )
}

#13 force of infection
corrective_factor <- function(sim) {
  if (is.null(sim$p$corrective_management))
    return(1)

  efficacy <- sim$p$corrective_management_efficacy[[sim$p$corrective_management]]
  max(0, 1 - efficacy)
}

foi <- function(sim, ac) {
  h <- sim$herd
  n <- h$n_age(ac)
  
  if (!n)
    return(0)
  
  corrective_factor(sim) * (
    sim$p$beta_L[[ac]] * h$n_state(ac, "I_L") +
      sim$p$beta_H[[ac]] * h$n_state(ac, "I_H")
  ) / n
}

#14 stochastic event rates
stochastic_rates <- function(sim) {
  p <- sim$p
  h <- sim$herd
  r <- list()
  
  add <- function(x, k, a) {
    if (x > 0) {
      list(
        rate = x,
        kind = k,
        ac = a
      )
    } else {
      NULL
    }
  }
  
  r <- c(r, list(add(
    h$n_age("calf") * p$calf_mortality_rate,
    "death_natural",
    "calf"
  )))
  
  r <- c(r, list(add(
    h$n_age("heifer") * p$heifer_mortality_rate,
    "death_natural",
    "heifer"
  )))
  
  r <- c(r, list(add(
    h$n_age("cow") * p$cow_culling_rate,
    "cull_background",
    "cow"
  )))
  
  for (ac in AGE_CLASSES) {
    r <- c(r, list(add(
      (
        h$n_state(ac, "I_L") +
          h$n_state(ac, "I_H")
      ) * p$blv_mortality_rate[[ac]],
      "death_blv",
      ac
    )))
  }
  
  r <- c(r, list(add(
    h$n_age("cow") * p$birth_rate,
    "birth",
    "cow"
  )))
  
  if (!is.null(p$maternal_waning_rate)) {
    r <- c(r, list(add(
      h$n_state("calf", "R") * p$maternal_waning_rate,
      "wane",
      "calf"
    )))
  }
  
  for (ac in AGE_CLASSES) {
    r <- c(r, list(add(
      h$n_state(ac, "S") * foi(sim, ac),
      "infect",
      ac
    )))
  }
  
  # Remove NULL entries
  Filter(Negate(is.null), r)
}

#15 random animal selection
random_member <- function(sim, ac, state = NULL) {
  s <- if (is.null(state)) {
    sim$herd$by_age[[ac]]
  } else {
    sim$herd$members[[ac]][[state]]
  }
  
  sim$herd$get_animal(
    s$random_one()
  )
}


random_infected <- function(sim, ac) {
  h <- sim$herd
  
  ids <- c(
    h$members[[ac]]$I_L$as_vector(),
    h$members[[ac]]$I_H$as_vector()
  )
  
  h$get_animal(
    ids[sample.int(length(ids), 1L)]
  )
}

#16 birth and vertical transmission
do_birth_from_dam <- function(sim, dam) {
  
  if (runif(1) >= sim$p$prop_female) {
    sim$acc$births_male <-
      sim$acc$births_male + 1L
    
    return(invisible(NULL))
  }
  
  sim$acc$births_female <-
    sim$acc$births_female + 1L
  
  if (dam$state %in% c("I_L", "I_H")) {
    
    if (runif(1) < sim$p$vertical_transmission_rate) {
      s <- dam$state
      src <- NULL
    } else {
      s <- "R"
      src <- "maternal"
    }
    
  } else {
    s <- "S"
    src <- NULL
  }
  
  calf <- sim$herd$new_animal(
    sim$t,
    "calf",
    s,
    src
  )
  
  schedule_age_up(sim, calf)
  
  if (s %in% c("I_L", "I_H")) {
    sim$acc$new_infected$calf <-
      sim$acc$new_infected$calf + 1L
    
    sim$ever_infected <- TRUE
  }
  
  invisible(NULL)
}

#17 horizontal infection
do_birth <- function(sim) {
  do_birth_from_dam(
    sim,
    random_member(sim, "cow")
  )
}


do_infection <- function(sim, ac) {
  h <- sim$herd
  a <- random_member(sim, ac, "S")
  
  nl <- h$n_state(ac, "I_L")
  nh <- h$n_state(ac, "I_H")
  
  h$set_state(
    a,
    if (runif(1) < nl / (nl + nh))
      "I_L"
    else
      "I_H"
  )
  
  sim$acc$new_infected[[ac]] <-
    sim$acc$new_infected[[ac]] + 1L
  
  sim$ever_infected <- TRUE
}

#18 preventive vaccination -------------------------------------------------
# Vaccine action is all-or-nothing.  Only S animals can become protected;
# vaccination never changes I_L or I_H and vaccine protection does not wane.
vaccinate_animal <- function(sim, a) {
  if (a$vaccinated)
    return(invisible(FALSE))

  a$vaccinated <- TRUE
  sim$acc$vaccinated_total <- sim$acc$vaccinated_total + 1L
  sim$acc$vaccinated_by_age[[a$age_class]] <-
    sim$acc$vaccinated_by_age[[a$age_class]] + 1L

  if (a$state == "S" && runif(1) < sim$p$vaccine_efficacy)
    sim$herd$set_state(a, "R", r_source = "vaccine")

  invisible(TRUE)
}

maybe_vaccinate <- function(sim, a) {
  p <- sim$p
  if (!isTRUE(p$vaccination_enabled) || a$vaccinated || !a$alive ||
      !(a$age_class %in% p$vaccination_groups))
    return(invisible(FALSE))

  # Coverage is the probability that an eligible, not-yet-vaccinated animal
  # receives its single lifetime dose at this opportunity.
  if (runif(1) >= p$vaccination_coverage)
    return(invisible(FALSE))

  vaccinate_animal(sim, a)
}

run_initial_vaccination_campaign <- function(sim) {
  for (ac in sim$p$vaccination_groups) {
    ids <- sim$herd$by_age[[ac]]$as_vector()
    for (id in ids) {
      a <- sim$herd$get_animal(id)
      if (!is.null(a)) maybe_vaccinate(sim, a)
    }
  }
  invisible(NULL)
}

#19 replacement cow
replace_cow <- function(sim, reason) {
  h <- sim$herd
  ids <- h$by_age$heifer$as_vector()
  
  if (length(ids)) {
    bt <- vapply(
      ids,
      function(id)
        h$get_animal(id)$birth_time,
      numeric(1)
    )
    
    a <- h$get_animal(ids[which.min(bt)])
    h$set_age_class(a, "cow")
    maybe_vaccinate(sim, a)
    
    return()
  }
  
  if (
    isTRUE(sim$p$external_replacements_tested) ||
    runif(1) >= sim$p$external_prevalence
  )
    s <- "S"
  else
    s <- if (
      runif(1) < sim$p$external_LH_ratio
    ) "I_L" else "I_H"
  
  a <- h$new_animal(
    sim$t - sim$p$heifer_to_cow_age,
    "cow",
    s
  )

  maybe_vaccinate(sim, a)
  
  sim$acc$external_purchases <-
    sim$acc$external_purchases + 1L
  
  if (s %in% c("I_L", "I_H"))
    sim$ever_infected <- TRUE
}

#20 execute a stochastic event
execute_stochastic <- function(sim, rates, total) {
  draw <- runif(1) * total
  cumulative <- 0
  event <- rates[[length(rates)]]
  
  for (x in rates) {
    cumulative <- cumulative + x$rate
    
    if (draw <= cumulative) {
      event <- x
      break
    }
  }
  
  h <- sim$herd
  
  if (event$kind == "death_natural") {
    
    h$remove_animal(
      random_member(sim, event$ac)
    )
    
    sim$acc$deaths_natural[[event$ac]] <-
      sim$acc$deaths_natural[[event$ac]] + 1L
  }
  
  else if (event$kind == "cull_background") {
    h$remove_animal(
      random_member(sim, "cow")
    )
    
    sim$acc$culled_cow <-
      sim$acc$culled_cow + 1L
    
    replace_cow(
      sim,
      "normal_cull"
    )
  }
  
  else if (event$kind == "death_blv") {
    h$remove_animal(
      random_infected(sim, event$ac)
    )
    
    sim$acc$dead_blv[[event$ac]] <-
      sim$acc$dead_blv[[event$ac]] + 1L
    
    if (event$ac == "cow")
      replace_cow(
        sim,
        "blv_death"
      )
  }
  
  else if (event$kind == "birth")
    do_birth(sim)
  
  else if (event$kind == "wane")
    h$set_state(
      random_member(sim, "calf", "R"),
      "S"
    )
  
  else if (event$kind == "infect")
    do_infection(
      sim,
      event$ac
    )
}

#21 surveillance
surveillance <- function(sim) {
  if (!sim$ever_infected)
    return()
  
  h <- sim$herd
  
  ids <- unlist(
    lapply(
      AGE_CLASSES,
      function(ac)
        c(
          h$members[[ac]]$I_L$as_vector(),
          h$members[[ac]]$I_H$as_vector()
        )
    )
  )
  
  found <- length(ids) > 0 &&
    any(
      runif(length(ids)) <
        sim$p$surveillance_sensitivity
    )
  
  if (found)
    sim$consecutive_clean_samples <- 0L
  else {
    sim$consecutive_clean_samples <-
      sim$consecutive_clean_samples + 1L
    
    if (
      !sim$eradicated &&
      sim$consecutive_clean_samples >=
      sim$p$n_clean_samples_required
    ) {
      sim$eradicated <- TRUE
      sim$eradication_time <- sim$t
    }
  }
}

#22 save one year's information
flush_year <- function(sim) {
  h <- sim$herd
  
  pa <- setNames(
    lapply(AGE_CLASSES, h$n_age),
    AGE_CLASSES
  )
  
  ps <- setNames(
    lapply(
      AGE_CLASSES,
      function(ac)
        setNames(
          lapply(
            STATES,
            function(s)
              h$n_state(ac, s)
          ),
          STATES
        )
    ),
    AGE_CLASSES
  )
  
  sim$year_records[[
    length(sim$year_records) + 1L
  ]] <- c(
    list(
      time = sim$t,
      year = floor(sim$t),
      population_by_age = pa,
      population_by_age_state = ps
    ),
    sim$acc
  )
  
  sim$acc <- blank_acc()
}

#23 scheduled events
handle_scheduled <- function(sim, kind, meta) {
  sim$processed_events[[
    length(sim$processed_events) + 1L
  ]] <- list(
    time = sim$t,
    kind = kind
  )
  
  h <- sim$herd
  
  if (kind == "age_to_heifer") {
    a <- h$get_animal(meta$id)
    
    if (
      !is.null(a) &&
      a$alive &&
      a$age_class == "calf"
    ) {
      h$set_age_class(
        a,
        "heifer"
      )
      
      if (identical(a$r_source, "maternal"))
        h$set_state(
          a,
          "S"
        )
      
      schedule_age_up(
        sim,
        a
      )

      maybe_vaccinate(sim, a)
    }
  }
  
  else if (kind == "age_to_cow") {
    a <- h$get_animal(meta$id)
    
    if (
      !is.null(a) &&
      a$alive &&
      a$age_class == "heifer"
    ) {
      total <- sum(
        vapply(
          AGE_CLASSES,
          h$n_age,
          integer(1)
        )
      )
      
      if (
        isTRUE(sim$p$cap_herd_at_target) &&
        total > sim$p$herd_size
      ) {
        h$remove_animal(a)
        
        sim$acc$heifers_sold_surplus <-
          sim$acc$heifers_sold_surplus + 1L
      } else {
        h$set_age_class(
          a,
          "cow"
        )

        maybe_vaccinate(sim, a)
      }
    }
  }

  else if (kind == "vaccination_campaign") {
    run_initial_vaccination_campaign(sim)
  }

  else if (kind == "control_test") {
    do_control_test(sim)

    sim$heap$push(
      sim$t + sim$p$testing_interval_years,
      "control_test"
    )
  }
  
  else if (kind == "snapshot") {
    flush_year(sim)
    
    sim$heap$push(
      sim$t + 1,
      "snapshot"
    )
  }
  
  else if (kind == "surveillance") {
    surveillance(sim)
    
    sim$heap$push(
      sim$t + sim$p$surveillance_interval_years,
      "surveillance"
    )
  }
}

# Test-and-remove intervention.  Vaccine-protected R animals are not
# considered infected and cannot be test-positive in this simplified model.
target_ids <- function(sim, target_group) {
  if (identical(target_group, "all")) {
    return(unlist(lapply(AGE_CLASSES, function(ac) {
      sim$herd$by_age[[ac]]$as_vector()
    }), use.names = FALSE))
  }

  sim$herd$by_age[[target_group]]$as_vector()
}

do_control_test <- function(sim) {
  p <- sim$p
  h <- sim$herd

  # Snapshot the IDs because removals mutate the herd indexes.
  for (id in target_ids(sim, p$testing_target_group)) {
    a <- h$get_animal(id)
    if (is.null(a) || !a$alive)
      next

    sim$acc$tested_by_age[[a$age_class]] <-
      sim$acc$tested_by_age[[a$age_class]] + 1L

    if (!(a$state %in% c("I_L", "I_H")) || runif(1) >= p$testing_sensitivity)
      next

    if (identical(p$testing_remove_policy, "high_only") && a$state != "I_H")
      next

    age_at_removal <- a$age_class
    h$remove_animal(a)
    sim$acc$test_positive_removed_by_age[[age_at_removal]] <-
      sim$acc$test_positive_removed_by_age[[age_at_removal]] + 1L

    if (identical(age_at_removal, "cow"))
      replace_cow(sim, "test_positive")
  }

  invisible(NULL)
}

#24 main SSA loop
run_ssa <- function(sim) {
  repeat {
    rates <- stochastic_rates(sim)
    
    total <- sum(
      vapply(
        rates,
        function(x) x$rate,
        numeric(1)
      )
    )
    
    ts <- if (total > 0)
      sim$t + rexp(1, total)
    else
      Inf
    
    if (sim$heap$peek_time() <= ts) {
      e <- sim$heap$pop()
      
      sim$t <- e$time
      
      handle_scheduled(
        sim,
        e$kind,
        e$meta
      )
      
      if (e$kind == "end")
        break
    } else {
      sim$t <- ts
      
      if (sim$t >= sim$p$horizon_years) {
        sim$t <- sim$p$horizon_years
        break
      }
      
      execute_stochastic(
        sim,
        rates,
        total
      )
    }
  }
  
  flush_year(sim)
  results(sim)
}

#25 return results
results <- function(sim) {
  list(
    seed = sim$seed,
    params_seed = sim$p$seed,
    year_records = sim$year_records,
    eradicated = sim$eradicated,
    eradication_time = sim$eradication_time,
    ever_infected = sim$ever_infected,
    final_population = setNames(
      lapply(
        AGE_CLASSES,
        sim$herd$n_age
      ),
      AGE_CLASSES
    )
  )
}

#26 run one replicate
run_one <- function(params, seed = NULL) {
  Simulation$new(
    params,
    seed
  )$run()
}

#27 seeds
run_strategy <- function(params, n_runs, master_seed = 12345L) {
  master_seed <- as.integer(master_seed)
  set.seed(master_seed)
  replicate_seeds <- sample.int(.Machine$integer.max, as.integer(n_runs))

  list(
    n_runs = n_runs,
    master_seed = master_seed,
    replicate_seeds = replicate_seeds,
    runs = lapply(replicate_seeds, function(seed) run_one(params, seed))
  )
}

#28 consistency check
herd_indexes_consistent <- function(h) {
  for (ac in AGE_CLASSES) {
    
    if (
      h$n_age(ac) !=
      sum(
        vapply(
          STATES,
          function(s)
            h$n_state(ac, s),
          integer(1)
        )
      )
    )
      return(FALSE)
    
    for (s in STATES)
      for (
        id in h$members[[ac]][[s]]$as_vector()
      ) {
        a <- h$get_animal(id)
        
        if (
          is.null(a) ||
          !a$alive ||
          a$age_class != ac ||
          a$state != s
        )
          return(FALSE)
      }
  }
  
  TRUE
}

