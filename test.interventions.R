# Focused checks for the three intervention mechanisms.
source("simulation.R")

# Corrective management changes only the force of infection.
p0 <- make_params(init_calves = 0, init_heifers = 1, init_cows = 0,
                  birth_rate = 0, seed = 1)
s0 <- Simulation$new(p0, seed = 1)
s0$herd$new_animal(0, "heifer", "I_L")
base_foi <- foi(s0, "heifer")

p1 <- p0
p1$corrective_management <- "high"
s1 <- Simulation$new(p1, seed = 1)
s1$herd$new_animal(0, "heifer", "I_L")
stopifnot(isTRUE(all.equal(foi(s1, "heifer"), base_foi * 0.40)))

# All eligible S animals receive one dose in the initial campaign.  With an
# efficacy of one, all of them become vaccine-protected R animals.
pv <- make_params(
  init_calves = 0, init_heifers = 2, init_cows = 2,
  init_prevalence = 0, birth_rate = 0,
  vaccination_enabled = TRUE,
  vaccination_groups = c("heifer", "cow"),
  vaccination_coverage = 1, vaccine_efficacy = 1, seed = 2
)
sv <- Simulation$new(pv, seed = 2)
handle_scheduled(sv, "vaccination_campaign", list())
stopifnot(
  sv$herd$n_state("heifer", "R") == 2L,
  sv$herd$n_state("cow", "R") == 2L,
  sv$acc$vaccinated_total == 4L
)

# A heifer vaccinated before ageing does not receive a second dose as a cow.
a <- sv$herd$get_animal(sv$herd$by_age$heifer$as_vector()[1])
handle_scheduled(sv, "age_to_cow", list(id = a$id))
stopifnot(a$vaccinated, sv$acc$vaccinated_total == 4L)

# A perfectly sensitive all-positive test removes an infected cow and fills
# the cow place through the existing replacement mechanism.
pt <- make_params(
  init_calves = 0, init_heifers = 0, init_cows = 1,
  init_prevalence = 1, init_LH_ratio = 1, birth_rate = 0,
  testing_enabled = TRUE, testing_sensitivity = 1,
  testing_target_group = "cow", testing_remove_policy = "all_positive",
  seed = 3
)
st <- Simulation$new(pt, seed = 3)
do_control_test(st)
stopifnot(
  st$herd$n_state("cow", "I_L") == 0L,
  st$herd$n_age("cow") == 1L,
  st$acc$test_positive_removed_by_age$cow == 1L,
  herd_indexes_consistent(st$herd)
)

# All interventions can run together in the SSA engine.
out <- run_one(make_params(
  herd_size = 50, horizon_years = 2,
  corrective_management = "medium",
  testing_enabled = TRUE, testing_interval_years = .5,
  vaccination_enabled = TRUE, vaccination_coverage = .7,
  seed = 4
), seed = 4)
stopifnot(length(out$year_records) >= 3L)

cat("All intervention tests passed.\n")
