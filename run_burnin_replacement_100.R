# Full burn-in: 18 scenarios, 100 runs each, using the revised cow-replacement
# policy. Complete annual age-by-disease-state and demographic records are saved.
# Existing compatible pilot runs are reused, so the first completed runs are not
# recalculated.
options(
  blv_burnin_n_runs=100L,
  blv_burnin_output_dir=file.path("results", "introduction_burnin_replacement_100"),
  blv_burnin_cache_dir=file.path("results", "introduction_burnin_replacement_pilot_20")
)
source("burnin_initialization.R")
