# BLV baseline cleanup
#- R/simulation.R: parameters, animals, herd, event queue, hybrid SSA engine, stochastic rates/events, scheduled events, and results
#- analysis.R: repeated-run result summary only.
#- test.baseline.R: base-R baseline tests.



#The core requires only `R6`; it does not load ggplot2. Intervention switches
#are retained as explicit inactive placeholders and fail fast if enabled.
#Surveillance is kept as an outcome protocol: its eradication flag means
#"not detected under the surveillance protocol," not true absence of BLV.
