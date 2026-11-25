# Current intervention-model parameters and assumptions

This checklist describes the configuration used by `simulation.R` and
`run_intervention_factorial_420.R` on 25 November 2026.  It is intended for
scientific checking before the 420-scenario batch is started.

## 1. Batch design

| Item | Value used | Notes |
|---|---:|---|
| Herd sizes | 100, 200, 800 | Three separate scenario strata. |
| Intervention scenarios | 420 | `3 herd sizes x 5 management levels x 7 testing settings x 4 vaccine coverages`. |
| Replicates per scenario | 100 | 42,000 stochastic runs in total. |
| Time horizon | 40 years | Annual records are saved for every replicate. |
| Master seed | 12345 | Generates the same 100 replicate seeds for every scenario, enabling matched comparisons. |
| Starting age structure | Original parameter structure | 12.5% calves, 37.5% heifers, 50% cows; it is *not* the burn-in demographic structure. |
| Starting counts, N = 100 | 12 calves, 38 heifers, 50 cows | R rounding is used. |
| Starting counts, N = 200 | 25 calves, 75 heifers, 100 cows |  |
| Starting counts, N = 800 | 100 calves, 300 heifers, 400 cows |  |

## 2. Infection state at intervention start

| Item | Value used | Interpretation |
|---|---:|---|
| Burn-in prevalence source | Median of the six scenario-level mean prevalences among Year-100 active (still infected) runs | Extinct burn-in runs are excluded when estimating the endemic level. |
| Initial prevalence, N = 100 | 0.8882 | Median across the six burn-in introduction scenarios for this herd size. |
| Initial prevalence, N = 200 | 0.8911 | Same calculation. |
| Initial prevalence, N = 800 | 0.8903 | Same calculation. |
| Prevalence scope | All age groups | The prevalence is independently assigned to calves, heifers, and cows at time zero. It is therefore intended as *whole-herd* prevalence, not cow-only prevalence. |
| Initial low/high split | 70% `I_L`, 30% `I_H` | Applied to initially infected animals after rounding. |
| Initial maternal protection | None explicitly assigned | Existing time-zero calves are initially susceptible or infected. Maternal protection applies to calves born during the simulation. |

## 3. Population and demographic processes

| Parameter | Value | How the code uses it |
|---|---:|---|
| Female-calf entry rate | 0.38 per cow-year | Each cow has a continuous-time birth event rate of 0.38. Each event adds one live female calf. |
| Natural calf mortality | 0.10 per calf-year | Continuous-time event rate. |
| Natural heifer mortality | 0.08 per heifer-year | Continuous-time event rate. |
| Background cow culling | 0.34 per cow-year | Continuous-time removal rate; every vacancy is replaced. |
| Calf-to-heifer age | 0.5 years | Scheduled transition at six months. |
| Heifer-to-cow age | 2 years | Scheduled transition at 24 months, conditional on an adult-cow vacancy. |
| Replacement-heifer window | 22–24 months | A cow vacancy first selects the oldest eligible heifer in this window. |
| Cow target | Initial cow count | 50, 100, and 400 for the three herd sizes. |
| No eligible heifer | Purchase one external cow | This keeps the adult-cow group at its target. |
| External replacement testing | `TRUE` | Purchased cows are susceptible; `external_prevalence` and `external_LH_ratio` have no effect while this remains true. |
| At 24 months with no cow vacancy | Sell the heifer as surplus | She does not automatically enter the adult-cow group. |
| BLV mortality | calf 0; heifer 0; cow 0.05 per infected cow-year | Extra mortality is applied only to infected cows; a resulting cow vacancy is replaced. |

## 4. Transmission and infection

| Parameter | Value used | How the code uses it |
|---|---:|---|
| Mixing | Herd-wide, frequency dependent | Age is retained for demography/testing, but not for a separate age-specific force of infection. |
| Force of infection | `management multiplier x (beta_L x I_L + beta_H x I_H) / N` | `I_L` and `I_H` are summed across calves, heifers, and cows; `N` is the current total herd population. |
| `beta_Snv_LInv` | `0.008722059 x 365 = 3.183552` per year | Low-infectious source coefficient; source values supplied per day and converted in code. |
| `beta_Snv_HInv` | `0.09607751 x 365 = 35.068791` per year | High-infectious source coefficient; source values supplied per day and converted in code. |
| Corrective-management multiplier | `1 - efficacy` | Multiplies the horizontal force of infection only. |
| New infection classification | 70% `I_L`, 30% `I_H` | Same `init_LH_ratio` is used for horizontal and perinatal infections. |
| Perinatal transmission | 0.186 | Applied only to calves born to an infected (`I_L` or `I_H`) dam. |
| Maternal immunity | Deterministic until 6 months | If the calf of an infected dam is not perinatally infected, it begins as maternal-protected `R` and becomes `S` at the calf-to-heifer transition. |

## 5. Surveillance and eradication definition

| Parameter | Value used | Interpretation |
|---|---:|---|
| Surveillance interval | 0.25 years | Every three months. |
| Surveillance sensitivity, calf | 0.91 | Applied to each infected calf at a surveillance round. |
| Surveillance sensitivity, heifer/cow | 0.95 | Applied to each infected heifer or cow at a surveillance round. |
| Clean rounds required | 3 | Three consecutive surveillance rounds with no detected infection. |
| Model eradication | Undetected under surveillance | It is not proof that no infected animal remains. |

## 6. Interventions varied in the factorial design

| Intervention | Levels |
|---|---|
| Corrective management | `CM_0`, `CM_25`, `CM_50`, `CM_90`, `CM_99`: 0%, 25%, 50%, 90%, 99% horizontal-FOI reduction. |
| Test-and-remove | None; high-only every 3, 6, or 12 months; all-positive every 3, 6, or 12 months. |
| Test target | All calves, heifers, and cows. |
| Test sensitivity | Calf 0.91; heifer/cow 0.95. |
| Removal policy | `high_only` removes detected `I_H`; `all_positive` removes detected `I_L` and `I_H`. |
| Replacement after removal | Removed cow is replaced under the replacement rule; removed calves/heifers are not directly replaced. |
| Vaccination coverage | 0%, 50%, 70%, 100%. The 70% level is provisional and should be changed if a final C50-derived level is chosen. |
| Vaccination target | Susceptible heifers and cows. |
| Initial campaign | Yes, at time zero. |
| Vaccine efficacy | 0.70 probability of take. |
| Vaccine action | All-or-nothing; successful doses move the animal to permanently protected `R`; failed doses leave it susceptible. |
| Repeat vaccination | No. Each animal may receive only one lifetime dose. |

## 7. Output and analysis quantities

| Output | Definition |
|---|---|
| Final prevalence | `(I_L + I_H) / current total herd` at year 40. |
| Eradication probability | Eradicated replicates / 100. |
| Eradication time | Time at the third consecutive clean surveillance round; non-eradicated runs remain missing for this outcome. |
| Full trajectories | Every annual record is retained inside each scenario `.rds` file. |
| Population and age structure | Calf, heifer, cow counts and infection states are saved annually. |

## 8. Items deliberately flagged for checking

1. **`cap_herd_at_target` is currently a legacy parameter and is not used by the code.** Adult-cow regulation is instead implemented through `target_cows`, the 22–24-month replacement window, and surplus-heifer sales.
2. **`beta_Sv_LInv`, `beta_Sv_HInv`, `beta_Snv_LIv`, and `beta_Sv_LIv` are stored but are not used by the current force-of-infection function.** This is consistent with the present all-or-nothing vaccine implementation: successful vaccinees are fully protected `R`; failed vaccinees remain ordinary susceptible animals.
3. **Original age composition is not the demographic equilibrium.** With these replacement rules, it will adjust during early simulation years. This is intentional here because the user selected original demographic parameters but burn-in-derived endemic prevalence.
4. **The initial infection assignment is a modelling assumption.** The burn-in gives a total endemic prevalence, but the model does not yet use age-specific endemic prevalence or maternal protection for calves already present at time zero.
5. **The 70% vaccination-coverage level is provisional.** Confirm or replace it before the long batch starts.
