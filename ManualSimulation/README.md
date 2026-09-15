# ManualSimulation

MATLAB simulation of an oscillating wave surge converter (OWSC) flap
coupled to one of several power take-off (PTO) drivetrains, used to compare
drivetrain/controller combinations against a theoretical optimal-power
benchmark across Humboldt Bay sea states.

## Drivetrains and controllers

| Drivetrain | Controller(s) | Notes |
|---|---|---|
| `PassivePump` | `CoulombDamping` | Fixed check valve, no real controller — a smoothed bang-bang torque. |
| `EHA` | `MPC_QP`, `PI` | Electro-hydrostatic actuator. Two speed/displacement strategies: fixed-speed/variable-displacement (default) or fixed-displacement/variable-speed (`runParams.ehaFixedDisplacement=true`, adds reflected inertia to the flap's own dynamics — see `parameters/getPhysical.m`). `considerLosses` selects mechanical- vs electrical-energy-optimized control. |
| `DHD` (Digital Hydraulic Drive) | `MPC_Astar`, `MPC_Astar_cont`, `MPC_QP`, `PI` | Discrete PTO force options from `pressure_rails` (2/3/4/5) independently-switched cap/rod sides. `MPC_Astar` is the validated combinatorial search; `MPC_QP` and `PI` rounds a continuous solution to the nearest discrete option. |
| — | `Sliding Mode` | **Known broken** — see `diagnostics/ReadMe.md`. Not part of the current comparison. |

## Directory structure

```
main_WEC_Simulation.m   Top-level script: one drivetrain/controller/sea-state case
Run_All_Cases.m         Runs main_WEC_Simulation.m over a table of cases

parameters/   The true parameter-builders: getParameters/getPhysical/getHydraulic/getSimulation
models/       EHA & DHD physics/loss models (makeEHALossMap*, makeSwitchLossMap)
  exploratory/  Comparison-only EHA loss models, not wired into any live run
optimization/ Pressure/rail/area optimization utilities (optimizePressure, sizeCylinderArea,
              evenlySpacedRails, buildDenseSwitchMap)
wave/         Wave spectrum, excitation torque time series, Humboldt sea-state bins
control/      Per-timestep controllers + the one-time ctrl-struct setup (getControl.m)
dynamics/     The time-domain integration loop and one forward-Euler plant step
evaluation/   Post-run loss/power accounting and RGP (relative generated power)
plotting/     Diagnostic plots for a single run (plotAll.m, plotEHAEfficiencyMap.m)
diagnostics/  Standalone validation scripts (run manually) + ReadMe.md of findings/fixes
HHEA_DP/      Older leftover results from a prior approach; not part of the active pipeline
archive/      Old one-off study scripts superseded by validated fixes (see diagnostics/ReadMe.md's
              "Fixed during the drivetrain comparison study" section) -- their pre-fix result sets
              and console logs live outside this repo instead, see archive/ARCHIVE_LOCATION.md,
              so they don't get swept up when syncing this directory to a compute cluster.
```

## Drivetrain comparison study -- current top-level scripts

The active (post-bugfix) pipeline for comparing drivetrain families across all 8 Humboldt sea
states, in the order they're meant to run:

1. **`rebuildSwitchMaps.m`** -- builds `results/denseSwitchMaps.mat` (200ms) and
   `results/denseSwitchMaps_100ms.mat` (100ms), ONE map each, valid for DHD2/DHD3/DHD4 alike (see
   `optimization/buildDenseSwitchMap.m`'s header). Rerun only if `models/makeSwitchLossMap.m`'s
   valve physics or `results/sizedAreas.mat`'s cylinder areas change.
2. **`Run_AstarConvergenceStudy_100ms.m`** -- A* search-horizon (`m_Astar`) convergence sweep for
   DHD2/DHD3/DHD4 at a single sea state (currently sea state 7, the most annual-energy bin by
   probability x theoretical power -- see its header), used to pick each family's `m_Astar` for
   the full sweep below without re-tuning it per sea state.
3. **`Run_FinalGridSearch_100ms.m`** -- the full pressure grid search, all 8 sea states x
   {PassivePump, DHD2, DHD3, DHD4, EHA_mech, EHA_elec}, using each DHD family's `m_Astar` chosen
   in step 2. Writes `results/gridSearchFinal/`.

Each writes its own results immediately per-task (survives an interrupted session) and has a
standalone, safely-rerunnable aggregator (`aggregateAstarConvergenceResults.m`,
`aggregateFinalGridSearchResults.m`) plus plotting script
(`plotAstarConvergenceResults.m`, `plotFinalPressureSweeps.m`).

## Data flow for one simulation run

```
runParams (drive, controller, pressure_rails, considerLosses, areas, ...)
        |
        v
getParameters.m  -->  getPhysical.m   (flap state-space, phys.sys)
        |         -->  getHydraulic.m  (cylinder, pressure rails, phys.hyd)
        |               |--> models/makeEHALossMap*.m, models/makeSwitchLossMap.m
        |               |--> optimization/evenlySpacedRails.m (3/4-rail DHD/PassivePump)
        |         -->  getSimulation.m (timing/sea-state defaults, phys.simu)
        v
generateExcitingTorque.m   (wave spectrum -> excitation torque time series)
        v
getControl.m                (theoretical-optimal trajectory + controller setup matrices)
        v
timeLoop.m  <-->  controlLaw.m  <-->  [PIcontrol / slidingMode / MPC_QP / MPC_Astar* / coulombDamping]
        |                                   (per-timestep control decision)
        v
advanceStep.m                (forward-Euler true-plant integration, one call per timestep)
        v
evaluate.m                   (getEHALoss / getMotorLoss / getValveLoss -> aveLoss, mechRGP, elecRGP)
        v
plotAll.m                     (diagnostic plots, opt-in via params.simu.makePlots)
```

`main_WEC_Simulation.m` runs exactly this sequence as a script. `Run_All_Cases.m`
sets `runParams` from a table and calls `main_WEC_Simulation.m` once per row.
`optimization/optimizePressure.m` and `optimization/sizeCylinderArea.m` instead call
the middle steps (`getParameters`→`generateExcitingTorque`→`getControl`→`timeLoop`)
directly in a loop, to grid-search a pressure or iteratively resize a cylinder.

## Entry points

- **Single case**: set `runParams` in the workspace, then run `main_WEC_Simulation.m`.
- **Table of cases**: `Run_All_Cases.m` (currently hardcoded to run only one row —
  see its header comment).
- **Pressure optimization for a sea state** (PassivePump/DHD): `optimization/optimizePressure.m`.
- **Cylinder-area sizing for a drivetrain family**: `optimization/sizeCylinderArea.m`.
- **Humboldt sea-state bins**: `wave/humboldtSeaStates.m` (8 occurrence-weighted
  representative sea states, from `wave/humboldtManualBinsCells.csv`; pass
  `true` to also plot the JPD heatmap and write `humboldtBinSummary.csv`).
- **Validate the math**: anything in `diagnostics/` — run directly, each prints
  `PASS`/`FAIL` per assertion. `diagnostics/ReadMe.md` has the full history of
  what's been checked and fixed.

## Known issues / where to look next

See `diagnostics/ReadMe.md` for the validated-fix history and the list of
issues found but intentionally left untouched (Sliding Mode's naming/setup
bugs). `control/MPC_Astar_cont2.m` (an unfinished DHD A* variant referencing
never-set `params.uInd`/`params.damping`, never wired into `controlLaw.m`)
has since been deleted as dead code. In brief:

- `models/exploratory/makeEHALossMap_v1p5.m` fits its quadratic loss
  surrogate against cylinder velocity, not `thetaDot` — a unit mismatch if
  used as a controller's cost directly. `models/makeEHALossMap_fixedDisp.m`
  is the corrected, controller-ready version of the same
  fixed-displacement/variable-speed physics.
- `models/exploratory/makeEHALossMap_v2.m` (both speed and displacement
  free) uses a different pump-oversizing margin (1.5x) than the other two
  files (1.2x) — not reconciled, matters if comparing its numbers directly
  to the others.
- `parameters/getHydraulic.m`'s DHD case loads a stale precomputed
  `SwitchMap.mat` by default rather than regenerating for the current
  pressure rails; callers that sweep pressure (`optimization/optimizePressure.m`,
  `optimization/sizeCylinderArea.m`) override `params.hyd.switchMap` themselves
  after the fact. `optimization/buildDenseSwitchMap.m` builds one reusable
  dense-grid map per cylinder-area family instead of regenerating per
  candidate pressure.
