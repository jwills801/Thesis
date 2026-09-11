# Diagnostics

Standalone MATLAB scripts (run each directly, e.g. `run('diagnostics/checkEnergyBalance.m')`
from the `ManualSimulation` folder) that verify the matrices/energy bookkeeping behind the
controllers, following up on the code review before the Humboldt sea-state comparison study.

Each script prints `PASS`/`FAIL` per assertion with the numeric mismatch.

- `checkPhysicalStateSpace.m` — `getPhysical.m`'s `A`,`B` vs. the main(1).tex "Full State Space" derivation.
- `checkTransitionMatrices.m` — `control/getTransition.m`'s `M`,`H` vs. brute-force ZOH propagation.
- `checkUtilityMatrices.m` — `control/getUtilityMatrices.m`'s `L`,`C` vs. direct construction.
- `checkMPC_EHA.m` — `getControl.m`'s `ctrl.MPC.{A,Bx,Bt}` vs. a finite-difference gradient of an
  independently-built `E_elec(u)`.
- `checkTerminalCost.m` — `ctrl.termCost(d).{P,rT,b}` vs. brute-force optimal-control energy.
- `checkAstarVsBruteForce.m` — confirms `MPC_Astar.m`'s branch-and-bound actually returns the true
  minimum-cost sequence over a short, exhaustively-enumerable horizon.
- `checkEnergyBalance.m` — `aveMechPow == aveElecPow + aveLoss` for one short run per drivetrain
  family, using `evaluate.m`'s own outputs.

## Fixed as part of this pass

- **`evaluation/evaluate.m`** (DHD/PassivePump branch): the mech-to-elec conversion loss was being
  computed from the *already-derated* `aveElecPow` instead of the pre-derated value, breaking the
  energy balance by ~2.25% of `aveElecPow` on every DHD/PassivePump run. `checkEnergyBalance.m`
  catches this (was FAIL, now PASS).
- **`evaluation/getMotorLoss.m`**: the unconditional `figure(...)` call is now gated behind
  `params.simu.makePlots` (default `true`, set `false` for batch/sweep runs) — needed before
  batch-running many sea states in phase 2.
- **`control/getTransition.m`, `control/getUtilityMatrices.m`**: extracted out of `getControl.m`'s
  local functions into their own files so they're directly testable, and the dead forward-Euler
  lines in `getTransition` (computed, then immediately overwritten by the ZOH version) were removed.
  Behavior is unchanged — `getControl.m` already used the ZOH result.
- **`Run_All_Cases.m` / `main_WEC_Simulation.m`**: added a `ConsiderLosses` column so EHA's
  mechanical- vs electrical-energy-optimized `MPC_QP` variants are both expressible (previously
  `considerLosses` was hardcoded to 0 for every EHA run).
- **`parameters/getHydraulic.m`**: `hyd.switchMap.valveConstant` (needed by `getValveLoss.m`'s
  `'PassivePump'` branch for its steady open-valve loss) was only being set for `case 'DHD'` --
  every `PassivePump`/`CoulombDamping` run (including the one actively selected in
  `Run_All_Cases.m`, row 4) errored inside `evaluate.m`. Fixing the crash by just reusing DHD's
  cached `SwitchMap.mat` (git history shows this flip-flopped between `case 'DHD'` and
  `case {'DHD','PassivePump'}` across several commits -- not a clean regression, an unresolved
  question) turned out to be physically wrong too: it borrows the sizing of a **DHD switching
  valve** ("2WRC-4x size 80") for PassivePump's **fixed check valve**, and measured `elecRGP`
  (0.17) came in far below the main(1).tex baseline (0.47) for the same nominal case. Per your
  direction, PassivePump now gets its own `valveConstant`, sized via the orifice equation
  `Q=valveConstant*sqrt(dP)` so it drops ~0.5MPa at the worst-case (cap-side) max flow, using the
  same `vMax=1 m/s` cylinder-velocity design point already used for EHA loss-map sizing.
- **`evaluation/evaluate.m`**: even with the check-valve fix above, `elecRGP` barely moved --
  because the dominant PassivePump loss wasn't coming from `getValveLoss` at all. `evaluate.m`
  called `getMotorLoss` (the "series hydraulic-to-electric converter" loss, main(1).tex's DHD-only
  architecture) for PassivePump too. `getMotorLoss` computes `u_elec = dyn.u - u_hyd` and runs it
  through the EHA loss function as if it were a real generator's torque; for PassivePump, that gap
  is just `coulombDamping.m`'s `tanh()`-smoothing of its bang-bang switch, not a real actuator
  (measured `rms(u_elec)` = 56% of `rms(u_hyd)` -- far too large to be numerical noise). Measured
  effect: `getValveLoss` alone gave `aveLoss=1.4kW`; adding `getMotorLoss` inflated it to `142kW`
  (100x). Fixed: PassivePump now only calls `getValveLoss`; the flat 85% main-motor conversion
  efficiency (the fixed accounting-order bug above) still applies to it via the second switch
  block. DHD is unaffected -- its discrete controllers already produce `u_elec~=0` with no
  smoothing, so `getMotorLoss` contributed negligibly there anyway.
  **Result after both PassivePump fixes**, known-good baseline (Hs=2.5m, Tp=8s, 20.6MPa):
  mechRGP=0.599, elecRGP=0.505 -- matches main(1).tex's reported 0.56/0.47 much more closely than
  the pre-fix 0.599/0.167. All sea states now give positive, sensible `elecRGP` (previously went
  as negative as -2.5 for a mild sea state at the default high pressure); the mild-sea/low-pressure
  case beats the mild-sea/high-pressure case (0.482 vs 0.314 elecRGP), which is exactly the effect
  phase 2's per-sea-state pressure optimization is meant to capture.
- **`evaluation/evaluate.m` / `plotting/plotAll.m`**: `mechRGP`/`elecRGP` were only ever computed
  inside `plotAll.m`, which also pops a full figure set -- meaning no automated sweep (pressure
  grid search, sea-state loop) could get an RGP number without spawning figures every iteration.
  Moved the two RGP lines into `evaluate.m` (now takes `ctrl` as a third argument, for
  `ctrl.optTraj.avePow`); `plotAll.m` just displays the already-computed values. Needed for
  phase 2's `optimizePressure.m` grid search.

## Additional fix, confirmed necessary by the phase-2 validation run

- **`getControl.m`'s `ctrl.numHorizons = 2.5*peakPeriod/timeHorizon` was never rounded.** For
  `MPC_QP` (`timeHorizon = 5*dt`) and `MPC_Astar`/`MPC_Astar_cont` (`timeHorizon = 0.2`), this is
  only an integer for specific `peakPeriod`/`dt` combinations (it happened to be exactly 100 and 400
  for the previous default `peakPeriod=8`, `dt=0.01`). First caught synthetically (a diagnostic
  using `peakPeriod=0.2` made `numHorizons=2.5`, crashing `getUtilityMatrices`'s `NaN(n*m,m)`), then
  confirmed as a real, immediate blocker: the very first placeholder Humboldt sea-state bin
  (`Tp=7s`) gives `numHorizons=87.5` and crashed `validatePhase2Subset.m`'s DHD/MPC_Astar case.
  Fixed by rounding `numHorizons` in both the `MPC_QP` and `MPC_Astar`/`MPC_Astar_cont` cases.

## EHA loss investigation (negative elecRGP in mild sea states)

Unlike PassivePump, this checked out as a **real physical effect**, not a bug:
`EHA.LossFunc(0,0)` (2000RPM pump idling, zero flow/pressure) = 21.3kW -- a genuine parasitic
floor from a fixed-speed pump (friction/windage terms scale with the fixed shaft speed, not
flow). For a mild sea state (`aveMechPow~50kW`), that floor alone is ~43% of the mechanical power,
and total losses exceed it. This matches main(1).tex's own documented tradeoff (constant-speed
chosen over variable-speed specifically to avoid shaft kinetic-energy costs, with an acknowledged
efficiency downside). Self-consistency check passed too: `considerLosses=1` produces lower loss
than `considerLosses=0` (64.2kW vs 75.7kW), i.e. the electrical-optimized controller is correctly
trading off mechanical capture to cut losses.

Added `parameters/plotEHAEfficiencyMap.m` (called from `getHydraulic.m` behind
`runParams.plotEHAEfficiency`, off by default) to visualize this: the real loss map has a kinked,
diamond-shaped minimum near zero torque that the quadratic fit used internally by the
electrical-optimized `MPC_QP` controller (`ctrl.MPC` in `getControl.m`) smooths into a round
elliptical bowl -- a fit-quality gap worth knowing about if `considerLosses=1` results look off in
a specific operating region, though not something fixed here.

## Known issues found but intentionally left untouched (per discussion)

- **Sliding Mode is broken.** `controlLaw.m` dispatches on `case 'SlidingMode'`, but
  `Run_All_Cases.m`/`plotAll.m` use `'Sliding Mode'` (with a space) — a real string mismatch, not a
  style nit: every 'Sliding Mode' case in `Run_All_Cases.m` falls through `controlLaw.m`'s switch
  with no matching case, so `out` is never assigned and it errors. Separately, `getControl.m` has
  no case at all for sliding mode, so `ctrl.lambda`, `ctrl.phi`, `ctrl.horizonInd`, and
  `ctrl.limitChoices` — all read by `slidingMode.m` — are never initialized. Also worth knowing:
  `slidingMode.m` has its own private copy of `getTransition` that still uses **forward Euler**,
  not the zero-order-hold version the rest of the codebase (`getControl.m`'s MPC path) now uses via
  `control/getTransition.m` — so even once the naming/`ctrl`-field bugs are fixed, Sliding Mode's
  internal prediction model would be discretized differently from MPC_QP/MPC_Astar's. Out of scope
  for this comparison study; not part of the Coulomb/EHA/DHD case matrix.
- **`control/MPC_Astar_cont2.m` looks unfinished.** It references `params.uInd` and
  `params.damping`, neither of which is ever set anywhere in `getParameters`/`getHydraulic`/
  `getPhysical`/`Run_All_Cases`. Not called from `Run_All_Cases.m`. Left as-is; flag if you want it
  either finished or deleted.
