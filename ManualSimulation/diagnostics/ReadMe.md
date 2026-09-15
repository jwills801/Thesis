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
  phase 2's `optimization/optimizePressure.m` grid search.

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

Added `plotting/plotEHAEfficiencyMap.m` (called from `getHydraulic.m` behind
`runParams.plotEHAEfficiency`, off by default) to visualize this: the real loss map has a kinked,
diamond-shaped minimum near zero torque that the quadratic fit used internally by the
electrical-optimized `MPC_QP` controller (`ctrl.MPC` in `getControl.m`) smooths into a round
elliptical bowl -- a fit-quality gap worth knowing about if `considerLosses=1` results look off in
a specific operating region, though not something fixed here.

## Fixed during the drivetrain comparison study (A* cap / 100ms / cylinder sizing pass)

- **`evaluation/getMotorLoss.m` was being charged unconditionally for every DHD run**, but it only
  models loss in a *continuous electric trim on top of the discrete rail choice* (`u_elec = dyn.u -
  u_hyd`) — real only for `MPC_Astar_cont`. For plain `MPC_Astar`, `u_elec` is exactly 0 (confirmed
  directly), but `params.hyd.EHA.LossFunc(Q, deltaP=0)` is *not* ~0 for nonzero flow `Q` — so every
  plain-`MPC_Astar` DHD run was being charged tens of kW for a conversion that wasn't happening.
  `evaluate.m` now only calls `getMotorLoss` when `params.runParams.controller` is
  `'MPC_Astar_cont'`. Measured impact: one DHD2 case swung from elecRGP=0.48 to ~0.69 once removed
  — every DHD number computed before this fix (including the published drivetrain-comparison
  artifact) is stale.
- **`models/makeSwitchLossMap.m`'s valve natural frequency didn't match its own datasheet comment.**
  The valve constant comment already said "26ms response time," but `valveNaturalFrequency` was set
  to give ~20.9ms (25 Hz). Corrected to 20.11 Hz (`(pi-acos(zeta))/(wn*sqrt(1-zeta^2))` solved for
  `wn` at `zeta=0.7`, `tr=26ms`) — verified numerically to reproduce exactly 26ms.
- **`control/MPC_Astar.m`'s `getSwitchingLoss` could return NaN mid-search** once `astarIterMax` was
  raised enough to let the search's own multi-step lookahead predict `velA`/`vol` outside
  `switchMap`'s built grid range (a real risk only once the cap stopped truncating the search early
  — see below). Fixed by passing an explicit large finite `interpn` extrapolation value (`1e12`)
  instead of the default NaN, so an out-of-range branch is heavily penalized but still sorts/sums
  correctly (NaN doesn't sort predictably; `getValveLoss.m`'s reporting-only equivalent call is left
  as-is, relying on its existing `omitnan` summation instead, since that's real trajectory data, not
  a search heuristic to discourage).
- **`astarIterMax` was masking both of the above.** The original sweep used `astarIterMax=20`
  (needed because DHD4 uncapped could take hours per sim) — tight enough that the search rarely if
  ever reached the regime where either bug mattered. Raising it to 100000 (confirmed to never
  actually bind, ~0.1-0.4 sec/simulated-sec depending on family/mAstar) surfaced both bugs; fixing
  them is what made the larger cap actually usable.
- **`control/getControl.m` now asserts `params.hyd.switchMap.finalTime == ctrl.timeHorizon`** for
  DHD runs — the coarse control step (`runParams.controlDT`, default 0.2s) and the switch-loss map's
  own simulated window (`hyd.switchTime` passed to `makeSwitchLossMap.m`, also default 0.2s) must
  match, or `getValveLoss.m` silently sizes the post-switch loss window wrong. Was previously
  enforced only by convention/comments; a mismatch now errors immediately instead of producing
  quietly-wrong numbers.
- **`optimization/buildDenseSwitchMap.m` need only be built ONCE per (capArea,rodArea,switchTime)
  combination, not once per DHD rail count** — DHD2/DHD3/DHD4 share identical capArea/rodArea (rail
  count is a control-layer choice, not a cylinder-sizing one), so one dense map is exactly valid for
  all three. The original overnight sweep built three identical maps (~56 wasted minutes); current
  scripts build one and reuse it.

## EHA cylinder resized 23.37in -> 20in bore (deliberate, not a bug fix)

EHA's cylinder was originally sized to its literal peak force demand (9.688 MN at the sizing sea
state, giving capArea=0.276822 m^2 via `sizeCylinderArea.m`). Checking the force trajectory showed
that peak is a rare spike, not a sustained requirement: only 0.093% of (post-ramp) time exceeds 90%
of it, RMS force is 2.595 MN, and the peak/RMS ratio is 3.73 (a pure sine wave would be ~1.41) --
this is a smooth baseline demand with a brief, extreme spike layered on top. Deliberately resized to
a 20in bore (capArea=rodArea=0.202683 m^2, 7.094 MN capacity at 35 MPa) to trade clipping that rare
spike for a much smaller cylinder. `control/MPC_QP.m`'s EHA case now saturates its unconstrained
continuous force at `capArea*hyd.maxPressure` (new field, `parameters/getHydraulic.m`, defaults to
35 MPa) instead of letting the QP solve imply an unphysical pressure -- this saturation is a no-op
whenever EHA is sized to its literal peak (never binds), so it's safe to leave in place regardless
of which sizing convention is used.

## EHA copper-loss coefficient recalibrated (5e-5 -> 1.868e-4)

`models/makeEHALossMap.m` and `makeEHALossMap_fixedDisp.m`'s copper-loss term
(`P_L_elect = coeff*T_Act^2`) used `coeff=5e-5`, with a comment claiming this gave "10kW at
1MNm." Checked directly: it doesn't -- 5e-5 gives 50MW at 1MNm of PUMP-SHAFT torque (T_Act), a
5000x mismatch with that comment. But T_Act is not flap torque (a common point of confusion
during this check -- they're related through the pump's own displacement/gearing, not the
flap's 2.7574 moment arm), and at 1MNm of FLAP torque specifically, the real T_Act is much
smaller, so the coefficient's actual behavior at realistic operating conditions is what
matters, not the stale comment. Checked that too: at real sea-state-8 flap torques (median
3.96MNm, 90th pct 11.0MNm, max 14.4MNm), a coefficient large enough to hit "10kW at 1MNm flap
torque" (2.316e-3) blows up to copper loss exceeding 100% of input power at the 90th-percentile
torque -- unusable. The original 5e-5 stayed physically sensible across that same real range,
suggesting the "10kW at 1MNm" comment described a smaller reference design from before the
flap/cylinder was resized to handle multi-MN forces, not a bug in the code.
Separately, the resulting overall EHA efficiency (85-88%) was judged too optimistic for a
real hydraulic+electric conversion chain running a WEC's mostly-partial-load, frequently-
reversing duty cycle (22-26% of the time spent with mechPower<=0 -- see the earlier
mean-vs-median instantaneous-efficiency discussion). Recalibrated to 5% copper loss at rated
power (vMax speed, 35MPa) instead: at the current 20in-bore EHA sizing that's 354.7kW at
7.094MW rated (real T_Act=43,576Nm there), giving coeff=1.868e-4 -- a modest 3.74x change from
the original 5e-5, not the 46x the stale comment would have implied.

## EHA shaft inertia and copper-loss coefficient rebuilt from real component data (1.868e-4 -> 5.0e-4 copperCoeff; shaftInertia 5 -> 56.2)

The 5% -of-rated-power derivation above was itself superseded after checking the main-motor
efficiency question: is the EHA's electric machine's efficiency comparable to the flat 85%
assumed for PassivePump/DHD's main motor? That led to sizing what the EHA's electric side would
actually look like, then pulling real component data instead of guessing at a target power
fraction.

**Hydraulic pump/motor side.** The EHA's fixed-displacement pump/motor, sized for the 20in-bore
cylinder's max flow (capArea*vMax) at the assumed 2000 RPM nominal design speed (20% margin),
works out to an effective displacement of ~7,297 cc/rev -- larger than any single real catalog
pump. Modeled instead as 14 Bosch Rexroth A4VSO500 axial piston pumps in parallel on a common
shaft (500 cc/rev each, real catalog part, close to Bosch Rexroth's largest standard size):
14 x 0.3325 kgm^2 (Rexroth RE 92050-01-X-B2/2019-08-23 datasheet, "Moment of inertia of the
rotary group J_TW") = 4.655 kgm^2, rounded to **5 kgm^2**.

**Electric motor/generator side.** Sized to match the aggregate pump power (14 x ~385-437kW =
5.4-6.1MW), modeled as 2 real ABB HXR 450LB-class 750kW/1492rpm induction motors (real catalog
row, ABB "High Voltage Induction Motors" technical catalogue EN 12-2007, p.25: 450LB frame,
T_N=4800 Nm, I_N=169A, full-load eff=97.0%, rotor inertia=25.6 kgm^2 each) instead of one
hypothetical giant motor: 2 x 25.6 = **51.2 kgm^2**. (One giant motor extrapolated to the same
aggregate power would come out to ~330-400 kgm^2 -- splitting power across several real,
catalog-sized units gives substantially lower combined inertia than one oversized unit, on both
the pump and motor sides, because rotating inertia scales roughly as size^(5/3) for a single
unit but only linearly with unit count for N identical parallel units.)

**Total: shaftInertia = 5 + 51.2 = 56.2 kgm^2** (was 5, an assumed placeholder with "no measured
value" per the old comment in getPhysical.m -- now built up from two real catalog components).

**Copper-loss coefficient**, same 750kW/1492rpm motor: total loss at rated load =
P_out*(1/eff-1) = 750000*(1/0.970-1) = 23.2kW, giving coeff_total = P_loss/T_N^2 =
23196/4800^2 = 1.007e-3 -- but that's *all* loss mechanisms (copper+iron+friction/windage+
stray), not copper alone. Using a typical ~50% copper-loss fraction for this class of large
induction motor: **coeff = 0.5 * 1.007e-3 ~= 5.0e-4**, the new default in both
`makeEHALossMap.m` and `makeEHALossMap_fixedDisp.m` (was 1.868e-4). Checked across the catalog's
full 110-750kW range: this total-loss-based coefficient shrinks steadily with motor size (from
1.24e-2 at 110kW down to 1.01e-3 at 750kW, since bigger motors are proportionally more
efficient) -- so 5.0e-4 is specifically calibrated around the 750kW size just adopted, not a
universal constant; it would need revisiting if the assumed motor size changes again.

Source for both the inertia and copper-loss numbers: ABB "High Voltage Induction Motors"
Technical Catalogue, EN 12-2007 (ABB/BU Machines) -- one single PDF covering both the 110-750kW
"Process Performance" cast-iron tables and the up-to-2800kW "Engineered motors" (HXR/AMA)
tables used to compare motor-count options.

## Electric motor count revised 2 -> 5 (shaftInertia 56.2 -> 133.0), motivated by peak vs. average power

The 2x750kW sizing above matched the aggregate PUMP power rating (~5.4-6.1MW aggregate from 14
pumps), but was never checked against the EHA's actual instantaneous power demand. It wasn't
close: SS8's average electrical power is only ~624kW (comfortably within one 750kW motor's
continuous rating), but instantaneous power in the same run peaks around 5,000-5,600kW, and the
absolute design-point peak (capArea*vMax*35MPa) is 7.094MW. Unlike HHEA (buffered by common
pressure rail accumulators, main pump/motor sized for mean power) or the check-valve PTO
(buffered by its accumulators too), a direct-drive EHA has no buffering at all -- every
instantaneous power swing passes straight through the motor (see the ASME/BATH 2021 HHEA-PTO
paper, WEC_PTO_comparison.pdf, Section 4.3: "the power is not smoothed in any way... the
electrical components need to be sized for the peak power").

Settled on 5 x 750kW motors (same ABB HXR 450LB-class unit as above, just more of them):
continuous rating 5*750kW=3750kW covers roughly half of the observed/design peak, with a 2x
short-term overload (7500kW) assumed to cover the rest -- a real, if aggressive, motor
overload assumption, not yet enforced anywhere in the simulation (see below).

New electric motor side: 5 x 25.6 kgm^2 = **128.0 kgm^2**.
**New total: shaftInertia = 5 + 128.0 = 133.0 kgm^2** (was 56.2).
New reflected inertia (same k~=481.3 gear ratio, unaffected by motor count): ~30.8 million
kgm^2, now 1.38x (phys.I+phys.Iinf) -- pushes the flap's natural period from ~12.44s (bare flap)
to ~19.18s, even further past every sea state's period (max 13.19s) than the 2-motor case's
~15.65s.

**Not yet implemented: an actual power cap reflecting the 5-motor (3750kW continuous / 7500kW
short-term) assumption.** The simulation currently lets the controller demand however much
electrical power it wants, unconstrained by any motor capacity -- only the CYLINDER's mechanical
force is capped (`MPC_QP.m`'s `capArea*hyd.maxPressure` saturation, parameters/getHydraulic.m).
Adding a genuine electrical-power cap to the control law is a materially harder problem than the
existing force cap: power is velocity*force (bilinear in state and control), not a simple bound
on the control alone, so it likely needs either a QP reformulation or a post-hoc saturation
applied consistently with an energy-balance check, not just clipped after the fact. Flagged as
follow-up work; the copperCoeff/shaftInertia numbers above are valid independent of whether/how
this cap gets added.

## Electric motor switched to one trend-extrapolated 4000kW unit (shaftInertia 133.0 -> 226.3)

Rather than summing 5 real 750kW catalog motors, refit the ABB HV catalog's rotor-inertia-vs-
power data (same 110-750kW dataset as above) as a power law: J = a*P^b, giving
a=2.4045e-3, b=1.3781 (log-log least-squares fit; see
results/motorInertiaTrend.png for the fitted curve against the real data points). Evaluated at
one 4000kW motor (close to the 3750kW continuous-rating target, rounded up): **J = 221.3 kgm^2**
-- notably larger than the 5x750kW real-component sum (128.0 kgm^2) would have given, consistent
with the pattern seen throughout this analysis (one big unit costs more combined inertia than
splitting the same total power across several real, smaller catalog units, since inertia scales
roughly as P^1.38 for a single unit but only linearly with count for N identical units). This is
a genuine extrapolation -- 750kW is the largest real datapoint in the catalog, so 4000kW is more
than 5x beyond it -- trusting the fitted exponent rather than a real catalog part, unlike the
5x750kW answer it replaces.

New electric motor side: **221.3 kgm^2** (one 4000kW motor, trend-extrapolated).
**New total: shaftInertia = 5 + 221.3 = 226.3 kgm^2** (was 133.0).
New reflected inertia (same k~=481.3 gear ratio): ~52.4 million kgm^2, now 2.34x
(phys.I+phys.Iinf) -- pushes the flap's natural period further still, to ~22.75s (from ~12.44s
bare-flap, ~15.65s at the original 2-motor sizing, ~19.18s at the 5-motor sizing), continuing the
same trend of moving resonance further from every sea state's period (max 13.19s) as the assumed
electric motor capacity grows to cover more of the peak.

The unenforced-power-cap caveat above still applies unchanged: nothing in the simulation yet
constrains electrical power to the assumed 3750-4000kW continuous / ~7500-8000kW short-term
capacity.

## Hydraulic pump switched to one trend-extrapolated ~7297cc/rev unit (shaftInertia 226.3 -> 249.9)

Same treatment as the electric motor above, applied to the pump side: rather than summing 14
real 500cc/rev Rexroth A4VSO pumps, refit that family's rotor-inertia-vs-displacement data (same
7-point dataset used originally, 40-500cc/rev) as a power law: J = a*D^b, giving a=9.60e-6,
b=1.676 (log-log least-squares fit; see results/ehaInertiaTrends_combined.png for the fitted
curve against the real data points, alongside the equivalent motor-side plot). Evaluated at the
pump's own actual required displacement (~7297cc/rev, from capArea*vMax at the 2000rpm nominal
design point -- unchanged from the original derivation): **J = 28.6 kgm^2** -- notably larger
than the 14x500cc/rev real-component sum (5 kgm^2) would have given, a bigger jump than the
motor side saw (1.7x, 128.0->221.3) because more real units are being replaced here (14 vs 5),
so the "sum-of-many-small" vs "one-big-extrapolated" gap is proportionally larger. Same
extrapolation caveat as before, more so: 500cc/rev is the largest real datapoint in the family,
so 7297cc/rev is nearly 15x beyond it (vs the motor's ~5x).

New hydraulic pump/motor side: **28.6 kgm^2** (one ~7297cc/rev pump, trend-extrapolated).
**New total: shaftInertia = 28.6 + 221.3 = 249.9 kgm^2** (was 226.3).
New reflected inertia (same k~=481.3 gear ratio): ~57.9 million kgm^2, now 2.59x
(phys.I+phys.Iinf) -- pushes the flap's natural period further still, to ~23.57s (from ~12.44s
bare-flap, ~15.65s/~19.18s/~22.75s at the previous three sizing iterations), continuing the same
trend: every step toward trusting the fitted trend over summing real discrete units, on either
side of the shaft, has pushed resonance further from every sea state's period (max 13.19s).

## DHD/PassivePump main-motor efficiency raised 85% -> 90% (evaluate.m)

Motivated by comparing DHD/PassivePump's flat main-motor efficiency against EHA's own copper-
loss physics: DHD's key architectural advantage is that its main motor runs at steady, buffered
conditions (via the pressure-rail system), unlike EHA's continuously-varying operating point --
so a flat 85% (the same number used regardless of that difference) was plausibly conservative
for DHD/PassivePump specifically. Changed `mechToElecLoss = 0.15*eval.aveElecPow` to
`0.10*eval.aveElecPow` in evaluate.m (both DHD and PassivePump cases share this line).

Since this is a flat multiplicative derate, its effect is exactly computable without rerunning
anything: every DHD/PassivePump elecRGP (and aveElecPow) scales by 0.90/0.85 = 1.0588 (a uniform
+5.88% relative increase), independent of sea state or family. Checked against the real
mAstarConvAndPressureOpt sweep results (partial coverage at the time: DHD2 7/8, DHD3 4/8, DHD4
1/8 sea states) and the current EHA_elec numbers (shaftInertia=249.9): DHD3 already exceeds
EHA_elec at SS2/3/4 and DHD4 at SS1 even before this change (at the old 85%); after the change,
DHD2 also edges ahead at SS8. PassivePump and DHD2 still trail EHA_elec at most sea states even
at 90% -- the gap there is too large for this one change to close. DHD3's SS3 result is likely
even better than the mAstar=5 sweep shows: at the more fully-converged mAstar=12 (own SS3-optimal
pressure, 23MPa), DHD3 already got elecRGP=0.691 at the old 85%, comfortably ahead of EHA_elec's
0.661 -- meaning mAstar=5 alone (used in the pressure sweep) was mildly under-converged for DHD3.

**Reverted back to 85%** after the EHA copper-loss `sign(T_Act)` fix below: the comparison that
motivated bumping DHD/PassivePump to 90% was against EHA_elec numbers that turned out to be
inflated ~20-22% (and EHA_mech far more, -55 to -71% at SS7/8) by that bug -- EHA was never as
strong a comparison point as it looked, so the case for treating DHD's buffered main motor as
meriting a better-than-EHA flat efficiency evaporates along with it. `mechToElecLoss` is back to
`0.15*eval.aveElecPow` in evaluate.m; every DHD/PassivePump elecRGP (and aveElecPow) scales back
down by 0.85/0.90 = 0.9444 (a uniform -5.56% relative decrease) from the 90%-derived numbers
above.

## EHA copper-loss sign(T_Act) removed (models/makeEHALossMap.m, makeEHALossMap_fixedDisp.m)

Root-caused via the hydraulic/electric stage efficiency maps (Run_EHAHydraulicElectricEfficiencyMaps.m):
`P_L_elect = copperCoeff*T_Act^2*sign(T_Act)` was wrong on its face -- I^2*R heat is never
negative (current appears squared), so the physical loss magnitude is `copperCoeff*T_Act^2`,
unsigned. The `sign(T_Act)` was apparently meant as a stand-in for motoring/generating direction
when combining the loss into `P_out = w*T_Act + P_L_elect`, but `sign(T_Act)` only tracks
`sign(deltaP)` (torque), not the real direction `sign(Q*deltaP) = sign(thetaDot)*sign(deltaP)`
(speed AND torque). Traced algebraically and confirmed at 4 concrete grid points (+-0.3 rad/s x
+-20MNm): because `sign(w)~=sign(thetaDot)` away from the low-speed leakage regime, and the real
direction is `sign(thetaDot)*sign(deltaP)`, the thetaDot-dependence cancelled out of the
self-consistency check entirely, leaving whether the old formula was "correct" or "backwards"
purely a function of `sign(deltaP)` (torque) alone -- explaining why the electric-efficiency map
showed the *entire* negative-torque half-plane broken (elecEff>1) regardless of speed or actual
direction: 1,856 of 3,600 grid points (51.6%).

Fix: removed `*sign(T_Act)` from `P_L_elect` in both files. After the fix, elecEff>1 points
dropped to 118/3600 (3.3%) -- the electric-efficiency map is now symmetric top-to-bottom
(matching the hydraulic map's shape) instead of having one whole half artificially pinned at 1.0.
The remaining 118 points are a *different*, still-open issue: the low-speed leakage regime where
`QLoss` can flip `sign(w)` away from `sign(thetaDot)` (confirmed unchanged by this fix at the
specific point traced earlier, since `T_Act` was already positive there, so `sign(T_Act)` was
already +1 -- this fix only changes points where `T_Act`'s sign actually flips something).

This changes the EHA's actual loss physics (not just a diagnostic script), so every EHA
result generated before this fix (SS8/SS3 spot checks, the full 8-sea-state
Run_EHA_AllSeaStates.m runs, etc.) reflects the old, direction-inconsistent formula and should
be treated as stale once this fix is exercised in a real simulation run (not yet re-run as of
this note -- only the standalone efficiency-map diagnostic has been checked against it so far).

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
- **`control/MPC_Astar_cont2.m` was unfinished** (referenced `params.uInd`/`params.damping`,
  neither ever set anywhere, and was never wired into `controlLaw.m`'s dispatch) — deleted as
  dead code during the directory restructuring.
