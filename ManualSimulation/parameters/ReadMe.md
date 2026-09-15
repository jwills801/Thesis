This folder holds the true parameter-builders: getParameters.m (calls the
other three), getPhysical.m, getHydraulic.m, getSimulation.m.

Some parameters change with the case that we are running.
These are held in the runParams structure.

getHydraulic.m also calls out to ../models/ (makeSwitchLossMap.m,
makeEHALossMap.m / makeEHALossMap_fixedDisp.m) and ../optimization/
(evenlySpacedRails.m) -- those physics/loss models and optimization
utilities live in their own top-level folders now, not here, since they're
a different kind of thing than parameter setup. See the repo-root
README.md for the full structure.

SwitchMap.mat (in this folder) is a precomputed cache getHydraulic.m loads
by default for DHD -- it's only valid for whatever pressure rails it was
built at, and goes stale silently if you change highPressure without
regenerating it. optimization/buildDenseSwitchMap.m is a better
alternative (builds one reusable dense-grid map instead).