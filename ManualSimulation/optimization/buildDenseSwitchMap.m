% buildDenseSwitchMap.m
% Builds a switch-loss map over a dense, fixed pressure grid (not the
% actual 2/3/4-rail set in use) so it can be reused across every sea
% state and candidate highPressure for a given (capArea,rodArea) family,
% instead of regenerating per candidate -- valid because consumers
% (getValveLoss.m, MPC_Astar.m, MPC_Astar_cont.m) already interpn against
% whatever switchMap.PR is given, not against the actual rails in use.
% Also reusable ACROSS DHD2/DHD3/DHD4 as long as they share the same
% capArea/rodArea (true in results/sizedAreas.mat -- rail count is a
% control-layer choice, not a cylinder-sizing one), since the map depends
% only on capArea/rodArea/stroke, not on how many rails will later be
% selected from its dense PR grid. Don't rebuild one per rail count.
% Calls: models/makeSwitchLossMap.m
% Called by: none yet within this repo (intended caller: build once per
%   drivetrain family, then pass the result into optimizePressure.m's
%   optional denseSwitchMap argument for its per-sea-state sweep)
function switchMap = buildDenseSwitchMap(hyd,pLow,pHigh,nPoints,switchTime)
% Builds a switch-loss map (makeSwitchLossMap.m) over a dense, fixed
% pressure grid spanning [pLow,pHigh], instead of the actual small rail
% set (2/3/4 values) a given DHD case uses. Callers (getValveLoss.m,
% MPC_Astar.m, MPC_Astar_cont.m) already query this map via interpn using
% whatever the ACTUAL rail pressures are at the time -- switchMap.PR is
% purely an interpolation grid and never has to equal those rails. So one
% dense map, built once per (capArea,rodArea) family, can be reused for
% every sea state and every candidate highPressure/rail-spacing for that
% family instead of regenerating per candidate.
%
% hyd must already have the family's correctly-sized capArea/rodArea/stroke
% (hyd.pressureRails is overwritten here, only used to size the map axis).
% switchTime (optional, default 0.2s): the coarse control step this map
% is valid for -- see makeSwitchLossMap.m's comment. Must match whatever
% control/getControl.m's ctrl.timeHorizon is set to for the run this map
% will be used in.
if nargin<4, nPoints=10; end
if nargin<5, switchTime=0.2; end
hyd.pressureRails = linspace(pLow,pHigh,nPoints);
hyd.switchTime = switchTime;
switchMap = makeSwitchLossMap(hyd);
end
