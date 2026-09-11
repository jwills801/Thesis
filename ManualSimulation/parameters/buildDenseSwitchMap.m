% buildDenseSwitchMap.m
% Builds a switch-loss map over a dense, fixed pressure grid (not the
% actual 2/3/4-rail set in use) so it can be reused across every sea
% state and candidate highPressure for a given (capArea,rodArea) family,
% instead of regenerating per candidate -- valid because consumers
% (getValveLoss.m, MPC_Astar.m, MPC_Astar_cont.m) already interpn against
% whatever switchMap.PR is given, not against the actual rails in use.
% Calls: parameters/makeSwitchLossMap.m
% Called by: none yet within this repo (intended caller: build once per
%   drivetrain family, then pass the result into optimizePressure.m's
%   optional denseSwitchMap argument for its per-sea-state sweep)
function switchMap = buildDenseSwitchMap(hyd,pLow,pHigh,nPoints)
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
if nargin<4, nPoints=10; end
hyd.pressureRails = linspace(pLow,pHigh,nPoints);
switchMap = makeSwitchLossMap(hyd);
end
