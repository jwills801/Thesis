% drivetrainStyle.m
% Central color/marker/line-style convention for every drivetrain family
% plotted in this repo, so plotting scripts stay visually consistent
% without each one re-picking its own palette.
%
% Color is fixed per top-level drivetrain type. DHD2/DHD3/DHD4 are a
% single sequential orange progression (more rails = darker), so they
% read as one related family; PassivePump and EHA each get their own
% hue. Chosen to stay distinguishable under red-green colorblindness (the
% most common form): blue/teal/orange span the blue-yellow axis, which
% colorblindness mostly preserves, and the DHD progression stays within
% orange -- never crossing into saturated red or purple -- so it never
% approaches EHA's teal.
%
% Marker shape is fixed per family too, so identification never depends
% on color alone (redundant channel, standard accessible-chart practice).
%
% Line style is NOT used for data completeness -- reserved instead for
% comparing different controllers/variants WITHIN the same drivetrain
% (e.g. EHA fixed- vs. variable-displacement, or elec- vs. mech-
% optimized) while keeping that drivetrain's color/marker fixed. Pass
% variant='secondary' for the alternate/de-emphasized variant (dashed);
% omit or pass 'primary' for the default (solid).
%
% Usage: s = drivetrainStyle('DHD3'); or s = drivetrainStyle('EHA','secondary');
%        s.color (1x3 RGB), s.marker (char), s.lineStyle ('-' or '--')
%
% Calls: none
% Called by: plotting scripts (plotMainResultsBySeaState.m,
%   plotMainResultsPressureSweeps.m, etc.)
function style = drivetrainStyle(family, variant)
if nargin < 2, variant = 'primary'; end

switch family
    case 'PassivePump'
        color = [0 114 178]/255; marker = 'o';
    case 'DHD2'
        color = [253 190 133]/255; marker = 's';
    case 'DHD3'
        color = [253 141 60]/255; marker = '^';
    case 'DHD4'
        color = [217 71 1]/255; marker = 'd';
    case 'EHA'
        color = [0 158 115]/255; marker = 'p';
    otherwise
        error('drivetrainStyle:unknownFamily', 'Unknown drivetrain family "%s"', family);
end

switch variant
    case 'primary'
        lineStyle = '-';
    case 'secondary'
        lineStyle = '--';
    otherwise
        error('drivetrainStyle:unknownVariant', 'variant must be "primary" or "secondary", got "%s"', variant);
end

style = struct('color',color,'marker',marker,'lineStyle',lineStyle);
end
