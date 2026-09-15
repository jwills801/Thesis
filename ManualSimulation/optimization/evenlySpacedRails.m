% evenlySpacedRails.m
% Places n=3 or n=4 pressure rails between fixed endpoints pLow/pHigh so
% the resulting n^2 pairwise DHD force options (cap-rail x rod-rail) are
% as evenly spaced as possible (grid search minimizing variance of sorted
% force gaps). n=2 has no interior rail to place.
% Calls: none
% Called by: parameters/getHydraulic.m
function pressureRails = evenlySpacedRails(pLow,pHigh,capArea,rodArea,n)
% Places n pressure rails between pLow and pHigh (both fixed endpoints, Pa)
% so that the resulting DHD force options -- all n^2 pairwise combinations
% of (cap rail, rod rail), since cap and rod sides switch independently
% among the same n rails, see getHydraulic.m's ptoForceOptions -- are as
% evenly spaced as possible. Only the interior (non-endpoint) rails are
% free; endpoints are exactly pLow and pHigh, unchanged from before.
%
% n=2 has no interior rails (nothing to place). n=3 has 1 free interior
% rail (search over a fine 1-D grid). n=4 has 2 free interior rails
% (search over a fine 2-D grid, parametrized as fractions of the
% remaining range so the ordering pLow<mid1<mid2<pHigh is automatic).
% Both grids are cheap -- no simulation involved, just arithmetic on a
% small force matrix -- so a plain grid search is used rather than a
% gradient-based optimizer, for robustness against local minima.

switch n
    case 2
        pressureRails = [pLow,pHigh];
    case 3
        candidates = linspace(pLow,pHigh,501);
        candidates = candidates(2:end-1);
        costs = arrayfun(@(pMid) evennessCost([pLow pMid pHigh],capArea,rodArea), candidates);
        [~,idx] = min(costs);
        pressureRails = [pLow, candidates(idx), pHigh];
    case 4
        ns = 80;
        s1 = linspace(0.01,0.99,ns);
        s2 = linspace(0.01,0.99,ns);
        bestCost = Inf; bestRails = [];
        for i = 1:ns
            pMid1 = pLow + s1(i)*(pHigh-pLow);
            for j = 1:ns
                pMid2 = pMid1 + s2(j)*(pHigh-pMid1);
                cost = evennessCost([pLow pMid1 pMid2 pHigh],capArea,rodArea);
                if cost < bestCost
                    bestCost = cost; bestRails = [pLow pMid1 pMid2 pHigh];
                end
            end
        end
        pressureRails = bestRails;
    otherwise
        error('evenlySpacedRails:unsupportedN','Only n=2,3,4 supported (n=%d requested).',n);
end
end

function cost = evennessCost(pressureRails,capArea,rodArea)
capForceOptions = pressureRails * capArea;
rodForceOptions = pressureRails * rodArea;
forces = capForceOptions' - rodForceOptions;
forces = sort(unique(forces(:)));
gaps = diff(forces);
cost = var(gaps);
end
