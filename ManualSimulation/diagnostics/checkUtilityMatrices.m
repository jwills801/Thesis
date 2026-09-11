% checkUtilityMatrices.m
% Verifies control/getUtilityMatrices.m's L (stretching matrix) and C
% (e1 block-diagonal) against direct construction from their definitions.
% Calls: control/getUtilityMatrices.m
% Called by: none (top-level diagnostic script, run manually)
clear; clc

here = fileparts(mfilename('fullpath')); root = fileparts(here);
addpath(fullfile(root,'control'));

m = 4; n = 3; % small coarse/fine horizon sizes
[L,C] = getUtilityMatrices(m,n);

% L should hold each coarse control constant for n fine steps
L_expected = zeros(n*m,m);
for col = 1:m
    L_expected((col-1)*n+1:col*n, col) = 1;
end

% C should place n copies of e1=[1;0;0;0] in the block for control col
e1 = [1;0;0;0];
C_expected = zeros(4*n*m,m);
for col = 1:m
    C_expected((col-1)*4*n+1:col*4*n, col) = repmat(e1,n,1);
end

errL = max(abs(L(:)-L_expected(:)));
errC = max(abs(C(:)-C_expected(:)));

fprintf('--- checkUtilityMatrices ---\n');
report('L matches direct stretching-matrix construction',errL,0);
report('C matches direct block-e1 construction',errC,0);

function report(name,err,tol)
if err <= tol
    fprintf('PASS  %-55s (err=%.3g)\n',name,err);
else
    fprintf('FAIL  %-55s (err=%.3g, tol=%.3g)\n',name,err,tol);
end
end
