% getUtilityMatrices.m
% Builds L (stretches an m-length coarse control onto the n*m fine grid)
% and C (selects the thetaDot channel at every fine step, block-arranged
% per coarse step) -- shared building blocks for the MPC_QP/MPC_Astar
% cost matrices.
% Calls: none
% Called by: control/getControl.m, diagnostics/checkMPC_EHA.m,
%   checkTerminalCost.m, checkUtilityMatrices.m
function [L,C] = getUtilityMatrices(m,n)
% L stretches an m-length coarse control vector onto the n*m fine time
% grid (each coarse control held for n fine steps).
% C picks out the theta-dot channel (e1=[1;0;0;0]) at every fine step, in
% the same block-diagonal layout used to weight u against X in E_mech/E_elec.
onesCol = ones(n,1);
zerosCol = zeros(n,1);

ei = [1;0;0;0];
c_block = repmat(ei,n,1);

L = NaN(n*m,m);
C = zeros(4*n*m,m);
for col = 1:m
    L(:,col) = [repmat(zerosCol,col-1,1);
        onesCol;
        repmat(zerosCol,m-col,1)];
    C(:,col) = [repmat(0*c_block,col-1,1);
        c_block;
        repmat(0*c_block,m-col,1)];
end
end
