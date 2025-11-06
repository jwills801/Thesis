clear
unsolvedNodes = {};
solvedNodes = {};
U = [1;3;5;7];
unsolvedNodes{1}.lb = [min(U) min(U)];
unsolvedNodes{1}.ub = [max(U) max(U)];
bestCost = 1e3;
flag = 0;
iter = 1;
iterMax = 10;
while flag == 0
% take next node to solve
    node = unsolvedNodes{1};

% solve node
    node = solveNode(node);

% Remove solved node
    unsolvedNodes(1) = [];
    solvedNodes{end+1} = node;

% At the beginning, round to the nearest node to get an integer solution
if iter == 1
    [M,I] = min(abs(solvedNodes{1}.u - [U U]));
    quickNode.lb = U(I)';
    quickNode.ub = U(I)';
    quickNode = solveNode(quickNode);
    bestCost = quickNode.cost
    solvedNodes{end+1} = quickNode;
    bestNodeInd = length(solvedNodes);
end

% Prune or branch
    if node.cost > bestCost % prune
        disp(['Pruned at node ', num2str(length(solvedNodes))])
    else % branch
        node.integer = min(abs(node.u - [U U])) < .001
        branchingVariable = find(~node.integer,1);
        if isempty(branchingVariable) % Integer solution
            % This is the best cost. If it wasn't we wouldn't be in the
            % branch part of this if statement
            bestCost = node.cost
            bestNodeInd = length(solvedNodes)
        else
        newNodeHigh.ub = node.ub;
        newNodeHigh.lb = node.lb;
        newNodeLow = newNodeHigh;
    % find which index to branch
        tmp = node.u(branchingVariable) - U;
        I = find(tmp>0);
        % I and (I+1) are the indices of the discrete input set that we're
        % bounding on
        newNodeHigh.lb(branchingVariable) = U(I+1);
        newNodeLow.ub(branchingVariable) = U(I);

        unsolvedNodes{end+1} = newNodeHigh;
        unsolvedNodes{end+1} = newNodeLow;
        end
    end

iter = iter + 1
% Check to see if we should stop
    if isempty(unsolvedNodes)
        flag = 1;
        disp('Solution found')
        bestNode = solvedNodes{bestNodeInd}
    elseif iter > iterMax
        flag = -1;
        disp('Max iterations')
    end
end

%%
uVals = 0:.5:10;
[uMat1, uMat2] = ndgrid(uVals,uVals);
costMat = NaN(size(uMat1));
for uInd = 1:length(uMat1(:))
costMat(uInd) = dynamics([uMat1(uInd),uMat2(uInd)],0);
end
figure, surf(uMat1,uMat2,costMat), xlabel('u_1'), ylabel('u_2'), zlabel('Cost')

function node = solveNode(node)
% options = optimoptions('fmincon','PlotFcn','optimplot');
options = optimoptions('fmincon');
A = []; b = [];
Aeq = []; beq = [];
lb = node.lb;
ub = node.ub;
nonlcon = []; params = 0;
u0=zeros(size(lb));
[uOpt,costOpt,exitflag,output] = fmincon(@(u) dynamics(u,params),u0,A,b,Aeq,beq,lb,ub,nonlcon,options);
node.u = uOpt;
node.cost = costOpt;
end

function cost = dynamics(u,params)
cost = 1;
for i = 1:length(u)
    cost = cost + (u(i)-2*i^2-.5).^2;
end
end