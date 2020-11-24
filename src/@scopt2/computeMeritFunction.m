function M = computeMeritFunction(obj,searchDirection,stepLength,Z_p)
% computes merit function on new point
X_n_p             = Z_p(obj.problem.phases.index.states);
U_n_p             = reshape(Z_p(obj.problem.phases.index.controls),size(obj.problem.phases.index.controls));
timeInitial_n_p   = Z_p(end-1);
timeFinal_n_p     = Z_p(end);
timeDisc_n_p      = (timeFinal_n_p - timeInitial_n_p)/(obj.problem.phases.Nodes-1);
time_n_p          = (timeInitial_n_p:timeDisc_n_p:timeFinal_n_p);

% Extract Points
Z             = computeNewPoint(obj,searchDirection,stepLength,Z_p);
X_n           = Z(obj.problem.phases.index.states);
U_n           = reshape(Z(obj.problem.phases.index.controls),size(obj.problem.phases.index.controls));
timeInitial_n = Z(end-1);
timeFinal_n   = Z(end);

timeDisc_n   = (timeFinal_n - timeInitial_n)/(obj.problem.phases.Nodes-1);
time_n       = (timeInitial_n:timeDisc_n:timeFinal_n);

% Compute Merit Function at Estimate Point
M          = computeAugmentedObjectiveFunctions(obj,time_n,X_n,U_n,time_n_p,X_n_p,U_n_p,0);

end

