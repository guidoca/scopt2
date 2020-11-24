function [J,L,Jreal,Lreal] = computeObjectiveFunctionsNormal(obj,time_n,X_n,U_n,time_n_p,X_n_p,U_n_p)
% Computes scaled and unscaled (real) nonlinear (J) and convexificed
% (L) objective functions
obj.problem.phases(1).BodyMap.evaluated = 0;

X            = obj.scaleStatesToReal(1,X_n);
U            = obj.scaleControlsToReal(1,U_n);
timeFinal    = obj.scaleTimeToReal(1,time_n(end));
timeInitial  = obj.scaleTimeToReal(1,time_n(1));

timeDisc      = (timeFinal - timeInitial)/(obj.problem.phases.Nodes-1);
time          = (timeInitial:timeDisc:timeFinal);

if nargin>4
    X_p            = obj.scaleStatesToReal(1,X_n_p);
    U_p            = obj.scaleControlsToReal(1,U_n_p);
    timeFinal_p    = obj.scaleTimeToReal(1,time_n_p(end));
    timeInitial_p  = obj.scaleTimeToReal(1,time_n_p(1));
    
    timeDisc_p      = (timeFinal_p - timeInitial_p)/(obj.problem.phases.Nodes-1);
    time_p          = (timeInitial_p:timeDisc_p:timeFinal_p);
    [Jreal,Lreal] = computeObjectiveFunctions(obj,time,X,U,time_p,X_p,U_p);
    J             = Jreal/obj.problem.objective.scale*obj.problem.objective.sign;
    L             = Lreal/obj.problem.objective.scale*obj.problem.objective.sign;
else
    [Jreal]       = computeObjectiveFunctions(obj,time,X,U);
    Lreal         = 0;
    J             = Jreal/obj.problem.objective.scale*obj.problem.objective.sign;
    L             = 0;
end
end

