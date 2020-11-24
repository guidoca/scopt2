function [Jaug,Laug] = computeAugmentedObjectiveFunctions(obj,time_n,X_n,U_n,time_n_p,X_n_p,U_n_p,linearToo)
if nargin<8
    linearToo = 1;
end
[penaltyDefectsNonLinear,penaltyBuffersNonLinear,penaltyTrustRegion,penaltyTime,penaltyDefectsLinear,penaltyBuffersLinear] = computeAugmentedObjectiveFunctionTerms(obj,time_n,X_n,U_n,time_n_p,X_n_p,U_n_p,linearToo);

Reg = computeRegularisation(obj,time_n,X_n,U_n) ;
if linearToo
    [J,L] = computeObjectiveFunctionsNormal(obj,time_n,X_n,U_n,time_n_p,X_n_p,U_n_p);
    Laug = L + Reg + penaltyDefectsLinear    + penaltyBuffersLinear + penaltyTrustRegion + penaltyTime ;
else
    [J] = computeObjectiveFunctionsNormal(obj,time_n,X_n,U_n);
end


Jaug = J + Reg + penaltyDefectsNonLinear + penaltyBuffersNonLinear + penaltyTrustRegion + penaltyTime ;
end

