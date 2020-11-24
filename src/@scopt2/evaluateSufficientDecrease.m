function [sufficientDecrease] = evaluateSufficientDecrease(obj,searchDirection,stepLength,directionalDerivative,Z_p,M_p)
% Evaluates the sufficient decrease condition for the
% contraction line search type

% Compute Merit Function at Estimate Point
M       = computeMeritFunction(obj,searchDirection,stepLength,Z_p);

% Compute Change
MChange    = (M_p - M)/stepLength; % should be positive for a decrease
% Compute Lower Bound
MChangeLow = -obj.algorithm.sequential.globalisation.lineSearch.k1*directionalDerivative;
% Compute Upper Bound
MChangeUpp = -obj.algorithm.sequential.globalisation.lineSearch.k2*directionalDerivative;

sufficientDecrease = (MChangeLow <= MChange) && (MChange <= MChangeUpp);
end

