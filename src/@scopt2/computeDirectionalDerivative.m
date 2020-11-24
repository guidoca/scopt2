function [directionalDerivative] = computeDirectionalDerivative(obj,searchDirection,Z_p,M_p)
% Computes directional derivatives for the contraction type
% line search routine
switch obj.algorithm.sequential.globalisation.lineSearch.directionalDerivativeType
    case 'numeric'
        epsilon     = obj.algorithm.sequential.globalisation.lineSearch.directionalDerivativeEpsilon;
        M_eps       = computeMeritFunction(obj,searchDirection,epsilon,Z_p);
        % Compute Directional Derivative
        directionalDerivative = (M_eps - M_p)/epsilon;
    case 'analytical'
        error('SCOPT Error, Analytical Directional Derivative not implemented \n')
    otherwise
        error('SCOPT Error, Directional Derivative type not defined \n')
end
end

