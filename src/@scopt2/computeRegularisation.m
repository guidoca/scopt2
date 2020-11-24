function [Reg] = computeRegularisation(obj,~,X_n,U_n)
% Regularisation Term for Final Time not Enabled.
% Does not descale the values and does not multiply by dt
Reg = obj.problem.objective.regularisation.weight*sum(sum(repmat(obj.problem.objective.regularisation.states,1,size(X_n,2)).*X_n,1) + sum(repmat(obj.problem.objective.regularisation.controls,1,size(X_n,2)).*U_n,1),2);
end


