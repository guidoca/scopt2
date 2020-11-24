function [radius,lambda,rejected] = adaptTrustRegion(obj,radius,lambda,rho,numericalIssues)
% Adapts trust region based on convexification errors
if nargin<5
    numericalIssues = 0;
end
radius_l  = obj.algorithm.sequential.trustRegion.radius_l;
radius_u  = obj.algorithm.sequential.trustRegion.radius_u;
lambda_l  = obj.algorithm.sequential.trustRegion.lambda_l;
lambda_u  = obj.algorithm.sequential.trustRegion.lambda_u;
rho0      = obj.algorithm.sequential.trustRegion.rho0;
rho1      = obj.algorithm.sequential.trustRegion.rho1;
rho2      = obj.algorithm.sequential.trustRegion.rho2;
alpha     = obj.algorithm.sequential.trustRegion.alpha;
beta      = obj.algorithm.sequential.trustRegion.beta;
rejected  = 0;
if strcmp(obj.algorithm.sequential.type,'trust-region') && obj.algorithm.sequential.trustRegion.adaptive
    %             Step 3
    if (rho <= rho0) || numericalIssues
        % Reject this step, contract trust region and go back to step 1
        radius = radius / alpha;
        radius = max(radius,radius_l);
        % Reject this step, increase trust region penalisation and go back to step 1
        lambda = lambda * alpha;
        lambda = min(lambda,lambda_u);
        
        if (radius>radius_l) || (lambda<lambda_u) || numericalIssues
            rejected = 1;
        end
        
        if obj.algorithm.sequential.trustRegion.phases.timeFinal.adaptive && obj.problem.phases(1).freeTimeFinal
            % Reject this step, contract trust region and go back to step 1
            obj.algorithm.sequential.trustRegion.phases(1).timeFinal.component_n = max(obj.algorithm.sequential.trustRegion.phases(1).timeFinal.component_n / alpha,radius_l);
            % Reject this step, increase trust region penalisation and go back to step 1
            obj.algorithm.sequential.trustRegion.phases(1).timeFinal.lambda      = min(obj.algorithm.sequential.trustRegion.phases(1).timeFinal.lambda * alpha,lambda_u);
        end
    else % Accept this step
        if rho < rho1
            % Contract Trust Region
            radius = radius / alpha;
            % Penalise Trust Region More
            lambda = lambda * alpha;
            
            if obj.algorithm.sequential.trustRegion.phases.timeFinal.adaptive && obj.problem.phases(1).freeTimeFinal
                % Reject this step, contract trust region and go back to step 1
                obj.algorithm.sequential.trustRegion.phases(1).timeFinal.component_n = max(obj.algorithm.sequential.trustRegion.phases(1).timeFinal.component_n / alpha,radius_l);
                % Reject this step, increase trust region penalisation and go back to step 1
                obj.algorithm.sequential.trustRegion.phases(1).timeFinal.lambda      = min(obj.algorithm.sequential.trustRegion.phases(1).timeFinal.lambda * alpha,lambda_u);
            end
        elseif rho >= rho2
            % Expand Trust Region
            radius = radius * beta;
            % Reduce Trust Region Penalisation
            lambda = lambda / beta;
            
            if obj.algorithm.sequential.trustRegion.phases.timeFinal.adaptive && obj.problem.phases(1).freeTimeFinal
                % Reject this step, contract trust region and go back to step 1
                obj.algorithm.sequential.trustRegion.phases(1).timeFinal.component_n = min(obj.algorithm.sequential.trustRegion.phases(1).timeFinal.component_n * beta,radius_u);
                % Reject this step, increase trust region penalisation and go back to step 1
                obj.algorithm.sequential.trustRegion.phases(1).timeFinal.lambda      = max(obj.algorithm.sequential.trustRegion.phases(1).timeFinal.lambda / beta,lambda_l);
            end
        end
        radius = min(max(radius,radius_l),radius_u);
        lambda = min(max(lambda,lambda_l),lambda_u);
    end
end
end

