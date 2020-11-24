function variablesNorm = computeConvergenceHistory(obj)
% Computes Convergence History as defined by the norm type
% variable tolerance
states_n_IG    = obj.problem.phases.initialGuess.states_n;
controls_n_IG  = obj.problem.phases.initialGuess.controls_n;
states_n_ALL   = cat(3,states_n_IG,obj.debugging.states_n_ALL);
controls_n_ALL = cat(3,controls_n_IG,obj.debugging.controls_n_ALL);
variables_ALL  = cat(1,states_n_ALL,controls_n_ALL);
trustRegionAux = diff(variables_ALL,1,3);
if obj.algorithm.sequential.variablesTol.include.timeFinal == 1
    timeFinal_n_IG  = obj.problem.phases.initialGuess.timeFinal_n;
    timeFinal_n_ALL = [timeFinal_n_IG,obj.debugging.timeFinal_n_ALL];
    trustRegionTimeFinalAux = diff(timeFinal_n_ALL,1,2);
end
for kk = 2:obj.solution.iterations+0
    if obj.debugging.rejected(kk-1)
        kkValid = find(obj.debugging.rejected(1:kk-1)==0,1,'last');
        trustRegionAux(:,:,kk) = variables_ALL(:,:,kk+1)-variables_ALL(:,:,kkValid+1);
        if obj.algorithm.sequential.variablesTol.include.timeFinal == 1
            trustRegionTimeFinalAux(1,kk) = timeFinal_n_ALL(1,kk+1) - timeFinal_n_ALL(1,kkValid+1);
        end
    end
end
trustRegionAux = reshape(trustRegionAux,(size(states_n_ALL,1)+size(controls_n_ALL,1))*size(states_n_ALL,2),obj.solution.iterations);
if obj.algorithm.sequential.variablesTol.include.timeFinal == 1
    trustRegionAux = [trustRegionAux;trustRegionTimeFinalAux];
end
variablesNorm  = vecnorm(trustRegionAux,obj.algorithm.sequential.variablesTol.norm,1);
end


