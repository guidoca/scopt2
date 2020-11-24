function varTol = terminationVariableConvergence(obj,timeInitial_n_opt,timeFinal_n_opt,X_n_opt,U_n_opt,timeInitial_n_prev,timeFinal_n_prev,X_n_prev,U_n_prev)
% Computes the termination convergence metrics
Nodes = obj.problem.phases.Nodes;
switch obj.algorithm.sequential.variablesTol.type
    case 'norm'
        terminationAuxVar = zeros((obj.problem.phases(1).n.states+obj.problem.phases(1).n.controls+2)*Nodes,1);
        terminationAuxVar(1:(obj.problem.phases(1).n.states+obj.problem.phases(1).n.controls)*Nodes) = reshape([X_n_opt - X_n_prev;U_n_opt - U_n_prev],[],1);
        if obj.problem.phases(1).freeTimeInitial
            terminationAuxVar(end-1) = timeInitial_n_opt-timeInitial_n_prev;
        end
        if obj.problem.phases(1).freeTimeFinal
            terminationAuxVar(end)   = timeFinal_n_opt-timeFinal_n_prev;
        end
        terminationAuxNorm= norm(terminationAuxVar,obj.algorithm.sequential.variablesTol.norm);
        if ~obj.quiet
            fprintf('Convergance Metric = %0.2f%% \n',terminationAuxNorm/obj.algorithm.sequential.variablesTol.etaTol*100)
        end
        varTol     = terminationAuxNorm<=obj.algorithm.sequential.variablesTol.etaTol;
    case 'component'
        if obj.algorithm.sequential.variablesTol.include.states
            convergedStates        = sum(reshape(abs(X_n_opt - X_n_prev),[],1)<=repmat(obj.algorithm.sequential.variablesTol.phases.states_n,obj.problem.phases.Nodes,1));
            requiredNumberStates   = obj.problem.phases(1).n.states*obj.problem.phases.Nodes;
            terminationAuxStates   = convergedStates==requiredNumberStates;
            if ~obj.quiet
                fprintf('Percentage of Converged States = %0.2f%% \n',convergedStates/requiredNumberStates*100)
            end
        else
            terminationAuxStates   = 0;
        end
        if obj.algorithm.sequential.variablesTol.include.controls
            convergedControls      = sum(reshape(abs(U_n_opt - U_n_prev),[],1)<=repmat(obj.algorithm.sequential.variablesTol.phases.controls_n,obj.problem.phases.Nodes,1));
            requiredNumberControls = obj.problem.phases(1).n.controls*obj.problem.phases.Nodes;
            terminationAuxControls = convergedControls==requiredNumberControls;
            if ~obj.quiet
                fprintf('Percentage of Converged Controls = %0.2f%% \n',convergedControls/requiredNumberControls*100)
            end
        else
            terminationAuxControls = 1;
        end
        if obj.problem.phases(1).freeTimeFinal && obj.algorithm.sequential.variablesTol.include.timeFinal
            terminationAuxTimeFinal = abs(timeFinal_n_opt - timeFinal_n_prev)<=obj.algorithm.sequential.variablesTol.phases.timeFinal_n;
            if ~obj.quiet
                fprintf('Converged Final Time = %i \n',terminationAuxTimeFinal)
            end
        else
            terminationAuxTimeFinal = 1;
        end
        if obj.problem.phases(1).freeTimeInitial && obj.algorithm.sequential.variablesTol.include.timeInitial
            terminationAuxTimeInitial = abs(timeInitial_n_opt - timeInitial_n_prev)<=obj.algorithm.sequential.variablesTol.phases.timeInitial_n;
            if ~obj.quiet
                fprintf('Converged Initial Time = %i \n',terminationAuxTimeFinal)
            end
        else
            terminationAuxTimeInitial = 1;
        end
        varTol     = terminationAuxStates && terminationAuxControls && terminationAuxTimeFinal && terminationAuxTimeInitial;
    case 'none'
        varTol  = 1;
    otherwise
        error('Invalid type %s for variable convergence termination setting',obj.algorithm.sequential.variablesTol.type)
end
end

