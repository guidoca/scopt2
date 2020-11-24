function solverOutput = solveSOCPProblem(obj)
%% Routine to solve convex problem. Requires the definition of the SOC and Linear Matrices
if ~obj.quiet,tic;end
%% Retrieve
maxIter = obj.algorithm.conicSolver.maxIter;
feasTol = obj.algorithm.conicSolver.feasTol;
relTol  = obj.algorithm.conicSolver.relTol;
absTol  = obj.algorithm.conicSolver.absTol;
Gaug    = obj.problem.transcription.Gaug;
haug    = obj.problem.transcription.haug;
MEQ     = obj.problem.transcription.MEQ;
pEQ     = obj.problem.transcription.pEQ;
f       = obj.problem.transcription.f;
dims    = obj.problem.transcription.dims;
solver_runTime    = nan;
solver_iterations = nan;
%% Run Convex Optimization  (internal)
switch obj.algorithm.conicSolver.type
    case {'cvx-sdpt3','cvx-sedumi','cvx-ecos'}
        % Runs CVX interface with the defined solver
        if obj.quiet
            cvx_begin quiet
        else
            cvx_begin
        end
        cvx_solver(obj.algorithm.conicSolver.type(5:end))
        %     cvx_solver
        if ~obj.quiet
            fprintf('Setting up Convex Solver \n')
        end
        Ntot = length(f);
        variables Z(Ntot)
        minimize( f'*Z )
        subject to
        %SOC
        if ~obj.quiet
            fprintf('Setting up SOC Constraints \n')
        end
        tic
        qprev = 0;
        for ii = 1:length(dims.q)
            rowInitial = dims.l +obj.MIOFF + qprev;
            rowEnd   = dims.l +obj.MIOFF + sum(dims.q(1:ii)) - 1;
            qprev    = sum(dims.q(1:ii));
            norm((-Gaug((rowInitial+1):rowEnd,:))*Z+haug((rowInitial+1):rowEnd))<=(-Gaug(rowInitial,:))*Z+haug(rowInitial);
        end
        if ~obj.quiet
            fprintf('SOC Constraints setted. Time elapsed: %0.2f [s] \n',toc)
        end
        % Inequalities
        Gaug(0+obj.MIOFF:dims.l - 1+obj.MIOFF,:)*Z <= haug(0+obj.MIOFF:dims.l - 1+obj.MIOFF);
        % Equalities
        MEQ*Z == pEQ;
        if ~obj.quiet
            fprintf('Calling solver \n')
        end
        tic
        cvx_end
        solver_status   = cvx_status;
        solver_optval   = cvx_optval;
        
        if strcmp(solver_status,'Failed') || strcmp(solver_status,'Infeasible')  || strcmp(solver_status,'Unbounded')
            solver_exitFlag = 1;
        else
            solver_exitFlag = 0;
        end
        Zopt = Z;
    case {'sdpt3'}
        %                     [obj,X,y,Z,info] = sdpt3(blk,At,C,b,OPTIONS);
        error('SCOPT Error: sdpt3 not implemented')
    case {'sedumi'}
        % Runs sedumi from Matlab
        pars.maxiter = maxIter;
        pars.eps     = feasTol;
        pars.numtol  = relTol;
        dims.l = dims.l + length(pEQ)*2;
        A  = [MEQ;-MEQ;Gaug];
        ct = [pEQ;-pEQ;haug];
        %                     ct = [pEQ-feasTol/2;-pEQ-feasTol/2;haug];
        [~,Z,info] = sedumi(A',-f,ct,dims,pars);
        Zopt = Z;
        solver_optval   = sum(f.*Z);
        pinf   = info.pinf;
        dinf   = info.dinf;
        if pinf
            solver_exitFlag  = 1;
            solver_status = 'Infeasible';
        elseif dinf
            solver_exitFlag  = 2;
            solver_status = 'Infeasible';
        else
            solver_exitFlag  = 0;
            solver_status = 'Feasible';
        end
        if  info.numerr
            solver_exitFlag = -2;
            solver_status = 'Numerical Problems';
        end
        
    case {'ecos'}
        % Runs ecos from Matlab
        if obj.quiet
            verbose = 0;
        else
            verbose = 2;
        end
        opts = ecosoptimset('verbose',verbose,'maxit',maxIter,'feastol',feasTol,'reltol',relTol,'abstol',absTol);
        [Zopt,~,info] = ecos(f,Gaug,haug,dims,MEQ,pEQ,opts);
        
        solver_status     = info.infostring;
        solver_optval     = info.pcost;
        solver_exitFlag   = info.exitflag;
        solver_runTime    = info.timing.tsolve;
        solver_iterations = info.iter;
        
    case {'dual-simplex','interior-point','interior-point-legacy'}
        % Runs matlab linear programming solvers
        if isempty(dims.q)
            options = optimoptions('linprog','Algorithm',obj.algorithm.conicSolver.type,'Display','iter','MaxIterations',maxIter,'OptimalityTolerance',feasTol);
            [Zopt,fval,exitflag,output,~] = linprog(f,Gaug,haug,MEQ,pEQ,[],[],options);
            
            switch exitflag
                case 0
                    solver_exitFlag =  10;
                case 1
                    solver_exitFlag =  0;
                case 3
                    solver_exitFlag = -2;
                case -2
                    solver_exitFlag = 1;
                case -5
                    solver_exitFlag = 2;
                case {-2,-3,-4,-6,-7,-8,-9}
                    solver_exitFlag = -3;
                otherwise
                    error('Check exit flags from linprog \n')
            end
            
            solver_optval     = fval;
            solver_status     = output.message;
        else
            error('The convex subproblem contains second order cone constraints and is not linear \n')
        end
    otherwise
        error('SCOPT Error: Invalid Solver. Time elapsed: %0.2f [s] \n Aborting Run \n',toc)
end
%% Store Output
obj.solution.Zopt     = Zopt;
solverOutput.Zopt     = Zopt;
solverOutput.Foptval  = solver_optval;
solverOutput.status   = solver_status;
solverOutput.exitFlag = solver_exitFlag;
solverOutput.runTime  = solver_runTime;
solverOutput.iterations  = solver_iterations;
if ~obj.quiet,fprintf('Solver Call Finished. Time elapsed: %0.2f [s] \n Retrieving State... \n',toc);end
end

