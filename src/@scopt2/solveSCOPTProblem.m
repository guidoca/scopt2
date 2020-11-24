function solveSCOPTProblem(obj)
%% Successive Convexification Routine.
% Solves Optimum Control Problem with Successive Convexification or with a single call to the Internal solver if sequential is inactive
% Requires level 1 and level 2 setup calls

%% Retrives some parameters
Nodes    = obj.problem.phases(1).Nodes;
index    = obj.problem.phases(1).index;

maxIter  = obj.algorithm.sequential.maxIter;
minIter  = obj.algorithm.sequential.minIter;
%% Initializes solution fields
X_n            = zeros(obj.problem.phases(1).n.states,Nodes,maxIter);
U_n            = zeros(obj.problem.phases(1).n.controls,Nodes,maxIter);
time_n         = zeros(1,Nodes,maxIter);
timeFinal_n    = zeros(1,maxIter);
timeInitial_n  = zeros(1,maxIter);

Z         = zeros(obj.problem.phases(1).Ntot,maxIter);
Foptval   = zeros(1,maxIter);
L         = zeros(1,maxIter);
Laug      = zeros(1,maxIter);
J         = zeros(1,maxIter);
Jreal     = zeros(1,maxIter);
Lreal     = zeros(1,maxIter);
Jaug      = zeros(1,maxIter);
rho       = zeros(1,maxIter);
radius    = zeros(1,maxIter);
lambda    = zeros(1,maxIter);
rejected  = zeros(1,maxIter);
exitFlag  = zeros(1,maxIter);
Reg       = zeros(1,maxIter);
runTimeSOCP    = zeros(1,maxIter);
iterationsSOCP   = zeros(1,maxIter);

%% Perform Successive Iteration Loop Until Convergence
radius(1) = obj.algorithm.sequential.trustRegion.radius;
lambda(1) = obj.algorithm.sequential.trustRegion.lambda;
kk = 1;
timeCall = tic;
while 1
    %% Initialise Loop with Convexification point
    if kk == 1 % Initial Guess
        ph = 1;
        
        X_n_prev            = obj.problem.phases(ph).initialGuess.states_n ;
        U_n_prev            = obj.problem.phases(ph).initialGuess.controls_n ;
        timeFinal_n_prev    = obj.problem.phases(ph).initialGuess.timeFinal_n ;
        timeInitial_n_prev  = obj.problem.phases(ph).initialGuess.timeInitial_n ;
        
        
        X_prev            = obj.problem.phases(ph).initialGuess.states      ;
        U_prev            = obj.problem.phases(ph).initialGuess.controls    ;
        timeFinal_prev    = obj.problem.phases(ph).initialGuess.timeFinal   ;
        timeInitial_prev  = obj.problem.phases(ph).initialGuess.timeInitial ;
        
        Jaug_prev           = 1e9;
        %                     Jaug_prev       = obj.problem.phases(ph).initialGuess.J;
        
        X_n(:,:,1)      = X_n_prev ; % Is this necessary?
        U_n(:,:,1)      = U_n_prev ;
        timeFinal_n(1)  = timeFinal_n_prev;
        timeInitial_n(1)= timeInitial_n_prev;
        
    else % Previous Trajectory (if successive pass) (internal)
        ph = 1;
        
        if rejected(kk-1) % If it was rejected in the previous iteration, retrieve the solution from the last valid iteration
            kkValid = find(rejected(1:kk-1)==0,1,'last');
            
            X_n_prev        = X_n(:,:,kkValid);
            U_n_prev        = U_n(:,:,kkValid);
            if obj.problem.phases(ph).freeTimeFinal
                % Free time of flight
                timeFinal_n_prev   = timeFinal_n(kkValid)   ;
            end
            if obj.problem.phases(ph).freeTimeInitial
                % Free time of flight
                timeInitial_n_prev = timeInitial_n(kkValid) ;
            end
            obj.problem.phases(ph).BodyMap.evaluated = 0;
            Jaug_prev       = Jaug(kkValid);
        else
            X_n_prev        = X_n(:,:,kk-1);
            U_n_prev        = U_n(:,:,kk-1);
            if obj.problem.phases(ph).freeTimeFinal
                % Free time of flight
                timeFinal_n_prev = timeFinal_n(kk-1) ;
            end
            if obj.problem.phases(ph).freeTimeInitial
                % Free time of flight
                timeInitial_n_prev = timeInitial_n(kk-1) ;
            end
            Jaug_prev       = Jaug(kk-1);
        end
        
        X_prev            = obj.scaleStatesToReal(1,X_n_prev);
        U_prev            = obj.scaleControlsToReal(1,U_n_prev);
        timeFinal_prev    = obj.scaleTimeToReal(1,timeFinal_n_prev);
        timeInitial_prev  = obj.scaleTimeToReal(1,timeInitial_n_prev);
    end
    
    if ~obj.quiet
        fprintf('Iteration: %0.0f \n',kk)
        fprintf('Sequential Type: %s \n',obj.algorithm.sequential.type)
        fprintf('Trust Region Radius = %0.16f [-] \n',radius(kk) )
        fprintf('Trust Region Weight = %0.16f [-] \n',lambda(kk) )
    end
    
    timeDisc_prev     = (timeFinal_prev - timeInitial_prev)/(Nodes-1);
    time_prev         = (timeInitial_prev:timeDisc_prev:timeFinal_prev);
    timeDisc_n_prev   = (timeFinal_n_prev - timeInitial_n_prev)/(Nodes-1);
    time_n_prev       = (timeInitial_n_prev:timeDisc_n_prev:timeFinal_n_prev);
    %% Transcription of SOCP Problem
    obj.problemTranscription(timeInitial_n_prev,timeFinal_n_prev,X_n_prev,U_n_prev,radius(kk),lambda(kk));
    
    %% Solve problem
    solverOutput = obj.solveSOCPProblem();
    
    Z_opt        = solverOutput.Zopt    ;
    Foptval(kk)  = solverOutput.Foptval ;
    status       = solverOutput.status  ;
    exitFlag(kk) = solverOutput.exitFlag;
    runTimeSOCP(kk)  = solverOutput.runTime;
    iterationsSOCP(kk)  = solverOutput.iterations;
    
    % Check if feasible
    if kk >1
        if ~(exitFlag(kk)==-2 || exitFlag(kk)==-1 || exitFlag(kk)==0 || exitFlag(kk)==10)
            failed          = true;
            numericalIssues = true;
            warning('Failed with exit flag %0.0f in iteration %0.0f \n Aborting Optimization and retrieving last solution',exitFlag(kk),kk )
        elseif (exitFlag(kk)==-2)
            numericalIssues = true ;
            if obj.algorithm.sequential.trustRegion.adaptive
                failed          = false;
                warning('SOCP subproblem having numerical problems in iteration %0.0f, \n Rejecting step if adaptive \n',kk )
            else
                failed          = true;
                warning('SOCP subproblem having numerical problems in iteration %0.0f \n Rejecting step and terminating with last accepeted step \n',kk )
            end
            
        elseif (exitFlag(kk)==-1)
            failed          = false;
            numericalIssues = true ;
            warning('SOCP subproblem reached maximum number of iterations in iteration %0.0f \n',kk )
        elseif (exitFlag(kk)==10)
            failed          = false;
            numericalIssues = true ;
            warning('SOCP subproblem returning an inaccurate result at iteration %0.0f \n',kk )
        else
            failed          = false ;
            numericalIssues = false ;
        end
    elseif (exitFlag(1)==-1)
        failed          = false;
        numericalIssues = true;
        warning('Maximum iterations in the first iteration with exit flag %0.0f \n ',exitFlag(1) )
    elseif ~(exitFlag(1)==0)
        failed          = true;
        numericalIssues = true;
        warning('Failed in the first iteration with exit flag %0.0f \n Aborting Optimization and retrieving initial guess',exitFlag(1) )
    else
        failed          = false ;
        numericalIssues = false ;
    end
    %% Retrieve States, Controls, and Time
    if ~failed
        X_n_opt = Z_opt(index.states);
        U_n_opt = reshape(Z_opt(index.controls),size(index.controls));
        
        if obj.problem.phases(1).freeTimeFinal
            timeFinal_n_opt = Z_opt(index.timeFinal);
        else
            timeFinal_n_opt = timeFinal_n_prev;
        end
        if obj.problem.phases(1).freeTimeInitial
            timeInitial_n_opt = Z_opt(index.timeInitial);
        else
            timeInitial_n_opt = timeInitial_n_prev;
        end
    else % Retrieve the previous optimum value
        X_n_opt  = X_n_prev;
        U_n_opt  = U_n_prev;
        timeFinal_n_opt   = timeFinal_n_prev;
        timeInitial_n_opt = timeInitial_n_prev;
    end
    
    X_opt            = obj.scaleStatesToReal(1,X_n_opt);
    U_opt            = obj.scaleControlsToReal(1,U_n_opt);
    timeFinal_opt    = obj.scaleTimeToReal(1,timeFinal_n_opt);
    timeInitial_opt  = obj.scaleTimeToReal(1,timeInitial_n_opt);
    
    timeDisc_n_opt = (timeFinal_n_opt - timeInitial_n_opt)/(Nodes-1);
    time_n_opt     = (timeInitial_n_opt:timeDisc_n_opt:timeFinal_n_opt);
    
    timeDisc_opt   = (timeFinal_opt - timeInitial_opt)/(Nodes-1);
    time_opt       = (timeInitial_opt:timeDisc_opt:timeFinal_opt);
    
    %% Interpolate Values for Next Iteration: Onluy if speedup.
    obj.problem.phases(1).BodyMap_prev = obj.problem.phases(1).BodyMap;
    if obj.problem.phases(1).interpolateParameters && obj.algorithm.speedUp
        [obj.problem.phases(1).BodyMap] = obj.problem.phases(1).interpolationsFunction(time_opt,X_opt,U_opt,obj.problem.phases(1).BodyMap);
        obj.problem.phases(1).BodyMap.evaluated = 1;
    else
        obj.problem.phases(1).BodyMap.evaluated = 0;
    end
    %% Compute Non-linear and Linear Objective Values
    [Jreal(kk),Lreal(kk)] = computeObjectiveFunctions(obj,time_opt,X_opt,U_opt,time_prev,X_prev,U_prev);
    J(kk)     = Jreal(kk)/obj.problem.objective.scale*obj.problem.objective.sign ;
    L(kk)     = Lreal(kk)/obj.problem.objective.scale*obj.problem.objective.sign ;
    
    Reg(kk)   = computeRegularisation(obj,time_n_opt,X_n_opt,U_n_opt);
    %% Adaptive Trust Region Algorithm
    if ~failed
        [penaltyDefectsNonLinear,penaltyBuffersNonLinear,penaltyTrustRegion,penaltyTime,penaltyDefectsLinear,penaltyBuffersLinear] = computeAugmentedObjectiveFunctionTerms(obj,time_n_opt,X_n_opt,U_n_opt,time_n_prev,X_n_prev,U_n_prev);
        %% Evaluation of Objective Functions
        % Step 1: Augment Objective Function
        Laug(kk)         = L(kk) + Reg(kk) + penaltyDefectsLinear    + penaltyBuffersLinear    + penaltyTrustRegion + penaltyTime;
        Jaug(kk)         = J(kk) + Reg(kk) + penaltyDefectsNonLinear + penaltyBuffersNonLinear + penaltyTrustRegion + penaltyTime;
        if kk > 1
            %% Step 2: Check Convexification Error
            ActualChange    = Jaug_prev - Jaug(kk)        ;
            PredictedChange = Jaug_prev - Laug(kk)        ;
            rho(kk)         = ActualChange/PredictedChange;
            RelativeChange  = ActualChange/Jaug_prev      ;
            %% Step 3: Adapt Trust Region based on Convexification Error
            [radius(kk+1),lambda(kk+1),rejected(kk)] = adaptTrustRegion(obj,radius(kk),lambda(kk),rho(kk),numericalIssues);
            %%  Globalisation Strategies (if enabled)
            if (~rejected(kk) || obj.algorithm.sequential.globalisation.lineSearch.forceIfRejected) && strcmp(obj.algorithm.sequential.globalisation.type,'line-search') && obj.algorithm.sequential.globalisation.activate && ~numericalIssues
                skipLineSearchStep = 0;
                % Compute Search Direction
                Z_opt_red       = [X_n_opt(:);U_n_opt(:);timeInitial_n_opt(:);timeFinal_n_opt(:)];
                Z_prev_red      = [X_n_prev(:);U_n_prev(:);timeInitial_n_prev(:);timeFinal_n_prev(:)];
                searchDirection = Z_opt_red - Z_prev_red;
                stepLength      = obj.algorithm.sequential.globalisation.lineSearch.stepLengthMax;
                flagDebugInternal = obj.algorithm.sequential.globalisation.lineSearch.debugging;
                switch obj.algorithm.sequential.globalisation.lineSearch.type
                    case 'contraction'
                        directionalDerivative = computeDirectionalDerivative(obj,searchDirection,Z_prev_red,Jaug(kk-1));
                        contractionFactor = obj.algorithm.sequential.globalisation.lineSearch.contractionFactor;
                        while 1
                            sufficientDecrease = evaluateSufficientDecrease(obj,searchDirection,stepLength,directionalDerivative,Z_prev_red,Jaug(kk-1));
                            if sufficientDecrease || directionalDerivative>=0
                                break
                            else
                                stepLength = contractionFactor*stepLength;
                            end
                        end
                    case 'polynomial-search'
                        bracket     = [0,stepLength];
                        Mbracket    = [Jaug_prev,J(kk)];
                        bracketTol  = obj.algorithm.sequential.globalisation.lineSearch.bracketTol;
                        while 1
                            x1      = bracket(1);
                            x3      = bracket(2);
                            M1      = Mbracket(1) ;
                            M3      = Mbracket(2) ;
                            
                            h       = (x3-x1)/2;
                            x2      = x1 + h;
                            
                            M2      = computeMeritFunction(obj,searchDirection,x2,Z_prev_red) ;
                            
                            if (M1 >= M2) && (M2 <= M3)
                                mu = M2 - M1;
                                eta= M3 - M1;
                                xNew = -(2*mu*h + (2*mu-eta)*(2*x1+h))/(2*(eta-2*mu));
                                if xNew>x2
                                    bracket(1) = x2;
                                    bracket(2) = x3;
                                    Mbracket(1) = M2;
                                    Mbracket(2) = M3;
                                else
                                    bracket(1) = x1;
                                    bracket(2) = x2;
                                    Mbracket(1) = M1;
                                    Mbracket(2) = M2;
                                end
                                if (abs(bracket(1)- bracket(2)) <= bracketTol)
                                    break
                                end
                            elseif (M1 > M2) && (M2 > M3)
                                stepLength = x3;
                                break
                            elseif (M1 < M2) && (M2 < M3)
                                stepLength = x1;
                                break
                            elseif (M1 < M2) && (M2 > M3)
                                stepLength = x3;
                                break
                            end
                        end
                    case 'golden-section'
                        [stepLength] = goldenSection (@(x) computeMeritFunction(obj,searchDirection,x,Z_prev_red),0,stepLength,obj.algorithm.sequential.globalisation.lineSearch.bracketTol,[],[Jaug_prev,Jaug(kk)]);
                    otherwise
                        error('SCOPT Error, line search algorithm not defined')
                end
                if stepLength == 0
                    warning('SCOPT Warning: Line search algorithm could not reduce objective along SOCP seach direction. If it was a rejected step, then do not update parameters \n')
                    if rejected(kk)
                        skipLineSearchStep = 1;
                    end
                end
                
                if ~skipLineSearchStep
                    Z_red      = computeNewPoint(obj,searchDirection,stepLength,Z_prev_red);
                    %                             Z_opt      = Z_red(index.states);
                    X_n_opt    = Z_red(index.states);
                    U_n_opt    = reshape(Z_red(index.controls),size(index.controls));
                    X_opt      = obj.scaleStatesToReal(1,X_n_opt);
                    U_opt      = obj.scaleControlsToReal(1,U_n_opt);
                    
                    timeFinal_n_opt    = Z_red(end);
                    
                    timeDisc_n_opt = (timeFinal_n_opt - timeInitial_n_opt)/(Nodes-1);
                    time_n_opt     = (timeInitial_n_opt:timeDisc_n_opt:timeFinal_n_opt);
                    timeFinal_opt    = obj.scaleTimeToReal(1,timeFinal_n_opt);
                    timeInitial_opt  = obj.scaleTimeToReal(1,timeInitial_n_opt);
                    timeDisc_opt   = (timeFinal_opt - timeInitial_opt)/(Nodes-1);
                    time_opt       = (timeInitial_opt:timeDisc_opt:timeFinal_opt);
                    
                    if obj.problem.phases(1).interpolateParameters && obj.algorithm.speedUp
                        [obj.problem.phases(1).BodyMap] = obj.problem.phases(1).interpolationsFunction(time_opt,X_opt,U_opt,obj.problem.phases(1).BodyMap);
                        obj.problem.phases(1).BodyMap.evaluated = 1;
                    else
                        obj.problem.phases(1).BodyMap.evaluated = 0;
                    end
                    
                    [Jreal(kk),Lreal(kk)] = computeObjectiveFunctions(obj,time_opt,X_opt,U_opt,time_prev,X_prev,U_prev);
                    %                             [Jreal(kk),Lreal(kk)] = obj.computeObjectiveFunctions(time_opt,X_opt,U_opt,time_prev,X_prev,U_prev);
                    
                    J(kk)     = Jreal(kk)/obj.problem.objective.scale*obj.problem.objective.sign;
                    L(kk)     = Lreal(kk)/obj.problem.objective.scale*obj.problem.objective.sign;
                    
                    Reg(kk)       = computeRegularisation(obj,time_n_opt,X_n_opt,U_n_opt);
                    [penaltyDefectsNonLinear,penaltyBuffersNonLinear,penaltyTrustRegion,penaltyTime,penaltyDefectsLinear,penaltyBuffersLinear] = computeAugmentedObjectiveFunctionTerms(obj,time_n_opt,X_n_opt,U_n_opt,time_n_prev,X_n_prev,U_n_prev);
                    Laug(kk)         = L(kk) + Reg(kk) + penaltyDefectsLinear    + penaltyBuffersLinear    + penaltyTrustRegion + penaltyTime;
                    Jaug(kk)         = J(kk) + Reg(kk) + penaltyDefectsNonLinear + penaltyBuffersNonLinear + penaltyTrustRegion + penaltyTime;
                    
                    
                    ActualChange    = Jaug_prev - Jaug(kk)        ;
                    PredictedChange = Jaug_prev - Laug(kk)        ;
                    rho(kk)         = ActualChange/PredictedChange;
                    RelativeChange  = ActualChange/Jaug_prev      ;
                    
                    if rejected(kk) && obj.algorithm.sequential.globalisation.lineSearch.forceIfRejected
                        % Reevaluate rejection criteria
                        [~,~,rejected(kk)] = adaptTrustRegion(obj,radius(kk),lambda(kk),rho(kk));
                    end
                end
            end
        else
            ActualChange    = Jaug_prev - J(kk)                ;
            PredictedChange = Jaug_prev - L(kk)                ;
            rho(kk)         = ActualChange/PredictedChange;
            RelativeChange  = ActualChange/Jaug_prev;
            radius(kk+1)    = radius(kk);
            lambda(kk+1)    = lambda(kk);
        end
        if ~obj.quiet
            fprintf('Change in Objective Function  = %0.2f [%%] \n',RelativeChange*100 )
        end
    else % If it is not the first iteration, store solutions
        Laug(kk)            = L(kk) + Reg(kk);
        Jaug(kk)            = J(kk) + Reg(kk);
        ActualChange        = Jaug_prev - Jaug(kk)  ;
        RelativeChange      = ActualChange/Jaug_prev;
        radius(kk+1)        = radius(kk);
        lambda(kk+1)        = lambda(kk);
    end
    
    %% Store Iteration History
    X_n(:,:,kk)      = X_n_opt           ;
    U_n(:,:,kk)      = U_n_opt           ;
    time_n(:,:,kk)   = time_n_opt        ;
    Z(:,kk)          = Z_opt             ;
    timeFinal_n(kk)  = timeFinal_n_opt   ;
    timeInitial_n(kk)= timeInitial_n_opt ;
    
    %% Monitor
    if ~obj.quiet
        fprintf('Real Optimal is Jopt  = %0.2f [] \n Convex Prediction is Lopt = %0.2f [] \n',Jreal(kk),Lreal(kk))
        fprintf('Debug:  Linear    Prediction diff is  = %0.8e [] \n',L(kk)-obj.problem.transcription.fnp'*Z_opt)
        fprintf('Debug:  Augmented Prediction diff is  = %0.8e [] \n',Laug(kk)-obj.problem.transcription.f'*Z_opt)
        fprintf('Final Time is t_f = %0.2f [-] \n',timeFinal_opt )
        fprintf('Initial Time is t_f = %0.2f [-] \n',timeInitial_opt )
    end
    
    %% Termination Settings Check and Loop ending
    % Compute Variable Convergence
    Termination.varTol = terminationVariableConvergence(obj,timeInitial_n_opt,timeFinal_n_opt,X_n_opt,U_n_opt,timeInitial_n_prev,timeFinal_n_prev,X_n_prev,U_n_prev);
    
    Termination.Jtol       = (abs(ActualChange)<=obj.algorithm.sequential.Jtol);
    Termination.JtolChange = (abs(RelativeChange)<=obj.algorithm.sequential.JtolChange);
    Termination.minIter    = (kk>=minIter);
    Termination.maxIter    = (kk>=maxIter);
    Termination.failed     = failed;
    
    Termination.converged     = (Termination.Jtol || Termination.JtolChange || Termination.varTol); % ~(rejected(kk) && ~numericalIssues) &&
    Termination.endSequential = (Termination.minIter && (Termination.converged || Termination.maxIter))  || Termination.failed || ~obj.algorithm.sequential.activate;
    % IF termination flag triggered, end convexification
    if Termination.endSequential
        obj.solution.sequential.termination = Termination;
        obj.solution.sequential.status      = status;
        obj.solution.sequential.exitFlag    = exitFlag(kk);
        
        totalIterations =kk;
        timeElapsed = toc(timeCall);
        X_n(:,:,totalIterations+1:maxIter)     = []   ;
        U_n(:,:,totalIterations+1:maxIter)     = []   ;
        time_n(:,:,totalIterations+1:maxIter)  = []   ;
        Z(:,totalIterations+1:maxIter)         = []   ;
        
        timeFinal_n(totalIterations+1:maxIter)    = []   ;
        timeInitial_n(totalIterations+1:maxIter)  = []   ;
        
        Foptval(totalIterations+1:maxIter) = []   ;
        Jreal(totalIterations+1:maxIter)   = []   ;
        Lreal(totalIterations+1:maxIter)   = []   ;
        L(totalIterations+1:maxIter)       = []   ;
        Laug(totalIterations+1:maxIter)    = []   ;
        J(totalIterations+1:maxIter)       = []   ;
        Jaug(totalIterations+1:maxIter)    = []   ;
        radius(totalIterations+1:end)      = []   ;
        lambda(totalIterations+1:end)      = []   ;
        rho(totalIterations+1:maxIter)     = []   ;
        rejected(totalIterations+1:maxIter)= []   ;
        exitFlag(totalIterations+1:maxIter)= []   ;
        runTimeSOCP(totalIterations+1:maxIter) = []   ;
        iterationsSOCP(totalIterations+1:maxIter) = []   ;
        if ~obj.quiet
            fprintf('Successive Convexification Finished at iteration %0.0f \n Elapsed Time: %0.2f \n',kk,timeElapsed);
            fprintf('Real Objectives are J= %0.2f [-] and L= %0.2f [-] \n',Jreal(totalIterations),Lreal(totalIterations));
        end
        break
    end
    kk = kk + 1;
end
%% Scale to dimensional values
X = zeros(size(X_n));
U = zeros(size(U_n));
for kk      = 1:totalIterations
    X(:,:,kk)           = obj.scaleStatesToReal(1,X_n(:,:,kk));
    U(:,:,kk)           = obj.scaleControlsToReal(1,U_n(:,:,kk));
end
time          = obj.scaleTimeToReal(1,time_n);
timeFinal     = obj.scaleTimeToReal(1,timeFinal_n);
timeInitial   = obj.scaleTimeToReal(1,timeInitial_n);
%% Retrieve Last Adequate Result
if rejected(totalIterations) % If failed also, return the last adequate estimate
    kkValid = find(rejected(1:totalIterations-1)==0,1,'last');
    kk_sol  = kkValid;
else
    kk_sol = totalIterations;
end
J_sol                  = J(kk_sol);
L_sol                  = L(kk_sol);
Jaug_sol               = Jaug(kk_sol);
Laug_sol               = Laug(kk_sol);
time_sol               = time(1,:,kk_sol);
timeFinal_sol          = timeFinal(kk_sol);
timeInitial_sol        = timeInitial(kk_sol);
X_sol                  = X(:,:,kk_sol);
U_sol                  = U(:,:,kk_sol);
Jreal_sol              = Jreal(kk_sol);

time_n_sol             = time_n(1,:,kk_sol);
timeFinal_n_sol        = timeFinal_n(kk_sol);
timeInitial_n_sol      = timeInitial_n(kk_sol);
X_n_sol                = X_n(:,:,kk_sol);
U_n_sol                = U_n(:,:,kk_sol);

Reg_sol                = Reg(kk_sol);
exitFlag_sol           = exitFlag(kk_sol);
%% Assess Solution Quality
%             obj.solution.objTol   = abs((Laug(Iter) - L(Iter))/L(Iter));
obj.solution.objTol   = abs(Laug_sol- Reg_sol - L_sol);
obj.solution.feasible = (obj.solution.objTol <= obj.algorithm.sequential.objTol) && (exitFlag_sol~=1) && (exitFlag_sol~=2) && (exitFlag_sol~=-3) ;
obj.solution.failed   = failed;
if ~failed
    if obj.solution.feasible
        if ~obj.quiet
            fprintf('Successive convexification ended successfully with a feasible augmented tolerance of = %0.2e [-] \n',obj.solution.objTol);
        end
        if rejected(totalIterations)
            if numericalIssues
                if ~obj.quiet
                    fprintf('Last iteration runned into numerical issues. Returning %i iteration [-] \n',kk_sol);
                end
            else
                fprintf('Last iteration performed an inadequate objective increase. Returning %i iteration [-] \n',kk_sol);
            end
        end
    else
        if ~obj.quiet
            warning('Successive convexification ended successfully but the feasible augmented tolernace is bigger than the tolerance = %0.2e > %0.2e \n Re-run scopt without virtual controls or buffers to ensure the trajectory is feasible',obj.solution.objTol,obj.algorithm.sequential.objTol);
        end
    end
else
    if ~obj.quiet
        fprintf('Solution is not feasible with a tolerance of = %0.2e [-] and certificate of infeasibility %0.0f \n Returning Previous feasible estimate',exitFlag_sol);
    end
end

obj.solution.numberOfWarnings = sum(exitFlag~=0);
if obj.solution.numberOfWarnings>0
    warning('There were %0.0f warnings in the SOCP subproblem \n',obj.solution.numberOfWarnings);
end
%% Store Solution Structure
obj.solution.Nodes              = Nodes;
obj.solution.iterations         = totalIterations;
obj.solution.timeElapsed        = timeElapsed;
obj.solution.J                  = J_sol;
obj.solution.L                  = L_sol;
obj.solution.Jaug               = Jaug_sol;
obj.solution.Laug               = Laug_sol;
obj.solution.time               = time_sol;
obj.solution.timeFinal          = timeFinal_sol;
obj.solution.timeInitial        = timeInitial_sol;
obj.solution.states             = X_sol;
obj.solution.controls           = U_sol;
obj.solution.Jreal              = Jreal_sol;

obj.solution.time_n             = time_n_sol;
obj.solution.timeFinal_n        = timeFinal_n_sol;
obj.solution.timeInitial_n      = timeInitial_n_sol;
obj.solution.states_n           = X_n_sol;
obj.solution.controls_n         = U_n_sol;
obj.solution.runTimeAverage     = mean(runTimeSOCP);
obj.solution.iterationsAverage  = mean(iterationsSOCP);

obj.debugging.Reg            = Reg;
obj.debugging.Foptval        = Foptval;
obj.debugging.time_n_ALL     = time_n;
obj.debugging.timeFinal_n_ALL= timeFinal_n;
obj.debugging.timeInitial_n_ALL= timeInitial_n;
obj.debugging.states_n_ALL   = X_n;
obj.debugging.controls_n_ALL = U_n;
obj.debugging.time_ALL       = time;
obj.debugging.timeFinal_ALL  = timeFinal;
obj.debugging.timeInitial_ALL  = timeInitial;
obj.debugging.states_ALL     = X;
obj.debugging.controls_ALL   = U;
obj.debugging.Z              = Z;
obj.debugging.Jreal          = Jreal;
obj.debugging.Lreal          = Lreal;
obj.debugging.J              = J;
obj.debugging.Jaug           = Jaug;
obj.debugging.L              = L;
obj.debugging.Laug           = Laug;
obj.debugging.rho            = rho;
obj.debugging.radius         = radius;
obj.debugging.lambda         = lambda;
obj.debugging.rejected       = rejected;
obj.debugging.exitFlag       = exitFlag;
if obj.debugging.saveVariableConvergence
    obj.debugging.variablesNorm = computeConvergenceHistory(obj);
    
end

if obj.debugging.saveVirtualControls && obj.algorithm.virtualControl.phases(1).include
    obj.debugging.virtualControls   = zeros(obj.problem.phases(1).n.virtualControls,Nodes-1,totalIterations);
    obj.debugging.virtualControlsP1 = zeros(obj.problem.phases(1).n.virtualControlsP1,Nodes-1,totalIterations);
    obj.debugging.virtualControlsP1s= zeros(obj.problem.phases(1).n.virtualControlsP1s,Nodes-1,totalIterations);
    obj.debugging.virtualControlsP2 = zeros(obj.problem.phases(1).n.virtualControlsP2,totalIterations);
    
    for kk = 1:totalIterations
        Zaux = Z(:,kk);
        obj.debugging.virtualControls(:,:,kk)    = Zaux(index.virtualControls);
        obj.debugging.virtualControlsP1(:,:,kk)  = Zaux(index.virtualControlsP1);
        obj.debugging.virtualControlsP1s(:,:,kk) = Zaux(index.virtualControlsP1s);
        obj.debugging.virtualControlsP2(:,kk)    = Zaux(index.virtualControlsP2);
    end
end
if obj.debugging.saveVirtualBuffers && obj.algorithm.virtualBuffer.phases(1).include
    obj.debugging.virtualBuffers.path    = nan(obj.problem.phases(1).n.virtualBuffers.path,Nodes,totalIterations);
    obj.debugging.virtualBuffers.pathPs  = nan(obj.problem.phases(1).n.virtualBuffers.pathPs,Nodes,totalIterations);
    obj.debugging.virtualBuffers.pathP   = nan(obj.problem.phases(1).n.virtualBuffers.pathP,totalIterations);
    obj.debugging.virtualBuffers.events  = nan(obj.problem.phases(1).n.virtualBuffers.events,totalIterations);
    obj.debugging.virtualBuffers.eventsP = nan(obj.problem.phases(1).n.virtualBuffers.eventsP,totalIterations);
    
    for kk = 1:totalIterations
        Zaux = Z(:,kk);
        obj.debugging.virtualBuffers.path(:,:,kk)   = Zaux(index.virtualBuffers.path);
        obj.debugging.virtualBuffers.pathPs(:,:,kk) = Zaux(index.virtualBuffers.pathPs);
        obj.debugging.virtualBuffers.pathP(:,kk)    = Zaux(index.virtualBuffers.pathP);
        if obj.problem.phases(1).n.virtualBuffers.events>0
            obj.debugging.virtualBuffers.events(:,kk)   = Zaux(index.virtualBuffers.events);
            obj.debugging.virtualBuffers.eventsP(:,kk)  = Zaux(index.virtualBuffers.eventsP);
        end
    end
end

if obj.debugging.saveTrustRegion && obj.algorithm.sequential.activate && ~strcmp(obj.algorithm.sequential.type,'none')
    obj.debugging.trustRegionP1     = zeros(1,Nodes,totalIterations);
    obj.debugging.trustRegionP2     = zeros(1,totalIterations);
    for kk = 1:totalIterations
        Zaux = Z(:,kk);
        obj.debugging.trustRegionP1(:,:,kk)  = Zaux(index.trustRegionP1);
        obj.debugging.trustRegionP2(:,kk)  = Zaux(index.trustRegionP2);
    end
end

%% PostProcessing if Enabled
if obj.solution.generateFigures
    obj.generateFigures();
end
end


