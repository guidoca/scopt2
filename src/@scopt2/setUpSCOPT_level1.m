function obj = setUpSCOPT_level1(obj)
%% Level 1 SCOPT Setup
% Setup Fields with Default Values based on the number of
% states, controls, parameters, free and final initial time
% fields, and phases
for ph = 0:obj.problem.nphases-1
    % Retrieve Phase Sepcific fields
    numberOfStates   = obj.problem.phases(ph+obj.MIOFF).n.states ;
    numberOfControls = obj.problem.phases(ph+obj.MIOFF).n.controls ;
    numberOfParameters =  obj.problem.phases(ph+obj.MIOFF).n.parameters ;
    numberOfNodes    = obj.problem.phases(ph+obj.MIOFF).Nodes ;
    
    % Only active in the phase if speedUp flag is enabled. Runs
    % some routines before. Not necessary anymroe with
    % griddedInterpolant fast use
    obj.problem.phases(ph+obj.MIOFF).interpolateParameters = 1;
    
    % Tolerances for each variable if specified as desired
    % tolernace criterion
    obj.algorithm.sequential.variablesTol.phases(ph+obj.MIOFF).states      = eps^(1/4)*ones(numberOfStates,1); %[-]
    obj.algorithm.sequential.variablesTol.phases(ph+obj.MIOFF).controls    = eps^(1/4)*ones(numberOfControls,1); %[-]
    obj.algorithm.sequential.variablesTol.phases(ph+obj.MIOFF).timeFinal   = eps^(1/4)*ones(1,1); %[-]
    obj.algorithm.sequential.variablesTol.phases(ph+obj.MIOFF).timeInitial = eps^(1/4)*ones(1,1); %[-]
    
    % Dynamics functions. Default uses sparse matrix
    % representation for state jacobians and control
    % jacobians. Dependency on ii (node number) has to be
    % specified by user
    obj.problem.phases(ph+obj.MIOFF).dynamics.stateDerivativeFunction        = @(t,X,U,BodyMap,ii) nan(numberOfStates,length(t));
    obj.problem.phases(ph+obj.MIOFF).dynamics.typeStateMatrix                = 1; % 0 for 3D, 1 for sparse representation
    obj.problem.phases(ph+obj.MIOFF).dynamics.stateMatrixFunction            = @(t,X,U,BodyMap,ii) nan(numberOfStates,numberOfStates,length(t));
    obj.problem.phases(ph+obj.MIOFF).dynamics.typeControlMatrix              = 1; % 0 for 3D, 1 for sparse representation
    obj.problem.phases(ph+obj.MIOFF).dynamics.controlMatrixFunction          = @(t,X,U,BodyMap,ii) nan(numberOfStates,numberOfControls,length(t));
    obj.problem.phases(ph+obj.MIOFF).dynamics.timeDerivativeFunction         = @(t,X,U,BodyMap,ii) zeros(numberOfStates,length(t));
    obj.problem.phases(ph+obj.MIOFF).dynamics.typeParameterMatrix            = 0; % 0 for 3D, 1 for sparse representation
    obj.problem.phases(ph+obj.MIOFF).dynamics.parametersMatrixFunction       = @(t,X,U,BodyMap,ii) nan(numberOfStates,numberOfParameters,length(t));
    
    % Initialize Scaling Terms
    obj.problem.phases(ph+obj.MIOFF).scale.states   = ones(numberOfStates,1);
    obj.problem.phases(ph+obj.MIOFF).scale.controls = ones(numberOfControls,1);
    obj.problem.phases(ph+obj.MIOFF).scale.time     = 1;
    obj.problem.phases(ph+obj.MIOFF).shift.states   = zeros(numberOfStates,1);
    obj.problem.phases(ph+obj.MIOFF).shift.controls = zeros(numberOfControls,1);
    obj.problem.phases(ph+obj.MIOFF).shift.time     = 0;
    
    % Initialzes Initial Conditions
    obj.problem.phases(ph+obj.MIOFF).initial.states.equal.index      = [];
    obj.problem.phases(ph+obj.MIOFF).initial.states.equal.value      = zeros(obj.problem.phases(ph+obj.MIOFF).initial.states.equal.index,1);
    obj.problem.phases(ph+obj.MIOFF).initial.controls.equal.index    = [];
    obj.problem.phases(ph+obj.MIOFF).initial.controls.equal.value    = zeros(obj.problem.phases(ph+obj.MIOFF).initial.controls.equal.index,1);
    obj.problem.phases(ph+obj.MIOFF).initial.states.lower.index      = [];
    obj.problem.phases(ph+obj.MIOFF).initial.states.lower.value      = zeros(obj.problem.phases(ph+obj.MIOFF).initial.states.lower.index,1);
    obj.problem.phases(ph+obj.MIOFF).initial.controls.lower.index    = [];
    obj.problem.phases(ph+obj.MIOFF).initial.controls.lower.value    = zeros(obj.problem.phases(ph+obj.MIOFF).initial.controls.lower.index,1);
    obj.problem.phases(ph+obj.MIOFF).initial.states.upper.index      = [];
    obj.problem.phases(ph+obj.MIOFF).initial.states.upper.value      = zeros(obj.problem.phases(ph+obj.MIOFF).initial.states.upper.index,1);
    obj.problem.phases(ph+obj.MIOFF).initial.controls.upper.index    = [];
    obj.problem.phases(ph+obj.MIOFF).initial.controls.upper.value    = zeros(obj.problem.phases(ph+obj.MIOFF).initial.controls.upper.index,1);
    
    % Initializes Final Conditions
    obj.problem.phases(ph+obj.MIOFF).final.states.equal.index      = [];
    obj.problem.phases(ph+obj.MIOFF).final.states.equal.value      = zeros(obj.problem.phases(ph+obj.MIOFF).final.states.equal.index,1);
    obj.problem.phases(ph+obj.MIOFF).final.controls.equal.index    = [];
    obj.problem.phases(ph+obj.MIOFF).final.controls.equal.value    = zeros(obj.problem.phases(ph+obj.MIOFF).final.controls.equal.index,1);
    obj.problem.phases(ph+obj.MIOFF).final.states.lower.index      = [];
    obj.problem.phases(ph+obj.MIOFF).final.states.lower.value      = zeros(obj.problem.phases(ph+obj.MIOFF).final.states.lower.index,1);
    obj.problem.phases(ph+obj.MIOFF).final.controls.lower.index    = [];
    obj.problem.phases(ph+obj.MIOFF).final.controls.lower.value    = zeros(obj.problem.phases(ph+obj.MIOFF).final.controls.lower.index,1);
    obj.problem.phases(ph+obj.MIOFF).final.states.upper.index      = [];
    obj.problem.phases(ph+obj.MIOFF).final.states.upper.value      = zeros(obj.problem.phases(ph+obj.MIOFF).final.states.upper.index,1);
    obj.problem.phases(ph+obj.MIOFF).final.controls.upper.index    = [];
    obj.problem.phases(ph+obj.MIOFF).final.controls.upper.value    = zeros(obj.problem.phases(ph+obj.MIOFF).final.controls.upper.index,1);
    
    % Initializes States and Controls Bounds
    obj.problem.phases(ph+obj.MIOFF).bounds.states.upper           =  zeros(numberOfStates,0);
    obj.problem.phases(ph+obj.MIOFF).bounds.states.lower           =  zeros(numberOfStates,0);
    obj.problem.phases(ph+obj.MIOFF).bounds.controls.upper         =  zeros(numberOfControls,0);
    obj.problem.phases(ph+obj.MIOFF).bounds.controls.lower         =  zeros(numberOfControls,0);
    obj.problem.phases(ph+obj.MIOFF).bounds.timeFinal.upper        =  [];
    obj.problem.phases(ph+obj.MIOFF).bounds.timeFinal.lower        =  [];
    obj.problem.phases(ph+obj.MIOFF).bounds.timeInitial.upper      =  [];
    obj.problem.phases(ph+obj.MIOFF).bounds.timeInitial.lower      =  [];
    
    
    
    % Initializes Event Constraints
    obj.problem.phases(ph+obj.MIOFF).events  = struct;
    for ii = 0:obj.problem.phases(ph+obj.MIOFF).n.events-1
        obj.problem.phases(ph+obj.MIOFF).events(ii+obj.MIOFF).limit             = 0;
        obj.problem.phases(ph+obj.MIOFF).events(ii+obj.MIOFF).funType           = 'undefined'; % default
        obj.problem.phases(ph+obj.MIOFF).events(ii+obj.MIOFF).type              = 'undefined'; % default
        obj.problem.phases(ph+obj.MIOFF).events(ii+obj.MIOFF).where             = 'undefined'; % default
        obj.problem.phases(ph+obj.MIOFF).events(ii+obj.MIOFF).indexWhere        = []; % default is empty, only if index where type
        obj.problem.phases(ph+obj.MIOFF).events(ii+obj.MIOFF).cone              = obj.initCone(3,obj.problem.phases(ph+obj.MIOFF).n);
        obj.problem.phases(ph+obj.MIOFF).events(ii+obj.MIOFF).states            = zeros(numberOfStates,1);
        obj.problem.phases(ph+obj.MIOFF).events(ii+obj.MIOFF).controls          = zeros(numberOfControls,1);
        obj.problem.phases(ph+obj.MIOFF).events(ii+obj.MIOFF).cons              = 0;
        obj.problem.phases(ph+obj.MIOFF).events(ii+obj.MIOFF).order             = 1; % or 2nd order
        obj.problem.phases(ph+obj.MIOFF).events(ii+obj.MIOFF).variableIndex     = nan;
        obj.problem.phases(ph+obj.MIOFF).events(ii+obj.MIOFF).statesIndex       = 1:numberOfStates;
        obj.problem.phases(ph+obj.MIOFF).events(ii+obj.MIOFF).controlsIndex     = 1:numberOfControls;
        obj.problem.phases(ph+obj.MIOFF).events(ii+obj.MIOFF).variableType      = 'undefined';
        obj.problem.phases(ph+obj.MIOFF).events(ii+obj.MIOFF).function          = @(t,X,U,BodyMap,ii) zeros(1,1);
        obj.problem.phases(ph+obj.MIOFF).events(ii+obj.MIOFF).jacobian.states   = @(t,X,U,BodyMap,ii) zeros(numberOfStates,1);
        obj.problem.phases(ph+obj.MIOFF).events(ii+obj.MIOFF).jacobian.controls = @(t,X,U,BodyMap,ii) zeros(numberOfControls,1);
        obj.problem.phases(ph+obj.MIOFF).events(ii+obj.MIOFF).jacobian.variable = @(t,X,U,BodyMap,ii) zeros(1,1); % this is for 2nd order 1var only
        obj.problem.phases(ph+obj.MIOFF).events(ii+obj.MIOFF).hessian.variable  = @(t,X,U,BodyMap,ii) zeros(1,1); % this is for 2nd order 1var only
        obj.problem.phases(ph+obj.MIOFF).events(ii+obj.MIOFF).scale             = 1;
        obj.problem.phases(ph+obj.MIOFF).events(ii+obj.MIOFF).buffer.include    = 0; % Only if Virtual Buffers are Enabled
        obj.problem.phases(ph+obj.MIOFF).events(ii+obj.MIOFF).buffer.lambda     = 1; % Buffer multiplier term. The higher the more penalty
        obj.problem.phases(ph+obj.MIOFF).events(ii+obj.MIOFF).buffer.penalty    = 1; % Buffer Penalty Type for violation (1 or 2)
        obj.problem.phases(ph+obj.MIOFF).events(ii+obj.MIOFF).buffer.penaltyAdaptive  = 1; % Buffer adaptive penalty type
    end
    
    % Initializes Path Constraints
    obj.problem.phases(ph+obj.MIOFF).path  = struct;
    for ii = 0:obj.problem.phases(ph+obj.MIOFF).n.path-1
        obj.problem.phases(ph+obj.MIOFF).path(ii+obj.MIOFF).limit             = 0;
        obj.problem.phases(ph+obj.MIOFF).path(ii+obj.MIOFF).funType           = 'undefined'; % default.
        obj.problem.phases(ph+obj.MIOFF).path(ii+obj.MIOFF).type              = 'undefined'; % default
        obj.problem.phases(ph+obj.MIOFF).path(ii+obj.MIOFF).indexWhere        = []; % default is empty, not implemented but if empty, should apply [1:Node]
        obj.problem.phases(ph+obj.MIOFF).path(ii+obj.MIOFF).integral          = 0;
        obj.problem.phases(ph+obj.MIOFF).path(ii+obj.MIOFF).derivative        = 0;
        obj.problem.phases(ph+obj.MIOFF).path(ii+obj.MIOFF).cone              = obj.initCone(3,obj.problem.phases(ph+obj.MIOFF).n);
        obj.problem.phases(ph+obj.MIOFF).path(ii+obj.MIOFF).states            = zeros(numberOfStates,1);
        obj.problem.phases(ph+obj.MIOFF).path(ii+obj.MIOFF).controls          = zeros(numberOfControls,1);
        obj.problem.phases(ph+obj.MIOFF).path(ii+obj.MIOFF).cons              = 0;
        obj.problem.phases(ph+obj.MIOFF).path(ii+obj.MIOFF).order             = 1; % or 2nd order
        obj.problem.phases(ph+obj.MIOFF).path(ii+obj.MIOFF).variableIndex     = nan;
        obj.problem.phases(ph+obj.MIOFF).path(ii+obj.MIOFF).statesIndex       = 1:numberOfStates;
        obj.problem.phases(ph+obj.MIOFF).path(ii+obj.MIOFF).controlsIndex     = 1:numberOfControls;
        obj.problem.phases(ph+obj.MIOFF).path(ii+obj.MIOFF).variableType      = 'undefined';
        obj.problem.phases(ph+obj.MIOFF).path(ii+obj.MIOFF).function          = @(t,X,U,BodyMap,ii) zeros(1,size(X,2));
        obj.problem.phases(ph+obj.MIOFF).path(ii+obj.MIOFF).jacobian.states   = @(t,X,U,BodyMap,ii) zeros(numberOfStates,size(X,2));
        obj.problem.phases(ph+obj.MIOFF).path(ii+obj.MIOFF).jacobian.controls = @(t,X,U,BodyMap,ii) zeros(numberOfControls,size(X,2));
        obj.problem.phases(ph+obj.MIOFF).path(ii+obj.MIOFF).jacobian.variable = @(t,X,U,BodyMap,ii) zeros(1,size(X,2));
        obj.problem.phases(ph+obj.MIOFF).path(ii+obj.MIOFF).hessian.issparse  = 0;
        obj.problem.phases(ph+obj.MIOFF).path(ii+obj.MIOFF).hessian.states    = @(t,X,U,BodyMap,ii) zeros(numberOfStates,numberOfStates,length(t));
        obj.problem.phases(ph+obj.MIOFF).path(ii+obj.MIOFF).hessian.controls  = @(t,X,U,BodyMap,ii) zeros(numberOfControls,numberOfControls,length(t));
        obj.problem.phases(ph+obj.MIOFF).path(ii+obj.MIOFF).hessian.crossXU   = @(t,X,U,BodyMap,ii) zeros(numberOfStates,numberOfControls,length(t));
        obj.problem.phases(ph+obj.MIOFF).path(ii+obj.MIOFF).hessian.crossUX   = @(t,X,U,BodyMap,ii) zeros(numberOfControls,numberOfStates,length(t));
        obj.problem.phases(ph+obj.MIOFF).path(ii+obj.MIOFF).hessian.variable  = @(t,X,U,BodyMap,ii) zeros(1,size(X,2));
        obj.problem.phases(ph+obj.MIOFF).path(ii+obj.MIOFF).scale             = 1;
        obj.problem.phases(ph+obj.MIOFF).path(ii+obj.MIOFF).buffer.include    = 0; % Only if Virtual Buffers are Enabled
        obj.problem.phases(ph+obj.MIOFF).path(ii+obj.MIOFF).buffer.lambda     = 1; % Buffer multiplier term. The higher the more penalty
        obj.problem.phases(ph+obj.MIOFF).path(ii+obj.MIOFF).buffer.penalty    = 1; % Buffer Penalty Type for violation (1 or 2)
        obj.problem.phases(ph+obj.MIOFF).path(ii+obj.MIOFF).buffer.penaltyAdaptive    = 1; % Buffer adaptive penalty type
        
    end
    
    % Initializes states and controls statewise trust regions
    % and final time trust regions
    obj.algorithm.sequential.trustRegion.phases(ph+obj.MIOFF).states.component    = ones(numberOfStates,1);
    obj.algorithm.sequential.trustRegion.phases(ph+obj.MIOFF).controls.component  = ones(numberOfControls,1);
    obj.algorithm.sequential.trustRegion.phases(ph+obj.MIOFF).timeFinal.radius    = 1.0;
    obj.algorithm.sequential.trustRegion.phases(ph+obj.MIOFF).timeFinal.lambda    = 1.0;
    obj.algorithm.sequential.trustRegion.phases(ph+obj.MIOFF).timeFinal.type      = 'soft';
    obj.algorithm.sequential.trustRegion.phases(ph+obj.MIOFF).timeFinal.adaptive  = 0;
    obj.algorithm.sequential.trustRegion.phases(ph+obj.MIOFF).controls.lambda     = 1.0e0;
    
    % Virtual Control fields.
    obj.algorithm.virtualControl.phases(ph+obj.MIOFF).include = 0; % Inactive by default
    obj.algorithm.virtualControl.phases(ph+obj.MIOFF).states  = 1:numberOfStates; % Index specifying states affected
    obj.algorithm.virtualControl.phases(ph+obj.MIOFF).lambda  = 1e0; % Augmentation term of the cost function
    obj.algorithm.virtualControl.phases(ph+obj.MIOFF).statePenalty = 1 ; % 1st or 2nd norm for the state-wise norm
    obj.algorithm.virtualControl.phases(ph+obj.MIOFF).nodePenalty  = 1 ; % 1st or 2nd norm for the node wise norm (applied to all state-wise norms)
    obj.algorithm.virtualControl.phases(ph+obj.MIOFF).scale   = nan; % Setted later
    
    % Virtual Buffer fields
    obj.algorithm.virtualBuffer.phases(ph+obj.MIOFF).include      = 0; % Inactive by default
    obj.algorithm.virtualBuffer.phases(ph+obj.MIOFF).nodePenalty  = 1 ; % 1st or 2nd norm for all buffers
    
    % Initial guess fields for all variables
    obj.problem.phases(ph+obj.MIOFF).initialGuess.timeFinal   = []; % empty by default, needs an initialization (if fixed, a value)
    obj.problem.phases(ph+obj.MIOFF).initialGuess.timeInitial = 0; % 0 by default
    obj.problem.phases(ph+obj.MIOFF).initialGuess.time      = zeros(1,numberOfNodes); % time is zero by default, for compatibility with SOCP solver
    obj.problem.phases(ph+obj.MIOFF).initialGuess.states    = zeros(numberOfStates,numberOfNodes); % time is zero by default, for compatibility with SOCP solver
    obj.problem.phases(ph+obj.MIOFF).initialGuess.controls  = zeros(numberOfControls,numberOfNodes); % time is zero by default, for compatibility with SOCP solver
    
end

% Objective
obj.problem.objective.type                       = 'minimize'; % minimize or maximize. Minimisation by default
obj.problem.objective.sign                       = [] ; % Initialized in level 2 with type
obj.problem.objective.scale                      = 1  ; % Scale term for the objective

% Objective Mayer term for a specific node
obj.problem.objective.mayer.funType              = 'undefined';
obj.problem.objective.mayer.where                = 'undefined';
obj.problem.objective.mayer.indexWhere           = [];
obj.problem.objective.mayer.absolute.term        = 0; % absolute of whole function
obj.problem.objective.mayer.absolute.states      = 0; % absolute function of states
obj.problem.objective.mayer.absolute.controls    = 0; % absolute function of controls
obj.problem.objective.mayer.function             = @(t,X,U,BodyMap,ii) 0;
obj.problem.objective.mayer.jacobian.states      = @(t,X,U,BodyMap,ii) zeros(numberOfStates,1);
obj.problem.objective.mayer.jacobian.controls    = @(t,X,U,BodyMap,ii) zeros(numberOfControls,1);
obj.problem.objective.mayer.states               = zeros(numberOfStates,1);
obj.problem.objective.mayer.controls             = zeros(numberOfControls,1);
obj.problem.objective.mayer.timeFinal            = 0;
obj.problem.objective.mayer.cone                 = obj.initCone(2,obj.problem.phases(1).n);

% Objective Lagrange term on all nodes
obj.problem.objective.lagrange.funType           = 'undefined';
obj.problem.objective.lagrange.absolute.term     = 0; % absolute of whole function
obj.problem.objective.lagrange.absolute.nodes    = 0; % absolute of all nodes
obj.problem.objective.lagrange.absolute.states   = 0; % absolute function of states
obj.problem.objective.lagrange.absolute.controls = 0; % absolute function of controls
obj.problem.objective.lagrange.function          = @(t,X,U,BodyMap,ii) zeros(1,size(X,2));
obj.problem.objective.lagrange.jacobian.states   = @(t,X,U,BodyMap,ii) zeros(numberOfStates,size(X,2));
obj.problem.objective.lagrange.jacobian.controls = @(t,X,U,BodyMap,ii) zeros(numberOfControls,size(X,2));
obj.problem.objective.lagrange.states            = zeros(numberOfStates,1);
obj.problem.objective.lagrange.controls          = zeros(numberOfControls,1);
obj.problem.objective.lagrange.timeFinal         = 0;
obj.problem.objective.lagrange.cone              = obj.initCone(2,obj.problem.phases(1).n);

% Objective minimisation(maximization) term on node with larger(lower) value
obj.problem.objective.minmax.funType           = 'undefined';
obj.problem.objective.minmax.absolute.term     = 0;
obj.problem.objective.minmax.absolute.nodes    = 0;
obj.problem.objective.minmax.absolute.states   = 0;
obj.problem.objective.minmax.absolute.controls = 0;
obj.problem.objective.minmax.function          = @(t,X,U,BodyMap,ii) zeros(1,size(X,2));
obj.problem.objective.minmax.jacobian.states   = @(t,X,U,BodyMap,ii) zeros(numberOfStates,size(X,2));
obj.problem.objective.minmax.jacobian.controls = @(t,X,U,BodyMap,ii) zeros(numberOfControls,size(X,2));
obj.problem.objective.minmax.states            = zeros(numberOfStates,1);
obj.problem.objective.minmax.controls          = zeros(numberOfControls,1);
obj.problem.objective.minmax.timeFinal         = 0;
obj.problem.objective.minmax.cone              = obj.initCone(2,obj.problem.phases(1).n);

% Objective regularisation term on whole function, can be
% inserted in lagrange, but included as additional type as
% dependency on time removed
obj.problem.objective.regularisation.weight            = 1;
obj.problem.objective.regularisation.states            = zeros(numberOfStates,1);
obj.problem.objective.regularisation.controls          = zeros(numberOfControls,1);
obj.problem.objective.regularisation.timeFinal         = 0;

end

