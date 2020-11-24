function obj = setUpSCOPT_level2(obj)
%% Level 2 SCOPT Setup
% Defines objective sign
% Defines number of internal slack variables and indexes
% Creates initial guess for states, controls and time based on
% mid values of bounds.
% Performs automatic scaling routines
if ~obj.algorithm.sequential.activate
    % If no successive convexificatiobn/linearisation, set
    % maximum iterations to 1 and minimum to 0
    obj.algorithm.sequential.maxIter = 1;
    obj.algorithm.sequential.minIter = 0;
end
%% Objective
switch obj.problem.objective.type
    case 'minimize'
        obj.problem.objective.sign = 1;
    case 'maximize'
        obj.problem.objective.sign = -1;
end
%% Definition of Number of Variables
Nodes       = obj.problem.phases.Nodes;
n           = obj.problem.phases.n;
if strcmp(obj.problem.objective.type,'feasibility')
    n.objectiveSlack = 0;
else
    n.objectiveSlack = 1;
end

if obj.problem.objective.lagrange.absolute.nodes || strcmp(obj.problem.objective.lagrange.funType,'convex') || strcmp(obj.problem.objective.lagrange.funType,'soc')
    n.objectiveLagrange = 1;
else
    n.objectiveLagrange = 0;
end

if obj.algorithm.virtualControl.phases(1).include
    n.virtualControls = length(obj.algorithm.virtualControl.phases(1).states);
    
    switch obj.algorithm.virtualControl.phases(1).statePenalty
        case 0
            n.virtualControlsP1      = 1;
            n.virtualControlsP1s = 0;
        case 1
            n.virtualControlsP1      = 1;
            n.virtualControlsP1s = length(obj.algorithm.virtualControl.phases(1).states);
        case 2
            n.virtualControlsP1      = 1;
            n.virtualControlsP1s = 0;
        case {inf,Inf}
            n.virtualControlsP1      = 1;
            n.virtualControlsP1s = 0;
        otherwise
            n.virtualControlsP1      = 0;
            n.virtualControlsP1s = 0;
    end
    n.virtualControlsP2 = 1;
else
    n.virtualControls         = 0;
    n.virtualControlsP1       = 0;
    n.virtualControlsP1s  = 0;
    n.virtualControlsP2       = 0;
end


if obj.problem.phases(1).freeTimeFinal
    % Free time of flight
    n.timeFinal = 1 ;
    switch obj.algorithm.sequential.trustRegion.phases(1).timeFinal.type
        case 'soft'
            n.timeFinalSlack = 1;
        case {'hard','none'}
            n.timeFinalSlack = 0;
    end
else
    % Fixed time of flight
    n.timeFinal = 0 ;
    n.timeFinalSlack = 0;
end


if obj.algorithm.virtualBuffer.phases(1).include ==1
    funIndex  = @(s) s.buffer.include;
    funLambda = @(s) s.buffer.lambda;
    if obj.problem.phases(1).n.path>0
        obj.algorithm.virtualBuffer.phases(1).path.index      = find(arrayfun(funIndex,obj.problem.phases(1).path));
        obj.algorithm.virtualBuffer.phases(1).path.lambda     = arrayfun(funLambda,obj.problem.phases(1).path(obj.algorithm.virtualBuffer.phases(1).path.index));
    else
        obj.algorithm.virtualBuffer.phases(1).path.index = [];
        obj.algorithm.virtualBuffer.phases(1).path.lambda = [];
    end
    if obj.problem.phases(1).n.events>0
        obj.algorithm.virtualBuffer.phases(1).events.index    = find(arrayfun(funIndex,obj.problem.phases(1).events));
        obj.algorithm.virtualBuffer.phases(1).events.lambda   = arrayfun(funLambda,obj.problem.phases(1).events(obj.algorithm.virtualBuffer.phases(1).events.index));
    else
        obj.algorithm.virtualBuffer.phases(1).events.index  = [];
        obj.algorithm.virtualBuffer.phases(1).events.lambda = [];
    end
    n.virtualBuffers.events       = length(obj.algorithm.virtualBuffer.phases(1).events.index);
    n.virtualBuffers.eventsP      = n.virtualBuffers.events ;
    
    n.virtualBuffers.path       = length(obj.algorithm.virtualBuffer.phases(1).path.index);
    n.virtualBuffers.pathP      = n.virtualBuffers.path;
    n.virtualBuffers.pathPs = n.virtualBuffers.path;
else
    obj.algorithm.virtualBuffer.phases(1).path.index = [];
    obj.algorithm.virtualBuffer.phases(1).path.lambda = [];
    obj.algorithm.virtualBuffer.phases(1).events.index  = [];
    obj.algorithm.virtualBuffer.phases(1).events.lambda = [];
    n.virtualBuffers.events       = 0;
    n.virtualBuffers.eventsP      = 0;
    n.virtualBuffers.path       = 0;
    n.virtualBuffers.pathP      = 0;
    n.virtualBuffers.pathPs = 0;
end


%             if strcmp(obj.algorithm.sequential.type,'trust-region-soft') || strcmp(obj.algorithm.sequential.type,'adaptive-trust-region-soft')
if strcmp(obj.algorithm.sequential.type,'trust-region') && obj.algorithm.sequential.activate
    n.trustRegionP1 = 1;
    
    n.trustRegionP1s = 0;
    if obj.algorithm.sequential.trustRegion.variablesPenalty == 1
        if obj.algorithm.sequential.trustRegion.include.controls
            n.trustRegionP1s = n.trustRegionP1s + n.controls;
        end
        if obj.algorithm.sequential.trustRegion.include.states
            n.trustRegionP1s = n.trustRegionP1s + n.states;
        end
    end
    
    n.trustRegionP2 = 1;
else
    n.trustRegionP1 = 0;
    n.trustRegionP1s = 0;
    n.trustRegionP2 = 0;
end

obj.problem.phases(1).n = n    ;

obj.problem.phases(1).Ntot = (n.states+n.controls)*Nodes ...
    + (Nodes-1)*(n.virtualControls + n.virtualControlsP1 + n.virtualControlsP1s) + n.virtualControlsP2...
    + n.virtualBuffers.events + n.virtualBuffers.eventsP + Nodes*(n.virtualBuffers.path + n.virtualBuffers.pathPs) + n.virtualBuffers.pathP...
    + n.objectiveLagrange*Nodes + n.objectiveSlack...
    + (n.trustRegionP1+n.trustRegionP1s)*Nodes + n.trustRegionP2  + n.timeFinalSlack + n.timeFinal;

%% Setup Indexing
index.nodes    = (0:(Nodes-1)) + 1;
index.states   = repmat(1 + (0:n.states:(n.states*(Nodes-1))),n.states,1) + repmat((0:(n.states-1))',1,Nodes);
index.controls = repmat(Nodes*n.states + 1 + (0:n.controls:((Nodes-1)*n.controls)),n.controls,1) + repmat((0:(n.controls-1))',1,Nodes);
if obj.algorithm.virtualControl.phases(1).include
    index.virtualControls    = repmat(Nodes*(n.states+n.controls) + obj.MIOFF + (0:n.virtualControls:((Nodes-2)*n.virtualControls)),n.virtualControls,1) + repmat((0:(n.virtualControls-1))',1,Nodes-1);
    index.virtualControlsP1  = index.virtualControls(end) + 1 + (0:((Nodes-2)))  ;
    if n.virtualControlsP1s>0
        index.virtualControlsP1s = repmat(index.virtualControlsP1(end) + 1 + (0:n.virtualControlsP1s:((Nodes-2)*n.virtualControlsP1s)),n.virtualControlsP1s,1) + repmat((0:(n.virtualControlsP1s-1))',1,Nodes-1);
    else
        index.virtualControlsP1s = [];
    end
    index.virtualControlsP2  = index.virtualControls(end) + 1 + (Nodes-1) + n.virtualControlsP1s*(Nodes-1);
else
    index.virtualControls    = [];
    index.virtualControlsP1  = [];
    index.virtualControlsP1s = [];
    index.virtualControlsP2  = [];
end
if obj.algorithm.virtualBuffer.phases(1).include ==1
    if n.virtualBuffers.events>0
        index.virtualBuffers.events  = Nodes*(n.states+n.controls)...
            + (Nodes-1)*(n.virtualControls + n.virtualControlsP1 + n.virtualControlsP1s) + n.virtualControlsP2 +obj.MIOFF...
            + (0:n.virtualBuffers.events-1)';
        index.virtualBuffers.eventsP = index.virtualBuffers.events(end) + 1  + (0:n.virtualBuffers.eventsP-1)';
    else
        index.virtualBuffers.events  = [];
        index.virtualBuffers.eventsP = [];
    end
    
    if n.virtualBuffers.path>0
        index.virtualBuffers.path    = repmat(...
            Nodes*(n.states+n.controls)...
            + (Nodes-1)*(n.virtualControls + n.virtualControlsP1 + n.virtualControlsP1s) + n.virtualControlsP2...
            + n.virtualBuffers.events + n.virtualBuffers.eventsP +obj.MIOFF...
            +(0:n.virtualBuffers.path:((Nodes-1)*n.virtualBuffers.path))...
            ,n.virtualBuffers.path,1) + repmat((0:(n.virtualBuffers.path-1))',1,Nodes);
        
        index.virtualBuffers.pathPs  = repmat(...
            index.virtualBuffers.path(end) + 1+(0:n.virtualBuffers.pathPs:((Nodes-1)*n.virtualBuffers.pathPs))...
            ,n.virtualBuffers.pathPs,1) + repmat((0:(n.virtualBuffers.pathPs-1))',1,Nodes);
        
        index.virtualBuffers.pathP   = index.virtualBuffers.pathPs(end) + 1 + (0:n.virtualBuffers.pathP-1)';
    else
        index.virtualBuffers.path   = [];
        index.virtualBuffers.pathPs = [];
        index.virtualBuffers.pathP  = [];
    end
else
    index.virtualBuffers.events  = [];
    index.virtualBuffers.eventsP = [];
    index.virtualBuffers.path    = [];
    index.virtualBuffers.pathPs  = [];
    index.virtualBuffers.pathP   = [];
end

if n.objectiveLagrange
    index.objectiveLagrange =  Nodes*(n.states+n.controls)...
        + (Nodes-1)*(n.virtualControls + n.virtualControlsP1 + n.virtualControlsP1s) + n.virtualControlsP2 ...
        + n.virtualBuffers.events + n.virtualBuffers.eventsP + Nodes*(n.virtualBuffers.path + n.virtualBuffers.pathPs) + n.virtualBuffers.pathP +obj.MIOFF...
        + (0:((Nodes-1)))  ;
else
    index.objectiveLagrange = [];
end

if n.objectiveSlack
    index.objectiveSlack  = Nodes*(n.states+n.controls)+ ...
        (Nodes-1)*(n.virtualControls + n.virtualControlsP1 + n.virtualControlsP1s) + n.virtualControlsP2  ...
        + n.virtualBuffers.events + n.virtualBuffers.eventsP + Nodes*(n.virtualBuffers.path + n.virtualBuffers.pathPs) + n.virtualBuffers.pathP ...
        + n.objectiveLagrange*Nodes +obj.MIOFF ;
else
    index.objectiveSlack  = [];
end

if obj.algorithm.sequential.activate && strcmp(obj.algorithm.sequential.type,'trust-region')
    index.trustRegionP1  = Nodes*(n.states+n.controls)+ (Nodes-1)*(n.virtualControls + n.virtualControlsP1 + n.virtualControlsP1s) + n.virtualControlsP2  + n.virtualBuffers.events + n.virtualBuffers.eventsP + Nodes*(n.virtualBuffers.path + n.virtualBuffers.pathPs) + n.virtualBuffers.pathP + n.objectiveLagrange*Nodes + n.objectiveSlack +obj.MIOFF + (0:((Nodes-1)))  ;
    if n.trustRegionP1s>0
        index.trustRegionP1s = repmat(index.trustRegionP1(end) + 1 + (0:n.trustRegionP1s:((Nodes-1)*n.trustRegionP1s))...
            ,n.trustRegionP1s,1) + repmat((0:(n.trustRegionP1s-1))',1,Nodes);
    else
        index.trustRegionP1s = [];
    end
    index.trustRegionP2  = Nodes*(n.states+n.controls)+ (Nodes-1)*(n.virtualControls + n.virtualControlsP1 + n.virtualControlsP1s) + n.virtualControlsP2  + n.virtualBuffers.events + n.virtualBuffers.eventsP + Nodes*(n.virtualBuffers.path + n.virtualBuffers.pathPs) + n.virtualBuffers.pathP + n.objectiveLagrange*Nodes + n.objectiveSlack + (n.trustRegionP1+n.trustRegionP1s)*Nodes +obj.MIOFF ;
else
    index.trustRegionP1  = [];
    index.trustRegionP1s = [];
    index.trustRegionP2  = [];
end

if obj.problem.phases(1).freeTimeFinal
    index.timeFinal      = obj.problem.phases(1).Ntot - 1 +obj.MIOFF;
    if strcmp(obj.algorithm.sequential.trustRegion.phases(1).timeFinal.type,'soft')
        index.timeFinalSlack = obj.problem.phases(1).Ntot - 2 +obj.MIOFF;
    else
        index.timeFinalSlack = [];
    end
else
    index.timeFinal      = [];
    index.timeFinalSlack = [];
end
obj.problem.phases(1).index = index;
 
%% Defects Penalty
% If not modified, the same virtual contorl type of penalties
% are used for teh adaptive trust region error metric
if isnan(obj.algorithm.sequential.trustRegion.defectsCompErrorPenalty)
    obj.algorithm.sequential.trustRegion.defectsCompErrorPenalty  = obj.algorithm.virtualControl.phases.statePenalty   ;
end
if isnan(obj.algorithm.sequential.trustRegion.defectsNodeErrorPenalty)
    obj.algorithm.sequential.trustRegion.defectsNodeErrorPenalty  = obj.algorithm.virtualControl.phases.nodePenalty   ;
end
%% Initial-Guess
switch obj.algorithm.initialGuessType
    case 'automatic' % Compute Initial Guess Based on Variable Bounds
        for ph=0:obj.problem.nphases-1
            if obj.problem.phases(ph+obj.MIOFF).freeTimeFinal && ~isempty(obj.problem.phases(ph+obj.MIOFF).initialGuess.timeFinal)
                obj.problem.phases(ph+obj.MIOFF).initialGuess.timeFinal = (obj.problem.phases(ph+obj.MIOFF).bounds.timeFinal.upper+obj.problem.phases(ph+obj.MIOFF).bounds.timeFinal.lower)/2;
            elseif isempty(obj.problem.phases(ph+obj.MIOFF).initialGuess.timeFinal) || obj.problem.phases(ph+obj.MIOFF).initialGuess.timeFinal==0
                error('SCOPT Error: Final time not specified or Zero \n Abort \n')
            end
            if obj.problem.phases(ph+obj.MIOFF).freeTimeInitial
                obj.problem.phases(ph+obj.MIOFF).initialGuess.timeInitial = (obj.problem.phases(ph+obj.MIOFF).bounds.timeInitial.upper+obj.problem.phases(ph+obj.MIOFF).bounds.timeInitial.lower)/2;
            end
            obj.problem.phases(ph+obj.MIOFF).initialGuess.time     = linspace(obj.problem.phases(ph+obj.MIOFF).initialGuess.timeInitial,obj.problem.phases(ph+obj.MIOFF).initialGuess.timeFinal,obj.problem.phases(ph+obj.MIOFF).Nodes);
            
            if ~(isempty(obj.problem.phases(ph+obj.MIOFF).bounds.states.upper) && isempty(obj.problem.phases(ph+obj.MIOFF).bounds.states.lower))
                obj.problem.phases(ph+obj.MIOFF).initialGuess.states   = repmat((obj.problem.phases(ph+obj.MIOFF).bounds.states.upper   + obj.problem.phases(ph+obj.MIOFF).bounds.states.lower  )/2,1,obj.problem.phases(ph+obj.MIOFF).Nodes);
            end
            if ~(isempty(obj.problem.phases(ph+obj.MIOFF).bounds.controls.upper) && isempty(obj.problem.phases(ph+obj.MIOFF).bounds.controls.lower))
                obj.problem.phases(ph+obj.MIOFF).initialGuess.controls   = repmat((obj.problem.phases(ph+obj.MIOFF).bounds.controls.upper   + obj.problem.phases(ph+obj.MIOFF).bounds.controls.lower  )/2,1,obj.problem.phases(ph+obj.MIOFF).Nodes);
            end
        end
    case 'user'
        for ph=0:obj.problem.nphases-1
            if ~isempty(obj.problem.phases(ph+obj.MIOFF).initialGuess.timeFinal) && ~isempty(obj.problem.phases(ph+obj.MIOFF).initialGuess.timeInitial)
                obj.problem.phases(ph+obj.MIOFF).initialGuess.time     = linspace(obj.problem.phases(ph+obj.MIOFF).initialGuess.timeInitial,obj.problem.phases(ph+obj.MIOFF).initialGuess.timeFinal,obj.problem.phases(ph+obj.MIOFF).Nodes);
            end
        end
    case 'none'
end

for ph=0:obj.problem.nphases-1
    obj.problem.phases(ph+obj.MIOFF).initialGuess.tauDisc            = 2/(obj.problem.phases(ph+obj.MIOFF).Nodes-1);
    obj.problem.phases(ph+obj.MIOFF).initialGuess.timeDisc           = obj.problem.phases(ph+obj.MIOFF).initialGuess.tauDisc/2*(obj.problem.phases(ph+obj.MIOFF).initialGuess.timeFinal-obj.problem.phases(ph+obj.MIOFF).initialGuess.timeInitial);
end
%% Quadrature Vector Initialisation
% dtau = 2*dt/(tf-t0)
for ph = 0:obj.problem.nphases-1
    switch obj.algorithm.quadrature.type
        case {'euler','forward'} % sum(dti*ui) = sum(1/2*dtaui*(tf-t0)*ui)
            obj.problem.phases(ph+obj.MIOFF).quadratureVector    = 1/2*obj.problem.phases(ph+obj.MIOFF).initialGuess.tauDisc*ones(1,obj.problem.phases(ph+obj.MIOFF).Nodes);
            obj.problem.phases(ph+obj.MIOFF).quadratureVector(obj.problem.phases(ph+obj.MIOFF).Nodes) = 0;
        case 'trapezoidal'
            obj.problem.phases(ph+obj.MIOFF).quadratureVector    = 1/2*obj.problem.phases(ph+obj.MIOFF).initialGuess.tauDisc*ones(1,obj.problem.phases(ph+obj.MIOFF).Nodes);
            obj.problem.phases(ph+obj.MIOFF).quadratureVector(1) = 1/4*obj.problem.phases(ph+obj.MIOFF).initialGuess.tauDisc;
            obj.problem.phases(ph+obj.MIOFF).quadratureVector(obj.problem.phases(ph+obj.MIOFF).Nodes) = 1/4*obj.problem.phases(ph+obj.MIOFF).initialGuess.tauDisc;
        otherwise
            error('SCOPT Error: Numerical integration quadrature type not defined \n')
    end
end
%% Interpolate Values for Initial Guess
if obj.problem.phases(1).interpolateParameters && obj.algorithm.speedUp
    [obj.problem.phases(1).BodyMap] = obj.problem.phases(1).interpolationsFunction(obj.problem.phases(1).initialGuess.time,obj.problem.phases(1).initialGuess.states,obj.problem.phases(1).initialGuess.controls,obj.problem.phases(1).BodyMap);
    obj.problem.phases(1).BodyMap.evaluated = 1;
else
    obj.problem.phases(1).BodyMap.evaluated = 0;
end
obj.problem.initialGuess.Jreal     = obj.computeObjectiveFunctions(obj.problem.phases.initialGuess.time,obj.problem.phases.initialGuess.states,obj.problem.phases.initialGuess.controls);

%% Scaling Settings

% Variables
switch obj.algorithm.scaling.variables
    case 'automatic'
        for ph = 0:obj.problem.nphases-1
            if ~isempty(obj.problem.phases(ph+obj.MIOFF).bounds.states.upper) && ~isempty(obj.problem.phases(ph+obj.MIOFF).bounds.states.lower)
                obj.problem.phases(ph+obj.MIOFF).scale.states   = (obj.problem.phases(ph+obj.MIOFF).bounds.states.upper(:) - obj.problem.phases(ph+obj.MIOFF).bounds.states.lower(:))/2;
                obj.problem.phases(ph+obj.MIOFF).shift.states   = (obj.problem.phases(ph+obj.MIOFF).bounds.states.upper(:) + obj.problem.phases(ph+obj.MIOFF).bounds.states.lower(:))/2;
                if obj.problem.phases(ph+obj.MIOFF).scale.states == 0
                    obj.problem.phases(ph+obj.MIOFF).scale.states = 1;
                    obj.problem.phases(ph+obj.MIOFF).shift.states = 0;
                end
                zeroScale = find(obj.problem.phases(ph+obj.MIOFF).scale.states==0);
                obj.problem.phases(ph+obj.MIOFF).scale.states(zeroScale) = 1;
                obj.problem.phases(ph+obj.MIOFF).shift.states(zeroScale) = 0;
            end
            if ~isempty(obj.problem.phases(ph+obj.MIOFF).bounds.controls.upper) && ~isempty(obj.problem.phases(ph+obj.MIOFF).bounds.controls.lower)
                for jj = 0:n.states-1
                    obj.problem.phases(ph+obj.MIOFF).scale.controls   = (obj.problem.phases(ph+obj.MIOFF).bounds.controls.upper(:) - obj.problem.phases(ph+obj.MIOFF).bounds.controls.lower(:))/2;
                    obj.problem.phases(ph+obj.MIOFF).shift.controls   = (obj.problem.phases(ph+obj.MIOFF).bounds.controls.upper(:) + obj.problem.phases(ph+obj.MIOFF).bounds.controls.lower(:))/2;
                end
                zeroScale = find(obj.problem.phases(ph+obj.MIOFF).scale.controls==0);
                obj.problem.phases(ph+obj.MIOFF).scale.controls(zeroScale) = 1;
                obj.problem.phases(ph+obj.MIOFF).shift.controls(zeroScale) = 0;
            end
            
            if obj.problem.phases(ph+obj.MIOFF).freeTimeFinal && ~obj.problem.phases(ph+obj.MIOFF).freeTimeInitial && ~isempty(obj.problem.phases(ph+obj.MIOFF).bounds.timeFinal.upper)
                obj.problem.phases(ph+obj.MIOFF).scale.time     = (obj.problem.phases(ph+obj.MIOFF).bounds.timeFinal.upper - obj.problem.phases(ph+obj.MIOFF).initialGuess.timeInitial)/2;
                obj.problem.phases(ph+obj.MIOFF).shift.time     = (obj.problem.phases(ph+obj.MIOFF).bounds.timeFinal.upper + obj.problem.phases(ph+obj.MIOFF).initialGuess.timeInitial)/2;
            elseif ~obj.problem.phases(ph+obj.MIOFF).freeTimeFinal && ~obj.problem.phases(ph+obj.MIOFF).freeTimeInitial
                obj.problem.phases(ph+obj.MIOFF).scale.time     = (obj.problem.phases(ph+obj.MIOFF).initialGuess.timeFinal - obj.problem.phases(ph+obj.MIOFF).initialGuess.timeInitial)/2;
                obj.problem.phases(ph+obj.MIOFF).shift.time     = (obj.problem.phases(ph+obj.MIOFF).initialGuess.timeFinal + obj.problem.phases(ph+obj.MIOFF).initialGuess.timeInitial)/2;
            end
        end
    case {'none','user'}
    otherwise
        error('SCOPT Error \nScaling Setting %s not defined \n Aborting \n',obj.algorithm.scaling.variables)
        
end

% Objective
switch obj.algorithm.scaling.objective
    case 'automatic-jacobian'
        % Mayer Term
        scaleTemp = 0;
        switch obj.problem.objective.mayer.funType
            case 'linear'
                scaleTemp = norm([obj.problem.phases(1).scale.time*obj.problem.objective.mayer.timeFinal;obj.problem.phases(1).scale.states(:).*obj.problem.objective.mayer.states(:);obj.problem.phases(1).scale.controls(:).*obj.problem.objective.mayer.controls(:)]);
            case {'convex','soc'}
                scaleTemp = norm([...
                    obj.problem.phases(1).scale.states.*obj.problem.objective.mayer.cone.right.states;...
                    obj.problem.phases(1).scale.controls.*obj.problem.objective.mayer.cone.right.controls;...
                    reshape(obj.problem.phases(1).scale.states.*obj.problem.objective.mayer.cone.norm.states,[],1);...
                    reshape(obj.problem.phases(1).scale.controls.*obj.problem.objective.mayer.cone.norm.controls,[],1)]);
                
            case 'non-linear'
                switch obj.problem.objective.mayer.where
                    case 'final'
                        indexWhere = Nodes - 1 +obj.MIOFF;
                    case 'initial'
                        indexWhere = 0 + obj.MIOFF;
                    case 'index'
                        indexWhere = obj.problem.objective.mayer.indexWhere;
                end
                dG_dX_ii = obj.problem.phases(1).scale.states.*obj.problem.objective.mayer.jacobian.states(obj.problem.phases(ph+obj.MIOFF).initialGuess.time(:,indexWhere),obj.problem.phases(ph+obj.MIOFF).initialGuess.states(:,indexWhere),obj.problem.phases(ph+obj.MIOFF).initialGuess.controls(:,indexWhere),obj.problem.phases(ph+obj.MIOFF).BodyMap,indexWhere);
                dG_dU_ii = obj.problem.phases(1).scale.controls.*obj.problem.objective.mayer.jacobian.controls(obj.problem.phases(ph+obj.MIOFF).initialGuess.time(:,indexWhere),obj.problem.phases(ph+obj.MIOFF).initialGuess.states(:,indexWhere),obj.problem.phases(ph+obj.MIOFF).initialGuess.controls(:,indexWhere),obj.problem.phases(ph+obj.MIOFF).BodyMap,indexWhere);
                scaleTemp = norm([dG_dX_ii;dG_dU_ii]);
        end
        if scaleTemp>0 && ~isnan(scaleTemp)
            obj.problem.objective.scale = scaleTemp;
        end
        
    case {'automatic-initialGuess'}
        scaleAuxObj = abs(obj.problem.initialGuess.Jreal);
        if scaleAuxObj>0
            obj.problem.objective.scale = scaleAuxObj;
        else
            obj.problem.objective.scale = 1;
        end
        
    case {'user'}
        
    case {'none'}
        obj.problem.objective.scale = 1;
    otherwise
        error('SCOPT Error \nScaling Setting %s not defined \n Aborting \n',obj.algorithm.scaling.objective)
        
end
obj.problem.initialGuess.J = obj.problem.initialGuess.Jreal/obj.problem.objective.scale*obj.problem.objective.sign;
% Nonlinear Events
switch obj.algorithm.scaling.events
    case 'automatic-limit'
        for ph = 0:obj.problem.nphases-1
            for ii = 0:obj.problem.phases(ph+obj.MIOFF).n.events-1
                limitTemp = obj.problem.phases(ph+obj.MIOFF).events(ii+obj.MIOFF).limit;
                if limitTemp~=0
                    obj.problem.phases(ph+obj.MIOFF).events(ii+obj.MIOFF).scale = abs(limitTemp);
                end
            end
        end
    case 'automatic-jacobian'
        for ph = 0:obj.problem.nphases-1
            for ii = 0:obj.problem.phases(ph+obj.MIOFF).n.events-1
                scaleTemp = 0;
                switch obj.problem.phases(ph+obj.MIOFF).events(ii+obj.MIOFF).funType
                    case 'linear'
                        scaleTemp = norm([obj.problem.phases(1).scale.states.*obj.problem.phases(ph+obj.MIOFF).events(ii+obj.MIOFF).states;obj.problem.phases(1).scale.controls.*obj.problem.phases(ph+obj.MIOFF).events(ii+obj.MIOFF).controls]);
                    case {'convex','soc'}
                        scaleTemp = norm([obj.problem.phases(1).scale.states.*obj.problem.phases(ph+obj.MIOFF).events(ii+obj.MIOFF).cone.right.states;...
                            obj.problem.phases(1).scale.controls.*obj.problem.phases(ph+obj.MIOFF).events(ii+obj.MIOFF).cone.right.controls;...
                            reshape(obj.problem.phases(1).scale.states.*obj.problem.phases(ph+obj.MIOFF).events(ii+obj.MIOFF).cone.norm.states,[],1);...
                            reshape(obj.problem.phases(1).scale.controls.*obj.problem.phases(ph+obj.MIOFF).events(ii+obj.MIOFF).cone.norm.controls,[],1)]);
                    case 'non-linear'
                        switch obj.problem.phases(ph+obj.MIOFF).events(ii+obj.MIOFF).where
                            case 'final'
                                indexWhere = Nodes - 1 +obj.MIOFF;
                            case 'initial'
                                indexWhere = 0 + obj.MIOFF;
                            case 'index'
                                indexWhere = obj.problem.phases(ph+obj.MIOFF).events(ii+obj.MIOFF).indexWhere;
                        end
                        dG_dX_ii = obj.problem.phases(1).scale.states.*obj.problem.phases(ph+obj.MIOFF).events(ii+obj.MIOFF).jacobian.states(obj.problem.phases(ph+obj.MIOFF).initialGuess.time(:,indexWhere),obj.problem.phases(ph+obj.MIOFF).initialGuess.states(:,indexWhere),obj.problem.phases(ph+obj.MIOFF).initialGuess.controls(:,indexWhere),obj.problem.phases(ph+obj.MIOFF).BodyMap,indexWhere);
                        dG_dU_ii = obj.problem.phases(1).scale.controls.*obj.problem.phases(ph+obj.MIOFF).events(ii+obj.MIOFF).jacobian.controls(obj.problem.phases(ph+obj.MIOFF).initialGuess.time(:,indexWhere),obj.problem.phases(ph+obj.MIOFF).initialGuess.states(:,indexWhere),obj.problem.phases(ph+obj.MIOFF).initialGuess.controls(:,indexWhere),obj.problem.phases(ph+obj.MIOFF).BodyMap,indexWhere);
                        scaleTemp = norm([dG_dX_ii;dG_dU_ii]);
                end
                if scaleTemp>0 || ~isnan(scaleTemp)
                    obj.problem.phases(ph+obj.MIOFF).events(ii+obj.MIOFF).scale = scaleTemp;
                end
            end
        end
    case 'user'
    case 'none'
        for ph = 0:obj.problem.nphases-1
            for ii = 0:obj.problem.phases(ph+obj.MIOFF).n.events-1
                obj.problem.phases(ph+obj.MIOFF).events(ii+obj.MIOFF).scale = 1;
            end
        end
    otherwise
        error('SCOPT Error \nScaling Setting %s not defined \n Aborting \n',obj.algorithm.scaling.events)
        
end
% Nonlinear Path
switch obj.algorithm.scaling.path
    case 'automatic-limit'
        for ph = 0:obj.problem.nphases-1
            for ii = 0:obj.problem.phases(ph+obj.MIOFF).n.path-1
                limitTemp = obj.problem.phases(ph+obj.MIOFF).path(ii+obj.MIOFF).limit;
                if limitTemp~=0
                    obj.problem.phases(ph+obj.MIOFF).path(ii+obj.MIOFF).scale = abs(limitTemp);
                end
            end
        end
    case 'automatic-jacobian'
        for ph = 0:obj.problem.nphases-1
            for ii = 0:obj.problem.phases(ph+obj.MIOFF).n.path-1
                scaleTemp = 0;
                switch obj.problem.phases(ph+obj.MIOFF).path(ii+obj.MIOFF).funType
                    case 'linear'
                        scaleTemp = norm([obj.problem.phases(1).scale.states.*obj.problem.phases(ph+obj.MIOFF).path(ii+obj.MIOFF).states;obj.problem.phases(1).scale.controls.*obj.problem.phases(ph+obj.MIOFF).path(ii+obj.MIOFF).controls]);
                        if ~obj.problem.phases(ph+obj.MIOFF).path(ii+obj.MIOFF).integral && obj.problem.phases(ph+obj.MIOFF).path(ii+obj.MIOFF).derivative
                            scaleTemp = (obj.problem.phases(ph+obj.MIOFF).initialGuess.timeFinal-obj.problem.phases(ph+obj.MIOFF).initialGuess.timeInitial)/(Nodes-1)*scaleTemp;
                        end
                    case {'convex','soc'}
                        scaleTemp = norm([obj.problem.phases(1).scale.states.*obj.problem.phases(ph+obj.MIOFF).path(ii+obj.MIOFF).cone.right.states...
                            ;obj.problem.phases(1).scale.controls.*obj.problem.phases(ph+obj.MIOFF).path(ii+obj.MIOFF).cone.right.controls...
                            ;reshape(repmat(obj.problem.phases(1).scale.states,1,obj.problem.phases(ph+obj.MIOFF).path(ii+obj.MIOFF).cone.dimensions-1).*obj.problem.phases(ph+obj.MIOFF).path(ii+obj.MIOFF).cone.norm.states,[],1)...
                            ;reshape(repmat(obj.problem.phases(1).scale.controls,1,obj.problem.phases(ph+obj.MIOFF).path(ii+obj.MIOFF).cone.dimensions-1).*obj.problem.phases(ph+obj.MIOFF).path(ii+obj.MIOFF).cone.norm.controls,[],1)]);
                    case 'non-linear'
                        dG_dX_ii = obj.problem.phases(1).scale.states.*obj.problem.phases(ph+obj.MIOFF).path(ii+obj.MIOFF).jacobian.states(obj.problem.phases(ph+obj.MIOFF).initialGuess.time,obj.problem.phases(ph+obj.MIOFF).initialGuess.states,obj.problem.phases(ph+obj.MIOFF).initialGuess.controls,obj.problem.phases(ph+obj.MIOFF).BodyMap,index.nodes);
                        dG_dU_ii = obj.problem.phases(1).scale.controls.*obj.problem.phases(ph+obj.MIOFF).path(ii+obj.MIOFF).jacobian.controls(obj.problem.phases(ph+obj.MIOFF).initialGuess.time,obj.problem.phases(ph+obj.MIOFF).initialGuess.states,obj.problem.phases(ph+obj.MIOFF).initialGuess.controls,obj.problem.phases(ph+obj.MIOFF).BodyMap,index.nodes);
                        
                        if obj.problem.phases(ph+obj.MIOFF).path(ii+obj.MIOFF).integral && ~obj.problem.phases(ph+obj.MIOFF).path(ii+obj.MIOFF).derivative
                            %scaleTemp = (obj.problem.phases(ph+obj.MIOFF).initialGuess.timeFinal-obj.problem.phases(ph+obj.MIOFF).initialGuess.timeInitial)/(Nodes-1)*norm([dG_dX_ii(:);dG_dU_ii(:)],2);
                            scaleTemp = abs(obj.problem.phases(ph+obj.MIOFF).path(ii+obj.MIOFF).limit); % JAcobian based not enabled for integral functions
                        elseif ~obj.problem.phases(ph+obj.MIOFF).path(ii+obj.MIOFF).integral && obj.problem.phases(ph+obj.MIOFF).path(ii+obj.MIOFF).derivative
                            scaleTemp = (obj.problem.phases(ph+obj.MIOFF).initialGuess.timeFinal-obj.problem.phases(ph+obj.MIOFF).initialGuess.timeInitial)/(Nodes-1)*...
                                vecnorm([dG_dX_ii(:,1:end-1);dG_dU_ii(:,1:end-1)],2,1);
                        elseif ~obj.problem.phases(ph+obj.MIOFF).path(ii+obj.MIOFF).integral && ~obj.problem.phases(ph+obj.MIOFF).path(ii+obj.MIOFF).derivative
                            scaleTemp = vecnorm([dG_dX_ii;dG_dU_ii],2,1);
                        end
                end
                scaleTemp(find((scaleTemp==0) | isnan(scaleTemp))) = 1;
                obj.problem.phases(ph+obj.MIOFF).path(ii+obj.MIOFF).scale = scaleTemp;
            end
        end
    case 'user'
    case 'none'
        for ph = 0:obj.problem.nphases-1
            for ii = 0:obj.problem.phases(ph+obj.MIOFF).n.path-1
                obj.problem.phases(ph+obj.MIOFF).path(ii+obj.MIOFF).scale = 1;
            end
        end
    otherwise
        error('SCOPT Error \nScaling Setting %s not defined \n Aborting \n',obj.algorithm.scaling.path)
        
end
%% Pre-scaling
% Prescales initial guess, tolerances, trust region, bounds and
% initial and final conditions
obj.algorithm.sequential.trustRegion.phases(1).states.component_n = (1./obj.problem.phases(1).scale.states).*obj.algorithm.sequential.trustRegion.phases(1).states.component;
obj.algorithm.sequential.trustRegion.phases(1).controls.component_n = (1./obj.problem.phases(1).scale.controls).*obj.algorithm.sequential.trustRegion.phases(1).controls.component;

if ~strcmp(obj.algorithm.sequential.trustRegion.phases(1).timeFinal.type,'none') && obj.problem.phases(1).freeTimeFinal
    obj.algorithm.sequential.trustRegion.phases(1).timeFinal.component_n = (1./obj.problem.phases(1).scale.time).*obj.algorithm.sequential.trustRegion.phases(1).timeFinal.component;
end

% Terminaion Settings

for ph = 0:obj.problem.nphases-1
    obj.algorithm.sequential.variablesTol.phases(ph+obj.MIOFF).states_n      = (1./obj.problem.phases(ph+obj.MIOFF).scale.states).*obj.algorithm.sequential.variablesTol.phases(ph+obj.MIOFF).states; %[-]
    obj.algorithm.sequential.variablesTol.phases(ph+obj.MIOFF).controls_n    = (1./obj.problem.phases(ph+obj.MIOFF).scale.controls).*obj.algorithm.sequential.variablesTol.phases(ph+obj.MIOFF).controls; %[-]
    obj.algorithm.sequential.variablesTol.phases(ph+obj.MIOFF).timeFinal_n   = (1./obj.problem.phases(ph+obj.MIOFF).scale.time).*obj.algorithm.sequential.variablesTol.phases(ph+obj.MIOFF).timeFinal; %[-]
    obj.algorithm.sequential.variablesTol.phases(ph+obj.MIOFF).timeInitial_n = (1./obj.problem.phases(ph+obj.MIOFF).scale.time).*obj.algorithm.sequential.variablesTol.phases(ph+obj.MIOFF).timeInitial; %[-]
end

% Variables Scaling
%   First routine if automatic scaling selected
%   Then, Routine to scale everything
obj.problem.phases(1).initialGuess.states_n           = obj.scaleStatesToNormalized(1,obj.problem.phases(1).initialGuess.states);
obj.problem.phases(1).initialGuess.controls_n         = obj.scaleControlsToNormalized(1,obj.problem.phases(1).initialGuess.controls);
obj.problem.phases(1).initialGuess.time_n             = obj.scaleTimeToNormalized(1,obj.problem.phases(1).initialGuess.time);
obj.problem.phases(1).initialGuess.timeFinal_n        = obj.scaleTimeToNormalized(1,obj.problem.phases(1).initialGuess.timeFinal);
obj.problem.phases(1).initialGuess.timeInitial_n      = obj.scaleTimeToNormalized(1,obj.problem.phases(1).initialGuess.timeInitial);


obj.problem.phases(1).bounds.states.upper_n           = obj.scaleStatesToNormalized(1,obj.problem.phases(1).bounds.states.upper(:));
obj.problem.phases(1).bounds.states.lower_n           = obj.scaleStatesToNormalized(1,obj.problem.phases(1).bounds.states.lower(:));
obj.problem.phases(1).bounds.controls.upper_n         = obj.scaleControlsToNormalized(1,obj.problem.phases(1).bounds.controls.upper);
obj.problem.phases(1).bounds.controls.lower_n         = obj.scaleControlsToNormalized(1,obj.problem.phases(1).bounds.controls.lower);
obj.problem.phases(1).bounds.timeFinal.upper_n        = obj.scaleTimeToNormalized(1,obj.problem.phases(1).bounds.timeFinal.upper);
obj.problem.phases(1).bounds.timeFinal.lower_n        = obj.scaleTimeToNormalized(1,obj.problem.phases(1).bounds.timeFinal.lower);
obj.problem.phases(1).bounds.timeInitial.upper_n      = obj.scaleTimeToNormalized(1,obj.problem.phases(1).bounds.timeInitial.upper);
obj.problem.phases(1).bounds.timeInitial.lower_n      = obj.scaleTimeToNormalized(1,obj.problem.phases(1).bounds.timeInitial.lower);


obj.problem.phases(1).initial.states.equal.value_n    = obj.scaleStatesToNormalized(1,obj.problem.phases(1).initial.states.equal.value,obj.problem.phases(1).initial.states.equal.index) ;
obj.problem.phases(1).initial.controls.equal.value_n  = obj.scaleControlsToNormalized(1,obj.problem.phases(1).initial.controls.equal.value,obj.problem.phases(1).initial.controls.equal.index) ;
obj.problem.phases(1).initial.states.upper.value_n    = obj.scaleStatesToNormalized(1,obj.problem.phases(1).initial.states.upper.value,obj.problem.phases(1).initial.states.upper.index) ;
obj.problem.phases(1).initial.controls.upper.value_n  = obj.scaleControlsToNormalized(1,obj.problem.phases(1).initial.controls.upper.value,obj.problem.phases(1).initial.controls.upper.index) ;
obj.problem.phases(1).initial.states.lower.value_n    = obj.scaleStatesToNormalized(1,obj.problem.phases(1).initial.states.lower.value,obj.problem.phases(1).initial.states.lower.index) ;
obj.problem.phases(1).initial.controls.lower.value_n  = obj.scaleControlsToNormalized(1,obj.problem.phases(1).initial.controls.lower.value,obj.problem.phases(1).initial.controls.lower.index) ;
obj.problem.phases(1).final.states.equal.value_n      = obj.scaleStatesToNormalized(1,obj.problem.phases(1).final.states.equal.value,obj.problem.phases(1).final.states.equal.index) ;
obj.problem.phases(1).final.controls.equal.value_n    = obj.scaleControlsToNormalized(1,obj.problem.phases(1).final.controls.equal.value,obj.problem.phases(1).final.controls.equal.index) ;
obj.problem.phases(1).final.states.upper.value_n      = obj.scaleStatesToNormalized(1,obj.problem.phases(1).final.states.upper.value,obj.problem.phases(1).final.states.upper.index) ;
obj.problem.phases(1).final.controls.upper.value_n    = obj.scaleControlsToNormalized(1,obj.problem.phases(1).final.controls.upper.value,obj.problem.phases(1).final.controls.upper.index) ;
obj.problem.phases(1).final.states.lower.value_n      = obj.scaleStatesToNormalized(1,obj.problem.phases(1).final.states.lower.value,obj.problem.phases(1).final.states.lower.index) ;
obj.problem.phases(1).final.controls.lower.value_n    = obj.scaleControlsToNormalized(1,obj.problem.phases(1).final.controls.lower.value,obj.problem.phases(1).final.controls.lower.index) ;


%% Sequential Algorithm Definition
if obj.algorithm.sequential.activate
    switch obj.algorithm.sequential.type
        case {'none','component-trust-region'}
            obj.algorithm.sequential.trustRegion.radius = 0;
            obj.algorithm.sequential.trustRegion.lambda = 0;
        case {'trust-region'}
            switch obj.algorithm.sequential.trustRegion.type
                case {'soft'}
                    obj.algorithm.sequential.trustRegion.radius = 0;
                    if isempty(obj.algorithm.sequential.trustRegion.lambda)
                        warning('SCOPT Warning, soft trust region weight not defined, using default')
                    end
                case {'hard'}
                    if isempty(obj.algorithm.sequential.trustRegion.radius)
                        auxArray = [];
                        if obj.algorithm.sequential.trustRegion.include.states
                            auxArray = [auxArray;obj.algorithm.sequential.trustRegion.phases(1).states.component_n];
                        end
                        if obj.algorithm.sequential.trustRegion.include.controls
                            auxArray = [auxArray;obj.algorithm.sequential.trustRegion.phases(1).controls.component_n];
                        end
                        
                        switch obj.algorithm.sequential.trustRegion.variablesPenalty
                            case 0
                                obj.algorithm.sequential.trustRegion.radius = sum(auxArray.*auxArray);
                            case 1
                                obj.algorithm.sequential.trustRegion.radius = norm(auxArray,1);
                            case 2
                                obj.algorithm.sequential.trustRegion.radius = norm(auxArray,2);
                            case inf
                                obj.algorithm.sequential.trustRegion.radius = norm(auxArray,inf);
                        end
                    else
                        if strcmp(obj.algorithm.scaling.trustRegion,'automatic')
                            % Adjust Radius with Scaling Term.
                            auxArray = [];
                            if obj.algorithm.sequential.trustRegion.include.states
                                auxArray = [auxArray;obj.problem.phases.scale.states];
                            end
                            if obj.algorithm.sequential.trustRegion.include.controls
                                auxArray = [auxArray;obj.problem.phases.scale.controls];
                            end
                            trustRegionPreScaling = mean(auxArray);
                        else
                            trustRegionPreScaling = 1;
                        end
                        
                        obj.algorithm.sequential.trustRegion.radius = obj.algorithm.sequential.trustRegion.radius/trustRegionPreScaling;
                    end
                    obj.algorithm.sequential.trustRegion.lambda = 0;
            end
        otherwise
            warning('Seqential obj.algorithm Not Defined \n')
    end
else
    obj.algorithm.sequential.trustRegion.radius = 0;
    obj.algorithm.sequential.trustRegion.lambda = 0;
end
for ph = 0:obj.problem.nphases-1
    if obj.problem.phases(ph+obj.MIOFF).freeTimeFinal
        switch obj.algorithm.sequential.trustRegion.phases(ph+obj.MIOFF).timeFinal.type
            case 'none'
            case 'soft'
            case 'hard'
                obj.algorithm.sequential.trustRegion.phases(ph+obj.MIOFF).timeFinal.radius = obj.algorithm.sequential.trustRegion.phases(ph+obj.MIOFF).timeFinal.component_n;
            otherwise
                warning('Free Time obj.algorithm Not Defined \n')
        end
    end
end

%% Other Virtual Control Settings
% Defines virtual control scaling term
for ph = 1:obj.problem.nphases
    if obj.algorithm.virtualControl.phases(ph).include % Is this right? How should this be scaled? (As it is basically changing with time in each iteration)
        if isnan(obj.algorithm.virtualControl.phases(ph).scale)
            obj.algorithm.virtualControl.phases(ph).scale = obj.problem.phases(ph).initialGuess.timeDisc;
        end
        
        obj.algorithm.virtualControl.phases(ph).rowE = obj.algorithm.virtualControl.phases(1).states(:);
        obj.algorithm.virtualControl.phases(ph).colE = (0:1:n.virtualControls-1)'+obj.MIOFF;
        obj.algorithm.virtualControl.phases(ph).valE = ones(n.virtualControls,1);
    else
        obj.algorithm.virtualControl.phases(ph).scale = 1;
    end
end
end

