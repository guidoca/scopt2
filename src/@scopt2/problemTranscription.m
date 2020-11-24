function problemTranscription(obj,timeInitial_n,timeFinal_n,X_n,U_n,radius,lambda)
%% Get Problem Fields
BodyMap             = obj.problem.phases(1).BodyMap;
Scale               = obj.problem.phases.scale;
Shift               = obj.problem.phases.shift;
Objective           = obj.problem.objective;
Events              = obj.problem.phases.events;
VirtualControl      = obj.algorithm.virtualControl;
VirtualBuffer       = obj.algorithm.virtualBuffer;
%% Retrieve Parameters from fields
% Previous solution (already processed for the current number of nodes, first one is the current state)
X         = obj.scaleStatesToReal(1,X_n);
U         = obj.scaleControlsToReal(1,U_n);
% Time and Nodes
timeFinal     = obj.scaleTimeToReal(1,timeFinal_n);
timeInitial   = obj.scaleTimeToReal(1,timeInitial_n);
Nodes         = obj.problem.phases(1).Nodes; % This sould account account for the initial node 0
% Convert to pseduo-time
tauDisc    = 2/(Nodes-1);
timeDisc   = tauDisc*(timeFinal - timeInitial)/2;
timeInterval= timeFinal-timeInitial;
time       = timeInitial:timeDisc:timeFinal;
%% Retrieve Auxiliary Structures
Ntot  = obj.problem.phases(1).Ntot; % Total number of variables
n     = obj.problem.phases(1).n;
index = obj.problem.phases(1).index;
%% Initial Conditions (Ineq & Eq)
% for ii = 1:obj.problem.nphases
Initial        = obj.problem.phases.initial;
% Equality
nEQ_states     = length(Initial.states.equal.index);
nEQ_controls   = length(Initial.controls.equal.index);
nEQ            = nEQ_states + nEQ_controls;
rowEQ_initial    = (0:nEQ-1)'+obj.MIOFF;
colEQ_initial    = [index.states(0+obj.MIOFF,1) + Initial.states.equal.index - 1;index.controls(0+obj.MIOFF,1) + Initial.controls.equal.index - 1];
valEQ_initial    = ones(nEQ,1);
MEQ_initial      = sparse(rowEQ_initial,colEQ_initial,valEQ_initial,nEQ,Ntot);
pEQ_initial      = [Initial.states.equal.value_n;Initial.controls.equal.value_n];
% Inequal
nINEQ_states_upp   = length(Initial.states.upper.index);
nINEQ_states_low   = length(Initial.states.lower.index);
nINEQ_controls_upp = length(Initial.controls.upper.index);
nINEQ_controls_low = length(Initial.controls.lower.index);
nINEQ              = nINEQ_states_upp + nINEQ_states_low + nINEQ_controls_upp + nINEQ_controls_low;
rowINEQ_initial      = (0:nINEQ-1)'+obj.MIOFF;
colINEQ_initial      = [index.states(0+obj.MIOFF,1) + Initial.states.upper.index - 1;index.states(0+obj.MIOFF,1) + Initial.states.lower.index - 1;...
    index.controls(0+obj.MIOFF,1)  + Initial.controls.upper.index - 1;index.controls(0+obj.MIOFF,1)  + Initial.controls.lower.index - 1];
valINEQ_initial      = [ones(nINEQ_states_upp,1);-ones(nINEQ_states_low,1);ones(nINEQ_controls_upp,1);-ones(nINEQ_controls_low,1)];
MINEQ_initial        = sparse(rowINEQ_initial,colINEQ_initial,valINEQ_initial,nINEQ,Ntot);
pINEQ_initial        = [Initial.states.upper.value_n;-Initial.states.lower.value_n;Initial.controls.upper.value_n;-Initial.controls.lower.value_n];
% end
obj.problem.phases(1).transcription.initial.MEQ   = MEQ_initial;
obj.problem.phases(1).transcription.initial.pEQ   = pEQ_initial;
obj.problem.phases(1).transcription.initial.MINEQ = MINEQ_initial;
obj.problem.phases(1).transcription.initial.pINEQ = pINEQ_initial;
%% Final Conditions  (Ineq & Eq)
Final              = obj.problem.phases.final;
% Equality
nEQ_states         = length(Final.states.equal.index);
nEQ_controls       = length(Final.controls.equal.index);
nEQ                = nEQ_states + nEQ_controls;
rowEQ_final        = (0:nEQ-1)'+obj.MIOFF;
colEQ_final        = [index.states(1,Nodes-1+obj.MIOFF) + Final.states.equal.index - 1;index.controls(0+obj.MIOFF,end)  + Final.controls.equal.index - 1];
valEQ_final        = ones(nEQ,1);
MEQ_final          = sparse(rowEQ_final,colEQ_final,valEQ_final,nEQ,Ntot);
pEQ_final          = [Final.states.equal.value_n;Final.controls.equal.value_n];
% Linear Inequal
nINEQ_states_upp   = length(Final.states.upper.index);
nINEQ_states_low   = length(Final.states.lower.index);
nINEQ_controls_upp = length(Final.controls.upper.index);
nINEQ_controls_low = length(Final.controls.lower.index);
nINEQ              = nINEQ_states_upp + nINEQ_states_low + nINEQ_controls_upp + nINEQ_controls_low;
rowINEQ_final      = (0:nINEQ-1)'+obj.MIOFF;
colINEQ_final      = [index.states(1,Nodes-1+obj.MIOFF) + Final.states.upper.index - 1;index.states(1,Nodes-1+obj.MIOFF) + Final.states.lower.index - 1;...
    index.controls(1,Nodes-1+obj.MIOFF) + Final.controls.upper.index - 1;index.controls(1,Nodes-1+obj.MIOFF) + Final.controls.lower.index - 1];
valINEQ_final      = [ones(nINEQ_states_upp,1);-ones(nINEQ_states_low,1);ones(nINEQ_controls_upp,1);-ones(nINEQ_controls_low,1)];
MINEQ_final        = sparse(rowINEQ_final,colINEQ_final,valINEQ_final,nINEQ,Ntot);
pINEQ_final        = [Final.states.upper.value_n;-Final.states.lower.value_n;Final.controls.upper.value_n;-Final.controls.lower.value_n];
% end
obj.problem.phases(1).transcription.final.MEQ = MEQ_final;
obj.problem.phases(1).transcription.final.pEQ = pEQ_final;
obj.problem.phases(1).transcription.final.MINEQ = MINEQ_final;
obj.problem.phases(1).transcription.final.pINEQ = pINEQ_final;
%% State and Control Constraints (Ineq)
Bounds            = obj.problem.phases.bounds;
if obj.problem.phases(1).freeTimeFinal && ~isempty(Bounds.timeFinal.upper_n) && ~isempty(Bounds.timeFinal.lower_n)
    % Free time of flight
    MINEQ_time = sparse((0:1)+obj.MIOFF,[index.timeFinal,index.timeFinal],[1 ,-1],2,Ntot);
    pINEQ_time = [Bounds.timeFinal.upper_n;-Bounds.timeFinal.lower_n];
else
    MINEQ_time = sparse(0,Ntot);
    pINEQ_time = zeros(0,1);
end
if obj.problem.phases(1).freeTimeInitial && ~isempty(Bounds.timeInitial.upper_n) && ~isempty(Bounds.timeInitial.lower_n)
    % Free time of flight
    MINEQ_time = [MINEQ_time ; sparse((0:1)+obj.MIOFF,[index.timeInitial,index.timeInitial],[1 ,-1],2,Ntot)];
    pINEQ_time = [pINEQ_time ; [Bounds.timeInitial.upper_n;-Bounds.timeInitial.lower_n]];
    error('SCOPT Error \n Free Initial Time not Implemented')
end
% if states upper and lower bounds present, create sparse matrix
if ~isempty(Bounds.states.upper_n) && ~isempty(Bounds.states.lower_n)
    rowINEQ_state     = reshape(repmat((0:(2*n.states):((2*n.states)*(Nodes-1))),2*n.states,1) + repmat((0:n.states-1)',2,Nodes) + repmat([zeros(n.states,1);n.states*ones(n.states,1)],1,Nodes) +obj.MIOFF,[],1);
    colINEQ_state     = reshape(repmat(index.states,2,1),[],1);
    valINEQ_state     = repmat([1*ones(n.states,1);-1*ones(n.states,1)],Nodes,1);
    MINEQ_states       = sparse(rowINEQ_state,colINEQ_state,valINEQ_state,2*Nodes*n.states,Ntot);
    pINEQ_states       = repmat([Bounds.states.upper_n;-Bounds.states.lower_n],Nodes,1);
else
    MINEQ_states       = sparse(0,Ntot);
    pINEQ_states       = zeros(0,1);
end
% if controls upper and lower bounds present, create sparse matrix
if ~isempty(Bounds.controls.upper_n) && ~isempty(Bounds.controls.lower_n)
    rowINEQ_control   = reshape(repmat((0:(2*n.controls):((2*n.controls)*(Nodes-1))),2*n.controls,1) + repmat((0:n.controls-1)',2,Nodes) + repmat([zeros(n.controls,1);n.controls*ones(n.controls,1)],1,Nodes) +obj.MIOFF,[],1);
    colINEQ_control   = reshape(repmat(index.controls,2,1),[],1);
    valINEQ_control   = repmat([1*ones(n.controls,1);-1*ones(n.controls,1)],Nodes,1);
    MINEQ_controls    = sparse(rowINEQ_control,colINEQ_control,valINEQ_control,2*Nodes*n.controls,Ntot);
    pINEQ_controls    = repmat([Bounds.controls.upper_n;-Bounds.controls.lower_n],Nodes,1);
else
    MINEQ_controls       = sparse(0,Ntot);
    pINEQ_controls       = zeros(0,1);
end
obj.problem.phases(1).transcription.bounds.states.MINEQ    = MINEQ_states;
obj.problem.phases(1).transcription.bounds.states.pINEQ    = pINEQ_states;
obj.problem.phases(1).transcription.bounds.controls.MINEQ  = MINEQ_controls;
obj.problem.phases(1).transcription.bounds.controls.pINEQ  = pINEQ_controls;
%% Additional Time Constraints
if obj.problem.phases(1).freeTimeFinal && ~obj.problem.phases(1).freeTimeInitial
    MINEQ_time = [MINEQ_time ; sparse(0+obj.MIOFF,[index.timeFinal],-1,1,Ntot)];
    pINEQ_time = [pINEQ_time ; -timeInitial_n];
elseif ~obj.problem.phases(1).freeTimeFinal && obj.problem.phases(1).freeTimeInitial
    MINEQ_time = [MINEQ_time ; sparse(0+obj.MIOFF,[index.timeInitial],1,1,Ntot)];
    pINEQ_time = [pINEQ_time ; timeFinal_n];
elseif obj.problem.phases(1).freeTimeFinal && obj.problem.phases(1).freeTimeInitial
    MINEQ_time = [MINEQ_time ; sparse([0,0]+obj.MIOFF,[index.timeFinal,index.timeInitial],[1 -1],1,Ntot)];
    pINEQ_time = [pINEQ_time ; 0];
end
%% Events
% Sets event constraints
if n.events
    pEQ_events      = [];
    rowEQ_events    = [];
    colEQ_events    = [];
    valEQ_events    = [];
    pINEQ_events    = [];
    rowINEQ_events  = [];
    colINEQ_events  = [];
    valINEQ_events  = [];
    rowSOC_events   = [];
    colSOC_events   = [];
    valSOC_events   = [];
    qSOC_events     = [];
    hSOC_events     = [];
    rowINEQ_Initial = 0;
    rowEQ_Initial   = 0;
    for ii = 0:n.events-1
        switch Events(ii+obj.MIOFF).where
            case 'final'
                indexWhere = Nodes - 1 +obj.MIOFF;
            case 'initial'
                indexWhere = 0 +obj.MIOFF;
            case 'index'
                indexWhere = Events(ii+obj.MIOFF).indexWhere;
                warning('SCOPT Error, Custom location for Event constraints not implemented \n Continuing without constraint \n')
                %                             continue
            otherwise
                error('SCOPT Warning, Event constraint at %s not implemented \n Continuing without constraint %i \n',Events(ii+obj.MIOFF).where,ii)
                %                             continue
        end
        scaleEvents    = Events(ii+obj.MIOFF).scale;
        limit          = Events(ii+obj.MIOFF).limit;
        switch Events(ii+obj.MIOFF).funType
            case 'linear'
                % Linear Event Constraints
                if strcmp(Events(ii+obj.MIOFF).type,'equal')
                    rowEQ_events     =  [rowEQ_events ; (zeros(n.states+n.controls,1) + rowEQ_Initial +obj.MIOFF)];
                    colEQ_events     =  [colEQ_events ; index.states(:,indexWhere);index.controls(:,indexWhere)];
                    valEQ_events     =  [valEQ_events ; 1/scaleEvents*Scale.states.*Events(ii+obj.MIOFF).states;1/scaleEvents*Scale.controls.*Events(ii+obj.MIOFF).controls];
                    pEQ_events       =  [pEQ_events   ; 1/scaleEvents*(limit - sum(Scale.states.*Events(ii+obj.MIOFF).cons,1) - sum(Shift.states.*Events(ii+obj.MIOFF).states,1)-  sum(Shift.controls.*Events(ii+obj.MIOFF).controls,1))];
                    
                else
                    switch Events(ii+obj.MIOFF).type % Does not affect Convex Constraint (can only have an upper constraint by definition)
                        case 'lower'
                            typeSign = -1;
                        case 'upper'
                            typeSign = 1;
                        otherwise
                            error('SCOPT Error, Events constraints of type %s not implemented \n Aborting',Events(ii+obj.MIOFF).type)
                    end
                    rowINEQ_events     =  [rowINEQ_events ; (zeros(n.states+n.controls,1) + rowINEQ_Initial +obj.MIOFF)];
                    colINEQ_events     =  [colINEQ_events ; index.states(:,indexWhere);index.controls(:,indexWhere)];
                    valINEQ_events     =  [valINEQ_events ; typeSign/scaleEvents*Scale.states.*Events(ii+obj.MIOFF).states;typeSign/scaleEvents*Scale.controls.*Events(ii+obj.MIOFF).controls];
                    pINEQ_events       =  [pINEQ_events   ; typeSign/scaleEvents*(limit - sum(Scale.states.*Events(ii+obj.MIOFF).cons,1) - sum(Shift.states.*Events(ii+obj.MIOFF).states,1)-  sum(Shift.controls.*Events(ii+obj.MIOFF).controls,1))];
                end
            case 'non-linear'
                % Event constraint defined by a nonlinear
                % function
                switch Events(ii+obj.MIOFF).type
                    case 'lower'
                        typeSign = -1;
                    case 'upper'
                        typeSign = 1;
                    otherwise
                        error('SCOPT Error, Events constraints of type %s not implemented \n Aborting',Events(ii+obj.MIOFF).type)
                        
                end
                f_n_ii         = Events(ii+obj.MIOFF).function(time(:,indexWhere),X(:,indexWhere),U(:,indexWhere),BodyMap,indexWhere)/scaleEvents;
                dG_dX_ii       = Events(ii+obj.MIOFF).jacobian.states(time(:,indexWhere),X(:,indexWhere),U(:,indexWhere),BodyMap,indexWhere);
                dG_dU_ii       = Events(ii+obj.MIOFF).jacobian.controls(time(:,indexWhere),X(:,indexWhere),U(:,indexWhere),BodyMap,indexWhere);
                dG_n_dX_n_ii   = Scale.states.*dG_dX_ii/scaleEvents;
                dG_n_dU_n_ii   = Scale.controls.*dG_dU_ii/scaleEvents;
                
                g_n_function      =  f_n_ii - limit/scaleEvents;
                dG_n_dX_n_X_n     =  dot(dG_n_dX_n_ii(:),X_n(:,indexWhere));
                dG_n_dU_n_U_n     =  dot(dG_n_dU_n_ii(:),U_n(:,indexWhere));
                
                rowINEQ_events     =  [rowINEQ_events ; (zeros(n.states+n.controls,1) + rowINEQ_Initial +obj.MIOFF)];
                colINEQ_events     =  [colINEQ_events ; index.states(:,indexWhere);index.controls(:,indexWhere)];
                valINEQ_events     =  [valINEQ_events ; typeSign*dG_n_dX_n_ii;typeSign*dG_n_dU_n_ii];
                
                if Events(ii+obj.MIOFF).buffer.include && VirtualBuffer.phases(1).include
                    rowINEQ_events  = [rowINEQ_events ; (0+ rowINEQ_Initial + obj.MIOFF)];
                    colINEQ_events  = [colINEQ_events ; index.virtualBuffers.events(ii+obj.MIOFF)];
                    valINEQ_events  = [valINEQ_events ; -1];
                end
                
                pINEQ_events        = [pINEQ_events ; - typeSign*(g_n_function - dG_n_dX_n_X_n - dG_n_dU_n_U_n)];
                
            case 'quasiconvex'
                % Event constraint defined with a nonlinear
                % quasiconvex function
                switch Events(ii+obj.MIOFF).type
                    case 'lower'
                        typeSign = -1;
                    case 'upper'
                        typeSign = 1;
                    otherwise
                        error('SCOPT Error, Event Constraints of constraint type %s not implemented. Use inequality constraints instead \n Aborting',Events(ii+obj.MIOFF).type)
                        
                end
                if (length(Events(ii+obj.MIOFF).statesIndex) + length(Events(ii+obj.MIOFF).controlsIndex) ==1)
                    % Single variable quasiconvex case
                    f_ii          = typeSign*Events(ii+obj.MIOFF).function(time(:,indexWhere),X(:,indexWhere),U(:,indexWhere),BodyMap,indexWhere);
                    if length(Events(ii+obj.MIOFF).statesIndex)==1 % variable is a state
                        df_dvar       = typeSign*Events(ii+obj.MIOFF).jacobian.states(time(:,indexWhere),X(:,indexWhere),U(:,indexWhere),BodyMap,indexWhere);
                        d2f_dvar2     = typeSign*Events(ii+obj.MIOFF).hessian.states(time(:,indexWhere),X(:,indexWhere),U(:,indexWhere),BodyMap,indexWhere);
                        var           = X(Events(ii+obj.MIOFF).variableIndex,indexWhere);
                        scaleVar      = Scale.states(Events(ii+obj.MIOFF).variableIndex);
                        shiftVar      = Shift.states(Events(ii+obj.MIOFF).variableIndex);
                        indexVector   = index.states(Events(ii+obj.MIOFF).variableIndex,indexWhere);
                    else % variable is a control
                        df_dvar       = typeSign*Events(ii+obj.MIOFF).jacobian.controls(time(:,indexWhere),X(:,indexWhere),U(:,indexWhere),BodyMap,indexWhere);
                        d2f_dvar2     = typeSign*Events(ii+obj.MIOFF).hessian.controls(time(:,indexWhere),X(:,indexWhere),U(:,indexWhere),BodyMap,indexWhere);
                        var           = U(Events(ii+obj.MIOFF).variableIndex,indexWhere);
                        scaleVar      = Scale.controls(Events(ii+obj.MIOFF).variableIndex);
                        shiftVar      = Shift.controls(Events(ii+obj.MIOFF).variableIndex);
                        indexVector   = index.controls(Events(ii+obj.MIOFF).variableIndex,indexWhere);
                    end
                    
                    LT  = sqrt(d2f_dvar2/2);
                    qT  = (df_dvar - var.*d2f_dvar2);
                    r   = f_ii - df_dvar.*var + 1/2*d2f_dvar2.*var.^2;
                    Ai  = [qT/2;LT];
                    bi  = [(1+ r)/2;0];
                    ci  = -1/2*qT;
                    di  = (1- r)/2;
                    
                    qSOC_cvx_ii      = 3;
                    qSOC_cvx         = qSOC_cvx_ii*ones(1,1);
                    rowSOC_events    = [rowSOC_events;(sum(qSOC_events) + reshape(repmat((0:qSOC_cvx_ii-1),1+n.states+n.controls,1),[],1) + obj.MIOFF)];
                    colSOC_events    = [colSOC_events;reshape(repmat([indexVector;index.states(:,indexWhere);index.controls(:,indexWhere)],qSOC_cvx_ii,1),[],1)];
                    valSOC_events    = [valSOC_events;reshape(-[scaleVar*ci;[typeSign*(-1/2)*Scale.states(:).*Events(ii+obj.MIOFF).states(:);typeSign*(-1/2)*Scale.controls(:).*Events(ii+obj.MIOFF).controls(:)];...
                        scaleVar*Ai(1);[typeSign*Scale.states(:).*Events(ii+obj.MIOFF).states(:)/2;typeSign*Scale.controls(:).*Events(ii+obj.MIOFF).controls(:)/2];scaleVar*Ai(2);zeros(n.states+n.controls,1)]./scaleEvents,[],1)];
                    
                    diShift          = -qT/2.*shiftVar - 1/2*sum(typeSign*Shift.states(:).*Events(ii+obj.MIOFF).states(:))- 1/2*sum(typeSign*Shift.controls(:).*Events(ii+obj.MIOFF).controls(:));
                    biShift          = [qT/2.*shiftVar + 1/2*sum(typeSign*Shift.states(:).*Events(ii+obj.MIOFF).states(:))+ 1/2*sum(typeSign*Shift.controls(:).*Events(ii+obj.MIOFF).controls(:));LT.*shiftVar] ;
                    
                    qSOC_events        = [qSOC_events,qSOC_cvx];
                    hSOC_events        = [hSOC_events;reshape([di + diShift;bi + biShift]./scaleEvents,[],1)];
                else % general case
                    error('SCOPT Error: Nonlinear General Quasiconvex Event Constraint not implemented yet \n')
                end
                
            case 'non-linear-1var'
                % Event constraint defined as quasiconvex with
                % a second derivative term, equivalent to the
                % quasiconvex with just one derivative variable
                switch Events(ii+obj.MIOFF).type
                    case 'lower'
                        typeSign = -1;
                    case 'upper'
                        typeSign = 1;
                    otherwise
                        error('SCOPT Error, Event Constraints of constraint type %s not implemented. Use inequality constraints instead \n Aborting',Events(ii+obj.MIOFF).type)
                        
                end
                f_ii          = typeSign*Events(ii+obj.MIOFF).function(time(:,indexWhere),X(:,indexWhere),U(:,indexWhere),BodyMap,indexWhere);
                df_dvar       = typeSign*Events(ii+obj.MIOFF).jacobian.variable(time(:,indexWhere),X(:,indexWhere),U(:,indexWhere),BodyMap,indexWhere);
                d2f_dvar2     = typeSign*Events(ii+obj.MIOFF).hessian.variable(time(:,indexWhere),X(:,indexWhere),U(:,indexWhere),BodyMap,indexWhere);
                switch Events(ii+obj.MIOFF).variableType
                    case 'state'
                        var      = X(Events(ii+obj.MIOFF).variableIndex,indexWhere);
                        scaleVar = Scale.states(Events(ii+obj.MIOFF).variableIndex);
                        shiftVar = Shift.states(Events(ii+obj.MIOFF).variableIndex);
                        indexVector = index.states(Events(ii+obj.MIOFF).variableIndex,indexWhere);
                    case 'control'
                        var      = U(Events(ii+obj.MIOFF).variableIndex,indexWhere);
                        scaleVar = Scale.controls(Events(ii+obj.MIOFF).variableIndex);
                        shiftVar = Shift.controls(Events(ii+obj.MIOFF).variableIndex);
                        indexVector = index.controls(Events(ii+obj.MIOFF).variableIndex,indexWhere);
                    case 'undefined'
                        error('Variable type not inserted')
                end
                switch Events(ii+obj.MIOFF).order
                    case 1
                    case 2
                        LT  = sqrt(d2f_dvar2/2);
                        qT  = (df_dvar - var.*d2f_dvar2);
                        r   = f_ii - df_dvar.*var + 1/2*d2f_dvar2.*var.^2;
                        Ai  = [qT/2;LT];
                        bi  = [(1+ r)/2;0];
                        ci  = -1/2*qT;
                        di  = (1- r)/2;
                        
                        qSOC_cvx_ii      = 3;
                        qSOC_cvx         = qSOC_cvx_ii*ones(1,1);
                        rowSOC_events    = [rowSOC_events;(sum(qSOC_events) + reshape(repmat((0:qSOC_cvx_ii-1),1+n.states+n.controls,1),[],1) + obj.MIOFF)];
                        colSOC_events    = [colSOC_events;reshape(repmat([indexVector;index.states(:,indexWhere);index.controls(:,indexWhere)],qSOC_cvx_ii,1),[],1)];
                        valSOC_events    = [valSOC_events;reshape(-[scaleVar*ci;[typeSign*(-1/2)*Scale.states(:).*Events(ii+obj.MIOFF).states(:);typeSign*(-1/2)*Scale.controls(:).*Events(ii+obj.MIOFF).controls(:)];...
                            scaleVar*Ai(1);[typeSign*Scale.states(:).*Events(ii+obj.MIOFF).states(:)/2;typeSign*Scale.controls(:).*Events(ii+obj.MIOFF).controls(:)/2];scaleVar*Ai(2);zeros(n.states+n.controls,1)]./scaleEvents,[],1)];
                        
                        diShift          = -qT/2.*shiftVar - 1/2*sum(typeSign*Shift.states(:).*Events(ii+obj.MIOFF).states(:))- 1/2*sum(typeSign*Shift.controls(:).*Events(ii+obj.MIOFF).controls(:));
                        biShift          = [qT/2.*shiftVar + 1/2*sum(typeSign*Shift.states(:).*Events(ii+obj.MIOFF).states(:))+ 1/2*sum(typeSign*Shift.controls(:).*Events(ii+obj.MIOFF).controls(:));LT.*shiftVar] ;
                        
                        qSOC_events        = [qSOC_events,qSOC_cvx];
                        hSOC_events        = [hSOC_events;reshape([di + diShift;bi + biShift]./scaleEvents,[],1)];
                    otherwise
                        error('SCOPT Error: Order not defined')
                end
            case {'convex','soc'}
                % Event constraint defined by a second ordeer
                % cone
                switch Events(ii+obj.MIOFF).type
                    case 'lower'
                        typeSign = -1;
                    case 'upper'
                        typeSign = 1;
                    otherwise
                        error('SCOPT Error, Event Constraints of constraint type %s not implemented. Use inequality constraints instead \n Aborting',Events(ii+obj.MIOFF).type)
                        
                end
                qSOC_cvx_ii      = Events(ii+obj.MIOFF).cone.dimensions;
                rowSOC_events    = [rowSOC_events; (sum(qSOC_events) + zeros(qSOC_cvx_ii*(n.states+n.controls),1) + reshape(repmat((0:qSOC_cvx_ii-1),n.states+n.controls,1),[],1) +obj.MIOFF)];
                colSOC_events    = [colSOC_events;reshape(repmat([index.states(:,indexWhere);index.controls(:,indexWhere)],qSOC_cvx_ii,1),[],1)];
                valSOC_events    = [valSOC_events;-[typeSign*Scale.states(:).*Events(ii+obj.MIOFF).cone.right.states(:);typeSign*Scale.controls(:).*Events(ii+obj.MIOFF).cone.right.controls(:);reshape([repmat(Scale.states(:),1,qSOC_cvx_ii-1).*Events(ii+obj.MIOFF).cone.norm.states;repmat(Scale.controls(:),1,qSOC_cvx_ii-1).*Events(ii+obj.MIOFF).cone.norm.controls],[],1)]/scaleEvents];
                
                qSOC_events      = [qSOC_events,qSOC_cvx_ii];
                hSOC_events      = [hSOC_events;[typeSign*(limit + Events(ii+obj.MIOFF).cone.right.cons + sum(Shift.states(:).*Events(ii+obj.MIOFF).cone.right.states,1)+ sum(Shift.controls(:).*Events(ii+obj.MIOFF).cone.right.controls,1)) ; (Events(ii+obj.MIOFF).cone.norm.cons(:) + reshape(sum(repmat(Shift.states(:),1,qSOC_cvx_ii-1).*Events(ii+obj.MIOFF).cone.norm.states,1) + sum(repmat(Shift.controls(:),1,qSOC_cvx_ii-1).*Events(ii+obj.MIOFF).cone.norm.controls,1),[],1)) ]/scaleEvents];
            otherwise
                error('SCOPT Error, Event Constraint of function type %s not implemented \n Aborting',Events(ii+obj.MIOFF).funType)
                
        end
        if ~isempty(rowINEQ_events)
            rowINEQ_Initial= max(rowINEQ_events);
        end
        if ~isempty(rowEQ_events)
            rowEQ_Initial= max(rowEQ_events);
        end
    end
    if isempty(rowINEQ_Initial)
        rowINEQ_Initial = 0;
    end
    if isempty(rowEQ_Initial)
        rowEQ_Initial = 0;
    end
    GSOC_events  = sparse(rowSOC_events,colSOC_events,valSOC_events,sum(qSOC_events),Ntot);
    MINEQ_events = sparse(rowINEQ_events,colINEQ_events,valINEQ_events,rowINEQ_Initial,Ntot);
    MEQ_events   = sparse(rowEQ_events,colEQ_events,valEQ_events,rowEQ_Initial,Ntot);
else
    GSOC_events  = sparse(0,Ntot);
    hSOC_events  = zeros(0,1);
    qSOC_events  = zeros(1,0);
    MINEQ_events = sparse(0,Ntot);
    pINEQ_events = zeros(0,1);
    MEQ_events   = sparse(0,Ntot);
    pEQ_events   = zeros(0,1);
end
obj.problem.phases(1).transcription.events.MINEQ = MINEQ_events;
obj.problem.phases(1).transcription.events.pINEQ = pINEQ_events;
obj.problem.phases(1).transcription.events.MEQ   = MEQ_events;
obj.problem.phases(1).transcription.events.pEQ   = pEQ_events;
obj.problem.phases(1).transcription.events.GSOC  = GSOC_events;
obj.problem.phases(1).transcription.events.hSOC  = hSOC_events;
obj.problem.phases(1).transcription.events.qSOC  = qSOC_events;
%% Path Constraints
% Sets Path Constraints
if n.path
    Path         = obj.problem.phases.path;
    rowEQ_path   = [];
    colEQ_path   = [];
    valEQ_path   = [];
    pEQ_path     = [];
    rowINEQ_path = [];
    colINEQ_path = [];
    valINEQ_path = [];
    pINEQ_path   = [];
    rowSOC_path  = [];
    colSOC_path  = [];
    valSOC_path  = [];
    hSOC_path    = [];
    qSOC_path    = [];
    jj = 0;
    rowINEQ_Initial = 0;
    rowEQ_Initial   = 0;
    for ii = 0:n.path-1
        scalePath   = Path(ii+obj.MIOFF).scale;
        if length(scalePath)==1
            Naux = Nodes;
        else
            Naux = 1;
        end
        limit       = Path(ii+obj.MIOFF).limit;
        switch Path(ii+obj.MIOFF).funType
            case 'linear'
                % Linear Path Constraints, may also be equal
                if strcmp(Path(ii+obj.MIOFF).type,'equal')
                    rowEQ_path     = [rowEQ_path ; (reshape(repmat((0:(Nodes-1)),n.states + n.controls,1),[],1) + rowEQ_Initial +obj.MIOFF) ];
                    colEQ_path     = [colEQ_path ; reshape([index.states;index.controls],[],1)];
                    valEQ_path     = [valEQ_path ; 1/scalePath*repmat([Scale.states.*Path(ii+obj.MIOFF).states;Scale.controls.*Path(ii+obj.MIOFF).controls],Nodes,1)];
                    pEQ_path       = [pEQ_path   ; (limit - Path(ii+obj.MIOFF).cons - sum(Shift.states.*Path(ii+obj.MIOFF).states,1)-  sum(Shift.controls.*Path(ii+obj.MIOFF).controls,1))/scalePath*ones(Nodes,1)];
                else
                    switch Path(ii+obj.MIOFF).type
                        case 'lower'
                            typeSign = -1;
                        case 'upper'
                            typeSign = 1;
                        otherwise
                            error('Path Constraint type not inserted')
                    end
                    if Path(ii+obj.MIOFF).derivative && ~Path(ii+obj.MIOFF).integral
                        % Derivative Linear Constraint
                        rowINEQ_path    = [rowINEQ_path; (reshape(repmat((0:(Nodes-2)),2*(n.states + n.controls),1),[],1) + rowINEQ_Initial +obj.MIOFF)];
                        colINEQ_path    = [colINEQ_path; reshape([index.states(:,0+obj.MIOFF:Nodes-2+obj.MIOFF);index.controls(:,0+obj.MIOFF:Nodes-2+obj.MIOFF);index.states(:,1+obj.MIOFF:Nodes-1+obj.MIOFF);index.controls(:,1+obj.MIOFF:Nodes-1+obj.MIOFF)],[],1)];
                        valINEQ_path    = [valINEQ_path; typeSign*reshape(repmat([-Scale.states.*Path(ii+obj.MIOFF).states;-Scale.controls.*Path(ii+obj.MIOFF).controls;Scale.states.*Path(ii+obj.MIOFF).states;Scale.controls.*Path(ii+obj.MIOFF).controls],1,Nodes-1),[],1)/scalePath];
                        
                        if obj.problem.phases(1).freeTimeFinal
                            dg_dtf_n      = tauDisc/2*Scale.time*limit;
                            
                            rowINEQ_path  = [rowINEQ_path; ((0:(Nodes-2))' + rowINEQ_Initial +obj.MIOFF)];
                            colINEQ_path  = [colINEQ_path; index.timeFinal*ones(Nodes-1,1)];
                            valINEQ_path  = [valINEQ_path; -typeSign*dg_dtf_n/scalePath*ones(Nodes-1,1)];
                            pINEQ_path    = [pINEQ_path ; typeSign*ones(Nodes-1,1)*(timeDisc*limit - dg_dtf_n*timeFinal_n)/scalePath];
                        else
                            pINEQ_path    = [pINEQ_path ; typeSign*ones(Nodes-1,1)*(timeDisc*limit)/scalePath];
                        end
                    elseif Path(ii+obj.MIOFF).integral && ~Path(ii+obj.MIOFF).derivative
                        % Integral Linear Constraint
                        rowINEQ_path    = [rowINEQ_path; (zeros(Nodes*(n.states + n.controls),1) + rowINEQ_Initial +obj.MIOFF)];
                        colINEQ_path    = [colINEQ_path; reshape([index.states;index.controls],[],1)];
                        valINEQ_path    = [valINEQ_path; typeSign*timeDisc/scalePath*reshape(repmat([Scale.states.*Path(ii+obj.MIOFF).states;Scale.controls.*Path(ii+obj.MIOFF).controls],1,Nodes),[],1)];
                        
                        if obj.problem.phases(1).freeTimeFinal
                            f_ii        = [Path(ii+obj.MIOFF).states;Path(ii+obj.MIOFF).controls]'*[X;U];
                            dI_dtf_n    = tauDisc/2*Scale.time*sum(f_ii);
                            
                            rowINEQ_path    = [rowINEQ_path; (0 + rowINEQ_Initial +obj.MIOFF)];
                            colINEQ_path    = [colINEQ_path; index.timeFinal];
                            valINEQ_path    = [valINEQ_path; typeSign*dI_dtf_n/scalePath];
                            pINEQ_path      = [pINEQ_path ; typeSign*(limit + dI_dtf_n*timeFinal_n + timeDisc*(- sum(Shift.states.*Path(ii+obj.MIOFF).states,1) -  sum(Shift.controls.*Path(ii+obj.MIOFF).controls,1)))/scalePath];
                        else
                            pINEQ_path      = [pINEQ_path ; typeSign*(limit + timeDisc*(- sum(Shift.states.*Path(ii+obj.MIOFF).states,1) -  sum(Shift.controls.*Path(ii+obj.MIOFF).controls,1)))/scalePath];
                        end
                    elseif ~Path(ii+obj.MIOFF).integral && ~Path(ii+obj.MIOFF).derivative
                        % Linear Constraint
                        rowINEQ_path     = [rowINEQ_path ; (reshape(repmat((0:(Nodes-1)),n.states + n.controls,1),[],1) + rowINEQ_Initial +obj.MIOFF) ];
                        colINEQ_path     = [colINEQ_path ; reshape([index.states;index.controls],[],1)];
                        valINEQ_path     = [valINEQ_path ; typeSign/scalePath*repmat([Scale.states.*Path(ii+obj.MIOFF).states;Scale.controls.*Path(ii+obj.MIOFF).controls],Nodes,1)];
                        pINEQ_path       = [pINEQ_path   ; typeSign*(limit - Path(ii+obj.MIOFF).cons - sum(Shift.states.*Path(ii+obj.MIOFF).states,1)-  sum(Shift.controls.*Path(ii+obj.MIOFF).controls,1))/scalePath*ones(Nodes,1)];
                    else
                        error('SCOPT Error, Path Constraint is not defined for linear integral and derivative terms \n Aborting')
                        
                    end
                end
            case 'non-linear'
                % Nonlinear path constraints
                switch Path(ii+obj.MIOFF).type
                    case 'lower'
                        typeSign = -1;
                    case 'upper'
                        typeSign = 1;
                    otherwise
                        error('SCOPT Error, Path Constraints of constraint type %s not implemented. Use inequality constraints instead \n Aborting',Path(ii+obj.MIOFF).type)
                        
                end
                f_ii          = Path(ii+obj.MIOFF).function(time,X,U,BodyMap,index.nodes);
                dG_dX_n_ii    = repmat(Scale.states,1,Nodes).*Path(ii+obj.MIOFF).jacobian.states(time,X,U,BodyMap,index.nodes);
                dG_dU_n_ii    = repmat(Scale.controls,1,Nodes).*Path(ii+obj.MIOFF).jacobian.controls(time,X,U,BodyMap,index.nodes);
                
                if Path(ii+obj.MIOFF).integral && ~Path(ii+obj.MIOFF).derivative
                    % Integral nonlinear path constraint
                    I_n          = timeInterval*sum(obj.problem.phases(1).quadratureVector.*f_ii)/scalePath;
                    dI_n_dX_n    = timeInterval*repmat(obj.problem.phases(1).quadratureVector,n.states,1).*dG_dX_n_ii/scalePath;
                    dI_n_dU_n    = timeInterval*repmat(obj.problem.phases(1).quadratureVector,n.controls,1).*dG_dU_n_ii/scalePath;
                    
                    rowINEQ_path    = [rowINEQ_path; (zeros(Nodes*(n.states + n.controls),1) + rowINEQ_Initial +obj.MIOFF)];
                    colINEQ_path    = [colINEQ_path; reshape([index.states;index.controls],[],1)];
                    valINEQ_path    = [valINEQ_path; reshape([typeSign*dI_n_dX_n;typeSign*dI_n_dU_n],[],1)];
                    
                    if obj.problem.phases(1).freeTimeFinal
                        dI_n_dtf_n   = Scale.time*sum(obj.problem.phases(1).quadratureVector.*f_ii)/scalePath;
                        rowINEQ_path = [rowINEQ_path; (0 + rowINEQ_Initial +obj.MIOFF)];
                        colINEQ_path = [colINEQ_path; index.timeFinal];
                        valINEQ_path = [valINEQ_path; typeSign*dI_n_dtf_n];
                        
                        pINEQ_path      = [pINEQ_path ; typeSign*(limit/scalePath - I_n + dI_n_dtf_n*timeFinal_n + dot(dI_n_dX_n(:),X_n(:)) + dot(dI_n_dU_n(:),U_n(:)))];
                    else
                        pINEQ_path      = [pINEQ_path ; typeSign*(limit/scalePath - I_n + dot(dI_n_dX_n(:),X_n(:)) + dot(dI_n_dU_n(:),U_n(:)))];
                    end
                    
                    if Path(ii+obj.MIOFF).buffer.include && VirtualBuffer.phases(1).include
                        rowINEQ_path  = [rowINEQ_path ; (0 + rowINEQ_Initial +obj.MIOFF)];
                        colINEQ_path  = [colINEQ_path ; index.virtualBuffers.integralPath];
                        valINEQ_path  = [valINEQ_path ; -1];
                        jj = jj +1;
                    end
                elseif ~Path(ii+obj.MIOFF).integral && Path(ii+obj.MIOFF).derivative
                    % NL derivative constraints may be setted the other way around, needs a check (sign switch for i and i+1 nodes)
                    % f(x)_(i+1)-f(x)_(i) <=
                    % tauDisc/2*Scale.time*(timeFinal_n-timeInitial_n)*limit
                    % fprev_(i+1)-fprev_(i) + dG_dX_n_ii(i+1)-dG_dX_n_ii(i)+ dG_dU_n_ii(i+1)-dG_dU_n_ii(i)<= tauDisc/2*Scale.time*(timeFinal_n-timeInitial_n)*limit
                    % dG_dX_n_ii(i+1)-dG_dX_n_ii(i)+ dG_dU_n_ii(i+1)-dG_dU_n_ii(i)<= - (fprev_(i+1) -fprev_(i)+ dG_dX_n_ii(i+1)*Xprev-dG_dX_n_ii(i)*Xprev+ dG_dU_n_ii(i+1)*Uprev-dG_dU_n_ii(i)Uprev + tauDisc/2*Scale.time*(timeFinal_n-timeInitial_n)*limit
                    dG_dX_n_X_n   =  accumarray(reshape(repmat((0:(Nodes-2)) +obj.MIOFF,n.states,1),[],1),reshape(dG_dX_n_ii(:,0+obj.MIOFF:Nodes-2+obj.MIOFF).*X_n(:,0+obj.MIOFF:Nodes-2+obj.MIOFF),[],1));
                    dG_dU_n_U_n   =  accumarray(reshape(repmat((0:(Nodes-2)) +obj.MIOFF,n.controls,1),[],1),reshape(dG_dU_n_ii(:,0+obj.MIOFF:Nodes-2+obj.MIOFF).*U_n(:,0+obj.MIOFF:Nodes-2+obj.MIOFF),[],1));
                    dG_dX_n_X_n1   =  accumarray(reshape(repmat((0:(Nodes-2)) +obj.MIOFF,n.states,1),[],1),reshape(dG_dX_n_ii(:,1+obj.MIOFF:Nodes-1+obj.MIOFF).*X_n(:,1+obj.MIOFF:Nodes-1+obj.MIOFF),[],1));
                    dG_dU_n_U_n1   =  accumarray(reshape(repmat((0:(Nodes-2)) +obj.MIOFF,n.controls,1),[],1),reshape(dG_dU_n_ii(:,1+obj.MIOFF:Nodes-1+obj.MIOFF).*U_n(:,1+obj.MIOFF:Nodes-1+obj.MIOFF),[],1));
                    
                    rowINEQ_path    = [rowINEQ_path; (reshape(repmat((0:(Nodes-2)),2*(n.states + n.controls),1),[],1) + rowINEQ_Initial +obj.MIOFF)];
                    colINEQ_path    = [colINEQ_path; reshape([index.states(:,0+obj.MIOFF:Nodes-2+obj.MIOFF);index.controls(:,0+obj.MIOFF:Nodes-2+obj.MIOFF);index.states(:,1+obj.MIOFF:Nodes-1+obj.MIOFF);index.controls(:,1+obj.MIOFF:Nodes-1+obj.MIOFF)],[],1)];
                    valINEQ_path    = [valINEQ_path; typeSign*reshape(...
                        [dG_dX_n_ii(:,0+obj.MIOFF:Nodes-2+obj.MIOFF);...
                        dG_dU_n_ii(:,0+obj.MIOFF:Nodes-2+obj.MIOFF);...
                        -dG_dX_n_ii(:,1+obj.MIOFF:Nodes-1+obj.MIOFF);...
                        -dG_dU_n_ii(:,1+obj.MIOFF:Nodes-1+obj.MIOFF)...
                        ],[],1)/scalePath];
                    
                    if obj.problem.phases(1).freeTimeFinal
                        dg_dtf_n      = tauDisc/2*Scale.time*limit;
                        
                        rowINEQ_path  = [rowINEQ_path; ((0:(Nodes-2))' + rowINEQ_Initial +obj.MIOFF)];
                        colINEQ_path  = [colINEQ_path; index.timeFinal*ones(Nodes-1,1)];
                        valINEQ_path  = [valINEQ_path; -typeSign*dg_dtf_n/scalePath*ones(Nodes-1,1)];
                        pINEQ_path    = [pINEQ_path ; typeSign*ones(Nodes-1,1)*(timeDisc*limit - dg_dtf_n*timeFinal_n) + typeSign*(reshape(-f_ii(1:Nodes-1)+f_ii(2:Nodes),[],1)-dG_dX_n_X_n-dG_dU_n_U_n + dG_dX_n_X_n1 + dG_dU_n_U_n1)/scalePath];
                    else
                        pINEQ_path    = [pINEQ_path ; typeSign*ones(Nodes-1,1)*(timeDisc*limit) + typeSign*(reshape(-f_ii(1:Nodes-1)+f_ii(2:Nodes),[],1)-dG_dX_n_X_n-dG_dU_n_U_n + dG_dX_n_X_n1 + dG_dU_n_U_n1)/scalePath];
                    end
                elseif ~Path(ii+obj.MIOFF).integral && ~Path(ii+obj.MIOFF).derivative
                    % Nonlinear path constraint
                    g_n_function    =  (f_ii - limit)./scalePath;
                    dG_n_dX_n_X_n   =  accumarray(reshape(repmat((0:(Nodes-1)) +obj.MIOFF,n.states,1),[],1),dG_dX_n_ii(:).*X_n(:))./scalePath(:);
                    dG_n_dU_n_U_n   =  accumarray(reshape(repmat((0:(Nodes-1)) +obj.MIOFF,n.controls,1),[],1),dG_dU_n_ii(:).*U_n(:))./scalePath(:);
                    
                    rowINEQ_path    = [rowINEQ_path; (reshape(repmat((0:(Nodes-1)),n.states + n.controls,1),[],1) + rowINEQ_Initial +obj.MIOFF) ];
                    colINEQ_path    = [colINEQ_path; reshape([index.states;index.controls],[],1)];
                    valINEQ_path    = [valINEQ_path; reshape([typeSign*dG_dX_n_ii;typeSign*dG_dU_n_ii]./repmat(scalePath,n.states+n.controls,Naux),[],1)];
                    
                    if Path(ii+obj.MIOFF).buffer.include && VirtualBuffer.phases(1).include
                        rowINEQ_path  = [rowINEQ_path ; ((0:(Nodes-1))' + rowINEQ_Initial +obj.MIOFF)];
                        colINEQ_path  = [colINEQ_path ; reshape(index.virtualBuffers.path(jj+obj.MIOFF,:),[],1)];
                        valINEQ_path  = [valINEQ_path ; -ones(Nodes,1)];
                        jj = jj +1;
                    end
                    pINEQ_path      = [pINEQ_path ; - typeSign*(g_n_function(:) - dG_n_dX_n_X_n - dG_n_dU_n_U_n)];
                else
                    error('SCOPT Error, Path Constraint is not defined for nonlinear integral and derivative, or derivative terms \n Aborting')
                    
                end
                
            case 'quasiconvex'
                % Nonlinear path constaints which are
                % quasiconvex. Requires the second order taylor
                % series expansion term as input
                switch Path(ii+obj.MIOFF).type
                    case 'lower'
                        typeSign = -1;
                    case 'upper'
                        typeSign = 1;
                    otherwise
                        error('SCOPT Error, Path Constraints of constraint type %s not implemented. Use inequality constraints instead \n Aborting',Path(ii+obj.MIOFF).type)
                        
                end
                % Temporary fix for indexWhere for path
                % constraints, only valid to be defined from
                % 1:Nodes
                indexWhere = 1:Nodes;
                
                f_ii          = typeSign*Path(ii+obj.MIOFF).function(time(:,indexWhere),X(:,indexWhere),U(:,indexWhere),BodyMap,indexWhere);
                if (length(Path(ii+obj.MIOFF).statesIndex) + length(Path(ii+obj.MIOFF).controlsIndex) ==1)
                    % Single variable quasiconvex case
                    if length(Path(ii+obj.MIOFF).statesIndex)==1 % variable is a state
                        df_dvar       = typeSign*Path(ii+obj.MIOFF).jacobian.states(time(:,indexWhere),X(:,indexWhere),U(:,indexWhere),BodyMap,indexWhere);
                        d2f_dvar2     = typeSign*Path(ii+obj.MIOFF).hessian.states(time(:,indexWhere),X(:,indexWhere),U(:,indexWhere),BodyMap,indexWhere);
                        var           = X(Path(ii+obj.MIOFF).statesIndex,indexWhere);
                        scaleVar      = Scale.states(Path(ii+obj.MIOFF).statesIndex);
                        shiftVar      = Shift.states(Path(ii+obj.MIOFF).statesIndex);
                        indexVector   = index.states(Path(ii+obj.MIOFF).statesIndex,indexWhere);
                    else % variable is a control
                        df_dvar       = typeSign*Path(ii+obj.MIOFF).jacobian.controls(time(:,indexWhere),X(:,indexWhere),U(:,indexWhere),BodyMap,indexWhere);
                        d2f_dvar2     = typeSign*Path(ii+obj.MIOFF).hessian.controls(time(:,indexWhere),X(:,indexWhere),U(:,indexWhere),BodyMap,indexWhere);
                        var           = U(Path(ii+obj.MIOFF).controlsIndex,indexWhere);
                        scaleVar      = Scale.controls(Path(ii+obj.MIOFF).controlsIndex);
                        shiftVar      = Shift.controls(Path(ii+obj.MIOFF).controlsIndex);
                        indexVector   = index.controls(Path(ii+obj.MIOFF).controlsIndex,indexWhere);
                    end
                    
                    LT  = sqrt(d2f_dvar2/2);
                    qT  = (df_dvar - var.*d2f_dvar2);
                    r   = f_ii - df_dvar.*var + 1/2*d2f_dvar2.*var.^2;
                    Ai  = [qT/2;LT];
                    bi  = [(1+ r)/2;zeros(1,Nodes)];
                    ci  = -1/2*qT;
                    di  = (1- r)/2;
                    
                    qSOC_cvx_ii      = 3;
                    qSOC_cvx         = qSOC_cvx_ii*ones(1,Nodes);
                    rowSOC_path      = [rowSOC_path;(sum(qSOC_path) + reshape(repmat(reshape(repmat((0:qSOC_cvx_ii-1),1+n.states+n.controls,1),[],1),1,Nodes) + repmat((0:qSOC_cvx_ii:(Nodes-1)*qSOC_cvx_ii),qSOC_cvx_ii*(1+n.states+n.controls),1),[],1) +obj.MIOFF)];
                    colSOC_path      = [colSOC_path;reshape(repmat([indexVector;index.states;index.controls],qSOC_cvx_ii,1),[],1)];
                    valSOC_path      = [valSOC_path;reshape(-[scaleVar*ci;repmat([typeSign*(-1/2)*Scale.states(:).*Path(ii+obj.MIOFF).states(:);typeSign*(-1/2)*Scale.controls(:).*Path(ii+obj.MIOFF).controls(:)],1,Nodes);...
                        scaleVar*Ai(1,:);repmat([typeSign*Scale.states(:).*Path(ii+obj.MIOFF).states(:)/2;typeSign*Scale.controls(:).*Path(ii+obj.MIOFF).controls(:)/2],1,Nodes);scaleVar*Ai(2,:);zeros(n.states+n.controls,Nodes)]./repmat(scalePath,qSOC_cvx_ii*(1+n.states+n.controls),Nodes),[],1)];
                    
                    diShift          = -qT/2.*shiftVar - 1/2*sum(typeSign*Shift.states(:).*Path(ii+obj.MIOFF).states(:))- 1/2*sum(typeSign*Shift.controls(:).*Path(ii+obj.MIOFF).controls(:));
                    biShift          = [qT/2.*shiftVar + 1/2*sum(typeSign*Shift.states(:).*Path(ii+obj.MIOFF).states(:))+ 1/2*sum(typeSign*Shift.controls(:).*Path(ii+obj.MIOFF).controls(:));LT.*shiftVar] ;
                    
                    qSOC_path        = [qSOC_path,qSOC_cvx];
                    hSOC_path        = [hSOC_path;reshape([di + diShift;bi + biShift]./repmat(scalePath,qSOC_cvx_ii,Nodes),[],1)];
                else % general case
                    df_dX         = typeSign*Path(ii+obj.MIOFF).jacobian.states(time,X,U,BodyMap,index.nodes);
                    df_dU         = typeSign*Path(ii+obj.MIOFF).jacobian.controls(time,X,U,BodyMap,index.nodes);
                    if  Path(ii+obj.MIOFF).hessian.issparse
                    else
                        df_dXX        = typeSign*Path(ii+obj.MIOFF).hessian.states(time,X,U,BodyMap,index.nodes);
                        df_dUU        = typeSign*Path(ii+obj.MIOFF).hessian.controls(time,X,U,BodyMap,index.nodes);
                        df_dXU        = typeSign*Path(ii+obj.MIOFF).hessian.crossXU(time,X,U,BodyMap,index.nodes);
                        df_dUX        = typeSign*Path(ii+obj.MIOFF).hessian.crossUX(time,X,U,BodyMap,index.nodes);
                        
                        H = cat(2,cat(1,df_dXX,df_dUX),cat(1,df_dXU,df_dUU));
                        
                        % routine to check values with zeros and build reduced jacobian and hessian
                        
                        % routine to check co
                    end
                    
                    
                    LT  = chol(d2f_dvar2/2);
                    qT  = (df_dvar - var.*d2f_dvar2);
                    r   = f_ii - sum(df_dvar.*var,1) + 1/2*d2f_dvar2.*var.^2;
                    Ai  = [qT/2;LT];
                    bi  = [(1+ r)/2;zeros(1,Nodes)];
                    ci  = -1/2*qT;
                    di  = (1- r)/2;
                    
                    qSOC_cvx_ii      = 3;
                    qSOC_cvx         = qSOC_cvx_ii*ones(1,Nodes);
                    rowSOC_path      = [rowSOC_path;(sum(qSOC_path) + reshape(repmat(reshape(repmat((0:qSOC_cvx_ii-1),1+n.states+n.controls,1),[],1),1,Nodes) + repmat((0:qSOC_cvx_ii:(Nodes-1)*qSOC_cvx_ii),qSOC_cvx_ii*(1+n.states+n.controls),1),[],1) +obj.MIOFF)];
                    colSOC_path      = [colSOC_path;reshape(repmat([indexVector;index.states;index.controls],qSOC_cvx_ii,1),[],1)];
                    valSOC_path      = [valSOC_path;reshape(-[scaleVar*ci;repmat([typeSign*(-1/2)*Scale.states(:).*Path(ii+obj.MIOFF).states(:);typeSign*(-1/2)*Scale.controls(:).*Path(ii+obj.MIOFF).controls(:)],1,Nodes);...
                        scaleVar*Ai(1,:);repmat([typeSign*Scale.states(:).*Path(ii+obj.MIOFF).states(:)/2;typeSign*Scale.controls(:).*Path(ii+obj.MIOFF).controls(:)/2],1,Nodes);scaleVar*Ai(2,:);zeros(n.states+n.controls,Nodes)]./repmat(scalePath,qSOC_cvx_ii*(1+n.states+n.controls),Nodes),[],1)];
                    
                    diShift          = -qT/2.*shiftVar - 1/2*sum(typeSign*Shift.states(:).*Path(ii+obj.MIOFF).states(:))- 1/2*sum(typeSign*Shift.controls(:).*Path(ii+obj.MIOFF).controls(:));
                    biShift          = [qT/2.*shiftVar + 1/2*sum(typeSign*Shift.states(:).*Path(ii+obj.MIOFF).states(:))+ 1/2*sum(typeSign*Shift.controls(:).*Path(ii+obj.MIOFF).controls(:));LT.*shiftVar] ;
                    
                    qSOC_path        = [qSOC_path,qSOC_cvx];
                    hSOC_path        = [hSOC_path;reshape([di + diShift;bi + biShift]./repmat(scalePath,qSOC_cvx_ii,Nodes),[],1)];
                    error('SCOPT Error: Nonlinear General Quasiconvex Path Constraint not implemented yet \n')
                end
            case 'non-linear-1var'
                switch Path(ii+obj.MIOFF).type
                    case 'lower'
                        typeSign = -1;
                    case 'upper'
                        typeSign = 1;
                    otherwise
                        error('SCOPT Error, Path Constraints of constraint type %s not implemented. Use inequality constraints instead \n Aborting',Path(ii+obj.MIOFF).type)
                        
                end
                f_ii          = typeSign*Path(ii+obj.MIOFF).function(time,X,U,BodyMap,index.nodes);
                df_dvar       = typeSign*Path(ii+obj.MIOFF).jacobian.variable(time,X,U,BodyMap,index.nodes);
                d2f_dvar2     = typeSign*Path(ii+obj.MIOFF).hessian.variable(time,X,U,BodyMap,index.nodes);
                switch Path(ii+obj.MIOFF).variableType
                    case 'state'
                        var      = X(Path(ii+obj.MIOFF).variableIndex,:);
                        scaleVar = Scale.states(Path(ii+obj.MIOFF).variableIndex);
                        shiftVar = Shift.states(Path(ii+obj.MIOFF).variableIndex);
                        indexVector = index.states(Path(ii+obj.MIOFF).variableIndex,:);
                    case 'control'
                        var      = U(Path(ii+obj.MIOFF).variableIndex,:);
                        scaleVar = Scale.controls(Path(ii+obj.MIOFF).variableIndex);
                        shiftVar = Shift.controls(Path(ii+obj.MIOFF).variableIndex);
                        indexVector = index.controls(Path(ii+obj.MIOFF).variableIndex,:);
                    case 'undefined'
                        error('variable type not inserted')
                end
                switch Path(ii+obj.MIOFF).order
                    case 1
                    case 2
                        LT  = sqrt(d2f_dvar2/2);
                        qT  = (df_dvar - var.*d2f_dvar2);
                        r   = f_ii - df_dvar.*var + 1/2*d2f_dvar2.*var.^2;
                        Ai  = [qT/2;LT];
                        bi  = [(1+ r)/2;zeros(1,Nodes)];
                        ci  = -1/2*qT;
                        di  = (1- r)/2;
                        
                        qSOC_cvx_ii      = 3;
                        qSOC_cvx         = qSOC_cvx_ii*ones(1,Nodes);
                        rowSOC_path      = [rowSOC_path;(sum(qSOC_path) + reshape(repmat(reshape(repmat((0:qSOC_cvx_ii-1),1+n.states+n.controls,1),[],1),1,Nodes) + repmat((0:qSOC_cvx_ii:(Nodes-1)*qSOC_cvx_ii),qSOC_cvx_ii*(1+n.states+n.controls),1),[],1) +obj.MIOFF)];
                        colSOC_path      = [colSOC_path;reshape(repmat([indexVector;index.states;index.controls],qSOC_cvx_ii,1),[],1)];
                        valSOC_path      = [valSOC_path;reshape(-[scaleVar*ci;repmat([typeSign*(-1/2)*Scale.states(:).*Path(ii+obj.MIOFF).states(:);typeSign*(-1/2)*Scale.controls(:).*Path(ii+obj.MIOFF).controls(:)],1,Nodes);...
                            scaleVar*Ai(1,:);repmat([typeSign*Scale.states(:).*Path(ii+obj.MIOFF).states(:)/2;typeSign*Scale.controls(:).*Path(ii+obj.MIOFF).controls(:)/2],1,Nodes);scaleVar*Ai(2,:);zeros(n.states+n.controls,Nodes)]./repmat(scalePath,qSOC_cvx_ii*(1+n.states+n.controls),Nodes),[],1)];
                        
                        diShift          = -qT/2.*shiftVar - 1/2*sum(typeSign*Shift.states(:).*Path(ii+obj.MIOFF).states(:))- 1/2*sum(typeSign*Shift.controls(:).*Path(ii+obj.MIOFF).controls(:));
                        biShift          = [qT/2.*shiftVar + 1/2*sum(typeSign*Shift.states(:).*Path(ii+obj.MIOFF).states(:))+ 1/2*sum(typeSign*Shift.controls(:).*Path(ii+obj.MIOFF).controls(:));LT.*shiftVar] ;
                        
                        qSOC_path        = [qSOC_path,qSOC_cvx];
                        hSOC_path        = [hSOC_path;reshape([di + diShift;bi + biShift]./repmat(scalePath,qSOC_cvx_ii,Nodes),[],1)];
                    otherwise
                        error('SCOPT Error: Order not defined')
                end
                
                
            case {'convex','soc'}
                % Convex SOC constraints enabled for time
                % variant and time invariant constraints
                switch Path(ii+obj.MIOFF).type
                    case 'lower'
                        typeSign = -1;
                    case 'upper'
                        typeSign = 1;
                    otherwise
                        error('SCOPT Error, Path Constraints of constraint type %s not implemented. Use inequality constraints instead \n Aborting',Path(ii+obj.MIOFF).type)
                        
                end
                
                if ~Path(ii+obj.MIOFF).integral && ~Path(ii+obj.MIOFF).derivative
                    qSOC_cvx_ii      = Path(ii+obj.MIOFF).cone.dimensions;
                    NodesAux         = Path(ii+obj.MIOFF).cone.nodes;
                    qSOC_cvx         = qSOC_cvx_ii*ones(1,Nodes);
                    if NodesAux==1 % constant SOCs
                        rowSOC_path      = [rowSOC_path;(sum(qSOC_path) + reshape(repmat(reshape(repmat((0:qSOC_cvx_ii-1),n.states+n.controls,1),[],1),1,Nodes) + repmat((0:qSOC_cvx_ii:(Nodes-1)*qSOC_cvx_ii),(qSOC_cvx_ii)*(n.states+n.controls),1),[],1) +obj.MIOFF)];
                        colSOC_path      = [colSOC_path;reshape(repmat([index.states;index.controls],qSOC_cvx_ii,1),[],1)];
                        
                        valSOC_path      = [valSOC_path;reshape(repmat(-[typeSign*Scale.states(:).*Path(ii+obj.MIOFF).cone.right.states(:);typeSign*Scale.controls(:).*Path(ii+obj.MIOFF).cone.right.controls(:);...
                            reshape([repmat(Scale.states(:),1,Path(ii+obj.MIOFF).cone.dimensions-1).*Path(ii+obj.MIOFF).cone.norm.states;repmat(Scale.controls(:),1,Path(ii+obj.MIOFF).cone.dimensions-1).*Path(ii+obj.MIOFF).cone.norm.controls],[],1)],1,Nodes)./repmat(scalePath,qSOC_cvx_ii*(n.states+n.controls),Naux),[],1)];
                        
                        qSOC_path        = [qSOC_path,qSOC_cvx];
                        hSOC_path        = [hSOC_path;reshape(repmat([typeSign*(limit + Path(ii+obj.MIOFF).cone.right.cons + sum(Shift.states.*Path(ii+obj.MIOFF).cone.right.states,1)+ sum(Shift.controls.*Path(ii+obj.MIOFF).cone.right.controls,1)) ;...
                            (Path(ii+obj.MIOFF).cone.norm.cons(:) + reshape(sum(repmat(Shift.states,1,Path(ii+obj.MIOFF).cone.dimensions-1).*Path(ii+obj.MIOFF).cone.norm.states,1) + sum(repmat(Shift.controls,1,Path(ii+obj.MIOFF).cone.dimensions-1).*Path(ii+obj.MIOFF).cone.norm.controls,1),[],1))],1,Nodes)./repmat(scalePath,qSOC_cvx_ii,Naux),[],1)];
                    else % node varying SOC
                        rowSOC_path      = [rowSOC_path;(sum(qSOC_path) + reshape(repmat(reshape(repmat((0:qSOC_cvx_ii-1),n.states+n.controls,1),[],1),1,Nodes) + repmat((0:qSOC_cvx_ii:(Nodes-1)*qSOC_cvx_ii),(qSOC_cvx_ii)*(n.states+n.controls),1),[],1) +obj.MIOFF)];
                        colSOC_path      = [colSOC_path;reshape(repmat([index.states;index.controls],qSOC_cvx_ii,1),[],1)];
                        
                        valSOC_path      = [valSOC_path;reshape(-[typeSign*repmat(Scale.states(:),1,NodesAux).*Path(ii+obj.MIOFF).cone.right.states;typeSign*repmat(Scale.controls(:),1,NodesAux).*Path(ii+obj.MIOFF).cone.right.controls;...
                            reshape([repmat(Scale.states(:),1,Path(ii+obj.MIOFF).cone.dimensions-1,NodesAux).*Path(ii+obj.MIOFF).cone.norm.states;repmat(Scale.controls(:),1,Path(ii+obj.MIOFF).cone.dimensions-1,NodesAux).*Path(ii+obj.MIOFF).cone.norm.controls],(n.states+n.controls)*(qSOC_cvx_ii-1),NodesAux)]./repmat(scalePath,qSOC_cvx_ii*(n.states+n.controls),Naux),[],1)];
                        
                        qSOC_path        = [qSOC_path,qSOC_cvx];
                        hSOC_path        = [hSOC_path;reshape([typeSign*(limit + Path(ii+obj.MIOFF).cone.right.cons + sum(repmat(Shift.states,1,NodesAux).*Path(ii+obj.MIOFF).cone.right.states,1)+ sum(repmat(Shift.controls,1,NodesAux).*Path(ii+obj.MIOFF).cone.right.controls,1)) ;...
                            (Path(ii+obj.MIOFF).cone.norm.cons + reshape(sum(repmat(Shift.states,1,qSOC_cvx_ii-1,NodesAux).*Path(ii+obj.MIOFF).cone.norm.states,1) + sum(repmat(Shift.controls,1,qSOC_cvx_ii-1,NodesAux).*Path(ii+obj.MIOFF).cone.norm.controls,1),(qSOC_cvx_ii-1),NodesAux))]./repmat(scalePath,qSOC_cvx_ii,Naux),[],1)];
                    end
                else
                    error('SCOPT Error, Path Constraint is not defined for soc integral and derivative, or derivative, or integral terms \n Aborting')
                    
                end
            otherwise
                warning('SCOPT Warning, Path function of type %s not implemented \n Continuing without constraint %i \n',Path(ii+obj.MIOFF).funType,ii)
                continue
        end
        if ~isempty(rowINEQ_path)
            rowINEQ_Initial= max(rowINEQ_path);
        end
        if ~isempty(rowEQ_path)
            rowEQ_Initial= max(rowEQ_path);
        end
    end
    if isempty(rowINEQ_Initial)
        rowINEQ_Initial = 0;
    end
    if isempty(rowEQ_Initial)
        rowEQ_Initial = 0;
    end
    
    GSOC_path  = sparse(rowSOC_path,colSOC_path,valSOC_path,sum(qSOC_path),Ntot);
    MINEQ_path = sparse(rowINEQ_path,colINEQ_path,valINEQ_path,rowINEQ_Initial,Ntot);
    MEQ_path   = sparse(rowEQ_path,colEQ_path,valEQ_path,rowEQ_Initial,Ntot);
else
    MEQ_path   = sparse(0,Ntot);
    pEQ_path   = zeros(0,1);
    MINEQ_path = sparse(0,Ntot);
    pINEQ_path = zeros(0,1);
    GSOC_path  = sparse(0,Ntot);
    hSOC_path  = zeros(0,1);
    qSOC_path  = zeros(1,0);
end
obj.problem.phases(1).transcription.path.MEQ   = MEQ_path;
obj.problem.phases(1).transcription.path.pEQ   = pEQ_path;
obj.problem.phases(1).transcription.path.MINEQ = MINEQ_path;
obj.problem.phases(1).transcription.path.pINEQ = pINEQ_path;
obj.problem.phases(1).transcription.path.GSOC  = GSOC_path;
obj.problem.phases(1).transcription.path.hSOC  = hSOC_path;
obj.problem.phases(1).transcription.path.qSOC  = qSOC_path;
%% Dynamics Matrix (Equality Constraints
% Sets dynamics matrices in sparse representation., with
% numerical or exact solutions depending on user input
rowcolEQ_dyn    = (1:((Nodes-1)*n.states ))';
valEQ_dynx_eye  = ones((Nodes-1)*n.states,1);

invScaleState = 1./Scale.states;
switch obj.algorithm.collocationMethod
    case {'euler','trapezoidal'}
        dX_n_dt_ii      = repmat((1./Scale.states(:)),1,Nodes) .* obj.problem.phases(1).dynamics.stateDerivativeFunction(time,X,U,BodyMap,index.nodes);
        dX_n_dt_col     = dX_n_dt_ii(:);
        switch obj.problem.phases(1).dynamics.typeStateMatrix
            case 0 % full matrix
                A = obj.problem.phases(1).dynamics.stateMatrixFunction(time,X,U,BodyMap,index.nodes);
                idxEQ_dyn_A_aux = find(A);
                valEQ_dyn_A_aux = A(idxEQ_dyn_A_aux);
                [rowEQ_dyn_A_aux,colEQ_dyn_A_aux,nodeEQ_dyn_A_aux]=ind2sub([n.states,n.states,Nodes],idxEQ_dyn_A_aux);
                rowEQ_dyn_A_n  = rowEQ_dyn_A_aux + (nodeEQ_dyn_A_aux-1)*n.states + index.states(1,1)-1;
                colEQ_dyn_A_n  = colEQ_dyn_A_aux + (nodeEQ_dyn_A_aux-1)*n.states + index.states(1,1)-1;
                valEQ_dyn_A_n  = invScaleState(rowEQ_dyn_A_aux).*valEQ_dyn_A_aux.*Scale.states(colEQ_dyn_A_aux);
                nEntriesA = histc(nodeEQ_dyn_A_aux, unique(nodeEQ_dyn_A_aux));
            case 1 % sparse matrix
                [rowEQ_dyn_A_aux,colEQ_dyn_A_aux,valEQ_dyn_A_aux,nEntriesA] = obj.problem.phases(1).dynamics.stateMatrixFunction(time,X,U,BodyMap,index.nodes);
                valEQ_dyn_A_n_aux = reshape(invScaleState(rowEQ_dyn_A_aux),size(rowEQ_dyn_A_aux)).*valEQ_dyn_A_aux.*reshape(Scale.states(colEQ_dyn_A_aux),size(colEQ_dyn_A_aux));
                rowEQ_dyn_A_n   = reshape(rowEQ_dyn_A_aux + repmat(index.states(1,:)-1,nEntriesA,1),[],1);
                colEQ_dyn_A_n   = reshape(colEQ_dyn_A_aux + repmat(index.states(1,:)-1,nEntriesA,1),[],1);
                valEQ_dyn_A_n   = reshape(valEQ_dyn_A_n_aux,[],1);
        end
        AX    = accumarray(rowEQ_dyn_A_n,valEQ_dyn_A_n.*reshape(X_n(colEQ_dyn_A_n),[],1),[n.states*Nodes,1]);
        switch obj.problem.phases(1).dynamics.typeControlMatrix
            case 0 % full matrix
                B = obj.problem.phases(1).dynamics.controlMatrixFunction(time,X,U,BodyMap,index.nodes);
                idxEQ_dyn_B_aux = find(B);
                valEQ_dyn_B_aux = B(idxEQ_dyn_B_aux);
                [rowEQ_dyn_B_aux,colEQ_dyn_B_aux,nodeEQ_dyn_B_aux]=ind2sub([n.states,n.controls,Nodes],idxEQ_dyn_B_aux);
                rowEQ_dyn_B_n  = rowEQ_dyn_B_aux + (nodeEQ_dyn_B_aux-1)*n.states + index.states(1,1)-1;
                colEQ_dyn_B_n  = colEQ_dyn_B_aux + (nodeEQ_dyn_B_aux-1)*n.controls + index.controls(1,1)-1;
                valEQ_dyn_B_n  = invScaleState(rowEQ_dyn_B_aux).*valEQ_dyn_B_aux.*Scale.controls(colEQ_dyn_B_aux);
                nEntriesB = histc(nodeEQ_dyn_B_aux, unique(nodeEQ_dyn_B_aux));
            case 1 % sparse matrix
                [rowEQ_dyn_B_aux,colEQ_dyn_B_aux,valEQ_dyn_B_aux,nEntriesB] = obj.problem.phases(1).dynamics.controlMatrixFunction(time,X,U,BodyMap,index.nodes);
                valEQ_dyn_B_n_aux = reshape(invScaleState(rowEQ_dyn_B_aux),size(rowEQ_dyn_B_aux)).*valEQ_dyn_B_aux.*reshape(Scale.controls(colEQ_dyn_B_aux),size(colEQ_dyn_B_aux));
                rowEQ_dyn_B_n   = reshape(rowEQ_dyn_B_aux + repmat(index.states(1,:)-1,nEntriesB,1),[],1);
                colEQ_dyn_B_n   = reshape(colEQ_dyn_B_aux + repmat(index.controls(1,:)-1,nEntriesB,1),[],1);
                valEQ_dyn_B_n   = reshape(valEQ_dyn_B_n_aux,[],1);
        end
        BU    = accumarray(rowEQ_dyn_B_n,valEQ_dyn_B_n.*reshape(U_n(colEQ_dyn_B_n- Nodes*n.states),[],1),[n.states*Nodes,1]);
    case {'exact'}
        f_n_ii      = repmat((1./Scale.states(:)),1,Nodes-1) .* obj.problem.phases(1).dynamics.stateDerivativeFunction(time,X,U,BodyMap,index.nodes);
        f_n_col     = f_n_ii(:);
        
        switch obj.problem.phases(1).dynamics.typeStateMatrix
            case 0 % full matrix
                [A,A1] = obj.problem.phases(1).dynamics.stateMatrixFunction(time,X,U,BodyMap,index.nodes);
                idxEQ_dyn_A_aux = find(A);
                valEQ_dyn_A_aux = A(idxEQ_dyn_A_aux);
                [rowEQ_dyn_A_aux,colEQ_dyn_A_aux,nodeEQ_dyn_A_aux]=ind2sub([n.states,n.states,Nodes-1],idxEQ_dyn_A_aux);
                rowEQ_dyn_A_n  = rowEQ_dyn_A_aux + (nodeEQ_dyn_A_aux-1)*n.states + index.states(1,1)-1;
                colEQ_dyn_A_n  = colEQ_dyn_A_aux + (nodeEQ_dyn_A_aux-1)*n.states + index.states(1,1)-1;
                valEQ_dyn_A_n  = invScaleState(rowEQ_dyn_A_aux).*valEQ_dyn_A_aux.*Scale.states(colEQ_dyn_A_aux);
                idxEQ_dyn_A_aux1 = find(A1);
                valEQ_dyn_A_aux1 = A1(idxEQ_dyn_A_aux1);
                [rowEQ_dyn_A_aux1,colEQ_dyn_A_aux1,nodeEQ_dyn_A_aux1]=ind2sub([n.states,n.states,Nodes-1],idxEQ_dyn_A_aux1);
                rowEQ_dyn_A_n1  = rowEQ_dyn_A_aux1 + (nodeEQ_dyn_A_aux1-1)*n.states + index.states(1,1)-1;
                colEQ_dyn_A_n1  = colEQ_dyn_A_aux1 + (nodeEQ_dyn_A_aux1  )*n.states + index.states(1,1)-1;
                valEQ_dyn_A_n1  = invScaleState(rowEQ_dyn_A_aux1).*valEQ_dyn_A_aux1.*Scale.states(colEQ_dyn_A_aux1);
                
            case 1 % sparse matrix
                [rowEQ_dyn_A_aux,colEQ_dyn_A_aux,valEQ_dyn_A_aux,nEntriesA,rowEQ_dyn_A_aux1,colEQ_dyn_A_aux1,valEQ_dyn_A_aux1,nEntriesA1] = obj.problem.phases(1).dynamics.stateMatrixFunction(time,X,U,BodyMap,index.nodes);
                valEQ_dyn_A_n_aux = reshape(invScaleState(rowEQ_dyn_A_aux),size(rowEQ_dyn_A_aux)).*valEQ_dyn_A_aux.*reshape(Scale.states(colEQ_dyn_A_aux),size(colEQ_dyn_A_aux));
                rowEQ_dyn_A_n   = reshape(rowEQ_dyn_A_aux + repmat(index.states(1,1:Nodes-1)-1,nEntriesA,1),[],1);
                colEQ_dyn_A_n   = reshape(colEQ_dyn_A_aux + repmat(index.states(1,1:Nodes-1)-1,nEntriesA,1),[],1);
                valEQ_dyn_A_n   = reshape(valEQ_dyn_A_n_aux,[],1);
                valEQ_dyn_A_n_aux1 = reshape(invScaleState(rowEQ_dyn_A_aux1),size(rowEQ_dyn_A_aux1)).*valEQ_dyn_A_aux1.*reshape(Scale.states(colEQ_dyn_A_aux1),size(colEQ_dyn_A_aux1));
                rowEQ_dyn_A_n1   = reshape(rowEQ_dyn_A_aux1 + repmat(index.states(1,1:Nodes-1)-1,nEntriesA1,1),[],1);
                colEQ_dyn_A_n1   = reshape(colEQ_dyn_A_aux1 + repmat(index.states(1,2:Nodes)-1  ,nEntriesA1,1),[],1);
                valEQ_dyn_A_n1   = reshape(valEQ_dyn_A_n_aux1,[],1);
        end
        AX     = accumarray(rowEQ_dyn_A_n ,valEQ_dyn_A_n .*reshape(X_n(colEQ_dyn_A_n) ,[],1),[n.states*(Nodes-1),1]);
        AX1    = accumarray(rowEQ_dyn_A_n1,valEQ_dyn_A_n1.*reshape(X_n(colEQ_dyn_A_n1),[],1),[n.states*(Nodes-1),1]);
        
        
        switch obj.problem.phases(1).dynamics.typeControlMatrix
            case 0 % full matrix
                [B,B1] = obj.problem.phases(1).dynamics.controlMatrixFunction(time,X,U,BodyMap,index.nodes);
                
                idxEQ_dyn_B_aux = find(B);
                valEQ_dyn_B_aux = B(idxEQ_dyn_B_aux);
                [rowEQ_dyn_B_aux,colEQ_dyn_B_aux,nodeEQ_dyn_B_aux]=ind2sub([n.states,n.controls,Nodes-1],idxEQ_dyn_B_aux);
                
                rowEQ_dyn_B_n  = rowEQ_dyn_B_aux + (nodeEQ_dyn_B_aux-1)*n.states + index.states(1,1)-1;
                colEQ_dyn_B_n  = colEQ_dyn_B_aux + (nodeEQ_dyn_B_aux-1)*n.controls + index.controls(1,1)-1;
                valEQ_dyn_B_n  = invScaleState(rowEQ_dyn_B_aux).*valEQ_dyn_B_aux.*Scale.controls(colEQ_dyn_B_aux);
                
                idxEQ_dyn_B_aux1 = find(B1);
                valEQ_dyn_B_aux1 = B1(idxEQ_dyn_B_aux1);
                [rowEQ_dyn_B_aux1,colEQ_dyn_B_aux1,nodeEQ_dyn_B_aux1]=ind2sub([n.states,n.controls,Nodes-1],idxEQ_dyn_B_aux1);
                
                rowEQ_dyn_B_n1  = rowEQ_dyn_B_aux1 + (nodeEQ_dyn_B_aux1-1)*n.states + index.states(1,1)-1;
                colEQ_dyn_B_n1  = colEQ_dyn_B_aux1 + (nodeEQ_dyn_B_aux1)*n.controls + index.controls(1,1)-1;
                valEQ_dyn_B_n1  = invScaleState(rowEQ_dyn_B_aux1).*valEQ_dyn_B_aux1.*Scale.controls(colEQ_dyn_B_aux1);
                
            case 1 % sparse matrix
                [rowEQ_dyn_B_aux,colEQ_dyn_B_aux,valEQ_dyn_B_aux,nEntriesB,rowEQ_dyn_B_aux1,colEQ_dyn_B_aux1,valEQ_dyn_B_aux1,nEntriesB1] = obj.problem.phases(1).dynamics.controlMatrixFunction(time,X,U,BodyMap,index.nodes);
                valEQ_dyn_B_n_aux = reshape(invScaleState(rowEQ_dyn_B_aux),size(rowEQ_dyn_B_aux)).*valEQ_dyn_B_aux.*reshape(Scale.controls(colEQ_dyn_B_aux),size(colEQ_dyn_B_aux));
                rowEQ_dyn_B_n   = reshape(rowEQ_dyn_B_aux  + repmat(index.states(1,1:Nodes-1)-1,nEntriesB,1),[],1);
                colEQ_dyn_B_n   = reshape(colEQ_dyn_B_aux  + repmat(index.controls(1,1:Nodes-1)-1,nEntriesB,1),[],1);
                valEQ_dyn_B_n   = reshape(valEQ_dyn_B_n_aux,[],1);
                valEQ_dyn_B_n_aux1 = reshape(invScaleState(rowEQ_dyn_B_aux1),size(rowEQ_dyn_B_aux1)).*valEQ_dyn_B_aux1.*reshape(Scale.controls(colEQ_dyn_B_aux1),size(colEQ_dyn_B_aux1));
                rowEQ_dyn_B_n1  = reshape(rowEQ_dyn_B_aux1 + repmat(index.states(1,1:Nodes-1)-1,nEntriesB1,1),[],1);
                colEQ_dyn_B_n1  = reshape(colEQ_dyn_B_aux1 + repmat(index.controls(1,2:Nodes)-1,nEntriesB1,1),[],1);
                valEQ_dyn_B_n1  = reshape(valEQ_dyn_B_n_aux1,[],1);
        end
        BU     = accumarray(rowEQ_dyn_B_n ,valEQ_dyn_B_n .*reshape(U_n(colEQ_dyn_B_n - Nodes*n.states),[],1),[n.states*(Nodes-1),1]);
        BU1    = accumarray(rowEQ_dyn_B_n1,valEQ_dyn_B_n1.*reshape(U_n(colEQ_dyn_B_n1- Nodes*n.states),[],1),[n.states*(Nodes-1),1]);
end

switch obj.algorithm.collocationMethod
    case 'euler' % Euler numerical derivative
        % xi+1 = xi + dt*fi
        % xi+1 = xi + dtau/2*(tfk-1-t0)*(Ai*(xi-xik-1) + Bi*(ui-uik-1) + fik-1) +
        % (dtau/2*fik-1 + dtk-1*dfdtk-1)*(tf-tfk-1)
        rowEQ_dyn   = [repmat(rowcolEQ_dyn,2,1) ; rowEQ_dyn_A_n(1:end-nEntriesA(end)) ; rowEQ_dyn_B_n(1:end-nEntriesB(end)) ];
        colEQ_dyn   = [rowcolEQ_dyn ;rowcolEQ_dyn + n.states ; colEQ_dyn_A_n(1:end-nEntriesA(end)) ; colEQ_dyn_B_n(1:end-nEntriesB(end))];
        valEQ_dyn   = [valEQ_dynx_eye ;-valEQ_dynx_eye ; timeDisc*valEQ_dyn_A_n(1:end-nEntriesA(end)) ; timeDisc*valEQ_dyn_B_n(1:end-nEntriesB(end)) ];
        f_col       = timeDisc*dX_n_dt_col(rowcolEQ_dyn);
        pEQ_dyn     = - f_col + timeDisc*(AX(rowcolEQ_dyn) + BU(rowcolEQ_dyn));
        if obj.problem.phases(1).freeTimeFinal % Free time of flight
            df_dtf_n_col = Scale.time*tauDisc*(1/2)*(dX_n_dt_col(rowcolEQ_dyn));
            rowEQ_dyn = [rowEQ_dyn ; rowcolEQ_dyn];
            colEQ_dyn = [colEQ_dyn ; index.timeFinal*valEQ_dynx_eye];
            valEQ_dyn = [valEQ_dyn ; df_dtf_n_col];
            pEQ_dyn   = pEQ_dyn + df_dtf_n_col*timeFinal_n;
        end
        rowEQnoVirtual_dyn = rowEQ_dyn;
        colEQnoVirtual_dyn = colEQ_dyn;
        valEQnoVirtual_dyn = valEQ_dyn;
        if VirtualControl.phases(1).include
            rowEQ_dyn    = [rowEQ_dyn ; reshape(repmat(n.states*(index.nodes(1:end-1)-1),n.virtualControls,1) + repmat(VirtualControl.phases(1).rowE(:),1,Nodes-1),[],1)];
            colEQ_dyn    = [colEQ_dyn ; reshape(index.virtualControls,[],1)];
            valEQ_dyn    = [valEQ_dyn ; reshape(repmat(VirtualControl.phases(1).scale*VirtualControl.phases(1).valE(:),1,Nodes-1),[],1)];
        end
        MEQ_dyn     = sparse(rowEQ_dyn,colEQ_dyn,valEQ_dyn,(Nodes-1)*n.states,Ntot);
    case 'trapezoidal' % Trapezoidal numerical derivative
        % xi+1 = xi + dt*(fi+fi+1)/2
        % xi+1 = xi + dtau/2*(tfk-1-t0)*1/2*(Ai*(xi-xik-1) + Bi*(ui-uik-1) + Ai+1*(xi+1-xi+1k-1) + Bi+1*(ui+1-ui+1k-1) + fik-1 + fi+1k-1) +
        % (dtau/2*(fik-1+fi+1k-1)/2 + dtk-1*(dfdtik-1+dfdti+1k-1)/2)*(tf-tfk-1)
        rowEQ_dyn   = [repmat(rowcolEQ_dyn,2,1) ; rowEQ_dyn_A_n(1:end-nEntriesA(end)) ; rowEQ_dyn_A_n(nEntriesA(1)+1:end)-n.states ; rowEQ_dyn_B_n(1:end-nEntriesB(end)) ; rowEQ_dyn_B_n(nEntriesB(1)+1:end)-n.states];
        colEQ_dyn   = [rowcolEQ_dyn ;rowcolEQ_dyn + n.states ; colEQ_dyn_A_n(1:end-nEntriesA(end)) ; colEQ_dyn_A_n(nEntriesA(1)+1:end) ; colEQ_dyn_B_n(1:end-nEntriesB(end)) ; colEQ_dyn_B_n(nEntriesB(1)+1:end)];
        valEQ_dyn   = [valEQ_dynx_eye ;-valEQ_dynx_eye ; timeDisc*(1/2)*valEQ_dyn_A_n(1:end-nEntriesA(end)) ; timeDisc*(1/2)*valEQ_dyn_A_n(nEntriesA(1)+1:end) ; timeDisc*(1/2)*valEQ_dyn_B_n(1:end-nEntriesB(end)) ; timeDisc*(1/2)*valEQ_dyn_B_n(nEntriesB(1)+1:end)];
        f_col       = timeDisc*(1/2)*(dX_n_dt_col(rowcolEQ_dyn) + dX_n_dt_col(rowcolEQ_dyn + n.states));
        pEQ_dyn     = -f_col + timeDisc*(1/2)*(AX(rowcolEQ_dyn) + BU(rowcolEQ_dyn)  + AX(rowcolEQ_dyn + n.states) + BU(rowcolEQ_dyn + n.states));
        if obj.problem.phases(1).freeTimeFinal % Free time of flight
            df_dtf_n_col = Scale.time*tauDisc*(1/2)*(1/2)*(dX_n_dt_col(rowcolEQ_dyn) + dX_n_dt_col(rowcolEQ_dyn + n.states));
            rowEQ_dyn = [rowEQ_dyn ; rowcolEQ_dyn];
            colEQ_dyn = [colEQ_dyn ; index.timeFinal*valEQ_dynx_eye];
            valEQ_dyn = [valEQ_dyn ; df_dtf_n_col];
            pEQ_dyn   = pEQ_dyn + df_dtf_n_col*timeFinal_n;
        end
        rowEQnoVirtual_dyn = rowEQ_dyn;
        colEQnoVirtual_dyn = colEQ_dyn;
        valEQnoVirtual_dyn = valEQ_dyn;
        if VirtualControl.phases(1).include
            rowEQ_dyn    = [rowEQ_dyn ; reshape(repmat(n.states*(index.nodes(1:end-1)-1),n.virtualControls,1) + repmat(VirtualControl.phases(1).rowE(:),1,Nodes-1),[],1)];
            colEQ_dyn    = [colEQ_dyn ; reshape(index.virtualControls,[],1)];
            valEQ_dyn    = [valEQ_dyn ; reshape(repmat(VirtualControl.phases(1).scale*VirtualControl.phases(1).valE(:),1,Nodes-1),[],1)];
        end
        MEQ_dyn     = sparse(rowEQ_dyn,colEQ_dyn,valEQ_dyn,(Nodes-1)*n.states,Ntot);
    case 'exact' % Exact ode solution
        % xi+1 = xi + fk-1 + Ai*(xi-xik-1) + Bi*(ui-uik-1) + Ani+1*(xi+1-xi+1k-1) + Bni+1*(ui+1-ui+1k-1) +
        % P*(p-pk-1)  + Tf*(tf-tfk-1) + T0*(t0-t0k-1) + Evi
        
        rowEQ_dyn   = [repmat(rowcolEQ_dyn,2,1) ; rowEQ_dyn_A_n ; rowEQ_dyn_A_n1 ; rowEQ_dyn_B_n ; rowEQ_dyn_B_n1];
        colEQ_dyn   = [rowcolEQ_dyn ;rowcolEQ_dyn + n.states ; colEQ_dyn_A_n ; colEQ_dyn_A_n1 ; colEQ_dyn_B_n ; colEQ_dyn_B_n1];
        valEQ_dyn   = [valEQ_dynx_eye ;-valEQ_dynx_eye ; valEQ_dyn_A_n ; valEQ_dyn_A_n1 ; valEQ_dyn_B_n ; valEQ_dyn_B_n1];
        pEQ_dyn     = -f_n_col(rowcolEQ_dyn) + AX + BU  + AX1 + BU1;
        if obj.problem.phases(1).freeTimeFinal || obj.problem.phases(1).freeTimeInitial % Free final or initial time
            df_n_dt_ii      = repmat((1./Scale.states(:)),1,Nodes-1) .* obj.problem.phases(1).dynamics.timeDerivativeFunction(time,X,U,BodyMap,index.nodes);
            df_n_dt_col     = df_n_dt_ii(:);
            df_n_dt_n_col   = Scale.time*tauDisc*(1/2)*(df_n_dt_col(rowcolEQ_dyn));
        end
        if obj.problem.phases(1).freeTimeFinal % Free final time
            rowEQ_dyn = [rowEQ_dyn ; rowcolEQ_dyn];
            colEQ_dyn = [colEQ_dyn ; index.timeFinal*valEQ_dynx_eye];
            valEQ_dyn = [valEQ_dyn ; df_n_dt_n_col];
            pEQ_dyn   = pEQ_dyn + df_n_dt_n_col*timeFinal_n;
        end
        if obj.problem.phases(1).freeTimeInitial % Free initial time
            rowEQ_dyn = [rowEQ_dyn ; rowcolEQ_dyn];
            colEQ_dyn = [colEQ_dyn ; index.timeInitial*valEQ_dynx_eye];
            valEQ_dyn = [valEQ_dyn ; -df_n_dt_n_col];
            pEQ_dyn   = pEQ_dyn - df_n_dt_n_col*timeInitial_n;
        end
        rowEQnoVirtual_dyn = rowEQ_dyn;
        colEQnoVirtual_dyn = colEQ_dyn;
        valEQnoVirtual_dyn = valEQ_dyn;
        if VirtualControl.phases(1).include
            rowEQ_dyn    = [rowEQ_dyn ; reshape(repmat(n.states*(index.nodes(1:end-1)-1),n.virtualControls,1) + repmat(VirtualControl.phases(1).rowE(:),1,Nodes-1),[],1)];
            colEQ_dyn    = [colEQ_dyn ; reshape(index.virtualControls,[],1)];
            valEQ_dyn    = [valEQ_dyn ; reshape(repmat(VirtualControl.phases(1).scale*VirtualControl.phases(1).valE(:),1,Nodes-1),[],1)];
        end
        MEQ_dyn     = sparse(rowEQ_dyn,colEQ_dyn,valEQ_dyn,(Nodes-1)*n.states,Ntot);
    otherwise
        error('SCOPT Error: Collocation method not defined')
end
obj.problem.phases(1).transcription.dynamics.MEQnoVirtual   = sparse(rowEQnoVirtual_dyn,colEQnoVirtual_dyn,valEQnoVirtual_dyn,(Nodes-1)*n.states,Ntot);
obj.problem.phases(1).transcription.dynamics.pEQ   = pEQ_dyn;
%% Sequential Algorithm Constraints
if obj.algorithm.sequential.activate
    switch obj.algorithm.sequential.type
        case 'component-trust-region'
            %% Inequalities (TODO: modify and adapt for inf input)
            rowINEQ_trustX     = reshape(repmat((0:(2*n.states):((2*n.states)*(Nodes-1))),2*n.states,1) + repmat((0:n.states-1)',2,1) + [zeros(n.states,1);n.states*ones(n.states,1)] +obj.MIOFF,[],1);
            colINEQ_trustX     = reshape(repmat(index.states,2,1),[],1);
            valINEQ_trustX     = repmat([1*ones(n.states,1);-1*ones(n.states,1)],Nodes,1);
            MINEQ_trustX       = sparse(rowINEQ_trustX,colINEQ_trustX,valINEQ_trustX,2*Nodes*n.states,Ntot);
            pINEQ_trustX       = reshape([obj.algorithm.sequential.trustRegion.phases(1).states.component_n + X_n;Tobj.algorithm.sequential.trustRegion.phases(1).states.component_n - X_n],[],1);
            
            if obj.problem.phases(1).freeTimeFinal
                rowINEQ_trustU     = reshape(repmat((0:(2*n.controls):((2*n.controls)*(Nodes-1))),2*n.controls,1) + repmat((0:n.controls-1)',2,1) + [zeros(n.controls,1);n.controls*ones(n.controls,1)] +obj.MIOFF,[],1);
                colINEQ_trustU     = reshape(repmat(index.controls,2,1),[],1);
                valINEQ_trustU     = repmat([1*ones(n.controls,1);-1*ones(n.controls,1)],Nodes,1);
                MINEQ_trustU       = sparse(rowINEQ_trustU,colINEQ_trustU,valINEQ_trustU,2*Nodes*n.controls,Ntot);
                pINEQ_trustU       = repmat([obj.algorithm.sequential.trustRegion.phases(1).controls.component_n + U_n;obj.algorithm.sequential.trustRegion.phases(1).controls.component_n - U_n],Nodes,1);
            else
                MINEQ_trustU = sparse(0,Ntot);
                pINEQ_trustU = zeros(0,1);
            end
            
            MINEQ_trust = [MINEQ_trustX;MINEQ_trustU];
            pINEQ_trust = [pINEQ_trustX;pINEQ_trustU];
            
            GSOC_trust = sparse(0,Ntot);
            hSOC_trust = zeros(0,1);
            qSOC_trust = 0*ones(1,0);
        case {'trust-region'}
            %% SOC trust region
            nVar      = 0;
            % gets variables to be considered for trust region
            % in statea and controls
            VarIndex  = zeros(0,Nodes);
            Var_n     = zeros(0,Nodes);
            if obj.algorithm.sequential.trustRegion.include.states
                nVar     = nVar + n.states;
                VarIndex = [VarIndex;index.states];
                Var_n    = [Var_n;X_n];
            end
            if obj.algorithm.sequential.trustRegion.include.controls
                nVar    = nVar + n.controls;
                VarIndex  = [VarIndex;index.controls];
                Var_n    = [Var_n;U_n];
            end
            switch obj.algorithm.sequential.trustRegion.variablesPenalty
                case 0 % quadratic trust region of variables per node
                    qSOC_trust_ii       = 2 + nVar;
                    r                   = sum(Var_n.*Var_n,1); %- obj.algorithm.sequential.trustRegion.radius;
                    qSOC_trustP1        = qSOC_trust_ii*ones(1,Nodes);
                    
                    rowSOC_trust        = reshape(repmat((0:qSOC_trust_ii:(qSOC_trust_ii*(Nodes-1))),nVar+1+nVar+1+nVar,1) + [zeros(nVar+1,1);ones(nVar+1,1);(2:qSOC_trust_ii-1)']+obj.MIOFF,[],1);
                    colSOC_trust        = reshape([VarIndex;index.trustRegionP1;VarIndex;index.trustRegionP1;VarIndex],[],1);
                    valSOC_trust        = -reshape([Var_n;1/2*ones(1,Nodes);-Var_n;-1/2*ones(1,Nodes);ones(nVar,Nodes)],[],1);
                    GSOC_trustP1        = sparse(rowSOC_trust,colSOC_trust,valSOC_trust,sum(qSOC_trustP1),Ntot);
                    
                    hSOC_trustP1        = reshape([(1- r)/2;(1+ r)/2;zeros(nVar,Nodes)],[],1);
                    
                    MINEQ_trustP1 = sparse(0,Ntot);
                    pINEQ_trustP1 = zeros(0,1);
                    MEQ_trustP1 = sparse(0,Ntot);
                    pEQ_trustP1 = zeros(0,1);
                case 1 % sum of absolutes trust region variables per node
                    qINEQ_trust_ii     = n.trustRegionP1s*2; % 3 for 1st and 3 for 2nd
                    
                    rowINEQ_trust_num  = (Nodes)*(qINEQ_trust_ii);
                    
                    rowINEQ_trust      = reshape(repmat((0:qINEQ_trust_ii:(qINEQ_trust_ii*(Nodes-1))),2*qINEQ_trust_ii,1) + repmat((0:n.trustRegionP1s-1)',4,Nodes) + repmat([zeros(qINEQ_trust_ii,1);n.trustRegionP1s*ones(qINEQ_trust_ii,1)],1,Nodes) +obj.MIOFF,[],1);
                    
                    colINEQ_trust      = reshape(repmat([VarIndex;index.trustRegionP1s],2,1),[],1);
                    valINEQ_trust      = reshape(repmat([ones(nVar,1);-ones(n.trustRegionP1s,1);-ones(nVar,1);-ones(n.trustRegionP1s,1)],1,Nodes),[],1);
                    MINEQ_trustP1        = sparse(rowINEQ_trust,colINEQ_trust,valINEQ_trust,rowINEQ_trust_num,Ntot);
                    pINEQ_trustP1        = reshape([Var_n;-Var_n],[],1);
                    
                    
                    rowEQ_trust = reshape(repmat((0:((Nodes-1))),n.trustRegionP1s+1,1)  +obj.MIOFF,[],1);
                    colEQ_trust = reshape([index.trustRegionP1s; index.trustRegionP1],[],1);
                    valEQ_trust = reshape(repmat([ones(n.trustRegionP1s,1); -1],1,Nodes),[],1);
                    
                    MEQ_trustP1   = sparse(rowEQ_trust,colEQ_trust,valEQ_trust,Nodes,Ntot);
                    pEQ_trustP1   = zeros(Nodes,1);
                    GSOC_trustP1 = sparse(0,Ntot);
                    hSOC_trustP1 = zeros(0,1);
                    qSOC_trustP1 = zeros(1,0);
                    
                case 2 % norm 2 of trust region variables per node
                    qSOC_trust_ii       = 1 + nVar;
                    qSOC_trustP1          = qSOC_trust_ii*ones(1,Nodes);
                    
                    rowSOC_trust        = reshape(repmat((0:qSOC_trust_ii:(qSOC_trust_ii*(Nodes-1))),qSOC_trust_ii,1) + repmat((0:(qSOC_trust_ii-1))',1,Nodes)+obj.MIOFF,[],1);
                    colSOC_trust        = reshape([index.trustRegionP1;VarIndex],[],1);
                    valSOC_trust        = -ones(Nodes*(qSOC_trust_ii),1);
                    GSOC_trustP1          = sparse(rowSOC_trust,colSOC_trust,valSOC_trust,sum(qSOC_trustP1),Ntot);
                    hSOC_trustP1          = reshape([zeros(1,Nodes);-Var_n],[],1);
                    
                    MINEQ_trustP1 = sparse(0,Ntot);
                    pINEQ_trustP1 = zeros(0,1);
                    MEQ_trustP1 = sparse(0,Ntot);
                    pEQ_trustP1 = zeros(0,1);
                case {inf,Inf} % maximum absolutes trust region variables per node
                    qINEQ_trust_ii     = n.trustRegionP1s*2;
                    
                    rowINEQ_trust_num  = (Nodes)*(qINEQ_trust_ii);
                    
                    rowINEQ_trust      = reshape(repmat((0:qINEQ_trust_ii:(qINEQ_trust_ii*(Nodes-1))),2*qINEQ_trust_ii,1) + repmat((0:n.trustRegionP1s-1)',4,1) + [zeros(qINEQ_trust_ii,1);n.trustRegionP1s*ones(qINEQ_trust_ii,1)] +obj.MIOFF,[],1);
                    colINEQ_trust      = reshape(repmat([index.VarIndex;index.trustRegionP1s],2,1),[],1);
                    valINEQ_trust      = reshape(repmat([ones(nVar,1);-ones(n.trustRegionP1s,1);-ones(nVar,1);-ones(n.trustRegionP1s,1)],1,Nodes),[],1);
                    
                    rowINEQ_trust      = [rowINEQ_trust;((Nodes)*(n.trustRegionP1s*2)+reshape(repmat((0:n.trustRegionP1s:(n.trustRegionP1s*(Nodes-1))),2*n.trustRegionP1s,1) + repmat((0:n.trustRegionP1s-1)',2,1) + obj.MIOFF,[],1))];
                    colINEQ_trust      = [colINEQ_trust;reshape([index.trustRegionP1s; repmat(index.trustRegionP1,n.trustRegionP1s,1)],[],1)];
                    valINEQ_trust      = [valINEQ_trust;reshape(repmat([ones(n.trustRegionP1s,1); -ones(n.trustRegionP1s,1)],1,Nodes),[],1)];
                    
                    MINEQ_trustP1        = sparse(rowINEQ_trust,colINEQ_trust,valINEQ_trust,(Nodes)*(n.trustRegionP1s*3),Ntot);
                    pINEQ_trustP1        = [reshape([Var_n;-Var_n],[],1);zeros(rowINEQ_trust_num,1)];
                    
                    MEQ_trustP1 = sparse(0,Ntot);
                    pEQ_trustP1 = zeros(0,1);
                    GSOC_trustP1 = sparse(0,Ntot);
                    hSOC_trustP1 = zeros(0,1);
                    qSOC_trustP1 = zeros(1,0);
            end
            
            switch obj.algorithm.sequential.trustRegion.nodePenalty
                case 0
                    error('SCOPT Error: Quadratic trust region for nodes not implemented \n')
                case 1 % Already Absolute values, so equality constraints only necessary. On all nodes
                    MEQ_trustP2      = sparse(1,[index.trustRegionP1(:);index.trustRegionP2],[ones(Nodes,1);-1],1,Ntot);
                    pEQ_trustP2      = zeros(1,1);
                    MINEQ_trustP2    = sparse(0,Ntot);
                    pINEQ_trustP2    = zeros(0,1);
                    GSOC_trustP2     = sparse(0,Ntot);
                    hSOC_trustP2     = zeros(0,1);
                    qSOC_trustP2     = zeros(1,0);
                case 2  % norm 2  trust region applied on all node
                    qSOC_trustP2       = Nodes+1;
                    
                    rowSOC_trust       = (0:Nodes)'+obj.MIOFF;
                    colSOC_trust       = [index.trustRegionP2;index.trustRegionP1'];
                    valSOC_trust       = -ones(Nodes+1,1);
                    GSOC_trustP2       = sparse(rowSOC_trust,colSOC_trust,valSOC_trust,Nodes+1,Ntot);
                    hSOC_trustP2       = zeros(Nodes+1,1);
                    MEQ_trustP2      = sparse(0,Ntot);
                    pEQ_trustP2      = zeros(0,1);
                    MINEQ_trustP2     = sparse(0,Ntot);
                    pINEQ_trustP2     = zeros(0,1);
                case {inf,Inf} % maximum absolutes trust region applied on all node
                    MINEQ_trustP2    = sparse(reshape(repmat((0:Nodes-1),2,1)+obj.MIOFF,[],1),reshape([index.trustRegionP1;repmat(index.trustRegionP2,1,Nodes)],[],1),repmat([1;-1],Nodes,1),Nodes,Ntot);
                    pINEQ_trustP2    = zeros(Nodes,1);
                    MEQ_trustP2      = sparse(0,Ntot);
                    pEQ_trustP2      = zeros(0,1);
                    GSOC_trustP2     = sparse(0,Ntot);
                    hSOC_trustP2     = zeros(0,1);
                    qSOC_trustP2     = zeros(1,0);
            end
            if strcmp(obj.algorithm.sequential.trustRegion.type,'hard')
                % if hard trust region constraints internal slack
                % variable to be lower than radius term
                MINEQ_trustType      = sparse(1,index.trustRegionP2,1,1,Ntot);
                pINEQ_trustType      = radius;
            else
                MINEQ_trustType      = sparse(0,Ntot);
                pINEQ_trustType      = zeros(0,1);
            end
            MEQ_trust = [MEQ_trustP1;MEQ_trustP2];
            pEQ_trust = [pEQ_trustP1;pEQ_trustP2];
            MINEQ_trust = [MINEQ_trustP1;MINEQ_trustP2;MINEQ_trustType];
            pINEQ_trust = [pINEQ_trustP1;pINEQ_trustP2;pINEQ_trustType];
            GSOC_trust = [GSOC_trustP1;GSOC_trustP2];
            hSOC_trust = [hSOC_trustP1;hSOC_trustP2];
            qSOC_trust = [qSOC_trustP1,qSOC_trustP2];
        case {'none'}
            MEQ_trust = sparse(0,Ntot);
            pEQ_trust = zeros(0,1);
            MINEQ_trust = sparse(0,Ntot);
            pINEQ_trust = zeros(0,1);
            GSOC_trust = sparse(0,Ntot);
            hSOC_trust = zeros(0,1);
            qSOC_trust = zeros(1,0);
        otherwise
            MEQ_trust = sparse(0,Ntot);
            pEQ_trust = zeros(0,1);
            MINEQ_trust = sparse(0,Ntot);
            pINEQ_trust = zeros(0,1);
            GSOC_trust = sparse(0,Ntot);
            hSOC_trust = zeros(0,1);
            qSOC_trust = zeros(1,0);
    end
    
    % Trust region for time of flight
    if  obj.problem.phases(1).freeTimeFinal
        switch obj.algorithm.sequential.trustRegion.phases(1).timeFinal.type
            case 'hard'
                MINEQ_trustTimeFinal  = sparse((0:1)+obj.MIOFF,[index.timeFinal,index.timeFinal],[1 ,-1],2,Ntot);
                pINEQ_trustTimeFinal  = [timeFinal_n;-timeFinal_n]+obj.algorithm.sequential.trustRegion.phases(1).timeFinal.radius;
            case 'soft'
                MINEQ_trustTimeFinal  = sparse([0,0,1,1]+obj.MIOFF,[index.timeFinal,index.timeFinalSlack,index.timeFinal,index.timeFinalSlack],[1 , -1,-1, -1],2,Ntot);
                pINEQ_trustTimeFinal  = [timeFinal_n;-timeFinal_n];
            case 'none'
                MINEQ_trustTimeFinal  = sparse(0,Ntot);
                pINEQ_trustTimeFinal  = zeros(0,1);
        end
    else
        MINEQ_trustTimeFinal  = sparse(0,Ntot);
        pINEQ_trustTimeFinal  = zeros(0,1);
    end
    
    MINEQ_trust = [MINEQ_trust ; MINEQ_trustTimeFinal];
    pINEQ_trust = [pINEQ_trust ; pINEQ_trustTimeFinal];
    
else
    MEQ_trust = sparse(0,Ntot);
    pEQ_trust = zeros(0,1);
    MINEQ_trust = sparse(0,Ntot);
    pINEQ_trust = zeros(0,1);
    GSOC_trust = sparse(0,Ntot);
    hSOC_trust = zeros(0,1);
    qSOC_trust = zeros(1,0);
end
obj.algorithm.sequential.trustRegion.phases(1).MEQ    = MEQ_trust;
obj.algorithm.sequential.trustRegion.phases(1).pEQ    = pEQ_trust;
obj.algorithm.sequential.trustRegion.phases(1).MINEQ  = MINEQ_trust;
obj.algorithm.sequential.trustRegion.phases(1).pINEQ  = pINEQ_trust;
obj.algorithm.sequential.trustRegion.phases(1).GSOC   = GSOC_trust;
obj.algorithm.sequential.trustRegion.phases(1).hSOC   = hSOC_trust;
obj.algorithm.sequential.trustRegion.phases(1).qSOC   = qSOC_trust;
%% Objective Function
f = zeros(Ntot,1);
if strcmp(Objective.type,'feasibility') % No objective, feasibility problem
    GSOC_obj    = sparse(0,Ntot);
    hSOC_obj    = zeros(0,1);
    qSOC_obj    = zeros(1,0);
    MEQ_obj     = sparse(0,Ntot);
    pEQ_obj     = zeros(0,1);
    MINEQ_obj   = sparse(0,Ntot);
    pINEQ_obj   = zeros(0,1);
else
    f(index.objectiveSlack) = Objective.sign;
    % Mayer Term
    switch Objective.mayer.where
        case 'final'
            indexWhere = Nodes - 1 +obj.MIOFF;
        case 'initial'
            indexWhere = 0 +obj.MIOFF;
        case 'index' % TODO
            indexWhere = Objective.mayer.indexWhere;
    end
    switch Objective.mayer.funType
        case 'linear'
            colObj    = [index.states(:,indexWhere);index.controls(:,indexWhere);index.objectiveSlack];
            valObj    = [Scale.states.*Objective.mayer.states/Objective.scale;Scale.controls.*Objective.mayer.controls/Objective.scale;-1];
            
            pEQ_obj   = ( - sum(Shift.states.*Objective.mayer.states,1) - sum(Shift.controls.*Objective.mayer.controls,1))/Objective.scale;
            
            if obj.problem.phases(1).freeTimeFinal
                colObj  = [colObj ; index.timeFinal];
                valObj  = [valObj ; Scale.time*Objective.mayer.timeFinal/Objective.scale];
                pEQ_obj = pEQ_obj - Shift.time*Objective.mayer.timeFinal/Objective.scale;
            end
            
            MEQ_obj   = sparse(0+obj.MIOFF,colObj,valObj,1,Ntot);
            MINEQ_obj = sparse(0,Ntot);
            pINEQ_obj = zeros(0,1);
            GSOC_obj  = sparse(0,Ntot);
            hSOC_obj  = zeros(0,1);
            qSOC_obj  = zeros(1,0);
        case {'convex','soc'}
            qSOC_obj_ii     = Objective.mayer.cone.dimensions;
            qSOC_obj        = qSOC_obj_ii*ones(1,1);
            rowSOC_obj      = [0;zeros(qSOC_obj_ii*(n.states+n.controls),1)] + [0;reshape(repmat((0:qSOC_obj_ii-1),n.states+n.controls,1),[],1)]+obj.MIOFF;
            colSOC_obj      = [index.objectiveSlack;repmat([index.states(:,indexWhere);index.controls(:,indexWhere)],qSOC_obj_ii,1)];
            valSOC_obj      = -[1;Scale.states(:).*Objective.mayer.cone.right.states(:)/Objective.scale;Scale.controls(:).*Objective.mayer.cone.right.controls(:)/Objective.scale;...
                reshape([Scale.states(:).*Objective.mayer.cone.norm.states;Scale.controls(:).*Objective.mayer.cone.norm.controls],[],1)/Objective.scale];
            GSOC_obj        = sparse(rowSOC_obj,colSOC_obj,valSOC_obj,sum(qSOC_obj),Ntot);
            hSOC_obj        = [(Objective.mayer.cone.right.cons + sum(Shift.controls.*Objective.mayer.cone.right.controls,1) + sum(Shift.states.*Objective.mayer.cone.right.states,1));...
                (Objective.mayer.cone.norm.cons(:) + reshape(sum(Shift.controls.*Objective.mayer.cone.norm.controls,1) + sum(Shift.states.*Objective.mayer.cone.norm.states,1),[],1))]/Objective.scale;
            
            MEQ_obj     = sparse(0,Ntot);
            pEQ_obj     = zeros(0,1);
            MINEQ_obj   = sparse(0,Ntot);
            pINEQ_obj   = zeros(0,1);
        case {'non-linear'}
            f_ii       = Objective.mayer.function(time(:,indexWhere),X(:,indexWhere),U(:,indexWhere),BodyMap,index.nodes(indexWhere));
            dF_dX_ii   = Objective.mayer.jacobian.states(time(:,indexWhere),X(:,indexWhere),U(:,indexWhere),BodyMap,index.nodes(indexWhere));
            dF_dU_ii   = Objective.mayer.jacobian.controls(time(:,indexWhere),X(:,indexWhere),U(:,indexWhere),BodyMap,index.nodes(indexWhere));
            
            dF_dX_n_ii = Scale.states(:).*dF_dX_ii;
            dF_dU_n_ii = Scale.controls(:).*dF_dU_ii;
            
            colObj     = [reshape([index.states(:,indexWhere);index.controls(:,indexWhere)],[],1);index.objectiveSlack];
            valObj     = [reshape([dF_dX_n_ii;dF_dU_n_ii],[],1)/Objective.scale;-1];
            
            MEQ_obj    = sparse(0+obj.MIOFF,colObj,valObj,1,Ntot);
            pEQ_obj    = (-f_ii + dot(dF_dX_ii(:),X(:,indexWhere)) + dot(dF_dU_ii(:),U(:,indexWhere)))/Objective.scale;
            
            MINEQ_obj  = sparse(0,Ntot);
            pINEQ_obj  = zeros(0,1);
            
            GSOC_obj = sparse(0,Ntot);
            hSOC_obj = zeros(0,1);
            qSOC_obj = zeros(1,0);
            obj.problem.objective.MEQnoSlack =  sparse(0+obj.MIOFF,colObj(1:end-1),valObj(1:end-1),1,Ntot);
    end
    
    % Lagrange
    switch Objective.lagrange.funType
        case 'linear'
            if Objective.lagrange.absolute.nodes
                dF_n_dX_n_ii      = Scale.states(:).*repmat(Objective.lagrange.states(:),1,Nodes)/Objective.scale;
                dF_n_dU_n_ii      = Scale.controls(:).*repmat(Objective.lagrange.controls(:),1,Nodes)/Objective.scale;
                
                rowINEQ_obj = reshape(repmat([zeros(n.states,1);zeros(n.controls,1);0;ones(n.states,1);ones(n.controls,1);1],1,Nodes) + repmat((0:2:(Nodes-1)*2),(n.states+n.controls+1)*2,1) +obj.MIOFF,[],1);
                colINEQ_obj = reshape(repmat([index.states;index.controls;index.objectiveLagrange],2,1),[],1);
                valINEQ_obj = reshape([dF_n_dX_n_ii;dF_n_dU_n_ii;-ones(1,Nodes);-dF_n_dX_n_ii;-dF_n_dU_n_ii;-ones(1,Nodes)],[],1);
                
                MINEQ_obj   = sparse(rowINEQ_obj,colINEQ_obj,valINEQ_obj,2*Nodes,Ntot);
                pINEQ_obj   = reshape([-sum(Shift.states.*repmat(Objective.lagrange.states(:),1,Nodes)/Objective.scale,1)-sum(Shift.controls.*repmat(Objective.lagrange.controls(:),1,Nodes)/Objective.scale,1);+sum(Shift.states.*repmat(Objective.lagrange.states(:),1,Nodes)/Objective.scale,1)+sum(Shift.controls.*repmat(Objective.lagrange.controls(:),1,Nodes)/Objective.scale,1)],[],1);
                
                colEQ_obj   =  index.objectiveLagrange(:);
                valEQ_obj   =  timeInterval*obj.problem.phases(1).quadratureVector(:);
                
                pEQ_obj = 0;
                if obj.problem.phases(1).freeTimeFinal
                    dI_n_dtf_n = Scale.time*sum(obj.problem.phases(1).quadratureVector.*abs(sum(Objective.lagrange.controls(:).*U,1) + sum(Objective.lagrange.states(:).*X,1)))/Objective.scale;
                    colEQ_obj  = [colEQ_obj;index.timeFinal];
                    valEQ_obj  = [valEQ_obj;dI_n_dtf_n];
                    pEQ_obj    = pEQ_obj + dI_n_dtf_n*timeFinal_n;
                end
                
                MEQ_obj     =  sparse(1,[colEQ_obj;index.objectiveSlack],[valEQ_obj;-1],1,Ntot);
                obj.problem.objective.MEQnoSlack =  sparse(1,colEQ_obj,valEQ_obj,1,Ntot);
            else
                dI_dX   = timeInterval*repmat(obj.problem.phases(1).quadratureVector,n.states,1).*repmat(Objective.lagrange.states(:),1,Nodes);
                dI_dU   = timeInterval*repmat(obj.problem.phases(1).quadratureVector,n.controls,1).*repmat(Objective.lagrange.controls(:),1,Nodes);
                dI_dX_n = repmat(Scale.states,1,Nodes).*dI_dX;
                dI_dU_n = repmat(Scale.controls,1,Nodes).*dI_dU;
                
                colObj    = reshape([index.states;index.controls],[],1);
                valObj    = reshape([dI_dX_n;dI_dU_n],[],1)/Objective.scale;
                
                pEQ_obj   = -(sum(sum(repmat(Shift.states,1,Nodes).*dI_dX,1)) + sum(sum(repmat(Shift.controls,1,Nodes).*dI_dU,1)))/Objective.scale;
                if obj.problem.phases(1).freeTimeFinal
                    dI_n_dtf_n  = Scale.time*(sum(obj.problem.phases(1).quadratureVector.*(sum(Objective.lagrange.controls(:).*U,1) + sum(Objective.lagrange.states(:).*X,1))))/Objective.scale;
                    colObj    = [colObj ; index.timeFinal];
                    valObj    = [valObj ; dI_n_dtf_n];
                    pEQ_obj   = pEQ_obj + (dI_n_dtf_n*timeFinal_n);
                end
                
                MEQ_obj   = sparse(0+obj.MIOFF,[colObj;index.objectiveSlack],[valObj;-1],1,Ntot);
                obj.problem.objective.MEQnoSlack =  sparse(0+obj.MIOFF,colObj,valObj,1,Ntot);
                
                MINEQ_obj = sparse(0,Ntot);
                pINEQ_obj = zeros(0,1);
            end
            GSOC_obj  = sparse(0,Ntot);
            hSOC_obj  = zeros(0,1);
            qSOC_obj  = zeros(1,0);
        case {'convex','soc'}
            qSOC_obj_ii    = Objective.lagrange.cone.dimensions;
            qSOC_obj       = qSOC_obj_ii*ones(1,Nodes);
            rowSOC_obj     = reshape(repmat([0;zeros(qSOC_obj_ii*(n.states+n.controls),1) + reshape(repmat((0:qSOC_obj_ii-1),n.states+n.controls,1),[],1)],1,Nodes) + repmat((0:qSOC_obj_ii:(Nodes-1)*qSOC_obj_ii),qSOC_obj_ii*(n.states+n.controls)+1,1),[],1) +obj.MIOFF;
            colSOC_obj     = reshape([index.objectiveLagrange;repmat([index.states;index.controls],qSOC_obj_ii,1)],[],1);
            valSOC_obj     = reshape([-ones(1,Nodes);repmat(-[Scale.states(:).*Objective.lagrange.cone.right.states(:);Scale.controls(:).*Objective.lagrange.cone.right.controls(:);reshape([repmat(Scale.states(:),1,Objective.lagrange.cone.dimensions-1).*Objective.lagrange.cone.norm.states;repmat(Scale.controls(:),1,Objective.lagrange.cone.dimensions-1).*Objective.lagrange.cone.norm.controls],[],1)]/Objective.scale,1,Nodes)],[],1);
            
            GSOC_obj       = sparse(rowSOC_obj,colSOC_obj,valSOC_obj,sum(qSOC_obj),Ntot);
            hSOC_obj       = repmat([(Objective.lagrange.cone.right.cons + sum(Shift.controls.*Objective.lagrange.cone.right.controls,1) + sum(Shift.states.*Objective.lagrange.cone.right.states,1)) ; (Objective.lagrange.cone.norm.cons(:) + reshape(sum(repmat(Shift.controls,1,Objective.lagrange.cone.dimensions-1).*Objective.lagrange.cone.norm.controls,1) + sum(repmat(Shift.states,1,Objective.lagrange.cone.dimensions-1).*Objective.lagrange.cone.norm.states,1),[],1))]/Objective.scale,Nodes,1);
            
            colEQ_obj      =  index.objectiveLagrange(:);
            valEQ_obj      =  timeInterval*obj.problem.phases(1).quadratureVector(:);
            
            pEQ_obj = 0;
            if obj.problem.phases(1).freeTimeFinal
                % Norm Computation
                NormTermDim = zeros(qSOC_obj_ii-1,Nodes);
                for dd = 0:qSOC_obj_ii-2
                    NormTermDim(dd+obj.MIOFF,:)=[Objective.lagrange.cone.norm.states(:,dd+obj.MIOFF);Objective.lagrange.cone.norm.controls(:,dd+obj.MIOFF)]'*[X;U] + Objective.lagrange.cone.norm.cons(dd+obj.MIOFF);
                end
                NormTerm = vecnorm(NormTermDim,2,1);
                f_n_ii   = (NormTerm  - Objective.lagrange.cone.right.cons-[Objective.lagrange.cone.right.states(:);Objective.lagrange.cone.right.controls(:)]'*[X;U])/Objective.scale;
                
                dI_n_dtf_n    = Scale.time*sum(obj.problem.phases(1).quadratureVector.*f_n_ii);
                
                colEQ_obj = [colEQ_obj;index.timeFinal];
                valEQ_obj = [valEQ_obj;dI_n_dtf_n];
                pEQ_obj   = pEQ_obj + dI_n_dtf_n*timeFinal_n;
            end
            
            MEQ_obj   = sparse(1,[colEQ_obj;index.objectiveSlack],[valEQ_obj;-1],1,Ntot);
            MINEQ_obj = sparse(0,Ntot);
            pINEQ_obj = zeros(0,1);
            obj.problem.objective.MEQnoSlack =  sparse(0+obj.MIOFF,colEQ_obj,valEQ_obj,1,Ntot);
        case {'non-linear'}
            f_ii           = Objective.lagrange.function(time,X,U,BodyMap,index.nodes);
            dF_dX_ii       = Objective.lagrange.jacobian.states(time,X,U,BodyMap,index.nodes);
            dF_dU_ii       = Objective.lagrange.jacobian.controls(time,X,U,BodyMap,index.nodes);
            
            if Objective.lagrange.absolute.nodes
                dF_n_dX_n_ii      = Scale.states(:).*dF_dX_ii/Objective.scale;
                dF_n_dU_n_ii      = Scale.controls(:).*dF_dU_ii/Objective.scale;
                %                             dF_n_dX_n_ii      = Scale.states(:).*dF_dX_ii./repmat(Objective.scaleAux,n.states,1);
                %                             dF_n_dU_n_ii      = Scale.controls(:).*dF_dU_ii./repmat(Objective.scaleAux,n.controls,1);
                %                             dF_n_dX_X         = accumarray(reshape(repmat((0:Nodes-1) +obj.MIOFF,n.states,1),[],1),dF_n_dX_n_ii(:).*X_n(:));
                %                             dF_n_dU_U         = accumarray(reshape(repmat((0:Nodes-1) +obj.MIOFF,n.controls,1),[],1),dF_n_dU_n_ii(:).*U_n(:));
                dF_n_dX_n_X_n         = sum(dF_n_dX_n_ii.*X_n,1);
                dF_n_dU_n_U_n         = sum(dF_n_dU_n_ii.*U_n,1);
                
                rowINEQ_obj = reshape(repmat([zeros(n.states,1);zeros(n.controls,1);0;ones(n.states,1);ones(n.controls,1);1],1,Nodes) + repmat((0:2:(Nodes-1)*2),(n.states+n.controls+1)*2,1) +obj.MIOFF,[],1);
                colINEQ_obj = reshape(repmat([index.states;index.controls;index.objectiveLagrange],2,1),[],1);
                valINEQ_obj = reshape([dF_n_dX_n_ii;dF_n_dU_n_ii;-ones(1,Nodes);-dF_n_dX_n_ii;-dF_n_dU_n_ii;-ones(1,Nodes)],[],1);
                
                MINEQ_obj   = sparse(rowINEQ_obj,colINEQ_obj,valINEQ_obj,2*Nodes,Ntot);
                pINEQ_obj   = reshape([-f_ii/Objective.scale + dF_n_dX_n_X_n + dF_n_dU_n_U_n;+f_ii/Objective.scale- dF_n_dX_n_X_n - dF_n_dU_n_U_n],[],1);
                %                             pINEQ_obj   = reshape([-f_ii./Objective.scaleAux + dF_n_dX_n_X_n + dF_n_dU_n_U_n;+f_ii./Objective.scaleAux- dF_n_dX_n_X_n - dF_n_dU_n_U_n],[],1);
                
                
                colEQ_obj   =  index.objectiveLagrange(:);
                valEQ_obj   =  timeInterval*obj.problem.phases(1).quadratureVector(:);
                %                             valEQ_obj   =  timeInterval*obj.problem.phases(1).quadratureVector(:).*Objective.scaleAux(:)/Objective.scale;
                
                pEQ_obj = 0;
                if obj.problem.phases(1).freeTimeFinal
                    dI_n_dtf_n = Scale.time*sum(obj.problem.phases(1).quadratureVector.*abs(f_ii))/Objective.scale;
                    colEQ_obj  = [colEQ_obj;index.timeFinal];
                    valEQ_obj  = [valEQ_obj;dI_n_dtf_n];
                    pEQ_obj    = pEQ_obj + dI_n_dtf_n*timeFinal_n;
                end
                
                MEQ_obj     =  sparse(1,[colEQ_obj;index.objectiveSlack],[valEQ_obj;-1],1,Ntot);
                obj.problem.objective.MEQnoSlack =  sparse(0+obj.MIOFF,colEQ_obj,valEQ_obj,1,Ntot);
            else
                I_n          = timeInterval*sum(obj.problem.phases(1).quadratureVector.*f_ii)/Objective.scale;
                dI_dX        = timeInterval*repmat(obj.problem.phases(1).quadratureVector,n.states,1).*(dF_dX_ii);
                dI_dU        = timeInterval*repmat(obj.problem.phases(1).quadratureVector,n.controls,1).*(dF_dU_ii);
                dI_n_dX_n    = Scale.states(:).*dI_dX/Objective.scale;
                dI_n_dU_n    = Scale.controls(:).*dI_dU/Objective.scale;
                
                colEQ_Obj       = reshape([index.states;index.controls],[],1);
                valEQ_obj       = reshape([dI_n_dX_n;dI_n_dU_n],[],1);
                
                %                             pEQ_obj      = -I_n + (dot(dI_dX(:),X(:)) + dot(dI_dU(:),U(:)))/Objective.scale;
                pEQ_obj      = -I_n + (dot(dI_n_dX_n(:),X_n(:)) + dot(dI_n_dU_n(:),U_n(:)));
                if obj.problem.phases(1).freeTimeFinal
                    dI_n_dtf_n = Scale.time*sum(obj.problem.phases(1).quadratureVector.*f_ii)/Objective.scale;
                    colEQ_Obj  = [colEQ_Obj;index.timeFinal];
                    valEQ_obj  = [valEQ_obj;dI_n_dtf_n];
                    pEQ_obj    = pEQ_obj + dI_n_dtf_n*timeFinal_n;
                end
                
                MEQ_obj   = sparse(1,[colEQ_Obj;index.objectiveSlack],[valEQ_obj;-1],1,Ntot);
                
                MINEQ_obj = sparse(0,Ntot);
                pINEQ_obj = zeros(0,1);
                obj.problem.objective.MEQnoSlack =  sparse(0+obj.MIOFF,colEQ_Obj,valEQ_obj,1,Ntot);
            end
            GSOC_obj  = sparse(0,Ntot);
            hSOC_obj  = zeros(0,1);
            qSOC_obj  = zeros(1,0);
    end
    
    
    % Minmax
    switch Objective.minmax.funType
        case 'linear'
            eqNodes   = Objective.minmax.absolute.nodes + 1;
            rowObj    = reshape(repmat((0:1:(Nodes-1)),(n.states+n.controls),1) +obj.MIOFF,[],1);
            colObj    = reshape([index.states;index.controls],[],1);
            valObj    = Objective.sign*repmat([Scale.states.*Objective.minmax.states(:);Scale.controls.*Objective.minmax.controls(:)],Nodes,1)/Objective.scale;
            
            pINEQ_obj = repmat(-Objective.sign*(sum(Shift.states.*Objective.minmax.states(:)) + sum(Shift.controls.*Objective.minmax.controls(:)))/Objective.scale,Nodes,1);
            
            if Objective.minmax.absolute.nodes
                rowObj    = [rowObj;(rowObj + Nodes)];
                colObj    = [colObj;colObj];
                valObj    = [valObj;-valObj];
                pINEQ_obj = [pINEQ_obj;-pINEQ_obj];
            end
            MINEQ_obj = sparse([rowObj;((0:(Nodes*eqNodes-1))' + obj.MIOFF)],[colObj;repmat(index.objectiveSlack,Nodes*eqNodes,1)],[valObj;-Objective.sign*ones(Nodes*eqNodes,1)],Nodes*eqNodes,Ntot);
            
            MEQ_obj   = sparse(0,Ntot);
            pEQ_obj   = zeros(0,1);
            GSOC_obj  = sparse(0,Ntot);
            hSOC_obj  = zeros(0,1);
            qSOC_obj  = zeros(1,0);
        case {'convex','soc'}
            %                         error('Convex minmax not implemented yet \n')
            qSOC_obj_ii    = Objective.minmax.cone.dimensions;
            qSOC_obj       = qSOC_obj_ii*ones(1,Nodes);
            rowSOC_obj     = reshape(repmat([0;zeros(qSOC_obj_ii*(n.states+n.controls),1) + reshape(repmat((0:qSOC_obj_ii-1),n.states+n.controls,1),[],1)],1,Nodes) + repmat((0:qSOC_obj_ii:(Nodes-1)*qSOC_obj_ii),(n.states+n.controls)*qSOC_obj_ii+1,1),[],1) +obj.MIOFF;
            colSOC_obj     = reshape([repmat(index.objectiveSlack,1,Nodes);repmat([index.states;index.controls],qSOC_obj_ii,1)],[],1);
            valSOC_obj     = reshape([-ones(1,Nodes);repmat(-[Scale.states(:).*Objective.minmax.cone.right.states(:);Scale.controls(:).*Objective.minmax.cone.right.controls(:);reshape([Scale.states(:).*Objective.minmax.cone.norm.states;Scale.controls(:).*Objective.minmax.cone.norm.controls],[],1)]/Objective.scale,1,Nodes)],[],1);
            
            GSOC_obj       = sparse(rowSOC_obj,colSOC_obj,valSOC_obj,sum(qSOC_obj),Ntot);
            hSOC_obj       = repmat([(Objective.minmax.cone.right.cons + sum(Shift.controls.*Objective.minmax.cone.right.controls,1) + sum(Shift.states.*Objective.minmax.cone.right.states,1)) ; (Objective.minmax.cone.norm.cons(:) + reshape(sum(Shift.controls.*Objective.minmax.cone.norm.controls,1) + sum(Shift.states.*Objective.minmax.cone.norm.states,1),[],1))]/Objective.scale,Nodes,1);
            
            MEQ_obj     = sparse(0,Ntot);
            pEQ_obj     = zeros(0,1);
            MINEQ_obj   = sparse(0,Ntot);
            pINEQ_obj   = zeros(0,1);
        case {'non-linear'}
            f_ii           = Objective.minmax.function(time,X,U,BodyMap,index.nodes);
            dF_dX_ii       = Objective.minmax.jacobian.states(time,X,U,BodyMap,index.nodes);
            dF_dU_ii       = Objective.minmax.jacobian.controls(time,X,U,BodyMap,index.nodes);
            dF_n_dX_n_ii   = Scale.states(:).*dF_dX_ii/Objective.scale;
            dF_n_dU_n_ii   = Scale.controls(:).*dF_dU_ii/Objective.scale;
            dF_n_dX_n_X_n  = sum(dF_n_dX_n_ii.*X_n,1);
            dF_n_dU_n_U_n  = sum(dF_n_dU_n_ii.*U_n,1);
            
            eqNodes   = Objective.minmax.absolute.nodes + 1;
            
            rowObj    = reshape(repmat((0:1:(Nodes-1)),(n.states+n.controls),1) +obj.MIOFF,[],1);
            colObj    = reshape([index.states;index.controls],[],1);
            valObj    = Objective.sign*reshape([dF_n_dX_n_ii;dF_n_dU_n_ii],[],1);
            
            pINEQ_obj = reshape(Objective.sign*(-f_ii/Objective.scale + dF_n_dX_n_X_n + dF_n_dU_n_U_n),[],1);
            
            %                         obj.problem.objective.MINEQnoSlack =  sparse(rowObj,colObj,valObj,1,Ntot);
            
            if Objective.minmax.absolute.nodes
                rowObj    = [rowObj;(rowObj + Nodes)];
                colObj    = [colObj;colObj];
                valObj    = [valObj;-valObj];
                pINEQ_obj = [pINEQ_obj;-pINEQ_obj];
            end
            MINEQ_obj = sparse([rowObj;((0:(Nodes*eqNodes-1))' + obj.MIOFF)],[colObj;repmat(index.objectiveSlack,Nodes*eqNodes,1)],[valObj;-Objective.sign*ones(Nodes*eqNodes,1)],Nodes*eqNodes,Ntot);
            
            MEQ_obj   = sparse(0,Ntot);
            pEQ_obj   = zeros(0,1);
            GSOC_obj  = sparse(0,Ntot);
            hSOC_obj  = zeros(0,1);
            qSOC_obj  = zeros(1,0);
    end
end
fnp = f; % Non-penalized Cost Function

% Regularisation terms
f(index.states)    = f(index.states) + Objective.regularisation.weight*repmat(Objective.regularisation.states(:),1,Nodes);
f(index.controls)  = reshape(f(index.controls),size(index.controls)) + Objective.regularisation.weight*repmat(Objective.regularisation.controls(:),1,Nodes);
f(index.timeFinal) = f(index.timeFinal) + Objective.regularisation.weight*Objective.regularisation.timeFinal;

% If soft trust region for final time, aguments trust region slack variable
if obj.problem.phases(1).freeTimeFinal && strcmp(obj.algorithm.sequential.trustRegion.phases(1).timeFinal.type,'soft')
    f(index.timeFinalSlack)    = obj.algorithm.sequential.trustRegion.phases(1).timeFinal.lambda;
end
% If virtual controls, augments penalization term
if VirtualControl.phases(1).include
    f(index.virtualControlsP2) = VirtualControl.phases(1).lambda;
end
% if virtual buffers, augments of all internal slack variables
if VirtualBuffer.phases(1).include
    f(index.virtualBuffers.eventsP) = VirtualBuffer.phases(1).events.lambda(:);
    f(index.virtualBuffers.pathP) = VirtualBuffer.phases(1).path.lambda(:);
end
% If soft trust region for states and controls, aguments trust region slack variable
if obj.algorithm.sequential.activate && (strcmp(obj.algorithm.sequential.type,'trust-region') && strcmp(obj.algorithm.sequential.trustRegion.type,'soft'))
    f(index.trustRegionP2)    = lambda;
end
obj.problem.objective.f              =  f;
obj.problem.objective.MEQ            =  MEQ_obj;
obj.problem.objective.pEQ            =  pEQ_obj;
obj.problem.objective.MINEQ          =  MINEQ_obj;
obj.problem.objective.pINEQ          =  pINEQ_obj;
obj.problem.objective.GSOC           =  GSOC_obj;
obj.problem.objective.hSOC           =  hSOC_obj;
obj.problem.objective.qSOC           =  qSOC_obj;
%% Virtual Control and Buffer Zones Slack Variables
% for component wise penalties
if VirtualControl.phases(1).include
    if VirtualControl.phases(1).statePenalty == 0
        % Uses quadratic virtual controls
        qSOC_vconP1_ii     = 2 + n.virtualControls;
        qSOC_vconP1        = qSOC_vconP1_ii*ones(1,Nodes-1);
        
        rowSOC_vconP1      = reshape(repmat((0:qSOC_vconP1_ii:(qSOC_vconP1_ii*(Nodes-2))),1+1+n.virtualControls,1) + [0;1;(2:qSOC_vconP1_ii-1)']+obj.MIOFF,[],1);
        colSOC_vconP1      = reshape([index.virtualControlsP1;index.virtualControlsP1;index.virtualControls],[],1);
        valSOC_vconP1      = -reshape([1/2*ones(1,Nodes-1);-1/2*ones(1,Nodes-1);ones(n.virtualControls,Nodes-1)],[],1);
        GSOC_vconP1        = sparse(rowSOC_vconP1,colSOC_vconP1,valSOC_vconP1,sum(qSOC_vconP1),Ntot);
        hSOC_vconP1        = reshape(zeros(qSOC_vconP1_ii,Nodes-1),[],1);
        
        MINEQ_vconP1 = sparse(0,Ntot);
        pINEQ_vconP1 = zeros(0,1);
        MEQ_vconP1 = sparse(0,Ntot);
        pEQ_vconP1 = zeros(0,1);
    elseif VirtualControl.phases(1).statePenalty == 1
        % Represents an Absolute Value
        % E*vc - P1slack <=0,-E*vc-P1slack<=0,I*P1slack-P1=0 for each node
        % E*vc <= P1slack and -E*vc <= P1slack
        
        qineq_vconP1_ii     = n.virtualControlsP1s*2; % 3 for 1st and 3 for 2nd
        
        rowINEQ_vconP1_num  = (Nodes-1)*(qineq_vconP1_ii);
        
        rowINEQ_vconP1      = reshape(repmat((0:qineq_vconP1_ii:(qineq_vconP1_ii*(Nodes-2))),2*qineq_vconP1_ii,1) + repmat((0:n.virtualControlsP1s-1)',4,Nodes-1) + repmat([zeros(qineq_vconP1_ii,1);n.virtualControlsP1s*ones(qineq_vconP1_ii,1)],1,Nodes-1) +obj.MIOFF,[],1);
        
        colINEQ_vconP1      = reshape(repmat([index.virtualControls;index.virtualControlsP1s],2,1),[],1);
        valINEQ_vconP1      = reshape(repmat([VirtualControl.phases(1).valE;-ones(n.virtualControlsP1s,1);-VirtualControl.phases(1).valE;-ones(n.virtualControlsP1s,1)],1,Nodes-1),[],1);
        MINEQ_vconP1        = sparse(rowINEQ_vconP1,colINEQ_vconP1,valINEQ_vconP1,rowINEQ_vconP1_num,Ntot);
        pINEQ_vconP1        = zeros(rowINEQ_vconP1_num,1);
        
        
        rowEQ_vconP1       = reshape(repmat((0:(1*(Nodes-2))),n.virtualControlsP1s+1,1)  +obj.MIOFF,[],1);
        colEQ_vconP1       = reshape([index.virtualControlsP1s; index.virtualControlsP1],[],1);
        valEQ_vconP1       = reshape(repmat([ones(n.virtualControlsP1s,1); -1],1,Nodes-1),[],1);
        
        MEQ_vconP1   = sparse(rowEQ_vconP1,colEQ_vconP1,valEQ_vconP1,Nodes-1,Ntot);
        pEQ_vconP1   = zeros(Nodes-1,1);
        
        GSOC_vconP1 = sparse(0,Ntot);
        hSOC_vconP1 = zeros(0,1);
        qSOC_vconP1 = zeros(1,0);
    elseif VirtualControl.phases(1).statePenalty == 2
        qSOC_vconP1_ii      = (n.virtualControls + 1);
        qSOC_vconP1         = qSOC_vconP1_ii*ones(1,Nodes-1);
        
        rowSOC_vconP1_num   = (Nodes-1)*(qSOC_vconP1_ii);
        
        rowSOC_vconP1       = reshape(repmat((0:qSOC_vconP1_ii:(qSOC_vconP1_ii*(Nodes-2))),qSOC_vconP1_ii,1) + repmat([0;1 + (0:n.virtualControls-1)'],1,Nodes-1)+obj.MIOFF,[],1);
        colSOC_vconP1       = reshape([index.virtualControlsP1; index.virtualControls],[],1);
        valSOC_vconP1       = -reshape(repmat([1; VirtualControl.phases(1).valE],1,Nodes-1),[],1);
        GSOC_vconP1         = sparse(rowSOC_vconP1,colSOC_vconP1,valSOC_vconP1,rowSOC_vconP1_num,Ntot);
        hSOC_vconP1         = zeros(rowSOC_vconP1_num,1);
        
        MINEQ_vconP1        = sparse(0,Ntot);
        pINEQ_vconP1        = zeros(0,1);
        MEQ_vconP1        = sparse(0,Ntot);
        pEQ_vconP1        = zeros(0,1);
    else
        error('SCOPT Error: not defined for Virtual Control State Penalty of norm-%0.0 \n',VirtualControl.phases(1).statePenalty)
        
    end
    
    % For nodewise penalties
    if VirtualControl.phases(1).nodePenalty == 1
        MEQ_vconP2    = sparse(1,[index.virtualControlsP1(:);index.virtualControlsP2],[ones(Nodes-1,1);-1],1,Ntot);
        pEQ_vconP2    = zeros(1,1);
        
        GSOC_vconP2   = sparse(0,Ntot);
        hSOC_vconP2   = zeros(0,1);
        qSOC_vconP2   = zeros(1,0);
        MINEQ_vconP2  = sparse(0,Ntot);
        pINEQ_vconP2  = zeros(0,1);
    elseif VirtualControl.phases(1).nodePenalty == 2
        
        qSOC_vconP2         = (Nodes -1 + 1);
        
        rowSOC_vconP2       = (0:(qSOC_vconP2-1))+obj.MIOFF;
        colSOC_vconP2       = [index.virtualControlsP2; index.virtualControlsP1(:)];
        valSOC_vconP2       = -ones(qSOC_vconP2,1);
        GSOC_vconP2         = sparse(rowSOC_vconP2,colSOC_vconP2,valSOC_vconP2,qSOC_vconP2,Ntot);
        hSOC_vconP2         = zeros(qSOC_vconP2,1);
        
        MEQ_vconP2  = sparse(0,Ntot);
        pEQ_vconP2  = zeros(0,1);
        MINEQ_vconP2  = sparse(0,Ntot);
        pINEQ_vconP2  = zeros(0,1);
    elseif VirtualControl.phases(1).nodePenalty == inf
        MINEQ_vconP2    = sparse(reshape(repmat((0:Nodes-2),2,1)+obj.MIOFF,[],1),reshape([index.virtualControlsP1;repmat(index.virtualControlsP2,1,Nodes-1)],[],1),repmat([1;-1],Nodes-1,1),Nodes-1,Ntot);
        pINEQ_vconP2    = zeros(Nodes-1,1);
        
        MEQ_vconP2    = sparse(0,Ntot);
        pEQ_vconP2    = zeros(0,1);
        GSOC_vconP2   = sparse(0,Ntot);
        hSOC_vconP2   = zeros(0,1);
        qSOC_vconP2   = zeros(1,0);
    else
        error('SCOPT Error: not defined for Virtual Control Node Penalty of norm-%0.0 \n',VirtualControl.nodePenalty)
        
    end
    
    GSOC_vcon   = [GSOC_vconP1;GSOC_vconP2];
    hSOC_vcon   = [hSOC_vconP1;hSOC_vconP2];
    qSOC_vcon   = [qSOC_vconP1,qSOC_vconP2];
    MINEQ_vcon  = [MINEQ_vconP1;MINEQ_vconP2];
    pINEQ_vcon  = [pINEQ_vconP1;pINEQ_vconP2];
    MEQ_vcon    = [MEQ_vconP1;MEQ_vconP2];
    pEQ_vcon    = [pEQ_vconP1;pEQ_vconP2];
else
    GSOC_vcon   = sparse(0,Ntot);
    hSOC_vcon   = zeros(0,1);
    qSOC_vcon   = zeros(1,0);
    MINEQ_vcon  = sparse(0,Ntot);
    pINEQ_vcon  = zeros(0,1);
    MEQ_vcon  = sparse(0,Ntot);
    pEQ_vcon  = zeros(0,1);
end

% for virtual buffer zones
if VirtualBuffer.phases(1).include
    % Only Node Penalty 1 enabled
    % Represents an Absolute Value
    % vb - vbPs <=0,-vbPs<=0,I*vbPs-vbP=0 for each node
    % vb - vbP <=0,-vbP<=0 for Vertical Speed Constraint
    
    % Node Penalty 2
    % vb - vbPs <=0,-vbPs<=0,norm(vbPs)<=vbP for each node
    % Events
    rowINEQ_vbuf_events_each = 2;
    rowINEQ_vbuf_events_num  = rowINEQ_vbuf_events_each*n.virtualBuffers.events;
    rowINEQ_vbuf_events      = reshape(repmat([0;0;1],1,n.virtualBuffers.events) + repmat((0:rowINEQ_vbuf_events_each:rowINEQ_vbuf_events_each*(n.virtualBuffers.events-1)),3,1) +obj.MIOFF,[],1);
    colINEQ_vbuf_events      = reshape([index.virtualBuffers.events(:)';repmat(index.virtualBuffers.eventsP(:)',2,1)],[],1);
    valINEQ_vbuf_events      = reshape(repmat([1;-1;-1],1,n.virtualBuffers.events),[],1);
    MINEQ_vbuf_events        = sparse(rowINEQ_vbuf_events,colINEQ_vbuf_events,valINEQ_vbuf_events,rowINEQ_vbuf_events_num,Ntot);
    pINEQ_vbuf_events        = zeros(rowINEQ_vbuf_events_num,1);
    
    % Path
    MINEQ_vbuf_path = [];
    pINEQ_vbuf_path = zeros(2*Nodes*n.virtualBuffers.path,1);
    for ii = 0:n.virtualBuffers.path-1
        rowINEQ_vbuf_path_each = 2;
        rowINEQ_vbuf_path_num  = rowINEQ_vbuf_path_each*Nodes;
        rowINEQ_vbuf_path      = reshape(repmat([0;0;1],1,Nodes) + repmat((0:rowINEQ_vbuf_events_each:rowINEQ_vbuf_events_each*(Nodes-1)),3,1) +obj.MIOFF,[],1);
        colINEQ_vbuf_path      = reshape([index.virtualBuffers.path(ii+obj.MIOFF,:);repmat(index.virtualBuffers.pathPs(ii+obj.MIOFF,:),2,1)],[],1);
        valINEQ_vbuf_path      = reshape(repmat([1;-1;-1],1,Nodes),[],1);
        MINEQ_vbuf_path        = [MINEQ_vbuf_path ; sparse(rowINEQ_vbuf_path,colINEQ_vbuf_path,valINEQ_vbuf_path,rowINEQ_vbuf_path_num,Ntot)];
    end
    
    % Implements Norm 1 over all Nodes
    rowEQ_vbuf        = reshape(repmat((0:n.virtualBuffers.path-1)' +obj.MIOFF,1,Nodes+1),[],1);
    colEQ_vbuf        = reshape([index.virtualBuffers.pathPs,index.virtualBuffers.pathP(:)],[],1);
    valEQ_vbuf        = reshape(repmat([ones(1,Nodes),-1],n.virtualBuffers.path,1),[],1);
    MEQ_vbuf          = sparse(rowEQ_vbuf,colEQ_vbuf,valEQ_vbuf,n.virtualBuffers.path,Ntot);
    pEQ_vbuf          = zeros(n.virtualBuffers.path,1);
    
    
    MINEQ_vbuf = [MINEQ_vbuf_events;MINEQ_vbuf_path];
    pINEQ_vbuf = [pINEQ_vbuf_events;pINEQ_vbuf_path];
else
    MINEQ_vbuf= sparse(0,Ntot);
    pINEQ_vbuf= zeros(0,1);
    MEQ_vbuf  = sparse(0,Ntot);
    pEQ_vbuf  = zeros(0,1);
end
%% Assemble all Matrices
MEQ   = [ MEQ_events    ; MEQ_initial   ; MEQ_path   ; MEQ_dyn        ; MEQ_obj          ; MEQ_final    ; MEQ_vcon       ; MEQ_vbuf   ; MEQ_trust   ];
pEQ   = [ pEQ_events    ; pEQ_initial   ; pEQ_path   ; pEQ_dyn        ; pEQ_obj          ; pEQ_final    ; pEQ_vcon       ; pEQ_vbuf   ; pEQ_trust   ];
MINEQ = [ MINEQ_initial ; MINEQ_path ; MINEQ_final    ; MINEQ_events     ; MINEQ_states ; MINEQ_controls ; MINEQ_time ; MINEQ_trust ; MINEQ_obj ; MINEQ_vcon ; MINEQ_vbuf ];
pINEQ = [ pINEQ_initial ; pINEQ_path ; pINEQ_final    ; pINEQ_events     ; pINEQ_states ; pINEQ_controls ; pINEQ_time ; pINEQ_trust ; pINEQ_obj ; pINEQ_vcon ; pINEQ_vbuf ];
GSOC  = [ GSOC_events   ; GSOC_path  ; GSOC_trust     ; GSOC_obj         ; GSOC_vcon    ];
hSOC  = [ hSOC_events   ; hSOC_path  ; hSOC_trust     ; hSOC_obj         ; hSOC_vcon    ];
qSOC  = [ qSOC_events   , qSOC_path  , qSOC_trust     , qSOC_obj         , qSOC_vcon    ];
% if sum([sum(sum(isnan(MEQ))) (sum(isnan(pEQ))) sum(sum(isnan(MINEQ))) (sum(isnan(pINEQ))) sum(sum(isnan(GSOC))) (sum(isnan(hSOC)))]) >= 1
%     fprintf('Error, Nan Value in Affine Matrices \n')
%     return
% end
%% Compute additional matrices
dims.l = size(MINEQ,1);
dims.q = qSOC;
Gaug = [MINEQ;GSOC];
haug = [pINEQ;hSOC];
%% Store Matrices
obj.problem.transcription.MEQ = MEQ;
obj.problem.transcription.pEQ = pEQ;
obj.problem.transcription.dims = dims;
obj.problem.transcription.Gaug = Gaug;
obj.problem.transcription.haug = haug;
obj.problem.transcription.f    = f;
obj.problem.transcription.fnp  = fnp;
end


