function [penaltyDefectsNonLinear,penaltyBuffersNonLinear,penaltyTrustRegion,penaltyTime,penaltyDefectsLinear,penaltyBuffersLinear] = computeAugmentedObjectiveFunctionTerms(obj,time_n,X_n,U_n,time_n_p,X_n_p,U_n_p,linearToo)
% Computes Augmeted Objective Function terms to be used to
% measure the convexification errors
timeFinal_n  = time_n(end);
timeInitial_n= time_n(1);
X            = obj.scaleStatesToReal(1,X_n);
U            = obj.scaleControlsToReal(1,U_n);
timeFinal    = obj.scaleTimeToReal(1,timeFinal_n);
timeInitial  = obj.scaleTimeToReal(1,timeInitial_n);

timeDisc      = (timeFinal - timeInitial)/(obj.problem.phases.Nodes-1);
time          = (timeInitial:timeDisc:timeFinal);

if nargin<8
    linearToo = 1;
end
timeFinal_n_p  = time_n_p(end);
timeInitial_n_p= time_n_p(1);
X_p            = obj.scaleStatesToReal(1,X_n_p);
U_p            = obj.scaleControlsToReal(1,U_n_p);
timeFinal_p    = obj.scaleTimeToReal(1,timeFinal_n_p);
timeInitial_p  = obj.scaleTimeToReal(1,timeInitial_n_p);

timeDisc_p      = (timeFinal_p - timeInitial_p)/(obj.problem.phases.Nodes-1);
time_p          = (timeInitial_p:timeDisc_p:timeFinal_p);
%
index = obj.problem.phases(1).index;
Nodes = obj.problem.phases(1).Nodes;
n     = obj.problem.phases(1).n;
%%  Virtual Control Dynamics Error Measure
if strcmp(obj.algorithm.sequential.type,'trust-region') %&& obj.algorithm.sequential.trustRegion.adaptive % IF SEQUENTIAL BETTER. THERE IS ALWAYS A LINEARIZATION ERROR
    stateDiff_n  = X_n(:,2:end)-X_n(:,1:end-1);
    
    %                 ScalingVirtual = ones(n.states,1);
    %                 ScalingVirtual(obj.algorithm.virtualControl.phases(1).rowE) = obj.algorithm.virtualControl.phases(1).scale;
    ScalingVirtual = obj.algorithm.virtualControl.phases(1).scale;
    switch obj.algorithm.collocationMethod
        case 'euler'
            dX_n_dt_ii   = repmat((1./obj.problem.phases(1).scale.states(:)),1,Nodes) .* obj.problem.phases(1).dynamics.stateDerivativeFunction(time,X,U,obj.problem.phases(1).BodyMap,index.nodes);
            errorDefectsNonLinear_aux0          = (stateDiff_n - timeDisc*dX_n_dt_ii(:,1:end-1))./ScalingVirtual;
        case 'trapezoidal'
            dX_n_dt_ii   = repmat((1./obj.problem.phases(1).scale.states(:)),1,Nodes) .* obj.problem.phases(1).dynamics.stateDerivativeFunction(time,X,U,obj.problem.phases(1).BodyMap,index.nodes);
            errorDefectsNonLinear_aux0          = (stateDiff_n - timeDisc/2*(dX_n_dt_ii(:,1:Nodes-1)+dX_n_dt_ii(:,2:Nodes)))./ScalingVirtual;
        case 'exact'
            f_n_ii   = repmat((1./obj.problem.phases(1).scale.states(:)),1,Nodes-1) .* obj.problem.phases(1).dynamics.stateDerivativeFunction(time,X,U,obj.problem.phases(1).BodyMap,index.nodes);
            errorDefectsNonLinear_aux0          = (stateDiff_n - f_n_ii)./ScalingVirtual;
        otherwise
            error('SCOPT Error: Collocation Method does not exist')
    end
    
    errorDefectsLinear_aux0 =  reshape(-(obj.problem.phases(1).transcription.dynamics.MEQnoVirtual*obj.solution.Zopt-obj.problem.phases(1).transcription.dynamics.pEQ),obj.problem.phases(1).n.states,obj.problem.phases(1).Nodes-1)./ScalingVirtual;
    
    switch obj.algorithm.sequential.trustRegion.defectsCompErrorPenalty
        case 0
            errorDefectsNonLinear_aux1 = sum(errorDefectsNonLinear_aux0.*errorDefectsNonLinear_aux0,1);
            if linearToo,errorDefectsLinear_aux1    = sum(errorDefectsLinear_aux0.*errorDefectsLinear_aux0,1);end
        case {1,2,inf}
            errorDefectsNonLinear_aux1 = vecnorm(errorDefectsNonLinear_aux0,obj.algorithm.sequential.trustRegion.defectsCompErrorPenalty,1);
            if linearToo,errorDefectsLinear_aux1    = vecnorm(errorDefectsLinear_aux0,obj.algorithm.sequential.trustRegion.defectsCompErrorPenalty,1);end
        otherwise
            error('Error, not defined for adaptive trust region dynamic state penalty of type %0.0 \n',obj.algorithm.sequential.trustRegion.defectsCompErrorPenalty)
    end
    
    switch obj.algorithm.sequential.trustRegion.defectsNodeErrorPenalty
        case {1,2,inf}
            errorDefectsNonLinear  = vecnorm(errorDefectsNonLinear_aux1,obj.algorithm.sequential.trustRegion.defectsNodeErrorPenalty,2);
            if linearToo,errorDefectsLinear     = vecnorm(errorDefectsLinear_aux1,obj.algorithm.sequential.trustRegion.defectsNodeErrorPenalty,2);end
        otherwise
            error('Error, not defined for adaptive trust region dynamics penalty of type %0.0 \n',obj.algorithm.sequential.trustRegion.defectsNodeErrorPenalty)
    end
    penaltyDefectsNonLinear = obj.algorithm.virtualControl.phases(1).lambda*errorDefectsNonLinear;
    if linearToo,penaltyDefectsLinear    = obj.algorithm.virtualControl.phases(1).lambda*errorDefectsLinear;else,penaltyDefectsLinear = 0;end
else
    penaltyDefectsNonLinear = 0;
    penaltyDefectsLinear    = 0;
end
%%  Virtual Buffer Zone Error Measure
if strcmp(obj.algorithm.sequential.type,'trust-region') && obj.algorithm.sequential.trustRegion.include.pathAndEventsErrors && obj.algorithm.sequential.trustRegion.adaptive
    penaltyBufferNonLinearPath_ii    = zeros(1,obj.problem.phases.n.path);
    penaltyBufferLinearPath_ii       = zeros(1,obj.problem.phases.n.path);
    
    for ii = 0:obj.problem.phases.n.path-1
        scalePath   = obj.problem.phases(1).path(ii+obj.MIOFF).scale;
        limit       = obj.problem.phases(1).path(ii+obj.MIOFF).limit;
        switch obj.problem.phases(1).path(ii+obj.MIOFF).funType
            case 'linear'
                if ~strcmp(obj.problem.phases(1).path(ii+obj.MIOFF).type,'equal')
                    if obj.problem.phases(1).path(ii+obj.MIOFF).integral && ~obj.problem.phases(1).path(ii+obj.MIOFF).derivative
                        switch obj.problem.phases(1).path(ii+obj.MIOFF).type
                            case 'lower'
                                typeSign = -1;
                            case 'upper'
                                typeSign = 1;
                            otherwise
                                error('SCOPT Error: type of path constraint does not exist')
                        end
                        penaltyBufferNonLinearPath_ii(ii+obj.MIOFF) = obj.problem.phases(1).path(ii+obj.MIOFF).buffer.lambda*max(0,typeSign*(timeDisc*sum(sum((obj.problem.phases(1).path(ii+obj.MIOFF).states.*X + obj.problem.phases(1).path(ii+obj.MIOFF).controls.*U),1),2) - limit)/scalePath);
                        if linearToo,penaltyBufferLinearPath_ii(ii+obj.MIOFF)    = obj.problem.phases(1).path(ii+obj.MIOFF).buffer.lambda*max(0,typeSign*(tauDisc/2*(timeFinal_p-timeInitial_p)*sum(sum((obj.problem.phases(1).path(ii+obj.MIOFF).states.*X + obj.problem.phases(1).path(ii+obj.MIOFF).controls.*U),1),2) + tauDisc/2*(timeFinal - timeFinal_p)*sum(sum((obj.problem.phases(1).path(ii+obj.MIOFF).states.*X_p + obj.problem.phases(1).path(ii+obj.MIOFF).controls.*U_p),1),2) - limit)/scalePath);end
                    else
                        penaltyBufferNonLinearPath_ii(ii+obj.MIOFF) = 0;
                        penaltyBufferLinearPath_ii(ii+obj.MIOFF)    = 0;
                    end
                end
            case 'non-linear'
                switch obj.problem.phases(1).path(ii+obj.MIOFF).type
                    case 'lower'
                        typeSign = -1;
                    case 'upper'
                        typeSign = 1;
                    otherwise
                        error('SCOPT Error, Path Constraints of constraint type %s not implemented. Use inequality constraints instead \n Aborting',Path(ii+obj.MIOFF).type)
                        
                end
                f       = obj.problem.phases(1).path(ii+obj.MIOFF).function(time,X,U,obj.problem.phases(1).BodyMap,index.nodes);
                if linearToo
                    f_ii        = obj.problem.phases(1).path(ii+obj.MIOFF).function(time_p,X_p,U_p,obj.problem.phases(1).BodyMap_prev,index.nodes);
                    dG_dX_ii    = obj.problem.phases(1).path(ii+obj.MIOFF).jacobian.states(time_p,X_p,U_p,obj.problem.phases(1).BodyMap_prev,index.nodes);
                    dG_dU_ii    = obj.problem.phases(1).path(ii+obj.MIOFF).jacobian.controls(time_p,X_p,U_p,obj.problem.phases(1).BodyMap_prev,index.nodes);
                end
                if obj.problem.phases(1).path(ii+obj.MIOFF).integral && ~obj.problem.phases(1).path(ii+obj.MIOFF).derivative
                    penaltyBufferNonLinearPath_ii(ii+obj.MIOFF) = obj.problem.phases(1).path(ii+obj.MIOFF).buffer.lambda*max(0,typeSign*(timeDisc*sum(f,2) - limit)./scalePath);
                    if linearToo,penaltyBufferLinearPath_ii(ii+obj.MIOFF)    = obj.problem.phases(1).path(ii+obj.MIOFF).buffer.lambda*max(0,typeSign*(tauDisc/2*(timeFinal_p-timeInitial_p)*sum(f_ii + sum(dG_dX_ii.*(X - X_p),1) + sum(dG_dU_ii.*(U - U_p),1),2) + tauDisc/2*(timeFinal - timeFinal_p)*sum(f_ii,2) - limit)./scalePath);end
                elseif ~obj.problem.phases(1).path(ii+obj.MIOFF).integral && obj.problem.phases(1).path(ii+obj.MIOFF).derivative
                    penaltyBufferPathNonLinear_max  = max(0,typeSign*(f(2:end) - f(1:end-1) - limit*timeDisc)./scalePath);
                    penaltyBufferPathNonLinear_node = norm(penaltyBufferPathNonLinear_max,obj.problem.phases(1).path(ii+obj.MIOFF).buffer.penaltyAdaptive);
                    if linearToo,penaltyBufferPathLinear_max     = max(0,typeSign*(f_ii(2:end) + sum(dG_dX_ii(:,2:end).*(X(:,2:end) - X_p(:,2:end)),1) + sum(dG_dU_ii(:,2:end).*(U(:,2:end) - U_p(:,2:end)),1) - f_ii(1:end-1) - sum(dG_dX_ii(:,1:end-1).*(X(:,1:end-1) - X_p(:,1:end-1)),1) - sum(dG_dU_ii(:,1:end-1).*(U(:,1:end-1) - U_p(:,1:end-1)),1) - limit*timeDisc)./scalePath);end
                    if linearToo,penaltyBufferPathLinear_node    = norm(penaltyBufferPathLinear_max,obj.problem.phases(1).path(ii+obj.MIOFF).buffer.penaltyAdaptive);end
                    
                    penaltyBufferNonLinearPath_ii(ii+obj.MIOFF) = obj.problem.phases(1).path(ii+obj.MIOFF).buffer.lambda*penaltyBufferPathNonLinear_node;
                    if linearToo,penaltyBufferLinearPath_ii(ii+obj.MIOFF)    = obj.problem.phases(1).path(ii+obj.MIOFF).buffer.lambda*penaltyBufferPathLinear_node;end
                elseif ~obj.problem.phases(1).path(ii+obj.MIOFF).integral && ~obj.problem.phases(1).path(ii+obj.MIOFF).derivative
                    penaltyBufferPathNonLinear_max  = max(0,typeSign*(f - limit)./scalePath);
                    penaltyBufferPathNonLinear_node = norm(penaltyBufferPathNonLinear_max,obj.problem.phases(1).path(ii+obj.MIOFF).buffer.penaltyAdaptive);
                    if linearToo,penaltyBufferPathLinear_max     = max(0,typeSign*(f_ii + sum(dG_dX_ii.*(X - X_p),1) + sum(dG_dU_ii.*(U - U_p),1) - limit)./scalePath);end
                    if linearToo,penaltyBufferPathLinear_node    = norm(penaltyBufferPathLinear_max,obj.problem.phases(1).path(ii+obj.MIOFF).buffer.penaltyAdaptive);end
                    
                    penaltyBufferNonLinearPath_ii(ii+obj.MIOFF) = obj.problem.phases(1).path(ii+obj.MIOFF).buffer.lambda*penaltyBufferPathNonLinear_node;
                    if linearToo,penaltyBufferLinearPath_ii(ii+obj.MIOFF)    = obj.problem.phases(1).path(ii+obj.MIOFF).buffer.lambda*penaltyBufferPathLinear_node;end
                else
                    error('SCOPT Error, Path Constraint is not defined for nonlinear integral and derivative, or derivative terms \n Aborting')
                end
            case 'quasiconvex' % Could be included in the trust region error computation, but omitted
                penaltyBufferNonLinearPath_ii(ii+obj.MIOFF) = 0;
                penaltyBufferLinearPath_ii(ii+obj.MIOFF)    = 0;
            case 'non-linear-1var' % Could be included in the trust region error computation, but omitted
                penaltyBufferNonLinearPath_ii(ii+obj.MIOFF) = 0;
                penaltyBufferLinearPath_ii(ii+obj.MIOFF)    = 0;
            case {'convex','soc'}
                penaltyBufferNonLinearPath_ii(ii+obj.MIOFF) = 0;
                penaltyBufferLinearPath_ii(ii+obj.MIOFF)    = 0;
            otherwise
                warning('SCOPT Warning, Path function of type %s not implemented \n Continuing without constraint %i \n',Path(ii+obj.MIOFF).funType,ii)
                continue
        end
    end
    penaltyBuffersNonLinearPath = sum(penaltyBufferNonLinearPath_ii);
    penaltyBuffersLinearPath    = sum(penaltyBufferLinearPath_ii);
    
    penaltyBufferNonLinearEvents_ii = zeros(1,obj.problem.phases.n.events);
    penaltyBufferLinearEvents_ii    = zeros(1,obj.problem.phases.n.events);
    for ii = 0:obj.problem.phases.n.events-1
        switch obj.problem.phases(1).events(ii+obj.MIOFF).funType
            case 'non-linear'
                scaleEvents   = obj.problem.phases(1).events(ii+obj.MIOFF).scale;
                limit       = obj.problem.phases(1).events(ii+obj.MIOFF).limit;
                switch obj.problem.phases(1).events(ii+obj.MIOFF).where
                    case 'final'
                        indexWhere = Nodes - 1 +obj.MIOFF;
                    case 'initial'
                        indexWhere = 0 + obj.MIOFF;
                    case 'index'
                        indexWhere = obj.problem.phases(1).events(ii+obj.MIOFF).indexWhere;
                end
                switch obj.problem.phases(1).events(ii+obj.MIOFF).type
                    case 'lower'
                        typeSign = -1;
                    case 'upper'
                        typeSign = 1;
                    otherwise
                        error('SCOPT Error, Events Constraints of constraint type %s not implemented. Use inequality constraints instead \n Aborting',Path(ii+obj.MIOFF).type)
                end
                
                f             = obj.problem.phases(1).events(ii+obj.MIOFF).function(time(indexWhere),X(:,indexWhere),U(:,indexWhere),indexWhere);
                penaltyBufferNonLinearEvents_ii(ii+obj.MIOFF) = obj.problem.phases(1).events(ii+obj.MIOFF).buffer.lambda*max(0,typeSign*(f - limit)/scaleEvents);
                if linearToo
                    f_ii          = obj.problem.phases(1).events(ii+obj.MIOFF).function(time_p(indexWhere),X_p(:,indexWhere),U_p(:,indexWhere),indexWhere);
                    dG_dX_ii      = obj.problem.phases(1).events(ii+obj.MIOFF).jacobian.states(time_p(indexWhere),X_p(:,indexWhere),U_p(:,indexWhere),indexWhere);
                    dG_dU_ii      = obj.problem.phases(1).events(ii+obj.MIOFF).jacobian.controls(time_p(indexWhere),X_p(:,indexWhere),U_p(:,indexWhere),indexWhere);
                    penaltyBufferLinearEvents_ii(ii+obj.MIOFF)    = obj.problem.phases(1).events(ii+obj.MIOFF).buffer.lambda*max(0,typeSign*(f_ii + sum(dG_dX_ii.*(X(:,indexWhere) - X_p(:,indexWhere)),1) + sum(dG_dU_ii.*(U(:,indexWhere) - U_p(:,indexWhere)),1) - limit)/scaleEvents);
                end
            case 'quasiconvex' % Could be included in the trust region error computation, but omitted
                penaltyBufferNonLinearEvents_ii(ii+obj.MIOFF) = 0;
                penaltyBufferLinearEvents_ii(ii+obj.MIOFF)    = 0;
            case 'non-linear-1var' % Could be included in the trust region error computation, but omitted
                penaltyBufferNonLinearEvents_ii(ii+obj.MIOFF) = 0;
                penaltyBufferLinearEvents_ii(ii+obj.MIOFF)    = 0;
            otherwise
                penaltyBufferNonLinearEvents_ii(ii+obj.MIOFF) = 0;
                penaltyBufferLinearEvents_ii(ii+obj.MIOFF)    = 0;
        end
    end
    penaltyBuffersNonLinearEvents = sum(penaltyBufferNonLinearEvents_ii);
    penaltyBuffersLinearEvents    = sum(penaltyBufferLinearEvents_ii);
    
    penaltyBuffersNonLinear       = penaltyBuffersNonLinearEvents + penaltyBuffersNonLinearPath;
    penaltyBuffersLinear          = penaltyBuffersLinearEvents    + penaltyBuffersLinearPath   ;
    if ~obj.quiet
        fprintf(['Non-linear Penalisation of Constraints [' repmat('%0.2e ',1,obj.problem.phases.n.path) repmat('%0.2e ',1,obj.problem.phases.n.events) '] \n'],[penaltyBufferNonLinearPath_ii,penaltyBufferNonLinearEvents_ii]);
        fprintf(['Linear Penalisation of Constraints [' repmat('%0.2e ',1,obj.problem.phases.n.path) repmat('%0.2e ',1,obj.problem.phases.n.events) '] \n'],[penaltyBufferLinearPath_ii,penaltyBufferLinearEvents_ii]);
    end
else
    penaltyBuffersNonLinear = 0;
    penaltyBuffersLinear    = 0;
end
%% Soft Time Trust Region Penalty
if obj.problem.phases(1).freeTimeFinal && strcmp(obj.algorithm.sequential.trustRegion.phases(1).timeFinal.type,'soft')
    penaltyTime  = obj.algorithm.sequential.trustRegion.phases(1).timeFinal.lambda*norm(timeFinal_n - timeFinal_n_p);
else
    penaltyTime  = 0 ;
end
%% Trust Region Computation and Penalty
if strcmp(obj.algorithm.sequential.trustRegion.type,'soft') && strcmp(obj.algorithm.sequential.type,'trust-region') && obj.algorithm.sequential.activate
    AuxTrust = zeros(0,Nodes);
    if obj.algorithm.sequential.trustRegion.include.states
        AuxTrust = [AuxTrust;X_n - X_n_p];
    end
    if obj.algorithm.sequential.trustRegion.include.controls
        AuxTrust = [AuxTrust;U_n - U_n_p];
    end
    
    switch obj.algorithm.sequential.trustRegion.variablesPenalty
        case 0
            trustRegionP1 = sum(AuxTrust.*AuxTrust,1);
        case {1,2,inf}
            trustRegionP1 = vecnorm(AuxTrust,obj.algorithm.sequential.trustRegion.variablesPenalty,1);
        otherwise
            error('SCOPT ERROR \n Trust Region with variable penalty type %i not implemented \n',obj.algorithm.sequential.trustRegion.variablesPenalty )
    end
    
    switch obj.algorithm.sequential.trustRegion.nodePenalty
        case {1,2,inf}
            trustRegionP2 = vecnorm(trustRegionP1,obj.algorithm.sequential.trustRegion.nodePenalty,2);
        otherwise
            error('SCOPT ERROR \n Trust Region with node penalty type %i not implemented \n',obj.algorithm.sequential.trustRegion.nodePenalty )
    end
    penaltyTrustRegion  = obj.algorithm.sequential.trustRegion.lambda * trustRegionP2;
else
    penaltyTrustRegion  = 0 ;
end
end


