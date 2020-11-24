function [Jreal,Lreal] = computeObjectiveFunctions(obj,time,X,U,time_p,X_p,U_p)
% Computes unscaled nonlinear (Jreal) and convexified (Lreal)
% objective functions
Objective = obj.problem.objective;

Nodes      = obj.problem.phases.Nodes;
timeInterval = time(end)-time(1);

if nargin>4
    timeInterval_p = time_p(end)-time_p(1);
    BodyMap_prev  = obj.problem.phases(1).BodyMap_prev;BodyMap_prev.evaluated = 0;
end
if ~strcmp(Objective.type,'feasibility')
    % Mayer Term
    switch Objective.mayer.where
        case 'final'
            indexWhere = obj.problem.phases.Nodes - 1 + obj.MIOFF;
        case 'initial'
            indexWhere = 0 + obj.MIOFF;
        case 'index'
            indexWhere = Objective.mayer.indexWhere;
        otherwise
    end
    switch Objective.mayer.funType
        case 'linear'
            Jreal =  (obj.problem.phases(1).freeTimeFinal*time(indexWhere)*Objective.mayer.timeFinal + dot(X(:,indexWhere),Objective.mayer.states) + dot(U(:,indexWhere),Objective.mayer.controls));
            if nargin>4
                Lreal =  Jreal;
            else
                Lreal = 0;
            end
        case {'convex','soc'}
            Jreal = (norm([Objective.mayer.cone.norm.states;Objective.mayer.cone.norm.controls]'*[X(:,indexWhere);U(:,indexWhere)] + Objective.mayer.cone.norm.cons(:))...
                - Objective.mayer.cone.right.cons-[Objective.mayer.cone.right.states(:);Objective.mayer.cone.right.controls(:)]'*[X(:,indexWhere);U(:,indexWhere)]);
            if nargin>4
                Lreal =  Jreal;
            else
                Lreal = 0;
            end
        case 'non-linear'
            Jreal =  1/Objective.scale*Objective.mayer.function(time(indexWhere),X(:,indexWhere),U(:,indexWhere),obj.problem.phases(1).BodyMap,obj.problem.phases(1).index.nodes(indexWhere));
            if nargin>4
                Lreal =  (Objective.mayer.function(time_p(indexWhere),X_p(:,indexWhere),U_p(:,indexWhere),BodyMap_prev,obj.problem.phases(1).index.nodes(indexWhere))+...
                    sum(Objective.mayer.jacobian.states(time_p(indexWhere),X_p(:,indexWhere),U_p(:,indexWhere),BodyMap_prev,obj.problem.phases(1).index.nodes(indexWhere)).*(X(:,indexWhere) - X_p(:,indexWhere)),1)+...
                    sum(Objective.mayer.jacobian.controls(time_p(indexWhere),X_p(:,indexWhere),U_p(:,indexWhere),BodyMap_prev,obj.problem.phases(1).index.nodes(indexWhere)).*(U(:,indexWhere) - U_p(:,indexWhere)),1));
            else
                Lreal = 0;
            end
    end
    
    % Lagrange Term
    switch obj.problem.objective.lagrange.funType
        case 'linear'
            f_ii = sum(repmat(Objective.lagrange.states(:),1,Nodes).*X,1) + sum(repmat(Objective.lagrange.controls(:),1,Nodes).*U,1);if obj.problem.objective.lagrange.absolute.nodes,f_ii = abs(f_ii);end
            
            I  =  sum(obj.problem.phases(1).quadratureVector.*f_ii);
            
            Jreal =  timeInterval*I;
            if nargin>4
                fp  = sum(repmat(Objective.lagrange.states(:),1,Nodes).*X_p,1) + sum(repmat(Objective.lagrange.controls(:),1,Nodes).*U_p,1);if obj.problem.objective.lagrange.absolute.nodes,fp = abs(fp);end
                Ip  =  sum(obj.problem.phases(1).quadratureVector.*fp);
                
                Lreal  =  timeInterval_p*I + Ip*(time(end)-time_p(end));
            else
                Lreal  =   0;
            end
        case {'convex','soc'}
            %                             J(kk) =  1/Objective.scale*timeDisc_opt*sum(Z_opt(index.objectiveLagrange)); %ATTENTION LAGRANGE TERM COMPUTE THE FUNCTION AGAIN
            NormTermDim = zeros(Objective.lagrange.cone.dimensions-1,Nodes);
            for dd = 0:Objective.lagrange.cone.dimensions-2
                NormTermDim(dd+obj.MIOFF,:)=[Objective.lagrange.cone.norm.states(:,dd+obj.MIOFF);Objective.lagrange.cone.norm.controls(:,dd+obj.MIOFF)]'*[X;U] + Objective.lagrange.cone.norm.cons(dd+obj.MIOFF);
            end
            NormTerm = vecnorm(NormTermDim,2,1);
            f_ii = (NormTerm  - Objective.lagrange.cone.right.cons-[Objective.lagrange.cone.right.states(:);Objective.lagrange.cone.right.controls(:)]'*[X;U]);
            I    = sum(obj.problem.phases(1).quadratureVector.*f_ii);
            %                         J    =  timeDisc*sum(f_ii);
            Jreal    =  timeInterval*I;
            
            if nargin>4
                for dd = 0:Objective.lagrange.cone.dimensions-2
                    NormTermDim(dd+obj.MIOFF,:)=[Objective.lagrange.cone.norm.states(:,dd+obj.MIOFF);Objective.lagrange.cone.norm.controls(:,dd+obj.MIOFF)]'*[X_p;U_p] + Objective.lagrange.cone.norm.cons(dd+obj.MIOFF);
                end
                NormTerm = vecnorm(NormTermDim,2,1);
                f_ii_p = (NormTerm  - Objective.lagrange.cone.right.cons-[Objective.lagrange.cone.right.states(:);Objective.lagrange.cone.right.controls(:)]'*[X_p;U_p]);
                Ip      = sum(obj.problem.phases(1).quadratureVector.*f_ii_p);
                %                             L      =  (timeDisc_p*sum(f_ii) +  tauDisc*sum(f_ii_p)/2*(time(end)-time_p(end)));
                Lreal  =  timeInterval_p*I + Ip*(time(end)-time_p(end));
            else
                Lreal  = 0;
            end
        case 'non-linear'
            f_ii  = Objective.lagrange.function(time,X,U,obj.problem.phases(1).BodyMap,obj.problem.phases(1).index.nodes);if obj.problem.objective.lagrange.absolute.nodes,f_ii=abs(f_ii);end
            
            I        =  sum(obj.problem.phases(1).quadratureVector.*f_ii);
            Jreal    =  timeInterval*I;
            if nargin>4
                f_ii_p    = Objective.lagrange.function(time_p,X_p,U_p,BodyMap_prev,obj.problem.phases(1).index.nodes);
                AX_p      = sum(Objective.lagrange.jacobian.states(time_p,X_p,U_p,BodyMap_prev,obj.problem.phases(1).index.nodes).*(X-X_p),1);
                BU_p      = sum(Objective.lagrange.jacobian.controls(time_p,X_p,U_p,BodyMap_prev,obj.problem.phases(1).index.nodes).*(U-U_p),1);
                if obj.problem.objective.lagrange.absolute.nodes,f_ii_pabs=abs(f_ii_p);else,f_ii_pabs=f_ii_p ;end
                f_ii_lin  = (f_ii_p + AX_p + BU_p);if obj.problem.objective.lagrange.absolute.nodes,f_ii_lin=abs(f_ii_lin);end
                
                I_p    = sum(obj.problem.phases(1).quadratureVector.*f_ii_pabs);
                I_aux  = sum(obj.problem.phases(1).quadratureVector.*f_ii_lin);
                Lreal  = timeInterval_p*I_aux + I_p*(time(end)-time_p(end));
            else
                Lreal = 0;
            end
    end
    
    % minmax Term
    switch obj.problem.objective.minmax.funType
        case 'linear'
            f_ii = sum(Objective.minmax.states(:).*X,1) + sum(Objective.minmax.controls(:).*U,1);if obj.problem.objective.minmax.absolute.nodes,f_ii = abs(f_ii);end
            switch obj.problem.objective.type
                case 'minimize'
                    Jreal =  max(f_ii);
                case 'maximize'
                    Jreal =  min(f_ii);
            end
            if nargin>4
                Lreal  =  Jreal;
            else
                Lreal  =   0;
            end
        case {'convex','soc'}
            NormTermDim = zeros(Objective.minmax.cone.dimensions-1,Nodes);
            for dd = 0:Objective.minmax.cone.dimensions-2
                NormTermDim(dd+obj.MIOFF,:)=[Objective.minmax.cone.norm.states(:,dd+obj.MIOFF);Objective.minmax.cone.norm.controls(:,dd+obj.MIOFF)]'*[X;U] + Objective.minmax.cone.norm.cons(dd+obj.MIOFF);
            end
            NormTerm = vecnorm(NormTermDim,2,1);
            f_ii = (NormTerm  - Objective.minmax.cone.right.cons-[Objective.minmax.cone.right.states(:);Objective.minmax.cone.right.controls(:)]'*[X;U]);
            switch obj.problem.objective.type
                case 'minimize'
                    Jreal =  max(f_ii);
                case 'maximize'
                    Jreal =  min(f_ii);
            end
            if nargin>4
                Lreal  =  Jreal;
            else
                Lreal  =   0;
            end
        case 'non-linear'
            f_ii  = Objective.minmax.function(time,X,U,obj.problem.phases(1).BodyMap,obj.problem.phases(1).index.nodes);if obj.problem.objective.minmax.absolute.nodes,f_ii=abs(f_ii);end
            
            switch obj.problem.objective.type
                case 'minimize'
                    Jreal =  max(f_ii);
                case 'maximize'
                    Jreal =  min(f_ii);
            end
            if nargin>4
                f_ii_p    = Objective.minmax.function(time_p,X_p,U_p,BodyMap_prev,obj.problem.phases(1).index.nodes);
                AX_p      = sum(Objective.minmax.jacobian.states(time_p,X_p,U_p,BodyMap_prev,obj.problem.phases(1).index.nodes).*(X-X_p),1);
                BU_p      = sum(Objective.minmax.jacobian.controls(time_p,X_p,U_p,BodyMap_prev,obj.problem.phases(1).index.nodes).*(U-U_p),1);
                f_ii_lin  = (f_ii_p + AX_p + BU_p);if obj.problem.objective.minmax.absolute.nodes,f_ii_lin=abs(f_ii_lin);end
                
                switch obj.problem.objective.type
                    case 'minimize'
                        Lreal =  max(f_ii_lin);
                    case 'maximize'
                        Lreal =  min(f_ii_lin);
                end
            else
                Lreal = 0;
            end
    end
else
    Jreal = 0;
    Lreal = 0;
end
end


