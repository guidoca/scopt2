function [virtualControlsP1,virtualControlsP2] =  computePenaltiesVirtualControls(obj,virtualControls)
% Computes penalties P1 and P2 produced by the use of virtual
% controls based on the state penalties (per node) and their
% global norm. Handles multiple iterations arrays.
EvirtualControls   = repmat(obj.algorithm.virtualControl.phases(1).valE,1,size(virtualControls,2),size(virtualControls,3)).*virtualControls;
switch obj.algorithm.virtualControl.phases(1).statePenalty
    case 0
        virtualControlsP1 = sum(EvirtualControls.*EvirtualControls,1);
    case {1,2,Inf}
        virtualControlsP1 = vecnorm(EvirtualControls,obj.algorithm.virtualControl.phases(1).statePenalty,1);
    otherwise
        error('SCOPT ERROR \n Virtual Controls with state penalty type %i not implemented \n',obj.algorithm.sequential.trustRegion.variablesPenalty )
end

switch obj.algorithm.virtualControl.phases(1).nodePenalty
    case {1,2,Inf}
        virtualControlsP2 = reshape(vecnorm(virtualControlsP1,obj.algorithm.virtualControl.phases(1).nodePenalty,2),1,[]);
    otherwise
        error('SCOPT ERROR \n Virtual Controls with node penalty type %i not implemented \n',obj.algorithm.sequential.trustRegion.nodePenalty )
end
end


