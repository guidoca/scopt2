function Us   = scaleControlsToReal(obj,ph,Us_n,index)
% Unscales controls and retrieves actual values
scale = obj.problem.phases(ph).scale.controls;
shift = obj.problem.phases(ph).shift.controls;
if nargin<4
    index = (1:obj.problem.phases(ph).n.controls)';
end
Nodes = size(Us_n,2);
Us         = repmat(scale(index),1,Nodes).*Us_n + repmat(shift(index),1,Nodes)   ;
end

