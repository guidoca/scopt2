function Xs   = scaleStatesToReal(obj,ph,Xs_n,index)
% Unscales states and retrieves actual values
scale = obj.problem.phases(ph).scale.states;
shift = obj.problem.phases(ph).shift.states;
if nargin<4
    index = (1:obj.problem.phases(ph).n.states)';
end
Nodes = size(Xs_n,2);
Xs    = repmat(scale(index),1,Nodes).*Xs_n + repmat(shift(index),1,Nodes)   ;
end

