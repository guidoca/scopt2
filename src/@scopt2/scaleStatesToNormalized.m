function Xs_n = scaleStatesToNormalized(obj,ph,Xs,index)
% Scales states
scale = obj.problem.phases(ph).scale.states;
shift = obj.problem.phases(ph).shift.states;
if nargin<4
    index = (1:obj.problem.phases(ph).n.states)';
end
if isempty(Xs)
    Xs_n = Xs;
else
    Nodes        = size(Xs,2);
    Xs_n         = repmat((1./scale(index)),1,Nodes).*(Xs - repmat(shift(index),1,Nodes))   ;
end
end

