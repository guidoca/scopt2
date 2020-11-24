function Us_n = scaleControlsToNormalized(obj,ph,Us,index)
% Scales controls
scale = obj.problem.phases(ph).scale.controls;
shift = obj.problem.phases(ph).shift.controls;
if nargin<4
    index = (1:obj.problem.phases(ph).n.controls)';
end
if isempty(Us)
    Us_n = Us;
else
    Nodes= size(Us,2);
    Us_n = repmat((1./scale(index)),1,Nodes).*(Us - repmat(shift(index),1,Nodes))   ;
end
end

