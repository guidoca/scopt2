function t_n  = scaleTimeToNormalized(obj,ph,t)
% Scales time
scale    = obj.problem.phases(ph).scale.time;
shift    = obj.problem.phases(ph).shift.time;
t_n      = (1./scale).*(t - shift)   ;
end

