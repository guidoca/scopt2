function t    = scaleTimeToReal(obj,ph,t_n)
% Unscales time and retrieves actual values
scale    = obj.problem.phases(ph).scale.time;
shift    = obj.problem.phases(ph).shift.time;
t        = scale.*(t_n) + shift   ;
end

