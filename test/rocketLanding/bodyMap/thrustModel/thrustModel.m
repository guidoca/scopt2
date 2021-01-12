classdef thrustModel
    %Thrust Model class with information on rocket
    % propulsion characteristics
    properties (Constant)
        g0 = 9.807  ; % Earth surface gravity constant [m/s^2] 
    end
    properties
        thrustVac             = 250e3+0.5*100.0e3    ; % Maximum Vaccum Thrust [N]
        throttleMax           = 1        ; % Maximum Throttle setting [-]
        throttleMin           = (100e3+0.5*100.0e3)/(250e3+0.5*100.0e3)   ; % Minimum Throttle setting [-]
        thrustRateMax         = 100e3     ; % [N/s]
        IspVac                = 300      ; % Effective vacuum exhaust speed [s]
        Ae                    = 0.5      ; % Nozzle Exit Area [m2]
        includeBackPressure   = 1        ;
		includeThrustRateLimit= 1        ;
        collocationMethod     = 'foh'    ; % thrust profile assumed collocation strategy (zoh/floor,foh/linear,spline,else)
    end
    properties (Dependent)
        ceffVac   ;
        mdotVac   ;
        mdotMax   ; 
        mdotMin   ; 
        thrustMax ;
        thrustMin ;
    end
    methods
        function thrustVac = computeVacuumThrust(obj,throttle)
            thrustVac = obj.thrustVac*throttle  ;
        end
        function thrust = computeThrust(obj,throttle,atmosPressure)
            thrust =  computeVacuumThrust(obj,throttle) ;
            if obj.includeBackPressure
                thrust =  thrust - obj.Ae*atmosPressure ;
            end
        end
        function Isp = computeEffectiveIsp(obj,atmosPressure)
            Isp = obj.IspVac - obj.Ae*atmosPressure./obj.mdotVac/obj.g0;
        end
        function ceff = computeEffectiveExhaustSpeed(obj,atmosPressure)
            ceff = obj.ceffVac - obj.Ae*atmosPressure./obj.mdotVac;
        end
        function ceff = computeEffectiveExhaustSpeedDerivative(obj)
            ceff = - obj.Ae./obj.mdotVac;
        end
        function mdot = computeMassFlow(obj,throttle)
            mdot = throttle*obj.mdotVac  ; 
        end
        function ceffVac = get.ceffVac(obj)
            ceffVac  = obj.IspVac*obj.g0 ;
        end
        function mdotVac = get.mdotVac(obj)
            mdotVac  = obj.thrustVac/obj.ceffVac ;
        end
        function thrustMax = get.thrustMax(obj)
            thrustMax  = computeVacuumThrust(obj,obj.throttleMax) ;
        end
        function thrustMin = get.thrustMin(obj)
            thrustMin  = computeVacuumThrust(obj,obj.throttleMin) ;
        end
        function mdotMax = get.mdotMax(obj)
            mdotMax  = computeMassFlow(obj,obj.throttleMax);
        end
        function mdotMin = get.mdotMin(obj)
            mdotMin  = computeMassFlow(obj,obj.throttleMin);
        end
    end
    
end

