function [ BodyMap ] = bodyMap( )
% Sets default bodyMap for problem
    BodyMap.GravityModel     = gravityModel ;
    BodyMap.ThrustModel      = thrustModel ;
    BodyMap.MassModel        = massModel ;
    BodyMap.OrientationModel = orientationModel ;
    BodyMap.AerodynamicModel = aerodynamicModel ;
    BodyMap.AtmosphereModel  = atmosphereModel ;
    BodyMap.DynamicModel.exactOde = 1;
    BodyMap.DynamicModel.controlVariable = 'thrust';

end

