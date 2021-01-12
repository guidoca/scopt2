function [Xdot] = stateDerivativeRocketLanding(t,X,U,BodyMap)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
Xdot    = zeros(size(X));
%% Retrieve State variables

R        = X(1:3,:);
V        = X(4:6,:);
switch BodyMap.DynamicModel.controlVariable
    case 'thrust' 
        M        = X(7,:);
        Minv     = 1./M;
        T        = U(1:3,:);
        Tmag     = U(4,:); 
        a        = T.*Minv; 
    case 'acceleration' 
        Z        = X(7,:);
        Minv     = exp(-Z);
        a        = U(1:3,:);
        amag     = U(4,:); 
end

%% Retrieve Environment

% Orientation Model
aR = BodyMap.OrientationModel.computeRotatingAcceleration(R,V);
siteNormal  = [1;0;0];

% Earth Shape Model
h = R(1,:);

% Earth Gravity
g        = BodyMap.GravityModel.g;
gVec     = repmat(-siteNormal*g,1,length(t));

% Atmosphere Model
if BodyMap.ThrustModel.includeBackPressure || BodyMap.AerodynamicModel.includeDrag
    [rho , pamb, ~ , soundSpeed] = BodyMap.AtmosphereModel.get_properties(h,[],[],[]);
end

% Thrust Model
cex         = BodyMap.ThrustModel.ceffVac; 
if BodyMap.ThrustModel.includeBackPressure  
	mdotbp       = -pamb*BodyMap.ThrustModel.Ae./cex;
else
	mdotbp       = 0;
end

% Aerodynamic
if BodyMap.AerodynamicModel.includeDrag
	Area    = BodyMap.AerodynamicModel.area;
    Vmag    = sqrt(V(1,:).^2+V(2,:).^2+V(3,:).^2);
    mach    = Vmag./soundSpeed;	
    cd      = BodyMap.AerodynamicModel.computeCd(mach);							  
    Drag    = - 0.5*Area*repmat(rho.*cd.*Vmag,3,1).*V;
    aD      = Drag.*repmat(Minv,3,1);
else
    aD      = zeros(3,length(t));
end

% Sum external accelerations
aE = aD + gVec + aR;
%% Compute State Derivatives
Xdot(1:3,:)      = V;
Xdot(4:6,:)      = a + aE; 
switch BodyMap.DynamicModel.controlVariable
    case 'thrust' 
        Xdot(7,:)        = -Tmag./cex + mdotbp; % mdot = -(Teff + Ae*(Pa))/c, z = log(m), dz = 1/m*dm = -T/m/c - Ae*Pa/c/m = -a/c - Ae*Pa/c*exp(-z)
    case 'acceleration' 
        Xdot(7,:)        = -amag./cex + mdotbp.*Minv; % mdot = -(Teff + Ae*(Pa))/c, z = log(m), dz = 1/m*dm = -T/m/c - Ae*Pa/c/m = -a/c - Ae*Pa/c*exp(-z)
end
end

