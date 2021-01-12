function [dfdt] = timeDerivativeSolutionRocketLanding(t,X,U,BodyMap)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
dfdt    = zeros(size(X,1),length(t)-1);
dt      = diff(t);
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
cex   = BodyMap.ThrustModel.ceffVac; 
    
if BodyMap.ThrustModel.includeBackPressure 
    mdotbp = -pamb.*BodyMap.ThrustModel.Ae./cex;
else
	mdotbp = zeros(1,length(t));
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
    aD = zeros(3,length(t));
end

% Sum external accelerations
aE = aD + gVec + aR;
%% Compute State Derivatives
if BodyMap.DynamicModel.exactOde==0
    dfdt(1:3,:)      = V(:,1:end-1)+ repmat(dt,3,1).*aE(:,1:end-1);
    dfdt(4:6,:)      = aE(:,1:end-1);
    switch BodyMap.DynamicModel.controlVariable
        case 'thrust' 
            dfdt(7,:) = mdotbp(1:end-1);
        case 'acceleration' 
            dfdt(7,:) = mdotbp(1:end-1).*Minv(1:end-1);
    end
elseif BodyMap.DynamicModel.exactOde==1
    dfdt(1:3,:)      = V(:,1:end-1)+ 2*repmat(dt,3,1)/3.*(aE(:,1:end-1)+aE(:,2:end)/2);
    dfdt(4:6,:)      = 1/2.*(aE(:,1:end-1)+aE(:,2:end));
    switch BodyMap.DynamicModel.controlVariable
        case 'thrust' 
    dfdt(7,:)        = 1/2.*(mdotbp(1:end-1)+mdotbp(2:end));
        case 'acceleration' 
    dfdt(7,:)        = 1/2.*(mdotbp(1:end-1).*Minv(1:end-1)+mdotbp(2:end).*Minv(2:end));
    end
else
    error('Dynamics Error: exactOde type not defined \n')
end

switch BodyMap.ThrustModel.collocationMethod
    case {'zoh','floor'}
        dfdt(1:3,:) = dfdt(1:3,:) + repmat(dt,3,1).*a(:,1:end-1);
        dfdt(4:6,:) = dfdt(4:6,:) + a(:,1:end-1);
        switch BodyMap.DynamicModel.controlVariable
            case 'thrust' 
        dfdt(7,:)   = dfdt(7,:)   - Tmag(1:end-1)./cex; 
            case 'acceleration' 
        dfdt(7,:)   = dfdt(7,:)   - amag(1:end-1)./cex; 
        end
    case {'foh','linear'}
        dfdt(1:3,:) = dfdt(1:3,:) + 2*repmat(dt,3,1).*(a(:,1:end-1)/3+a(:,2:end)/6);
        dfdt(4:6,:) = dfdt(4:6,:) + a(:,1:end-1)/2+a(:,2:end)/2;
        switch BodyMap.DynamicModel.controlVariable
            case 'thrust' 
        dfdt(7,:)   = dfdt(7,:)   - 1/2*(Tmag(1:end-1)+Tmag(2:end))/cex; 
            case 'acceleration' 
        dfdt(7,:)   = dfdt(7,:)   - 1/2*(amag(1:end-1)+amag(2:end))/cex; 
        end
end
end

