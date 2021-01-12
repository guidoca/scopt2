function [f] = stateDerivativeSolutionRocketLanding(t,X,U,BodyMap)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
f    = zeros(size(X,1),length(t)-1);
dt   = diff(t);
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
    f(1:3,:)      = repmat(dt,3,1).*V(:,1:end-1)+ repmat(dt.^2,3,1).*aE(:,1:end-1)/2;
    f(4:6,:)      = repmat(dt,3,1).*(aE(:,1:end-1));
    switch BodyMap.DynamicModel.controlVariable
        case 'thrust' 
            f(7,:) = dt.*mdotbp(1:end-1);
        case 'acceleration' 
            f(7,:) = dt.*mdotbp(1:end-1).*Minv(1:end-1);
    end
elseif BodyMap.DynamicModel.exactOde==1
    f(1:3,:)      = repmat(dt,3,1).*V(:,1:end-1)+ repmat(dt.^2,3,1).*(aE(:,1:end-1)/3+aE(:,2:end)/6);
    f(4:6,:)      = repmat(dt,3,1).*(aE(:,1:end-1)/2+aE(:,2:end)/2);
    switch BodyMap.DynamicModel.controlVariable
        case 'thrust' 
            f(7,:) = dt*1/2.*(mdotbp(1:end-1)+mdotbp(2:end));
        case 'acceleration' 
            f(7,:) = dt*1/2.*(mdotbp(1:end-1).*Minv(1:end-1)+mdotbp(2:end).*Minv(2:end));
    end
else
    error('Dynamics Error: exactOde type not defined \n')
end

switch BodyMap.ThrustModel.collocationMethod
    case {'zoh','floor'}
        f(1:3,:) = f(1:3,:) + repmat(dt.^2,3,1).*a(:,1:end-1)/2;
        f(4:6,:) = f(4:6,:) + repmat(dt,3,1).*a(:,1:end-1);
        switch BodyMap.DynamicModel.controlVariable
            case 'thrust' 
                f(7,:) = f(7,:)   - dt.*Tmag(1:end-1)./cex; 
            case 'acceleration' 
                f(7,:) = f(7,:)   - dt.*amag(1:end-1)./cex; 
        end
    case {'foh','linear'}
        f(1:3,:) = f(1:3,:) + repmat(dt.^2,3,1).*(a(:,1:end-1)/3+a(:,2:end)/6);
        f(4:6,:) = f(4:6,:) + repmat(dt,3,1).*(a(:,1:end-1)/2+a(:,2:end)/2);
        switch BodyMap.DynamicModel.controlVariable
            case 'thrust' 
                f(7,:) = f(7,:)   - dt.*(Tmag(1:end-1)./cex/2+Tmag(2:end)./cex/2); 
            case 'acceleration' 
                f(7,:) = f(7,:)   - dt.*(amag(1:end-1)./cex/2+amag(2:end)./cex/2); 
        end
end
    
end

