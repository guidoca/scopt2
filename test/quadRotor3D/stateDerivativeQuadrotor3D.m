function [Xdot] = stateDerivativeQuadrotor3D(t,X,U,BodyMap,ii)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
planet  = BodyMap.CentralBody      ;
vehicle = BodyMap.ControlledBody  ;
Xdot     = zeros(size(X));
%% Retrieve State variables

V        = X(4:6,:);
T        = U(1:3,:);

%% Retrieve Environment

% Earth Gravity
g           = BodyMap.(planet).Gravity.g;

% Mass Model
Mass        = BodyMap.(vehicle).MassModel.mass ;

% Aerodynamic
cd          = BodyMap.(vehicle).Aerodynamic.cd;
includeAero = BodyMap.(vehicle).Aerodynamic.include;
%% Compute State Derivatives
Xdot(1:3,:)      = V;
Xdot(4:6,:)      = T/Mass - includeAero*cd*repmat(vecnorm(V,2,1),3,1).*V + repmat(g,1,size(V,2)) ; 

end

