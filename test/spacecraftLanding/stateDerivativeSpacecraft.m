function [Xdot] = stateDerivativeSpacecraft(t,X,U,BodyMap,ii)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
planet  = BodyMap.CentralBody      ;
vehicle = BodyMap.ControlledBody  ;
Xdot     = zeros(size(X));
%% Retrieve State variables

P        = X(1:3,:);
V        = X(4:6,:);
T        = U(1:3,:);

%% Retrieve Environment

% Earth Gravity
g       = BodyMap.(planet).Gravity.g;
vector  = BodyMap.(planet).Gravity.vector;
g_vec   = vector(:)*g;

% Mass Model
Mass     = BodyMap.(vehicle).MassModel.mass ;
%% Compute State Derivatives
Xdot(1:3,:)      = V;
Xdot(4:6,:)      = T/Mass + repmat(g_vec,1,size(T,2)); 

end

