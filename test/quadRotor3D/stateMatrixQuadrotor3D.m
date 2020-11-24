function [row,col,val,nEntries] = stateMatrixQuadrotor3D(t,X,U,BodyMap,ii)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
planet  = BodyMap.CentralBody      ;
vehicle = BodyMap.ControlledBody  ;
nodes   = size(X,2);
V        = X(4:6,:);

Vmag      = sqrt(V(1,:).^2+V(2,:).^2+V(3,:).^2);

IndexZero = find(Vmag~=0);

invVmag = zeros(1,nodes);
invVmag(IndexZero)   = 1./Vmag(IndexZero);

%% Retrieve Environment

% Earth Gravity
g           = BodyMap.(planet).Gravity.g;

% Mass Model
Mass        = BodyMap.(vehicle).MassModel.mass ;

% Aerodynamic
cd          = BodyMap.(vehicle).Aerodynamic.cd;
includeAero = BodyMap.(vehicle).Aerodynamic.include;
%% Compute State Derivatives
row     = repmat([(1:6)'],1,nodes);
col     = repmat([(4:6)';(4:6)'],1,nodes);

nEntries= size(row,1);

val = zeros(nEntries,nodes);

val( 1:3,:) = ones(3,nodes);

val( 4,:)   = -includeAero*cd*(V(1,:).^2.*invVmag + Vmag); % norm(V)*V = 1/sqrt(vx^2+vy^2)*[Vx;Vy]*(Vx+Vy) + sqrt(vx^2+vy^2)*[1;0] + sqrt(vx^2+vy^2)*[0;1]
val( 5,:)   = -includeAero*cd*(V(2,:).^2.*invVmag + Vmag); % norm(V)*V = 1/sqrt(vx^2+vy^2)*[Vx;Vy]*(Vx+Vy) + sqrt(vx^2+vy^2)*[1;0] + sqrt(vx^2+vy^2)*[0;1]
val( 6,:)   = -includeAero*cd*(V(3,:).^2.*invVmag + Vmag); % norm(V)*V = 1/sqrt(vx^2+vy^2)*[Vx;Vy]*(Vx+Vy) + sqrt(vx^2+vy^2)*[1;0] + sqrt(vx^2+vy^2)*[0;1]

end

