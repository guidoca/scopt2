% Rocket Landing problem verifying correct implementation of variable mass
% and free final time successive convexification for powered descent.

clc; clear all; close all;

addpath(genpath('../../src'));
addpath(genpath(pwd()));
FiguresFolder = 'Figures';
if ~isfolder(FiguresFolder)
    mkdir(FiguresFolder)
end

min2sec    = 60;
dbstop if error
%%
discretisationType = 'exact';
objectiveType = 1 ;
freeTime      = 1 ;
solveSCVx     = 1 ;

updateInitialGuessProblem2 = 1;
%% Environment And Vehicle
BodyMapRocket = bodyMap();
 
% Atmosphere
BodyMapRocket.AtmosphereModel.rho0  = 1.0; % kg/m3 
BodyMapRocket.AtmosphereModel.Hs    = 100.0e3 / BodyMapRocket.AtmosphereModel.g0/BodyMapRocket.AtmosphereModel.rho0;

% Dynamics Representation
BodyMapRocket.DynamicModel.controlVariable = 'acceleration';

% Gravity Model
BodyMapRocket.GravityModel.g  = 3.7114;
BodyMapRocket.GravityModel.g0 = 3.7114;

% Mass Model
BodyMapRocket.MassModel.massWet   = 1905;
BodyMapRocket.MassModel.massInert = 1205;

% Aerodynamic Model
BodyMapRocket.AerodynamicModel.dragMethod = 'none';
BodyMapRocket.AerodynamicModel.area = 1.0;

% Thrust Model
BodyMapRocket.ThrustModel.thrustVac             = 6*3.1e3*cosd(27)  ; % Maximum Vaccum Thrust [N]
BodyMapRocket.ThrustModel.throttleMax           = 0.8      ; % Maximum Throttle setting [-]
BodyMapRocket.ThrustModel.throttleMin           = 0.3      ; % Minimum Throttle setting [-]
BodyMapRocket.ThrustModel.thrustRateMax         = 10.0e3   ; % [N/s]
BodyMapRocket.ThrustModel.IspVac                = 225*cosd(27)      ; % Effective vacuum exhaust speed [s]
BodyMapRocket.ThrustModel.Ae                    = 0;0.5      ; % Nozzle Exit Area [m2]
BodyMapRocket.ThrustModel.includeBackPressure   = 0        ;
BodyMapRocket.ThrustModel.includeThrustRateLimit= 0        ;
BodyMapRocket.ThrustModel.collocationMethod     = 'foh'    ; % thrust profile assumed collocation strategy (zoh/floor,foh/linear,spline,else)
   
BodyMapRocket.DynamicModel.exactOde = 1;
% Retrieve Models
GravityModel     = BodyMapRocket.GravityModel;
AtmosphereModel  = BodyMapRocket.AtmosphereModel;
AerodynamicModel = BodyMapRocket.AerodynamicModel;
MassModel        = BodyMapRocket.MassModel;
ThrustModel      = BodyMapRocket.ThrustModel;

FmaxThrottle = ThrustModel.computeThrust(1,AtmosphereModel.constantPressure);
Fmax = ThrustModel.computeThrust(ThrustModel.throttleMax,AtmosphereModel.constantPressure);
Fmin = ThrustModel.computeThrust(ThrustModel.throttleMin,AtmosphereModel.constantPressure);

amax = Fmax./MassModel.massInert;
amin = Fmin./MassModel.massWet;

%% Initial, final conditions and constraints

Nodes      = 82;
% Nodes      = 70;

n0 = [1,0,0]';

Tmag0 = nan;
R0    = [1.5e3,0,2e3]';
V0    = [-75,0,100]';
V0    = [-75,0,50]';
M0    = MassModel.massWet;
z0    = log(M0);
t0    = 0.0;
amag0 = Tmag0/M0;

Rmag0 = norm(R0);
Vmag0 = norm(V0);

zU    = log(M0);
zL    = log(MassModel.massInert);

site.n = [1,0,0]';
site.R = [0,0,0]';
nf  = [1,0,0]';
Rf  = site.R;
Vfm = 0;
tfGuess = 60;81;15.0;

tfmax = 90;
tfmin = 40;

% Constraints

includePointingP = 0 ;
includePointing0 = 0 ;
includePointingF = 1 ;
includeNoFlyZone = 0 ;
nobs    = 1;
Robs1   = 200;
Robs2   = 200;
Pobs1   = [170,400,200]';
Pobs2   = [0,300,0]';
glideSlope  = deg2rad(86);
tiltAngleP  = deg2rad(15);
tiltAngle0  = deg2rad(1);
tiltAngleF  = deg2rad(0);

%% Scaling terms

R_sc = Rmag0*0.5;
V_sc = Vmag0*0.5;

T_sc = (Fmax + Fmin)/2;
a_sc = GravityModel.g0;%(amax + amin)/2;

t_sc = sqrt(norm(R0)*2/GravityModel.g0)*2;


M_sc = (MassModel.massInert+M0)/2;
z_sc = 1;(zL+zU)/2;
%% User Initial Guess


t0_IG   = t0;
tf_IG   = tfGuess;
time_IG = linspace(t0_IG,tf_IG,Nodes);
timeDisc_IG = (tf_IG- t0_IG)/(Nodes-1);


% Acikmese mass profile initial guess
% M_maxThrust = max(M0 - Fmax/ThrustModel.ceffVac*time_IG,MassModel.massInert);
M_maxThrust = M0 - Fmax/ThrustModel.ceffVac*time_IG;
z_maxThrust = log(M_maxThrust);

switch BodyMapRocket.DynamicModel.controlVariable
    case 'thrust'
X_0  = [R0;V0;M0];
X_f  = [Rf;0;0;0;MassModel.massInert];
    case 'acceleration'
X_0  = [R0;V0;z0];
X_f  = [Rf;0;0;0;zL];
end
X_IG = zeros(length(X_0),Nodes);
for ii = 1:7
    X_IG(ii,:)    = linspace(X_0(ii),X_f(ii),Nodes);
end
switch BodyMapRocket.DynamicModel.controlVariable
    case 'thrust'
X_IG(end,:) = M_maxThrust;
    case 'acceleration'
X_IG(end,:) = z_maxThrust;
end

switch BodyMapRocket.DynamicModel.controlVariable
    case 'thrust'
U_IG      = 0.8*repmat(Fmax,4,Nodes).*ones(4,Nodes);
U_IG(1:3,1:end) = + repmat(site.n(:)*(Fmin+Fmax)/2,1,Nodes);
U_IG(1:3,1:end) = + repmat(site.n(:)*Fmax,1,Nodes);
    case 'acceleration'
U_IG      = 0.8*repmat(Fmax,4,Nodes)./exp(X_IG(7,:)).*ones(4,Nodes);
U_IG(1:3,1:end) = + repmat(site.n(:)*(Fmin+Fmax)/2,1,Nodes)./exp(X_IG(7,:));
U_IG(1:3,1:end) = + repmat(site.n(:)*Fmax,1,Nodes)./exp(X_IG(7,:));
end

%% Setup SCOPT

Scopt     = scopt2;

Scopt.problem.nphases                      = 1;
Scopt.problem.phases(1).Nodes              = Nodes;
Scopt.problem.phases(1).n.states           = 7;
Scopt.problem.phases(1).n.controls         = 4;
Scopt.problem.phases(1).n.parameters       = 0;
Scopt.problem.phases(1).n.events           = 1 + includePointing0 + includePointingF;
Scopt.problem.phases(1).n.path             = 2 + includePointingP + ThrustModel.includeThrustRateLimit*2 + includeNoFlyZone*nobs + strcmp(BodyMapRocket.DynamicModel.controlVariable,'acceleration')*2;
Scopt.problem.phases(1).freeTimeFinal      = freeTime*solveSCVx;
Scopt.problem.phases(1).freeTimeInitial    = 0; % Not Implemented

Scopt.setUpSCOPT_level1();
%%  Interface Main Input
index.states.position               = (1:3)' ;
index.states.velocity               = (4:6)' ;
index.states.mass                   = 7 ;
index.controls.T                    = (1:3)' ;
index.controls.TMag                 = 4 ;

% Initial Guess (if user supplied)

Scopt.problem.phases(1).initialGuess.states      = X_IG;
Scopt.problem.phases(1).initialGuess.controls    = U_IG;
Scopt.problem.phases(1).initialGuess.timeInitial = t0_IG;
Scopt.problem.phases(1).initialGuess.timeFinal   = tf_IG;
Scopt.problem.phases(1).initialGuess.time        = time_IG;

% Dynamics
%   this is to use full matrix representation instead of sparse in state and control matrix definition
Scopt.problem.phases(1).dynamics.typeStateMatrix   = 0; 
Scopt.problem.phases(1).dynamics.typeControlMatrix = 0;
%   choice between numerical trapezoidal/euler, or exact dynamics solution
Scopt.algorithm.collocationMethod         = discretisationType;
switch discretisationType
    case {'trapezoidal','euler'} 
    Scopt.problem.phases(1).dynamics.stateDerivativeFunction = @(t,X,U,BodyMap,ii) stateDerivativeRocketLanding(t,X,U,BodyMapRocket);
    Scopt.problem.phases(1).dynamics.stateMatrixFunction     = @(t,X,U,BodyMap,ii) stateMatrixRocketLandingFull(t,X,U,BodyMapRocket);
    Scopt.problem.phases(1).dynamics.controlMatrixFunction   = @(t,X,U,BodyMap,ii) controlMatrixRocketLandingFull(t,X,U,BodyMapRocket);
    case 'exact'
    Scopt.problem.phases(1).dynamics.stateDerivativeFunction = @(t,X,U,BodyMap,ii) stateDerivativeSolutionRocketLanding(t,X,U,BodyMapRocket);
    Scopt.problem.phases(1).dynamics.stateMatrixFunction     = @(t,X,U,BodyMap,ii) stateMatrixSolutionRocketLandingFull(t,X,U,BodyMapRocket);
    Scopt.problem.phases(1).dynamics.controlMatrixFunction   = @(t,X,U,BodyMap,ii) controlMatrixSolutionRocketLandingFull(t,X,U,BodyMapRocket);
    Scopt.problem.phases(1).dynamics.timeDerivativeFunction  = @(t,X,U,BodyMap,ii) timeDerivativeSolutionRocketLanding(t,X,U,BodyMapRocket);
end
% User Scaling
switch BodyMapRocket.DynamicModel.controlVariable
    case 'thrust'
    Scopt.problem.phases(1).scale.states   = [R_sc;R_sc;R_sc;V_sc;V_sc;V_sc;M_sc];
    Scopt.problem.phases(1).scale.controls = [T_sc;T_sc;T_sc;T_sc];
    Scopt.problem.phases(1).scale.time     = t_sc;
    case 'acceleration'
    Scopt.problem.phases(1).scale.states   = [R_sc;R_sc;R_sc;V_sc;V_sc;V_sc;z_sc];
    Scopt.problem.phases(1).scale.controls = [a_sc;a_sc;a_sc;a_sc];
    Scopt.problem.phases(1).scale.time     = t_sc;
end
% Bounds

switch BodyMapRocket.DynamicModel.controlVariable
    case 'thrust'
    Scopt.problem.phases(1).bounds.states.upper    = [Rmag0*2.0;Rmag0*2.0;Rmag0*2.0;Vmag0*1.5;Vmag0*1.5;Vmag0*1.5;MassModel.massWet] ;
    Scopt.problem.phases(1).bounds.states.lower    = [0; -Rmag0*2.0 ;-Rmag0*2.0;-Vmag0*1.5;-Vmag0*1.5;-Vmag0*1.5;MassModel.massInert] ;
    Scopt.problem.phases(1).bounds.controls.upper  = [+Fmax;+Fmax;+Fmax;+Fmax] ;
    Scopt.problem.phases(1).bounds.controls.lower  = [-Fmax;-Fmax;-Fmax;+Fmin] ;
    case 'acceleration'
    Scopt.problem.phases(1).bounds.states.upper    = [Rmag0*2.0;Rmag0*2.0;R0(3)*2.0;Vmag0*1.5;Vmag0*1.5;Vmag0*1.5;zU] ;
    Scopt.problem.phases(1).bounds.states.lower    = [0; -Rmag0*2.0 ;-Rmag0*2.0;-Vmag0*1.5;-Vmag0*1.5;-Vmag0*1.5;zL] ;
    Scopt.problem.phases(1).bounds.controls.upper  = [+amax;+amax;+amax;+amax] ;
    Scopt.problem.phases(1).bounds.controls.lower  = [-amax;-amax;-amax;+amin] ;
end
if freeTime*solveSCVx
Scopt.problem.phases(1).bounds.timeFinal.upper = tfmax ;
Scopt.problem.phases(1).bounds.timeFinal.lower = tfmin ;
end
% Event Constraints

naux = 1;
Scopt.problem.phases(1).events(naux).cone    = Scopt.initCone(4,Scopt.problem.phases(1).n);
Scopt.problem.phases(1).events(naux).type    = 'upper';
Scopt.problem.phases(1).events(naux).funType = 'convex';
Scopt.problem.phases(1).events(naux).where   = 'final';
Scopt.problem.phases(1).events(naux).cone.norm.states(4,1) = 1;
Scopt.problem.phases(1).events(naux).cone.norm.states(5,2) = 1;
Scopt.problem.phases(1).events(naux).cone.norm.states(6,3) = 1;
Scopt.problem.phases(1).events(naux).limit = Vfm;
if includePointing0>=1
    naux = naux+1;
    Scopt.problem.phases(1).events(naux).type        = 'upper';
    Scopt.problem.phases(1).events(naux).funType     = 'linear';
    Scopt.problem.phases(1).events(naux).where       = 'initial';
    Scopt.problem.phases(1).events(naux).controls(4) = n0(1)*cos(tiltAngle0);
    Scopt.problem.phases(1).events(naux).controls(1) = -1;
    Scopt.problem.phases(1).events(naux).limit       = 0;
    % naux = naux+1;
    % Scopt.problem.phases(1).events(naux).type        = 'equal';
    % Scopt.problem.phases(1).events(naux).funType     = 'linear';
    % Scopt.problem.phases(1).events(naux).where       = 'initial';
    % Scopt.problem.phases(1).events(naux).controls(4) = n0(2);
    % Scopt.problem.phases(1).events(naux).controls(2) = -1;
    % Scopt.problem.phases(1).events(naux).limit       = 0;
    % naux = naux+1;
    % Scopt.problem.phases(1).events(naux).type        = 'equal';
    % Scopt.problem.phases(1).events(naux).funType     = 'linear';
    % Scopt.problem.phases(1).events(naux).where       = 'initial';
    % Scopt.problem.phases(1).events(naux).controls(4) = n0(3);
    % Scopt.problem.phases(1).events(naux).controls(3) = -1;
    % Scopt.problem.phases(1).events(naux).limit       = 0;
end
if includePointingF>=1
    naux = naux+1;
    Scopt.problem.phases(1).events(naux).type        = 'upper';
    Scopt.problem.phases(1).events(naux).funType     = 'linear';
    Scopt.problem.phases(1).events(naux).where       = 'final';
    Scopt.problem.phases(1).events(naux).controls(4) = nf(1)*cos(tiltAngleF);
    Scopt.problem.phases(1).events(naux).controls(1) = -1;
    Scopt.problem.phases(1).events(naux).limit       = 0;
    % naux = naux+1;
    % Scopt.problem.phases(1).events(naux).type        = 'equal';
    % Scopt.problem.phases(1).events(naux).funType     = 'linear';
    % Scopt.problem.phases(1).events(naux).where       = 'final';
    % Scopt.problem.phases(1).events(naux).controls(4) = nf(2);
    % Scopt.problem.phases(1).events(naux).controls(2) = -1;
    % Scopt.problem.phases(1).events(naux).limit       = 0;
    % naux = naux+1;
    % Scopt.problem.phases(1).events(naux).type        = 'equal';
    % Scopt.problem.phases(1).events(naux).funType     = 'linear';
    % Scopt.problem.phases(1).events(naux).where       = 'final';
    % Scopt.problem.phases(1).events(naux).controls(4) = nf(3);
    % Scopt.problem.phases(1).events(naux).controls(3) = -1;
    % Scopt.problem.phases(1).events(naux).limit       = 0;
end
% Path Constraints

naux = 1; %norm(A*X+B*U+d)<=CX+R*U+g
Scopt.problem.phases(1).path(naux).cone    = Scopt.initCone(4,Scopt.problem.phases(1).n);
Scopt.problem.phases(1).path(naux).type    = 'upper';
Scopt.problem.phases(1).path(naux).funType = 'convex';
Scopt.problem.phases(1).path(naux).cone.norm.controls(1,1) = 1;
Scopt.problem.phases(1).path(naux).cone.norm.controls(2,2) = 1;
Scopt.problem.phases(1).path(naux).cone.norm.controls(3,3) = 1;
Scopt.problem.phases(1).path(naux).cone.right.controls(4)  = 1;
naux = naux+1;
%
if includePointingP
    Scopt.problem.phases(1).path(naux).type    = 'upper';
    Scopt.problem.phases(1).path(naux).funType = 'linear';
    Scopt.problem.phases(1).path(naux).controls(4)  = cos(tiltAngleP);
    Scopt.problem.phases(1).path(naux).controls(1)  = -1;
    naux = naux+1;
end

Scopt.problem.phases(1).path(naux).cone    = Scopt.initCone(4,Scopt.problem.phases(1).n);
Scopt.problem.phases(1).path(naux).type    = 'upper';
Scopt.problem.phases(1).path(naux).funType = 'convex';
Scopt.problem.phases(1).path(naux).cone.norm.states(1,1) = 1;
Scopt.problem.phases(1).path(naux).cone.norm.states(2,2) = 1;
Scopt.problem.phases(1).path(naux).cone.norm.states(3,3) = 1;
Scopt.problem.phases(1).path(naux).cone.norm.cons        = -site.R(:);
Scopt.problem.phases(1).path(naux).cone.right.states(1:3)= site.n/cos(glideSlope);
Scopt.problem.phases(1).path(naux).cone.right.cons       = -site.R'*site.n/cos(glideSlope);
naux = naux+1;

if ThrustModel.includeThrustRateLimit
    Scopt.problem.phases(1).path(naux).type              = 'upper';
    Scopt.problem.phases(1).path(naux).funType           = 'linear';
    Scopt.problem.phases(1).path(naux).derivative        = 1;
    Scopt.problem.phases(1).path(naux).controls(4)       = 1;  
switch BodyMapRocket.DynamicModel.controlVariable
    case 'thrust'
    Scopt.problem.phases(1).path(naux).limit             = ThrustModel.thrustRateMax;
    case 'acceleration'
    Scopt.problem.phases(1).path(naux).limit             = ThrustModel.thrustRateMax/MassModel.massInert;
end
    naux = naux+1;
    
    Scopt.problem.phases(1).path(naux).type              = 'upper';
    Scopt.problem.phases(1).path(naux).funType           = 'linear';
    Scopt.problem.phases(1).path(naux).derivative        = 1;
    Scopt.problem.phases(1).path(naux).controls(4)       = -1;
switch BodyMapRocket.DynamicModel.controlVariable
    case 'thrust'
    Scopt.problem.phases(1).path(naux).limit             = ThrustModel.thrustRateMax;
    case 'acceleration'
    Scopt.problem.phases(1).path(naux).limit             = ThrustModel.thrustRateMax/MassModel.massInert;
end
    naux = naux+1;
end

if includeNoFlyZone >0
    if nobs>=1
        Scopt.problem.phases(1).path(naux).type              = 'lower';
        Scopt.problem.phases(1).path(naux).funType           = 'non-linear';
        Scopt.problem.phases(1).path(naux).function          = @(t,X,U,BodyMap,ii) vecnorm((X(1:3,:)-repmat(Pobs1,1,size(X,2))),2,1);
        Scopt.problem.phases(1).path(naux).jacobian.states   = @(t,X,U,BodyMap,ii) [(repmat(1./vecnorm((X(1:3,:)-repmat(Pobs1,1,size(X,2))),2,1),3,1)).*(X(1:3,:)-repmat(Pobs1,1,size(X,2)));zeros(4,length(t))];
        Scopt.problem.phases(1).path(naux).limit             = Robs1;
        Scopt.problem.phases(1).path(naux).scale             = 1;
        Scopt.problem.phases(1).path(naux).buffer.include    = 1;
        Scopt.problem.phases(1).path(naux).buffer.penaly     = 1;
        Scopt.problem.phases(1).path(naux).buffer.lambda     = 1e5;
        naux = naux+1;
    end
    if nobs>=2
        Scopt.problem.phases(1).path(naux).type              = 'lower';
        Scopt.problem.phases(1).path(naux).funType           = 'non-linear';
        Scopt.problem.phases(1).path(naux).function          = @(t,X,U,BodyMap,ii) vecnorm((X(1:3,:)-repmat(Pobs2,1,size(X,2))),2,1);
        Scopt.problem.phases(1).path(naux).jacobian.states   = @(t,X,U,BodyMap,ii) [(repmat(1./vecnorm((X(1:3,:)-repmat(Pobs2,1,size(X,2))),2,1),3,1)).*(X(1:3,:)-repmat(Pobs2,1,size(X,2)));zeros(4,length(t))];
        Scopt.problem.phases(1).path(naux).limit             = Robs2;
        Scopt.problem.phases(1).path(naux).scale             = 1;
        Scopt.problem.phases(1).path(naux).buffer.include    = 1;
        Scopt.problem.phases(1).path(naux).buffer.penaly     = 1;
        Scopt.problem.phases(1).path(naux).buffer.lambda     = 1e5;
        naux = naux+1;
    end
end

if strcmp(BodyMapRocket.DynamicModel.controlVariable,'acceleration')
    Scopt.problem.phases(1).path(naux).type                 = 'upper';
    Scopt.problem.phases(1).path(naux).funType              = 'non-linear';
    Scopt.problem.phases(1).path(naux).function             = @(t,X,U,BodyMap,ii) U(4,:) - Fmax*exp(-X(7,:));
    Scopt.problem.phases(1).path(naux).jacobian.states      = @(t,X,U,BodyMap,ii) [zeros(6,length(t));Fmax*exp(-X(7,:))];
    Scopt.problem.phases(1).path(naux).jacobian.controls    = @(t,X,U,BodyMap,ii) [zeros(3,length(t));ones(1,length(t))];
    Scopt.problem.phases(1).path(naux).limit                = 0;
    naux = naux+1;
    Scopt.problem.phases(1).path(naux).limit              = 0;
    Scopt.problem.phases(1).path(naux).type               = 'upper';
    Scopt.problem.phases(1).path(naux).funType            = 'quasiconvex'; 
    Scopt.problem.phases(1).path(naux).controls(4)        = -1; 
    Scopt.problem.phases(1).path(naux).statesIndex        = 7;
    Scopt.problem.phases(1).path(naux).controlsIndex      = [];
    Scopt.problem.phases(1).path(naux).function           = @(t,X,U,BodyMap,ii)   Fmin*exp(-X(7,:));
    Scopt.problem.phases(1).path(naux).jacobian.states    = @(t,X,U,BodyMap,ii) - Fmin*exp(-X(7,:));
    Scopt.problem.phases(1).path(naux).hessian.states     = @(t,X,U,BodyMap,ii)   Fmin*exp(-X(7,:));
    naux = naux+1;
end
% Initial Conditions

Scopt.problem.phases(1).initial.states.equal.index = (1:7)';
if ~isnan(Tmag0)
    Scopt.problem.phases(1).initial.controls.equal.index = 4;
end
Scopt.problem.phases(1).final.states.equal.index   = (1:3)';

switch BodyMapRocket.DynamicModel.controlVariable
    case 'thrust'
        Scopt.problem.phases(1).initial.states.equal.value = [R0;V0;M0];
        if ~isnan(Tmag0)
            Scopt.problem.phases(1).initial.controls.equal.value = Tmag0;
        end
    case 'acceleration'
        Scopt.problem.phases(1).initial.states.equal.value = [R0;V0;z0];
        if ~isnan(Tmag0)
            Scopt.problem.phases(1).initial.controls.equal.value = amag0;
        end
end
Scopt.problem.phases(1).final.states.equal.value   = Rf;

switch objectiveType
    case 0
        Scopt.problem.objective.type          = 'feasibility';
    case 1
        Scopt.problem.objective.type               = 'minimize'; 
        Scopt.problem.objective.mayer.funType      = 'linear';
        Scopt.problem.objective.mayer.where        = 'final';
        Scopt.problem.objective.mayer.states(7)    = -1;
        switch BodyMapRocket.DynamicModel.controlVariable
            case 'thrust'
                Scopt.problem.objective.scale              = (0.5*(MassModel.massWet - MassModel.massInert) + MassModel.massInert);
            case 'acceleration' 
                Scopt.problem.objective.scale              = 1;
        end
    case 2
        Scopt.problem.objective.type               = 'minimize';
        Scopt.problem.objective.scale              = 10*min2sec;
        Scopt.problem.objective.mayer.funType      = 'linear'; % 1 for linear, 2 for convex, 3 for non-convex
        Scopt.problem.objective.mayer.where        = 'final';
        Scopt.problem.objective.mayer.timeFinal    = 1;
end

Scopt.algorithm.initialGuessType          = 'user';
Scopt.algorithm.scaling.variables         = 'automatic';
Scopt.algorithm.scaling.events            = 'automatic-limit'; % Alternative is automatic-limit automatic-jacobian
Scopt.algorithm.scaling.path              = 'automatic-limit'; % Alternative is automatic-limit automatic-jacobian
Scopt.algorithm.scaling.objective         = 'user'; % Alternative is automatic-jacobian
Scopt.algorithm.meshRefinement            = 'none'; % No Mesh Refinement Implemented (YET)

% Scopt.algorithm
if solveSCVx
    Scopt.algorithm.sequential.activate      = 1;
    switch BodyMapRocket.DynamicModel.controlVariable
        case 'thrust'  
            Scopt.algorithm.sequential.activate      = 1;
            Scopt.algorithm.sequential.type          = 'trust-region'; % trust region or line-search
%             Scopt.algorithm.sequential.type          = 'none'; % trust region or line-search
            Scopt.algorithm.sequential.maxIter       = 50;
            Scopt.algorithm.sequential.minIter       = 0;
            Scopt.algorithm.sequential.Jtol          = eps^(1/3);
            Scopt.algorithm.sequential.JtolChange    = 1e-6;
            Scopt.algorithm.sequential.trustRegion.adaptive = 0;
            % Scopt.algorithm.sequential.trustRegion.type     = 'none';
            Scopt.algorithm.sequential.trustRegion.type     = 'hard';
        %     Scopt.algorithm.sequential.trustRegion.type     = 'soft';
            Scopt.algorithm.sequential.trustRegion.include.states   = 1;
            Scopt.algorithm.sequential.trustRegion.include.controls = 1;
            Scopt.algorithm.sequential.trustRegion.variablesPenalty = 2; % 2nd norm, 0 means quadratic norm
            Scopt.algorithm.sequential.trustRegion.nodePenalty      = 1; % 2nd norm
            Scopt.algorithm.sequential.trustRegion.rho0     = 0.0;
            Scopt.algorithm.sequential.trustRegion.rho1     = 0.25;
            Scopt.algorithm.sequential.trustRegion.rho2     = 0.7;
            Scopt.algorithm.sequential.trustRegion.alpha    = 1.2;
            Scopt.algorithm.sequential.trustRegion.beta     = 2.0;
            Scopt.algorithm.sequential.trustRegion.radius     = 1e-2;
        %     Scopt.algorithm.sequential.trustRegion.lambda     = 1e-4;


            Scopt.algorithm.virtualControl.phases(1).states       = 4:6; %
            Scopt.algorithm.virtualControl.phases(1).states       = 1:7; %
        %     Scopt.algorithm.virtualControl.phases(1).states       = 1:7; %
            % Scopt.algorithm.virtualControl.phases(1).scale   = 2/(Nodes-1)*Scopt.problem.phases(1).initialGuess.timeFinal;
            Scopt.algorithm.virtualControl.phases(1).lambda       = 1.0e4;
            Scopt.algorithm.virtualControl.phases(1).lambda       = 1.0e7;
            Scopt.algorithm.virtualControl.phases(1).statePenalty = 2 ; % 1st or 2nd norm, 0 for quadratic
            Scopt.algorithm.virtualControl.phases(1).nodePenalty  = 1 ; % 1st or 2nd norm

            Scopt.algorithm.virtualControl.phases(1).include  = 1 ;
            Scopt.algorithm.virtualBuffer.phases(1).include   = 0 ;


            Scopt.algorithm.sequential.globalisation.type                       = 'line-search'; % line-search
            Scopt.algorithm.sequential.globalisation.activate                   = 0; % line-search
            Scopt.algorithm.sequential.globalisation.lineSearch.type            = 'golden-section'; % line-searchcontraction
            Scopt.algorithm.sequential.globalisation.lineSearch.bracketTol      = 1e-2; % line-search

            if freeTime
            Scopt.algorithm.sequential.trustRegion.phases(1).timeFinal.component = 2;
            Scopt.algorithm.sequential.trustRegion.phases(1).timeFinal.lambda    = 1e-2; 
            Scopt.algorithm.sequential.trustRegion.phases(1).timeFinal.type      = 'hard';
            end
        case 'acceleration' 
            Scopt.algorithm.sequential.type          = 'none';
            if freeTime
                Scopt.algorithm.sequential.trustRegion.phases(1).timeFinal.component = 2;
                Scopt.algorithm.sequential.trustRegion.phases(1).timeFinal.lambda    = 1e-2; 
                Scopt.algorithm.sequential.trustRegion.phases(1).timeFinal.type      = 'soft';
%                 Scopt.algorithm.sequential.trustRegion.phases(1).timeFinal.type      = 'hard';
            end
    end
else
    Scopt.algorithm.sequential.activate      = 0;
end
Scopt2=copy(Scopt);
%% Set Up Level 2
Scopt.setUpSCOPT_level2();

%% Solve Scopt.problem
Scopt.solveSCOPTProblem();
 
%% setUp With Drag

if solveSCVx
    Scopt2.algorithm.sequential.activate      = 1;
%     Scopt2.algorithm.sequential.type          = 'trust-region'; % trust region or line-search
    Scopt2.algorithm.sequential.type          = 'trust-region'; % trust region or line-search
    Scopt2.algorithm.sequential.maxIter       = 50;
    Scopt2.algorithm.sequential.minIter       = 0;
    Scopt2.algorithm.sequential.Jtol          = eps^(1/3);
    Scopt2.algorithm.sequential.JtolChange    = 1e-6;
    Scopt2.algorithm.sequential.trustRegion.adaptive = 1;
    % Scopt2.algorithm.sequential.trustRegion.type     = 'none';
    Scopt2.algorithm.sequential.trustRegion.type     = 'hard';
%     Scopt2.algorithm.sequential.trustRegion.type     = 'soft';
    Scopt2.algorithm.sequential.trustRegion.include.states   = 1;
    Scopt2.algorithm.sequential.trustRegion.include.controls = 1;
    Scopt2.algorithm.sequential.trustRegion.variablesPenalty = 2; % 2nd norm, 0 means quadratic norm
    Scopt2.algorithm.sequential.trustRegion.nodePenalty      = 1; % 2nd norm
    Scopt2.algorithm.sequential.trustRegion.rho0     = 0.0;
    Scopt2.algorithm.sequential.trustRegion.rho1     = 0.25;
    Scopt2.algorithm.sequential.trustRegion.rho2     = 0.7;
    Scopt2.algorithm.sequential.trustRegion.alpha    = 1.2;
    Scopt2.algorithm.sequential.trustRegion.beta     = 2.0;
    Scopt2.algorithm.sequential.trustRegion.radius     = 1e1;
%     Scopt2.algorithm.sequential.trustRegion.lambda     = 1e-4;
    
    
    Scopt2.algorithm.virtualControl.phases(1).states       = 4:6; %
    Scopt2.algorithm.virtualControl.phases(1).states       = 1:7; %
%     Scopt2.algorithm.virtualControl.phases(1).states       = 1:7; %
    % Scopt2.algorithm.virtualControl.phases(1).scale   = 2/(Nodes-1)*Scopt2.problem.phases(1).initialGuess.timeFinal;
    Scopt2.algorithm.virtualControl.phases(1).lambda       = 1.0e4;
    Scopt2.algorithm.virtualControl.phases(1).lambda       = 1.0e2;
    Scopt2.algorithm.virtualControl.phases(1).statePenalty = 2 ; % 1st or 2nd norm, 0 for quadratic
    Scopt2.algorithm.virtualControl.phases(1).nodePenalty  = 1 ; % 1st or 2nd norm
    
    Scopt2.algorithm.virtualControl.phases(1).include  = 1 ;
    Scopt2.algorithm.virtualBuffer.phases(1).include   = 0 ;
    
    
    Scopt2.algorithm.sequential.globalisation.type                       = 'line-search'; % line-search
    Scopt2.algorithm.sequential.globalisation.activate                   = 0; % line-search
    Scopt2.algorithm.sequential.globalisation.lineSearch.type            = 'golden-section'; % line-searchcontraction
    Scopt2.algorithm.sequential.globalisation.lineSearch.bracketTol      = 1e-2; % line-search
    
    if freeTime
    Scopt2.algorithm.sequential.trustRegion.phases(1).timeFinal.component = 2;
    Scopt2.algorithm.sequential.trustRegion.phases(1).timeFinal.lambda    = 1e-2; 
    Scopt2.algorithm.sequential.trustRegion.phases(1).timeFinal.type      = 'soft';
    end
else
    Scopt2.algorithm.sequential.activate      = 0;    
end
Scopt2.solution.generateFigures = 1;
if updateInitialGuessProblem2 
Scopt2.problem.phases(1).initialGuess.states      = Scopt.solution.states;
Scopt2.problem.phases(1).initialGuess.controls    = Scopt.solution.controls;
Scopt2.problem.phases(1).initialGuess.timeInitial = Scopt.solution.timeInitial;
Scopt2.problem.phases(1).initialGuess.timeFinal   = Scopt.solution.timeFinal;
Scopt2.problem.phases(1).initialGuess.time        = Scopt.solution.time;
end
BodyMapRocketDrag = BodyMapRocket;
BodyMapRocketDrag.AerodynamicModel.dragMethod = 'constant';

Scopt2.algorithm.collocationMethod         = discretisationType;
switch discretisationType
    case {'trapezoidal','euler'}
Scopt2.problem.phases(1).dynamics.stateDerivativeFunction = @(t,X,U,BodyMap,ii) stateDerivativeRocketLanding(t,X,U,BodyMapRocketDrag);
Scopt2.problem.phases(1).dynamics.stateMatrixFunction     = @(t,X,U,BodyMap,ii) stateMatrixRocketLandingFull(t,X,U,BodyMapRocketDrag);
Scopt2.problem.phases(1).dynamics.controlMatrixFunction   = @(t,X,U,BodyMap,ii) controlMatrixRocketLandingFull(t,X,U,BodyMapRocketDrag);
    case 'exact' 
Scopt2.problem.phases(1).dynamics.stateDerivativeFunction = @(t,X,U,BodyMap,ii) stateDerivativeSolutionRocketLanding(t,X,U,BodyMapRocketDrag);
Scopt2.problem.phases(1).dynamics.stateMatrixFunction     = @(t,X,U,BodyMap,ii) stateMatrixSolutionRocketLandingFull(t,X,U,BodyMapRocketDrag);
Scopt2.problem.phases(1).dynamics.controlMatrixFunction   = @(t,X,U,BodyMap,ii) controlMatrixSolutionRocketLandingFull(t,X,U,BodyMapRocketDrag);
Scopt2.problem.phases(1).dynamics.timeDerivativeFunction  = @(t,X,U,BodyMap,ii) timeDerivativeSolutionRocketLanding(t,X,U,BodyMapRocketDrag);
end


Scopt2.setUpSCOPT_level2();
Scopt2.solveSCOPTProblem();
%% Retrieve Scopt.solutions

Solution = Scopt.solution;
Problem  = Scopt.problem;
Algorithm = Scopt.algorithm;
Debugging = Scopt.debugging;


Solution2 = Scopt2.solution;

Iter = Solution.iterations;

time_SOL         = Solution.time; 
X_SOL            = Solution.states;
U_SOL            = Solution.controls;
time_SOL2        = Solution2.time; 
X_SOL2           = Solution2.states;
U_SOL2           = Solution2.controls;


time_ALL         = Scopt.debugging.time_ALL; 
X_ALL            = Scopt.debugging.states_ALL;
U_ALL            = Scopt.debugging.controls_ALL;

U_IG               = Scopt.problem.phases.initialGuess.controls  ;
X_IG               = Scopt.problem.phases.initialGuess.states  ;
time_IG            = Scopt.problem.phases.initialGuess.time  ;
timeFinal_IG       = Scopt.problem.phases.initialGuess.timeFinal  ;

%% External Post Processing

lineOpacity_vector = linspace(1,0,Iter);
figure;defaultColor = get(gca,'colororder');close(gcf)

timeInitial_SOL= time_SOL(1);
timeFinal_SOL= time_SOL(end);
timeInitial_SOL2= time_SOL2(1);
timeFinal_SOL2= time_SOL2(end);


R_ALL   = X_ALL(1:3,:,:);
V_ALL   = X_ALL(4:6,:,:);
V_norm_ALL = squeeze(vecnorm(V_ALL,2,1));
R_IG    = X_IG(1:3,:);
V_IG   = X_IG(4:6,:);
V_norm_IG = vecnorm(V_IG,2,1);
R_SOL   = X_SOL(1:3,:);
V_SOL   = X_SOL(4:6,:);
V_norm_SOL = vecnorm(V_SOL,2,1);
R_SOL2   = X_SOL2(1:3,:);
V_SOL2   = X_SOL2(4:6,:);
V_norm_SOL2 = vecnorm(V_SOL2,2,1);

switch BodyMapRocket.DynamicModel.controlVariable
    case 'thrust'
M_ALL   = X_ALL(7,:,:);
Ft_ALL  = U_ALL(1:3,:,:);
Ftm_ALL = U_ALL(4,:,:);
M_IG   = X_IG(7,:);
Ft_IG  = U_IG(1:3,:);
Ftm_IG = U_IG(4,:);
M_SOL   = X_SOL(7,:);

Ft_SOL  = U_SOL(1:3,:);
Ftm_SOL = U_SOL(4,:);

M_SOL2   = X_SOL2(7,:);

Ft_SOL2  = U_SOL2(1:3,:);
Ftm_SOL2 = U_SOL2(4,:);
    case 'acceleration'
z_ALL   = X_ALL(7,:,:);
M_ALL   = exp(z_ALL);
at_ALL  = U_ALL(1:3,:,:);
atm_ALL = U_ALL(4,:,:);
Ft_ALL  = at_ALL.*repmat(M_ALL,3,1,1);
Ftm_ALL = atm_ALL.*M_ALL;
z_IG   = X_IG(7,:);
M_IG   = exp(z_IG);
at_IG  = U_IG(1:3,:);
atm_IG = U_IG(4,:);
Ft_IG  = at_IG.*repmat(M_IG,3,1);
Ftm_IG = atm_IG.*M_IG;

z_SOL   = X_SOL(7,:);
M_SOL   = exp(z_SOL);

at_SOL  = U_SOL(1:3,:);
atm_SOL = U_SOL(4,:);
Ft_SOL   = at_SOL.*repmat(M_SOL,3,1);
Ftm_SOL  = atm_SOL.*M_SOL;

z_SOL2   = X_SOL2(7,:);
M_SOL2   = exp(z_SOL2);
at_SOL2  = U_SOL2(1:3,:);
atm_SOL2 = U_SOL2(4,:);
Ft_SOL2   = at_SOL2.*repmat(M_SOL2,3,1);
Ftm_SOL2  = atm_SOL2.*M_SOL2;
end
F_norm_IG = vecnorm(Ft_IG,2,1);
Ft_norm_SOL = vecnorm(Ft_SOL,2,1);
Ft_norm_SOL2 = vecnorm(Ft_SOL2,2,1);

Ft_norm_ALL = zeros(1,Nodes,Solution.iterations);
for kk  = 1:Solution.iterations
    for ii = 1:Nodes
        Ft_norm_ALL(1,ii,kk) = vecnorm(Ft_ALL(:,ii,kk),2,1);
    end
end
Ft_norm_ALL = squeeze(Ft_norm_ALL);

fprintf('No Drag Case\n' )
fprintf('Mass Consumption %0.5f [kg]\n',M_SOL(1)-M_SOL(end))
fprintf('Time to Go       %0.5f [s]\n',timeFinal_SOL)
fprintf('Drag Case\n' )
fprintf('Mass Consumption %0.5f [kg]\n',M_SOL2(1)-M_SOL2(end))
fprintf('Time to Go       %0.5f [s]\n',timeFinal_SOL2)
%% Plotting Options
lineOpacity_vector = linspace(1,0,Iter+1);
lineOpacity_vector(end) = [];
defaultColor = get(gca,'colororder');close(gcf)
%% Plot 3D PLot
figure

t = [0;tan(glideSlope)];
[X_n,Y_n,Z_n] = cylinder(t);

% Rgo_Vertical=dot(site.n,Rgo_norm_SOL(round(Nodes/2))-site.R);
Rgo_Vertical  = 50;

X   = X_n*Rgo_Vertical;
Y   = Y_n*Rgo_Vertical;
Z   = Z_n*Rgo_Vertical;

% s = surf(X,Y,Z);
% s.FaceAlpha = 0.1;
% vec_s = [0 ; 0 ; 1];
% direction = (cross(vec_s,site.n)); direction = direction/norm(direction);
% angle = acosd(dot(vec_s,site.n)/norm(vec_s)/norm(site.n));
% rotate(s,direction,angle,[0 0 0])
% s.XData = s.XData + site.R(1);
% s.YData = s.YData + site.R(2);
% s.ZData = s.ZData + site.R(3);


hold on


h3dinitial = plot3(R0(1),R0(2),R0(3),'ro','MarkerSize',10); hold on
h3dlanding = plot3(Rf(1),Rf(2),Rf(3),'r*','MarkerSize',10); hold on


h3obs1 = plot3(Pobs1(1),Pobs1(2),Pobs1(3),'g*','MarkerSize',10); hold on
h3obs2 = plot3(Pobs2(1),Pobs2(2),Pobs2(3),'g*','MarkerSize',10); hold on

if includeNoFlyZone
    [x,y,z] = sphere;
    surf(x*Robs1 + Pobs1(1),y*Robs1 + Pobs1(2),z*Robs1 + Pobs1(3), 'FaceColor','none','EdgeColor','k');hold on
    surf(x*Robs2 + Pobs2(1),y*Robs2 + Pobs2(2),z*Robs2 + Pobs2(3), 'FaceColor','none','EdgeColor','k');hold on
end


plot3(R_IG(1,:),R_IG(2,:),R_IG(3,:),'b--','linewidth',1); hold on
for kk = 1:Iter
    plot3(R_ALL(1,:,kk),R_ALL(2,:,kk),R_ALL(3,:,kk),'Color',[lineOpacity_vector(kk) lineOpacity_vector(kk) 1]); hold on
end

h3dtraject = plot3(R_SOL(1,:),R_SOL(2,:),R_SOL(3,:),'b-','linewidth',3); hold on
h3dtraject2 = plot3(R_SOL2(1,:),R_SOL2(2,:),R_SOL2(3,:),'r-','linewidth',3); hold on
axis equal
xlabel('$P_{Up}$ [m]','interpreter','latex')
ylabel('$P_{East}$ [m]','interpreter','latex')
zlabel('$P_{North}$ [m]','interpreter','latex')
grid on
legend([h3dinitial,h3dlanding,h3dtraject,h3dtraject2],{'Initial Condition','Landing Site','No Drag','With Drag'})


quiver3(R_SOL(1,:),R_SOL(2,:),R_SOL(3,:),...
    Ft_SOL(1,:),Ft_SOL(2,:),Ft_SOL(3,:),1,'k','linewidth',1);

view(-45,-45)

saveas(gcf,fullfile(FiguresFolder, '3D View.fig'))
saveas(gcf,fullfile(FiguresFolder, '3D View.jpg'))

%% Plot States
figure
set(gcf,'units','normalized','position',[0 0 0.7 1])
subplot 211 % Positions

plot(time_SOL(1),R0(1),'ro','MarkerSize',5); hold on
plot(time_SOL(1),R0(2),'bo','MarkerSize',5); hold on
plot(time_SOL(1),R0(3),'go','MarkerSize',5); hold on
plot(time_SOL(end),Rf(1),'r*','MarkerSize',5); hold on
plot(time_SOL(end),Rf(2),'b*','MarkerSize',5); hold on
plot(time_SOL(end),Rf(3),'g*','MarkerSize',5); hold on
% plot([time_SOL(1) time_SOL(end)],r_max*R_sc*ones(1,2),'k--','linewidth',1); hold on

hsubplotrx = plot(time_SOL,R_SOL(1,:),'r-','linewidth',2); hold on
hsubplotry = plot(time_SOL,R_SOL(2,:),'b-','linewidth',2); hold on
hsubplotrz = plot(time_SOL,R_SOL(3,:),'g-','linewidth',2); hold on
xlabel('time [s]','interpreter','latex')
ylabel('Position [m]','interpreter','latex')
grid on
legend([hsubplotrx,hsubplotry,hsubplotrz],{'x','y','z'},'location','best')

subplot 212 % Velocities
plot(time_SOL(1),V0(1),'ro','MarkerSize',5); hold on
plot(time_SOL(1),V0(2),'bo','MarkerSize',5); hold on
plot(time_SOL(1),V0(3),'go','MarkerSize',5); hold on
plot(time_SOL(1),norm(V0),'ko','MarkerSize',10); hold on
plot(time_SOL(end),0,'k*','MarkerSize',10); hold on
% plot([time_SOL(1) time_SOL(end)],v_max*Vc_sc*ones(1,2),'k--','linewidth',1); hold on
% plot([time_SOL(1) time_SOL(end)],-v_max*Vc_sc*ones(1,2),'k--','linewidth',1); hold on
% for ii = 1:kk
%     plot(time_ALL(:,ii),R_ALL_norm(:,ii),'Color',lineOpacity_vector(ii)*[1 1 1]); hold on
% end
hsubplotrnorm = plot(time_SOL,V_norm_SOL,'k-','linewidth',3); hold on
hsubplotrx = plot(time_SOL,V_SOL(1,:),'r-','linewidth',2); hold on
hsubplotry = plot(time_SOL,V_SOL(2,:),'b-','linewidth',2); hold on
hsubplotrz = plot(time_SOL,V_SOL(3,:),'g-','linewidth',2); hold on

ylim([-min(max(V_norm_SOL)) , min(max(V_norm_SOL))]*1.2);
xlabel('time [s]','interpreter','latex')
ylabel('Relative Speed [m/s]','interpreter','latex')
grid on
legend([hsubplotrnorm,hsubplotrx,hsubplotry,hsubplotrz],{'norm','x','y','z'},'location','best')

grid on

saveas(gcf,fullfile(FiguresFolder ,'State Profile.jpg'))
saveas(gcf,fullfile(FiguresFolder ,'State Profile.fig'))


figure

hsubplotminT = plot([time_SOL(1) time_SOL(end)],Fmin*ones(1,2)/1000,'m-.','linewidth',1); hold on
hsubplotminT = plot([time_SOL(1) time_SOL(end)],-Fmin*ones(1,2)/1000,'m-.','linewidth',1); hold on
hsubplotmaxT = plot([time_SOL(1) time_SOL(end)],Fmax*ones(1,2)/1000,'m--','linewidth',1); hold on
hsubplotmaxT = plot([time_SOL(1) time_SOL(end)],-Fmax*ones(1,2)/1000,'m--','linewidth',1); hold on
% plot(time_IG,Ft_norm_IG,'k--','linewidth',1); hold on
hsubplotrxIG = plot(time_IG,Ft_IG(1,:)/1000,'r--','linewidth',1); hold on
hsubplotryIG = plot(time_IG,Ft_IG(2,:)/1000,'b--','linewidth',1); hold on
hsubplotrzIG = plot(time_IG,Ft_IG(3,:)/1000,'g--','linewidth',1); hold on
for ii = 1:kk
    plot(time_ALL(:,ii),Ft_norm_ALL(:,ii)/1000,'Color',lineOpacity_vector(ii)*[1 1 1],'linewidth',1); hold on
    plot(time_SOL,Ft_ALL(1,:,ii)/1000,'-','Color',lineOpacity_vector(ii)*[1 0 0],'linewidth',1); hold on
    plot(time_SOL,Ft_ALL(2,:,ii)/1000,'-','Color',lineOpacity_vector(ii)*[0 0 1],'linewidth',1); hold on
    plot(time_SOL,Ft_ALL(3,:,ii)/1000,'-','Color',lineOpacity_vector(ii)*[0 1 0],'linewidth',1); hold on
end
hsubplotrnorm = plot(time_SOL,Ft_norm_SOL/1000,'k-','linewidth',3); hold on
hsubplotrmag = plot(time_SOL,Ftm_SOL/1000,'y--','linewidth',1); hold on
hsubplotrx = plot(time_SOL,Ft_SOL(1,:)/1000,'r-','linewidth',2); hold on
hsubplotry = plot(time_SOL,Ft_SOL(2,:)/1000,'b-','linewidth',2); hold on
hsubplotrz = plot(time_SOL,Ft_SOL(3,:)/1000,'g-','linewidth',2); hold on
ylim([-Fmax Fmax]*1.2/1000)
grid on
xlabel('time [s]','interpreter','latex')
ylabel('Thrust [kN]','interpreter','latex')
legend([hsubplotrnorm,hsubplotrmag,hsubplotrx,hsubplotry,hsubplotrz,hsubplotmaxT,hsubplotminT],{'norm','mag','x','y','z','Max Thrust','Min Thrust'},'location','best')

saveas(gcf,fullfile(FiguresFolder ,'Thrust Profile.jpg'))
saveas(gcf,fullfile(FiguresFolder ,'Thrust Profile.fig'))


figure

hsubplotminT = plot([time_SOL(1) time_SOL(end)],Fmin/FmaxThrottle*ones(1,2),'m-.','linewidth',1); hold on
hsubplotminT = plot([time_SOL(1) time_SOL(end)],-Fmin/FmaxThrottle*ones(1,2),'m-.','linewidth',1); hold on
hsubplotmaxT = plot([time_SOL(1) time_SOL(end)],Fmax/FmaxThrottle*ones(1,2),'m--','linewidth',1); hold on
hsubplotmaxT = plot([time_SOL(1) time_SOL(end)],-Fmax/FmaxThrottle*ones(1,2),'m--','linewidth',1); hold on
% plot(time_IG,Ft_norm_IG,'k--','linewidth',1); hold on
for ii = 1:kk
    plot(time_ALL(:,ii),Ft_norm_ALL(:,ii)/FmaxThrottle,'Color',lineOpacity_vector(ii)*[1 1 1],'linewidth',1); hold on
end
hsubplotrnorm = plot(time_SOL,Ft_norm_SOL/FmaxThrottle,'k-','linewidth',3); hold on
hsubplotrmag = plot(time_SOL,Ftm_SOL/FmaxThrottle,'y--','linewidth',1); hold on
ylim([-Fmax Fmax]*1.2/FmaxThrottle)
grid on
xlabel('time [s]','interpreter','latex')
ylabel('Throttle [-]','interpreter','latex')
legend([hsubplotrnorm,hsubplotrmag,hsubplotmaxT,hsubplotminT],{'norm','mag','Max Thrust','Min Thrust'},'location','best')

saveas(gcf,fullfile(FiguresFolder ,'Throttle Profile.jpg'))
saveas(gcf,fullfile(FiguresFolder ,'Throttle Profile.fig'))