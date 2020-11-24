% Quadrotor path planning unit test
% Verifies correct implementation of Successive convexification,
% particularly
%    - trust region
%    - virtual controls
%    - virtual buffers on path constraints
clc; clear all; close all;

addpath(genpath('../../src'));
%%
min2sec    = 60;

vehicle       = 'vehicle';
planet        = 'planet';

dbstop if error


ProblemType = 1;

if ProblemType==2 || ProblemType==3
    freeTime = 1;
else
    freeTime = 0;
end

%% Environment And Vehicle



% Gravity
Gravity.type    = 0;
Gravity.g       = [-9.81;0;0];

% Drag
Aerodynamic.include  = 1;
Aerodynamic.cd       = 0.5;

ThrustModel.Fmax        = 4.0;
ThrustModel.Fmin        = 1.0;

MassModel.type = 0;
MassModel.mass = 0.3;

Nodes      = 31;
%% Define Body Map

BodyMap.(planet).Gravity                = Gravity;
BodyMap.(vehicle).MassModel             = MassModel;
BodyMap.(vehicle).ThrustModel           = ThrustModel;
BodyMap.(vehicle).Aerodynamic           = Aerodynamic;

BodyMap.CentralBody    = planet  ;
BodyMap.ControlledBody = vehicle ;
BodyMapQuad = BodyMap;
%% Initial, final conditions and constraints

T0 = -Gravity.g*MassModel.mass;
P0 = [0,0,0]';
V0 = [0,0.5,0]';
t0 = 0.0;

Tf = -Gravity.g*MassModel.mass;
Pf = [0,10,0]';
Vf = [0,0.5,0]';
tf = 3.0;
% tf = 4;

% Constraints
nobs    = 2;
Robs1   = 1;
Robs2   = 1;
Pobs1   = [0,3,0.45]';
Pobs2   = [0,7,-0.45]';
tiltAngle = deg2rad(45);
%% User Initial Guess
X_0  = [P0;V0];
X_f  = [Pf;Vf];
X_IG = zeros(length(X_0),Nodes);
for ii = 1:6
X_IG(ii,:)    = linspace(X_0(ii),X_f(ii),Nodes);
end
% for ii = 4:6
% X_IG(ii,2:end)    = (X_IG(ii-3,2:end) - X_IG(ii-3,1:end-1))/((tf-t0)/(Nodes-1));
% end
% 
unitVectorUp     = -Gravity.g/norm(Gravity.g);
unitVectorTarget = (Pf-P0)/norm((Pf-P0));
rot = deg2rad(45);
% U_IG       = zeros(4,Nodes);
U_IG       = ThrustModel.Fmax*ones(4,Nodes);
U_IG(1:3,:)=  ThrustModel.Fmax*repmat(cos(rot)*unitVectorUp + sin(rot)*norm(unitVectorUp)/norm(unitVectorTarget)*unitVectorTarget,1,Nodes);
U_IG(1:3,1)  =T0;
U_IG(1:3,end)=Tf;
% U_IG(3,:)  = linspace(ThrustModel.Fmax,0,Nodes);
% U_IG(3,:)  = ThrustModel.Fmax;
% 
t0_IG   = t0;
tf_IG   = tf;
time_IG = linspace(t0_IG,tf_IG,Nodes);

%% Setup SCOPT

Scopt     = scopt2;

Scopt.problem.nphases                      = 1;
Scopt.problem.phases(1).Nodes              = Nodes;
Scopt.problem.phases(1).n.states           = 6;
Scopt.problem.phases(1).n.controls         = 4;
Scopt.problem.phases(1).n.parameters       = 0;
Scopt.problem.phases(1).n.events           = 0;
Scopt.problem.phases(1).n.path             = 3+nobs;
Scopt.problem.phases(1).freeTimeFinal      = freeTime;
Scopt.problem.phases(1).freeTimeInitial    = 0; % Not Implemented

Scopt.setUpSCOPT_level1();
%%  Interface Main Input
index.states.position             = (1:3)' ;
index.states.velocity             = (4:6)' ;
index.controls.T                  = (1:3)' ;
index.controls.TMag               = 4 ;

% Initial Guess (if user supplied)

Scopt.problem.phases(1).initialGuess.states      = X_IG;
Scopt.problem.phases(1).initialGuess.controls    = U_IG;
Scopt.problem.phases(1).initialGuess.timeInitial = t0_IG;
Scopt.problem.phases(1).initialGuess.timeFinal   = tf_IG;
Scopt.problem.phases(1).initialGuess.time        = time_IG;

% Dynamics
Scopt.problem.phases(1).dynamics.stateDerivativeFunction = @(t,X,U,BodyMap,ii) stateDerivativeQuadrotor3D(t,X,U,BodyMapQuad,ii);
Scopt.problem.phases(1).dynamics.stateMatrixFunction     = @(t,X,U,BodyMap,ii) stateMatrixQuadrotor3D(t,X,U,BodyMapQuad,ii);
Scopt.problem.phases(1).dynamics.controlMatrixFunction   = @(t,X,U,BodyMap,ii) controlMatrixQuadrotor3D(t,X,U,BodyMapQuad,ii);

%Bounds
Scopt.problem.phases(1).bounds.states.upper    = [0;+10;+5;0;10;10] ;
Scopt.problem.phases(1).bounds.states.lower    = [0; 0 ;-5;-0;-10;-10] ;
Scopt.problem.phases(1).bounds.controls.upper  = [+ThrustModel.Fmax;+ThrustModel.Fmax;+ThrustModel.Fmax;+ThrustModel.Fmax] ;
Scopt.problem.phases(1).bounds.controls.lower  = [-ThrustModel.Fmax;-ThrustModel.Fmax;-ThrustModel.Fmax;+ThrustModel.Fmin] ;
% Scopt.problem.phases(1).bounds.timeFinal.upper = 10 ;
% Scopt.problem.phases(1).bounds.timeFinal.lower = 1 ;

% Path Constraints

naux = 1;
Scopt.problem.phases(1).path(naux).cone    = Scopt.initCone(4,Scopt.problem.phases(1).n);
Scopt.problem.phases(1).path(naux).type    = 'upper';
Scopt.problem.phases(1).path(naux).funType = 'convex';
Scopt.problem.phases(1).path(naux).cone.norm.controls(1,1) = 1;
Scopt.problem.phases(1).path(naux).cone.norm.controls(2,2) = 1;
Scopt.problem.phases(1).path(naux).cone.norm.controls(3,3) = 1;
Scopt.problem.phases(1).path(naux).cone.right.controls(4)  = 1;
naux = naux+1;
% 
Scopt.problem.phases(1).path(naux).type    = 'upper';
Scopt.problem.phases(1).path(naux).funType = 'linear';
Scopt.problem.phases(1).path(naux).controls(4)  = cos(tiltAngle);
Scopt.problem.phases(1).path(naux).controls(1)  = -1;
naux = naux+1;


Scopt.problem.phases(1).path(naux).type    = 'equal';
Scopt.problem.phases(1).path(naux).funType = 'linear';
Scopt.problem.phases(1).path(naux).states(1)  = 1;
naux = naux+1;

if nobs >0
Scopt.problem.phases(1).path(naux).type              = 'upper';
Scopt.problem.phases(1).path(naux).funType           = 'non-linear';
Scopt.problem.phases(1).path(naux).function          = @(t,X,U,BodyMap,ii) -vecnorm((X(1:3,:)-repmat(Pobs1,1,size(X,2))),2,1)+Robs1;
Scopt.problem.phases(1).path(naux).jacobian.states   = @(t,X,U,BodyMap,ii) -[(repmat(1./vecnorm((X(1:3,:)-repmat(Pobs1,1,size(X,2))),2,1),3,1)).*(X(1:3,:)-repmat(Pobs1,1,size(X,2)));zeros(3,length(t))];
Scopt.problem.phases(1).path(naux).scale             = 1;
Scopt.problem.phases(1).path(naux).buffer.include    = 1;
Scopt.problem.phases(1).path(naux).buffer.penalty     = 1;
Scopt.problem.phases(1).path(naux).buffer.lambda     = 1e3;
naux = naux+1;
Scopt.problem.phases(1).path(naux).type              = 'upper';
Scopt.problem.phases(1).path(naux).funType           = 'non-linear';
Scopt.problem.phases(1).path(naux).function          = @(t,X,U,BodyMap,ii) -vecnorm((X(1:3,:)-repmat(Pobs2,1,size(X,2))),2,1)+Robs2;
Scopt.problem.phases(1).path(naux).jacobian.states   = @(t,X,U,BodyMap,ii) -[(repmat(1./vecnorm((X(1:3,:)-repmat(Pobs2,1,size(X,2))),2,1),3,1)).*(X(1:3,:)-repmat(Pobs2,1,size(X,2)));zeros(3,length(t))];
Scopt.problem.phases(1).path(naux).scale             = 1;
Scopt.problem.phases(1).path(naux).buffer.include    = 1;
Scopt.problem.phases(1).path(naux).buffer.penalty     = 1;
Scopt.problem.phases(1).path(naux).buffer.lambda     = 1e3;
naux = naux+1;
end

% d(atm*exp(q))dt <= rate
% Scopt.problem.phases(1).path(naux).type              = 'upper';
% Scopt.problem.phases(1).path(naux).funType           = 'linear';
% Scopt.problem.phases(1).path(naux).derivative        = 1;
% Scopt.problem.phases(1).path(naux).controls(4)       = 1;
% Scopt.problem.phases(1).path(naux).limit             = Orientation.maxAccelerationRate;
% Scopt.problem.phases(1).path(naux).scale             = Orientation.maxAccelerationRate*timeDisc_IG;
% % Scopt.problem.phases(1).path(naux).function          = @(t,X,U,BodyMap,ii) U(4,:).*exp(X(7,:));
% % Scopt.problem.phases(1).path(naux).jacobian.states   = @(t,X,U,BodyMap,ii) [zeros(6,length(t));U(4,:).*exp(X(7,:))];
% % Scopt.problem.phases(1).path(naux).jacobian.controls = @(t,X,U,BodyMap,ii) [zeros(3,length(t));exp(X(7,:))];
% % Scopt.problem.phases(1).path(naux).limit             = Thrusters.maxThrustRate;
% % Scopt.problem.phases(1).path(naux).scale             = Thrusters.maxThrustRate*timeDisc_IG;
% naux = naux+1;

% Scopt.problem.phases(1).path(naux).type              = 'lower';
% Scopt.problem.phases(1).path(naux).funType           = 'linear';
% Scopt.problem.phases(1).path(naux).derivative        = 1;
% Scopt.problem.phases(1).path(naux).controls(4)       = 1;
% Scopt.problem.phases(1).path(naux).limit             = -Orientation.maxAccelerationRate;
% Scopt.problem.phases(1).path(naux).scale             = Orientation.maxAccelerationRate*timeDisc_IG;
% % Scopt.problem.phases(1).path(naux).function          = @(t,X,U,BodyMap,ii) U(4,:).*exp(X(7,:));
% % Scopt.problem.phases(1).path(naux).jacobian.states   = @(t,X,U,BodyMap,ii) [zeros(6,length(t));U(4,:).*exp(X(7,:))];
% % Scopt.problem.phases(1).path(naux).jacobian.controls = @(t,X,U,BodyMap,ii) [zeros(3,length(t));exp(X(7,:))];
% % Scopt.problem.phases(1).path(naux).limit             = -Thrusters.maxThrustRate;
% % Scopt.problem.phases(1).path(naux).scale             = Thrusters.maxThrustRate*timeDisc_IG;
% naux = naux+1;

% Initial Conditions

Scopt.problem.phases(1).initial.states.equal.index = (1:6)';
Scopt.problem.phases(1).initial.states.equal.value = [P0;V0];
Scopt.problem.phases(1).initial.controls.equal.index = (1:3)';
Scopt.problem.phases(1).initial.controls.equal.value = T0;
Scopt.problem.phases(1).final.states.equal.index   = (1:6)';
Scopt.problem.phases(1).final.states.equal.value   = [Pf;Vf];
Scopt.problem.phases(1).final.controls.equal.index = (1:3)';
Scopt.problem.phases(1).final.controls.equal.value = Tf;

switch ProblemType
    case 0
        Scopt.problem.objective.type          = 'feasibility';
    case {1,3}
        Scopt.problem.objective.type                    = 'minimize';
        Scopt.problem.objective.scale                   = 1 ;
        Scopt.problem.objective.lagrange.funType        = 'linear';
        Scopt.problem.objective.lagrange.controls(4)    = 1;
    case 2
        Scopt.problem.objective.type                    = 'minimize';
        Scopt.problem.objective.scale                   = 3*min2sec;
        Scopt.problem.objective.mayer.funType           = 'linear'; % 1 for linear, 2 for convex, 3 for non-convex
        Scopt.problem.objective.mayer.where             = 'final';
        Scopt.problem.objective.mayer.timeFinal         = 1;
end

% Scopt.algorithm

Scopt.algorithm.conicSolver.type         = 'ecos';
Scopt.algorithm.conicSolver.maxIter      = 100;
Scopt.algorithm.conicSolver.feasTol       = eps^(1.0/2);
Scopt.algorithm.conicSolver.relTol        = eps^(1.0/2);
Scopt.algorithm.conicSolver.absTol        = eps^(1.0/2);

Scopt.algorithm.sequential.activate      = 1;
Scopt.algorithm.sequential.type          = 'trust-region'; % trust region or line-search
% Scopt.algorithm.sequential.type          = 'none'; % trust region or line-search
Scopt.algorithm.sequential.maxIter       = 20; 
Scopt.algorithm.sequential.minIter       = 0;
Scopt.algorithm.sequential.Jtol          = 0;
Scopt.algorithm.sequential.JtolChange    = 1e-3;


Scopt.algorithm.sequential.trustRegion.adaptive = 1;
Scopt.algorithm.sequential.trustRegion.type     = 'hard';
Scopt.algorithm.sequential.trustRegion.include.states   = 1;
Scopt.algorithm.sequential.trustRegion.include.controls = 1;
Scopt.algorithm.sequential.trustRegion.variablesPenalty = 1; % Problem uses 1st norm. 2nd norm, 0 means quadratic norm
Scopt.algorithm.sequential.trustRegion.nodePenalty      = inf; % 2nd norm
Scopt.algorithm.sequential.trustRegion.defectsCompErrorPenalty  = 1   ; % 2nd norm, 0 means quadratic norm This is for the computation of the adaptive trust region error
Scopt.algorithm.sequential.trustRegion.defectsNodeErrorPenalty  = 1   ; % 2nd norm This is for the computation of the adaptive trust region error
Scopt.algorithm.sequential.trustRegion.include.pathAndEventsErrors = 1;0;
Scopt.algorithm.sequential.trustRegion.rho0     = 0.0;
Scopt.algorithm.sequential.trustRegion.rho1     = 0.25;
Scopt.algorithm.sequential.trustRegion.rho2     = 0.7;
Scopt.algorithm.sequential.trustRegion.alpha    = 2.0;
Scopt.algorithm.sequential.trustRegion.beta     = 3.2;
Scopt.algorithm.initialGuessType          = 'user';
Scopt.algorithm.scaling.variables         = 'user';
Scopt.algorithm.collocationMethod         = 'trapezoidal';
Scopt.algorithm.scaling.events            = 'user'; % Alternative is automatic-limit automatic-jacobian
Scopt.algorithm.scaling.path              = 'user'; % Alternative is automatic-limit automatic-jacobian
Scopt.algorithm.scaling.objective         = 'user'; % Alternative is automatic-jacobian
Scopt.algorithm.scaling.trustRegion       = 'none'; % Alternative is automatic-jacobian
Scopt.algorithm.meshRefinement            = 'none'; % No Mesh Refinement Implemented (YET)

Scopt.algorithm.virtualControl.phases(1).states       = 1:6; % 
Scopt.algorithm.virtualControl.phases(1).lambda       = 1.0e2; % Problem uses 1e5
Scopt.algorithm.virtualControl.phases(1).statePenalty = 1 ; % 1st or 2nd norm, 0 for quadratic
Scopt.algorithm.virtualControl.phases(1).nodePenalty  = 1 ; % 1st or 2nd norm

Scopt.algorithm.virtualControl.phases(1).include  = 1 ;
Scopt.algorithm.virtualBuffer.phases(1).include   = 1 ;

Scopt.algorithm.sequential.trustRegion.radius     = 1.0e0; % % ADAPT WHEN SCALING, THIS SHOULD BE THE UNSCALED RADIUS
Scopt.algorithm.sequential.trustRegion.lambda     = 1e2;


Scopt.algorithm.sequential.globalisation.type                       = 'line-search'; % line-search
Scopt.algorithm.sequential.globalisation.activate                   = 1; % line-search
Scopt.algorithm.sequential.globalisation.lineSearch.type            = 'golden-section'; % line-searchcontraction
Scopt.algorithm.sequential.globalisation.lineSearch.bracketTol      = 1e-2; % line-search
Scopt.algorithm.sequential.globalisation.lineSearch.stepLengthMax   = 1.0; % line-search
Scopt.algorithm.sequential.globalisation.lineSearch.forceIfRejected = 1; % line-search


Scopt.setUpSCOPT_level2();
%% Solve Scopt.problem
Scopt.solveSCOPTProblem();
%% Retrieve Scopt.solutions

Solution = Scopt.solution;
Problem  = Scopt.problem;
Algorithm = Scopt.algorithm;
Debugging = Scopt.debugging;

Iter = Solution.iterations;

time_SOL         = Solution.time;
timeFinal_SOL    = Solution.timeFinal;
X_SOL            = Solution.states;
U_SOL            = Solution.controls;
time_n_SOL       = Solution.time_n;
timeFinal_n_SOL  = Solution.timeFinal_n;
X_n_SOL          = Solution.states_n;
U_n_SOL          = Solution.controls_n;

time_ALL         = Scopt.debugging.time_ALL;
timeFinal_ALL    = Scopt.debugging.timeFinal_ALL;
X_ALL            = Scopt.debugging.states_ALL;
U_ALL            = Scopt.debugging.controls_ALL;
time_n_ALL       = Scopt.debugging.time_n_ALL;
timeFinal_n_ALL  = Scopt.debugging.timeFinal_n_ALL;
X_n_ALL          = Scopt.debugging.states_n_ALL;
U_n_ALL          = Scopt.debugging.controls_n_ALL;

VC_ALL           = Scopt.debugging.virtualControls   ;
VCP1_ALL         = Scopt.debugging.virtualControlsP1 ;
VCP1s_ALL        = Scopt.debugging.virtualControlsP1s;
VCP2_ALL         = Scopt.debugging.virtualControlsP2 ;
trustRegionP1_ALL= Scopt.debugging.trustRegionP1 ;
trustRegionP2_ALL= Scopt.debugging.trustRegionP2 ;


VC_norm_ALL = zeros(1,Nodes-1,Solution.iterations);
for kk  = 1:Solution.iterations
    for ii = 1:Nodes-1
        VC_norm_ALL(1,ii,kk) = vecnorm(VC_ALL(:,ii,kk),2,1);
    end
end
VC_norm_ALL = squeeze(VC_norm_ALL);

VC_SOL           = VC_ALL(:,:,end);
VCP1_SOL         = VCP1_ALL(:,:,end);
VCP1s_SOL        = VCP1s_ALL(:,:,end);
VCP2_SOL         = VCP2_ALL(:,end);

U_IG             = Scopt.problem.phases.initialGuess.controls  ;
X_IG             = Scopt.problem.phases.initialGuess.states  ;
time_IG          = Scopt.problem.phases.initialGuess.time  ;
timeFinal_IG     = Scopt.problem.phases.initialGuess.timeFinal  ;
U_n_IG           = Scopt.problem.phases.initialGuess.controls_n  ;
X_n_IG           = Scopt.problem.phases.initialGuess.states_n  ;
time_n_IG        = Scopt.problem.phases.initialGuess.time_n  ;
timeFinal_n_IG   = Scopt.problem.phases.initialGuess.timeFinal_n  ;

%% External Post Processing

lineOpacity_vector = linspace(1,0,kk);
figure;defaultColor = get(gca,'colororder');close(gcf)

timeInitial_SOL= time_SOL(1);
timeFinal_SOL= time_SOL(end);


r_ALL   = X_n_ALL(1:3,:,:); 
v_ALL   = X_n_ALL(4:6,:,:); 
r_SOL   = X_n_SOL(1:3,:); 
v_SOL   = X_n_SOL(4:6,:); 
r_IG    = X_n_IG(1:3,:); 
v_IG    = X_n_IG(4:6,:); 

P_ALL   = X_ALL(1:3,:,:); 
V_ALL   = X_ALL(4:6,:,:); 
V_norm_ALL = squeeze(vecnorm(V_ALL,2,1));
Ft_ALL  = U_ALL(1:3,:,:); 
Ftm_ALL = U_ALL(4,:,:); 
P_IG    = X_IG(1:3,:); 
V_IG   = X_IG(4:6,:); 
V_norm_IG = vecnorm(V_IG,2,1);
Ft_IG  = U_IG(1:3,:); 
Ftm_IG = U_IG(4,:); 
F_norm_IG = vecnorm(Ft_IG,2,1);
P_SOL   = X_SOL(1:3,:); 
V_SOL   = X_SOL(4:6,:); 
V_norm_SOL = vecnorm(V_SOL,2,1);

Ft_SOL  = U_SOL(1:3,:); 
Ftm_SOL = U_SOL(4,:); 
Ft_norm_SOL = vecnorm(Ft_SOL,2,1);

Ft_norm_ALL = zeros(1,Nodes,Solution.iterations);
for kk  = 1:Solution.iterations
    for ii = 1:Nodes
        Ft_norm_ALL(1,ii,kk) = vecnorm(Ft_ALL(:,ii,kk),2,1);
    end
end
Ft_norm_ALL = squeeze(Ft_norm_ALL);

%% Plotting Options
lineOpacity_vector = linspace(0,1,kk+1);
lineOpacity_vector(1) = [];
defaultColor = get(gca,'colororder');close(gcf)
%% Plot 3D PLot
figure(5)


h3dinitial = plot3(P0(1),P0(2),P0(3),'ro','MarkerSize',10); hold on
h3dlanding = plot3(Pf(1),Pf(2),Pf(3),'r*','MarkerSize',10); hold on


plot3(Pobs1(1),Pobs1(2),Pobs1(3),'g*','MarkerSize',10); hold on
plot3(Pobs2(1),Pobs2(2),Pobs2(3),'g*','MarkerSize',10); hold on


[x,y,z] = sphere;
h3obs1 = surf(x*Robs1 + Pobs1(1),y*Robs1 + Pobs1(2),z*Robs1 + Pobs1(3), 'FaceColor','r','EdgeColor','k');hold on
h3obs1.FaceAlpha = 0.1;
h3obs2 = surf(x*Robs2 + Pobs2(1),y*Robs2 + Pobs2(2),z*Robs2 + Pobs2(3), 'FaceColor','y','EdgeColor','k');hold on
h3obs2.FaceAlpha = 0.1;


h3dtrajectIG = plot3(P_IG(1,:),P_IG(2,:),P_IG(3,:),'y--','linewidth',2); hold on
for kk = 1:Iter
    plot3(P_ALL(1,:,kk),P_ALL(2,:,kk),P_ALL(3,:,kk),'linewidth',2.0,'Color',[1 1 1]*(1-lineOpacity_vector(kk)) + lineOpacity_vector(kk)*[0 0 1]); hold on
end

h3dtraject = plot3(P_SOL(1,:),P_SOL(2,:),P_SOL(3,:),'b-','linewidth',3); hold on
axis equal
xlabel('$P_x$ [m]','interpreter','latex')
ylabel('$P_y$ [m]','interpreter','latex')
zlabel('$P_z$ [m]','interpreter','latex')
grid on

quiver3(P_SOL(1,:),P_SOL(2,:),P_SOL(3,:),...
        Ft_SOL(1,:),Ft_SOL(2,:),Ft_SOL(3,:),0.5,'k','linewidth',1.5);

if ~isa(h3dinitial,'double')
    legend([h3dinitial,h3dlanding,h3dtraject,h3dtrajectIG,h3obs1,h3obs2],{'Initial Position','Final Position','Trajectory','Initial Guess','$R_{obs,1}$','$R_{obs,2}$'},'interpreter','latex')
end
view(17.5923,6.2920);

% upFig(gcf,['verification_quadRotorPlanning3D_trajectory'],FiguresFolder,'.pdf','wide')
% upFig(gcf,['verification_quadRotorPlanning3D_trajectory'],FiguresFolder,'.pdf','wide-geo',1)
%% Plot States
figure(6)
set(gcf,'units','normalized','position',[0 0 0.7 1])
subplot 211 % Positions

plot(time_SOL(1),P0(1),'ro','MarkerSize',5); hold on
plot(time_SOL(1),P0(2),'bo','MarkerSize',5); hold on
plot(time_SOL(1),P0(3),'go','MarkerSize',5); hold on
plot(time_SOL(end),Pf(1),'r*','MarkerSize',5); hold on
plot(time_SOL(end),Pf(2),'b*','MarkerSize',5); hold on
plot(time_SOL(end),Pf(3),'g*','MarkerSize',5); hold on
% plot([time_SOL(1) time_SOL(end)],r_max*R_sc*ones(1,2),'k--','linewidth',1); hold on

hsubplotrx = plot(time_SOL,P_SOL(1,:),'r-','linewidth',2); hold on
hsubplotry = plot(time_SOL,P_SOL(2,:),'b-','linewidth',2); hold on
hsubplotrz = plot(time_SOL,P_SOL(3,:),'g-','linewidth',2); hold on
xlabel('time [s]','interpreter','latex')
ylabel('Position [m]','interpreter','latex')
grid on
legend([hsubplotrx,hsubplotry,hsubplotrz],{'x','y','z'},'location','best')

subplot 212 % Velocities
plot(time_SOL(1),V0(1),'ro','MarkerSize',5); hold on
plot(time_SOL(1),V0(2),'bo','MarkerSize',5); hold on
plot(time_SOL(1),V0(3),'go','MarkerSize',5); hold on
plot(time_SOL(1),norm(V0),'ko','MarkerSize',10); hold on
plot(time_SOL(end),norm(Vf),'k*','MarkerSize',10); hold on
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

% saveas(gcf,[FiguresFolder 'State Profile.fig'])


figure(7)

hsubplotminT = plot([time_SOL(1) time_SOL(end)],ThrustModel.Fmin*ones(1,2),'m-.','linewidth',1); hold on
hsubplotminT = plot([time_SOL(1) time_SOL(end)],-ThrustModel.Fmin*ones(1,2),'m-.','linewidth',1); hold on
hsubplotmaxT = plot([time_SOL(1) time_SOL(end)],ThrustModel.Fmax*ones(1,2),'m--','linewidth',1); hold on
hsubplotmaxT = plot([time_SOL(1) time_SOL(end)],-ThrustModel.Fmax*ones(1,2),'m--','linewidth',1); hold on
% plot(time_IG,Ft_norm_IG,'k--','linewidth',1); hold on
hsubplotrxIG = plot(time_IG,Ft_IG(1,:),'r--','linewidth',1); hold on
hsubplotryIG = plot(time_IG,Ft_IG(2,:),'b--','linewidth',1); hold on
hsubplotrzIG = plot(time_IG,Ft_IG(3,:),'g--','linewidth',1); hold on
% for ii = 1:kk
%     plot(time_ALL(1,:,ii),Ft_norm_ALL(:,ii),'Color',[1 1 1]*(1-lineOpacity_vector(kk)) + lineOpacity_vector(kk)*[0 0 0],'linewidth',1); hold on
%     plot(time_ALL(1,:,ii),Ft_ALL(1,:,ii),'-','Color',[1 1 1]*(1-lineOpacity_vector(kk)) + lineOpacity_vector(kk)*[1 0 0],'linewidth',1); hold on
%     plot(time_ALL(1,:,ii),Ft_ALL(2,:,ii),'-','Color',[1 1 1]*(1-lineOpacity_vector(kk)) + lineOpacity_vector(kk)*[0 0 1],'linewidth',1); hold on
%     plot(time_ALL(1,:,ii),Ft_ALL(3,:,ii),'-','Color',[1 1 1]*(1-lineOpacity_vector(kk)) + lineOpacity_vector(kk)*[0 1 0],'linewidth',1); hold on
% end
hsubplotrnorm = plot(time_SOL,Ft_norm_SOL,'k-','linewidth',3); hold on
% hsubplotrmag = plot(time_SOL,Ftm_SOL,'y--','linewidth',1); hold on
hsubplotrx = plot(time_SOL,Ft_SOL(1,:),'r-','linewidth',2); hold on
hsubplotry = plot(time_SOL,Ft_SOL(2,:),'b-','linewidth',2); hold on
hsubplotrz = plot(time_SOL,Ft_SOL(3,:),'g-','linewidth',2); hold on
ylim([-ThrustModel.Fmax ThrustModel.Fmax]*1.2)
grid on
xlabel('time [s]','interpreter','latex')
ylabel('Thrust [N]','interpreter','latex')
legend([hsubplotrnorm,hsubplotrx,hsubplotry,hsubplotrz,hsubplotmaxT,hsubplotminT],{'$\Gamma$','$F_1$','$F_2$','$F_3$','Max Thrust','Min Thrust'},'location','best','interpreter','latex')

% upFig(gcf,['verification_quadRotorPlanning3D_thrustProfile'],FiguresFolder,'.pdf')

%% Internal PostProcessing Tools

% Scopt.generateFigures();
%%
Scopt.plotConvergenceHistory();
% upFig(gcf,['verification_quadRotorPlanning3D_convergence'],FiguresFolder,'.pdf')


[~,~,~,hfig2] =  Scopt.plotVirtualControlHistory(40,41);
% upFig(hfig2,['verification_quadRotorPlanning3D_virtualControlsPn'],FiguresFolder,'.pdf')

[~,~,hfig1,hfig2] =  Scopt.plotVirtualControlErrorHistory(44,45);
figure(hfig1)
% upFig(hfig1,['verification_quadRotorPlanning3D_virtualControlsPcError'],FiguresFolder,'.pdf')
figure(hfig2)
% upFig(hfig2,['verification_quadRotorPlanning3D_virtualControlsPnError'],FiguresFolder,'.pdf')

[~,~,hfig1,hfig2] =  Scopt.plotTrustRegionHistory(48,49);
figure(hfig1)
% upFig(hfig1,['verification_quadRotorPlanning3D_trustRegionPc'],FiguresFolder,'.pdf')
figure(hfig2)
% upFig(hfig2,['verification_quadRotorPlanning3D_trustRegionPn'],FiguresFolder,'.pdf')

[~,~,hfig1,hfig2] =  Scopt.plotTrustRegionErrorHistory(50,51);
figure(hfig1)
% upFig(hfig1,['verification_quadRotorPlanning3D_trustRegionPcError'],FiguresFolder,'.pdf')
figure(hfig2)
% upFig(hfig2,['verification_quadRotorPlanning3D_trustRegionPnError'],FiguresFolder,'.pdf')


[hfig1,~,hfig3] =  Scopt.plotVerificationDefectConstraints(60,61,62);
figure(hfig1)
% upFig(hfig1,['verification_quadRotorPlanning3D_defectij'],FiguresFolder,'.pdf')
figure(hfig3)
% upFig(hfig3,['verification_quadRotorPlanning3D_totalDefect'],FiguresFolder,'.pdf')

[~,~,~,hfig1,hfig2] =  Scopt.plotVirtualBuffersHistory(70,71);
figure(hfig1)
% upFig(hfig1,['verification_quadRotorPlanning3D_virtualBuffersPc'],FiguresFolder,'.pdf','narrow')
figure(hfig2)
% upFig(hfig2,['verification_quadRotorPlanning3D_virtualBuffersPn'],FiguresFolder,'.pdf')

[~,~,~,hfig1,hfig2] =  Scopt.plotVirtualBuffersErrorHistory(72,73);
figure(hfig1)
% upFig(hfig1,['verification_quadRotorPlanning3D_virtualBuffersPcError'],FiguresFolder,'.pdf')
figure(hfig2)
% upFig(hfig2,['verification_quadRotorPlanning3D_virtualBuffersPnError'],FiguresFolder,'.pdf')


hfig1 =  Scopt.plotAdaptiveTrustRegionMonitor(80);
set(gcf,'units','normalized','position',[0 0 0.35 1])
% upFig(hfig1,['verification_quadRotorPlanning3D_adaptiveTrustRegionMonitor'],FiguresFolder,'.pdf','3wide')