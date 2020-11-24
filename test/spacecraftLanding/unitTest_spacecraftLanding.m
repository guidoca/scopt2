% Spacecraft Landing with fixed final mass unit test. problem is fully
% convex and can be solved with one iteration. Solved with free final time
% to check for differences
% A full Description is available in TU Delft Repository:
%  'Optimum On-board Abort Guidance based on Successive Convexification
%  for Atmospheric Re-Entry Vehicles', Master Thesis, Guillermo
%  Joaquin Dominguez Calabuig. 
% <https://repository.tudelft.nl/islandora/object/uuid%3A3a33c54b-ed54-407a-88ce-c9d0b3d8f738?collection=education>
 
clc; clear all; close all;
addpath(genpath('../../src'));


vehicle       = 'Spacecraft';
planet        = 'Mars';

dbstop if error

%% Environment And Vehicle


P0 = [50,50,100]';
V0 = [-10,0,-10]';

t0_IG   = 0;

tanGlide = 0.5;

% Gravity
Gravity.type    = 0;
Gravity.g       = 0.1;
Gravity.vector  = [0;0;-1];

ThrustModel.Fmax   = 10;
ThrustModel.gamma  = 1;

MassModel.type = 0;
MassModel.mass = 10;

Nodes      = 36;
%% Define Body Map
BodyMap.(planet).Gravity              = Gravity;
BodyMap.(vehicle).ThrustModel           = ThrustModel;
BodyMap.(vehicle).MassModel           = MassModel;
BodyMap.CentralBody    = planet  ;
BodyMap.ControlledBody = vehicle ;
%% Initial Guess
X_0  = [P0;V0];
X_IG = zeros(length(X_0),Nodes);
for ii = 1:length(X_0)
X_IG(ii,:)    = linspace(X_0(ii),0,Nodes);
end
% 
unitVectorTarget = ([0;0;0]-P0)/norm(([0;0;0]-P0));
U_IG    = -repmat(unitVectorTarget*ThrustModel.Fmax*3/4,1,Nodes);

%%
problems = [1 2 3 3];

tf_IG   = [35 35 35 35];

NumberProblems = length(problems);



X            = zeros(6,Nodes,NumberProblems);
U            = zeros(3,Nodes,NumberProblems);
Jreal        = zeros(1,NumberProblems);
timeFinal    = zeros(1,NumberProblems);
timeElapsed  = zeros(1,NumberProblems);

for ii = 1:NumberProblems

ProblemType = problems(ii);

if ProblemType==2 || ProblemType==3
    freeTime = 1;
else
    freeTime = 0;
end
%% Setup SCOPT

Scopt     = scopt2;

Scopt.problem.nphases                      = 1;
Scopt.problem.phases(1).Nodes              = Nodes;
Scopt.problem.phases(1).n.states           = 6;
Scopt.problem.phases(1).n.controls         = 3;
Scopt.problem.phases(1).n.path             = 2;
Scopt.problem.phases(1).freeTimeFinal      = freeTime;

Scopt.setUpSCOPT_level1();
%%  Interface Main Input
index.states.position               = (1:3)' ;
index.states.velocity               = (4:6)' ;
index.controls.thrust               = (1:3)' ;

Scopt.problem.phases(1).initialGuess.states      = X_IG;
Scopt.problem.phases(1).initialGuess.controls    = U_IG;
Scopt.problem.phases(1).initialGuess.timeInitial = t0_IG;
Scopt.problem.phases(1).initialGuess.timeFinal   = tf_IG(ii);

% Internal Body Model
Scopt.problem.phases(1).BodyMap       = BodyMap;

% Dynamics
Scopt.problem.phases(1).dynamics.stateDerivativeFunction = @(t,X,U,BodyMap,ii) stateDerivativeSpacecraft(t,X,U,BodyMap,ii);
Scopt.problem.phases(1).dynamics.stateMatrixFunction     = @(t,X,U,BodyMap,ii) stateMatrixSpacecraft(t,X,U,BodyMap,ii);
Scopt.problem.phases(1).dynamics.controlMatrixFunction   = @(t,X,U,BodyMap,ii) controlMatrixSpacecraft(t,X,U,BodyMap,ii);

% Successive Convexification Speed Up
Scopt.problem.phases(1).interpolateParameters    = 0;

% Bounds
Scopt.problem.phases(1).bounds.timeFinal.upper = 50 ;
Scopt.problem.phases(1).bounds.timeFinal.lower = 20 ;

% Convex Path Constraints (usually resulting from lossless convexification)
Scopt.problem.phases(1).path(1).cone  = Scopt.initCone(4,Scopt.problem.phases(1).n);
Scopt.problem.phases(1).path(1).type  = 'upper';
Scopt.problem.phases(1).path(1).funType = 'convex';
Scopt.problem.phases(1).path(1).cone.norm.controls(1,1) = 1;
Scopt.problem.phases(1).path(1).cone.norm.controls(2,2) = 1;
Scopt.problem.phases(1).path(1).cone.norm.controls(3,3) = 1;
Scopt.problem.phases(1).path(1).cone.right.cons         = ThrustModel.Fmax;

Scopt.problem.phases(1).path(2).cone  = Scopt.initCone(3,Scopt.problem.phases(1).n);
Scopt.problem.phases(1).path(2).type  = 'upper';
Scopt.problem.phases(1).path(2).funType = 'convex';
Scopt.problem.phases(1).path(2).cone.norm.states(1,1) = 1;
Scopt.problem.phases(1).path(2).cone.norm.states(2,2) = 1;
Scopt.problem.phases(1).path(2).cone.right.states(3)  = 1/tanGlide;

% Initialing Conditions

Scopt.problem.phases(1).initial.states.equal.index = ...
    [index.states.position;...
    index.states.velocity];
Scopt.problem.phases(1).initial.states.equal.value = ...
    [P0;...
    V0];

% Final Conditions

Scopt.problem.phases(1).final.states.equal.index = [index.states.position;index.states.velocity];
Scopt.problem.phases(1).final.states.equal.value = zeros(6,1);

switch ProblemType
    case 0
        Scopt.problem.objective.type           = 'feasibility';
    case 1
        Scopt.problem.objective.type          = 'minimize';
        Scopt.problem.objective.scale                            = 200 ;
        Scopt.problem.objective.lagrange.funType           = 'convex';
        Scopt.problem.objective.lagrange.cone              = Scopt.initCone(4,Scopt.problem.phases(1).n);
        Scopt.problem.objective.lagrange.cone.norm.controls(1,1) = ThrustModel.gamma;
        Scopt.problem.objective.lagrange.cone.norm.controls(2,2) = ThrustModel.gamma;
        Scopt.problem.objective.lagrange.cone.norm.controls(3,3) = ThrustModel.gamma;
        
        Scopt.algorithm.sequential.activate      = 0;
    case 2
        Scopt.problem.objective.type          = 'minimize';
        Scopt.problem.objective.scale                   = tf_IG(ii);
        Scopt.problem.objective.mayer.funType           = 'linear'; % 1 for linear, 2 for convex, 3 for non-convex
        Scopt.problem.objective.mayer.where             = 'final';
        Scopt.problem.objective.mayer.timeFinal         = 1;
        Scopt.algorithm.sequential.activate      = 1;
    case 3
        Scopt.problem.objective.type          = 'minimize';
        Scopt.problem.objective.scale                            = 200 ;
        Scopt.problem.objective.lagrange.funType           = 'convex';
        Scopt.problem.objective.lagrange.cone              = Scopt.initCone(4,Scopt.problem.phases(1).n);
        Scopt.problem.objective.lagrange.cone.norm.controls(1,1) = ThrustModel.gamma;
        Scopt.problem.objective.lagrange.cone.norm.controls(2,2) = ThrustModel.gamma;
        Scopt.problem.objective.lagrange.cone.norm.controls(3,3) = ThrustModel.gamma;
        Scopt.algorithm.sequential.activate      = 1;
end

% Scopt.algorithm
Scopt.algorithm.conicSolver.feasTol      = eps^(1/2);
Scopt.algorithm.conicSolver.relTol       = eps^(1/2);
Scopt.algorithm.conicSolver.absTol       = eps^(1/2);
Scopt.algorithm.conicSolver.maxIter      = 100;

Scopt.algorithm.sequential.type          = 'none';
Scopt.algorithm.sequential.maxIter       = 20; % trust region or line-search
Scopt.algorithm.sequential.minIter       = 0;
% Scopt.algorithm.sequential.Jtol          = 0.0001;
% Scopt.algorithm.sequential.JtolChange    = 0.0001;

Scopt.algorithm.initialGuessType          = 'user';
Scopt.algorithm.scaling.variables         = 'automatic';
Scopt.algorithm.collocationMethod         = 'trapezoidal';
Scopt.algorithm.scaling.events            = 'automatic-limit'; % Alternative is automatic-limit automatic-jacobian
Scopt.algorithm.scaling.path              = 'automatic-jacobian'; % Alternative is automatic-limit automatic-jacobian
Scopt.algorithm.scaling.objective         = 'automatic-jacobian'; % Alternative is automatic-jacobian
Scopt.algorithm.meshRefinement            = 'none'; % No Mesh Refinement Implemented (YET)

Scopt.algorithm.sequential.trustRegion.phases(1).timeFinal.component = 5;
Scopt.algorithm.sequential.trustRegion.phases(1).timeFinal.lambda    = 1.0e0;
Scopt.algorithm.sequential.trustRegion.phases(1).timeFinal.type      = 'hard';

if ii==4 
    Scopt.algorithm.sequential.trustRegion.phases(1).timeFinal.type      = 'soft';
    Scopt.algorithm.sequential.trustRegion.phases(1).timeFinal.lambda    = 3e-1;
end

Scopt.setUpSCOPT_level2();
%% Solve Scopt.problem
Scopt.solveSCOPTProblem();
%% Retrieve Scopt.solutions

Solution = Scopt.solution;

X_IG                 = Scopt.problem.phases.initialGuess.states  ;
X(:,:,ii)            = Solution.states;
U(:,:,ii)            = Solution.controls;
Jreal(ii)            = Solution.Jreal;
timeFinal(ii)        = Solution.timeFinal;
timeElapsed(ii)      = Solution.timeElapsed;
end


%% External Post Processing


% figure;defaultColor = get(gca,'colororder');close(gcf)
defaultColor = [1 0 0;0 1 0;0 0 1;0.5 1 1];
linestyle    = {'-','-','-','-.'};

p_IG  = X_IG(1:3,:); 
p = X(1:3,:,:); 
f = U(1:3,:,:); 

%% Plot 3D
% Plot Glide Cone
x = linspace(-40,55,30); y = linspace(0,55,30);
[X,Y] = meshgrid(x,y);
Z = tanGlide*sqrt(X.^2+Y.^2);
figure(1); colormap autumn ; s=surf(X,Y,Z);s.FaceAlpha = 0.1;
axis([-40,55,0,55,0,105]);
grid on; hold on;

% INSERT YOUR VARIABLES HERE:
% -------------------------------------------------------

hplotIG=plot3(p_IG(1,:),p_IG(2,:),p_IG(3,:),'k--','linewidth',1);
for kk = 1:NumberProblems
    hplot(kk)=plot3(p(1,:,kk),p(2,:,kk),p(3,:,kk),'linestyle',linestyle{kk},'Color',defaultColor(kk,:),'linewidth',1.5);
    quiver3(p(1,:,kk),p(2,:,kk),p(3,:,kk),...
        f(1,:,kk),f(2,:,kk),f(3,:,kk),0.3,'k','linewidth',1.5);
end

xlabel('$P_x$ [-]','interpreter','latex')
ylabel('$P_y$ [-]','interpreter','latex')
zlabel('$P_z$ [-]','interpreter','latex')
legend([hplotIG,hplot(1),hplot(2),hplot(3),hplot(4)],{'Initial Guess','Problem 1','Problem 2 (hard)','Problem 3 (hard)','Problem 3 (soft)'},'interpreter','latex','location','east')


% upFig(gcf,['verification_spacecraftLandingBoyd_trajectory'],'Guidance/Convex/Verification','.pdf','wide')


