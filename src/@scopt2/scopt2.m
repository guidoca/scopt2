%% *SCOPT2* 
% Successive Convexification OPtimal ConTrol 2.0
% Optimal control tool to solve single phase autonomous problems
% using successive convexification.
% Copyright (C) 2020  Guillermo J. Dominguez Calabuig

% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.

%% Description
% Handle class object defining and solving an optimal control problem
% based on successive convexification and linearisation. Can also solve
% basic SOCP/Linear problems with one iteration.
% It containts three main structures:
%
% * algorithm: Defines the settings necessary for the sequential socp
%       solver, as the socp solver settings, virtual controls, buffer
%       zones, trust regions, scaling, line search strategies, and
%       collocation method
% * problem: Defines the Optimum Control Problem per phase. Currently
%   it only supports single phase problems.
% * solution: Stores the solution of the Optimum Control Problem
%
% A full Description is available in the Git Repository:
% https://github.com/guidoca/scopt2.git

% and in a Master Thesis from the TU Delft Repository:
%  'Optimum On-board Abort Guidance based on Successive Convexification
%  for Atmospheric Re-Entry Vehicles', Master Thesis, Guillermo
%  Joaquin Dominguez Calabuig. 
% <https://repository.tudelft.nl/islandora/object/uuid%3A3a33c54b-ed54-407a-88ce-c9d0b3d8f738?collection=education>
 
classdef scopt2 < matlab.mixin.Copyable
    
    properties (Constant)
        infCOPT = 1e19;
        MIOFF   = 1   ; % Matlab Index Offset
    end
    properties
        algorithm = struct;
        problem   = struct;
        solution  = struct;
        debugging = struct;
        quiet     = 0;
    end
    
    methods
        %% SCOPT Constructor
        function obj = scopt2()
            %scopt Construct an instance of the scopt class loading the
            %default fields for the problem and algorithm
            addpath(genpath('toolbox'));
            warning('on','all');
            %% problem Default Settings 
            obj.problem.nphases                                 = 1 ; % number of phases
            for ph = 1:obj.problem.nphases
                obj.problem.phases(ph).n.states              = []; % number of states, requires initialisation
                obj.problem.phases(ph).n.controls            = []; % number of controls, requires initialisation
                obj.problem.phases(ph).n.parameters          = 0 ; % number of parameters, 0 by default
                obj.problem.phases(ph).n.events              = 0 ; % number of events, 0 by default
                obj.problem.phases(ph).n.path                = 0 ; % number of path constraints, 0 by default
                obj.problem.phases(ph).freeTimeFinal         = 0 ; % indicated if free final   time phase, false by default
                obj.problem.phases(ph).freeTimeInitial       = 0 ; % indicated if free initial time phase, false by default
            end
            %% algorithm Default Settings
            % SOCP Solver settings. Defines the maximum number of SOCP
            % subiterations, solver type, and feasibility, relative and absolute
            % tolerances
            obj.algorithm.conicSolver.type          = 'ecos';
            obj.algorithm.conicSolver.maxIter       = 100;
            obj.algorithm.conicSolver.feasTol       = eps^(1/2);
            obj.algorithm.conicSolver.relTol        = eps^(1/2);
            obj.algorithm.conicSolver.absTol        = eps^(1/2);
            
            % Sequential algorithm settings
            obj.algorithm.sequential.activate       = 1; % successive convexification is active by default
            obj.algorithm.sequential.type           = 'trust-region'; % trust region
            obj.algorithm.sequential.maxIter        = 100; % Maximum number of successive convexification iterations
            obj.algorithm.sequential.minIter        = 3; % Minimum number of successive convexification iterations
            obj.algorithm.sequential.objTol         = 0.05;  % Maximum percentage of augmentation to consider a solution feasible, heuristic term.
            obj.algorithm.sequential.Jtol           = 1e-4; % Stopping condition tolerance for the objective function (if difference is lower than this number, stops)
            obj.algorithm.sequential.JtolChange     = 1e-4; % Stopping condition tolerance for the objective function (if relative change is lower than this number, stops)
            obj.algorithm.sequential.variablesTol.etaTol = 1.0e-4; % Stopping condition tolerance the change in variables 
            obj.algorithm.sequential.variablesTol.type   = 'norm'; % Type of norm for variables and states, norm or componentwise, default is norm
            obj.algorithm.sequential.variablesTol.norm   = 2; % Type of norm used
            obj.algorithm.sequential.variablesTol.include.states        = 1; % flag to include states
            obj.algorithm.sequential.variablesTol.include.controls      = 1; % flag to include controls
            obj.algorithm.sequential.variablesTol.include.timeFinal     = 1; % flag to include final time
            obj.algorithm.sequential.variablesTol.include.timeInitial   = 1; % flag to include initial time
            
            % Sequential algorithm line search settings
            obj.algorithm.sequential.globalisation.type                       = 'line-search'; % line-search
            obj.algorithm.sequential.globalisation.activate                   = 0; % inactive by default
            obj.algorithm.sequential.globalisation.lineSearch.type            = 'golden-section'; % line-search type, options are contraction, polynomial, and golden section
            obj.algorithm.sequential.globalisation.lineSearch.contractionFactor       = 0.5; % contraction factor for the contraction line search type
            obj.algorithm.sequential.globalisation.lineSearch.directionalDerivativeType = 'numeric';% Numerical directional derivative
            obj.algorithm.sequential.globalisation.lineSearch.directionalDerivativeEpsilon  = 1e-4; % Used for numerical directional derivative steps, only for contraint type
            obj.algorithm.sequential.globalisation.lineSearch.k1              = 1e-4; % Used for contraction line search type
            obj.algorithm.sequential.globalisation.lineSearch.k2              = 0.9 ; % Used for contraction line search type
            obj.algorithm.sequential.globalisation.lineSearch.bracketTol      = 1e-2; % Used for bracketing methods as golden section
            obj.algorithm.sequential.globalisation.lineSearch.golden          = (sqrt(5)-1)/2; % Golden ratio
            obj.algorithm.sequential.globalisation.lineSearch.stepLengthMax   = 1.0; % Maximum step length along direction and magnitude defined by SOCP solver
            obj.algorithm.sequential.globalisation.lineSearch.forceIfRejected = 1 ; % Flag specifiying if the line search should still be called if the adaptive trust reegion algorithmn has rejected the line search step
            obj.algorithm.sequential.globalisation.lineSearch.debugging       = 0 ; % Flag to call debugging routines in line search
            
            % Sequential algorithm trust region settings
            obj.algorithm.sequential.trustRegion.defectsCompErrorPenalty  = nan ; % 2nd norm, 0 means quadratic norm This is for the computation of the adaptive trust region error
            obj.algorithm.sequential.trustRegion.defectsNodeErrorPenalty  = nan ; % 2nd norm This is for the computation of the adaptive trust region error
            obj.algorithm.sequential.trustRegion.radius_l = eps^(1/2)   ; % lower boundary for adaptive hard trust region radius
            obj.algorithm.sequential.trustRegion.radius_u = obj.infCOPT ; % upper boundary for adaptive hard trust region radius
            obj.algorithm.sequential.trustRegion.lambda_l = eps^(1/2)   ; % lower boundary for adaptive soft trust region multiplier
            obj.algorithm.sequential.trustRegion.lambda_u = obj.infCOPT ; % upper boundary for adaptive soft trust region multiplier
            obj.algorithm.sequential.trustRegion.rho0     = 0.0 ; % Adaptive trust region metric. Lower means it is a bad appropximation and should be rejected
            obj.algorithm.sequential.trustRegion.rho1     = 0.2 ; % Adaptive trust region metric. Lower means it is a bad appropximation and hard(soft) trust region should be reduced(increase)
            obj.algorithm.sequential.trustRegion.rho2     = 0.6 ; % Adaptive trust region metric. Lower Means approximation is OK. Higher means it is a good approximation and hard(soft) trust region should be increased(reduced)
            obj.algorithm.sequential.trustRegion.alpha    = 2.0 ; % Adaptive trust region increasing factor for good approximations
            obj.algorithm.sequential.trustRegion.beta     = 2.0 ; % Adaptive trust region decreasing factor for bad approximations
            obj.algorithm.sequential.trustRegion.type                 = 'hard'; % 'hard' or 'soft' trust region type
            obj.algorithm.sequential.trustRegion.adaptive             = 0     ; % 0 or 1. Default is a rigid trust region
            obj.algorithm.sequential.trustRegion.include.states       = 0     ; % 0 or 1. Default does not apply trust region to controls
            obj.algorithm.sequential.trustRegion.include.controls     = 1     ; % 0 or 1. Default does applies trust region to controls
            obj.algorithm.sequential.trustRegion.include.pathAndEventsErrors     = 1     ; % 0 or 1 if path and event convexification errors should be included in metric
            obj.algorithm.sequential.trustRegion.variablesPenalty     = 2     ; % 2nd norm, 0 means quadratic norm for trust region component terms
            obj.algorithm.sequential.trustRegion.nodePenalty          = inf   ; % 2nd norm, norm for trust region node wise terms (on components)
            obj.algorithm.sequential.trustRegion.radius               = []    ; % Shall be initialised
            obj.algorithm.sequential.trustRegion.lambda               = []    ; % Shall be initialised
            
            % Other Algorithm Settings
            obj.algorithm.initialGuessType          = 'automatic'; % Initial guess type. 'automatic', 'user', or 'none'
            obj.algorithm.collocationMethod         = 'trapezoidal'; % Collocation method. Numerical types: 'trapezoidal' or 'euler'. Or exact ode solution: 'exact'.
            obj.algorithm.scaling.variables         = 'automatic'; % Scaling for states and controls variable options.  'automatic', 'user', or 'none'
            obj.algorithm.scaling.events            = 'automatic-jacobian'; % Alternative is 'none' 'user' 'automatic-limit' 'automatic-jacobian'
            obj.algorithm.scaling.path              = 'automatic-jacobian'; % Alternative is 'none' 'user' 'automatic-limit' 'automatic-jacobian'
            obj.algorithm.scaling.objective         = 'automatic-jacobian'; % Alternative is 'automatic-jacobian'
            obj.algorithm.scaling.trustRegion       = 'none'; % Alternative is 'automatic' or 'none'
            obj.algorithm.meshRefinement            = 'none'; % No Mesh Refinement Implemented (YET)
            obj.algorithm.quadrature.type           = 'trapezoidal'; % Quadrature types for integral terms. Alternative is 'euler'
             
            
            obj.algorithm.speedUp      = 0; % To avoid double calculations in objective, dynamics, states and control variables. Default is inactive
            
            %% Default Post Processing Options in Solution Struct
            obj.solution.outputFolder   = [pwd '\scoptOutput']; % outpud folder directory
            obj.solution.figuresFolder  = 'Figures'; % Postprocessing Figures Folder Directioy
            obj.solution.dataFolder     = 'Data'; % Postprocessing Data Folder Directory
            obj.solution.saveFigures    = 1; % Flag to save figures (enabled by default)
            obj.solution.generateFigures= 0; % Flag to generate figures (disabled by default)
            obj.solution.generateOutput = 1; % Flag to generate output (enabled by default)
            obj.solution.figuresQuality = 600; % quality of figures in pdf format
            %% Debugging options 
            obj.debugging.saveVirtualControls    = 1; % flag to save virtual control verification monitors
            obj.debugging.saveTrustRegion        = 1; % flag to save trust region verification monitors
            obj.debugging.saveVirtualBuffers     = 1; % flag to save virtual buffers verification monitors
            obj.debugging.saveVariableConvergence= 1; % flag to save variable convergence monitor
            obj.debugging.generateFigureDefects  = 1; % flag to generate figures for the defect constraints (dynamics)
            obj.debugging.internalVerifications  = 1; % flag to do all of the above
        end
        %% Setup Methods
        obj = setUpSCOPT_level1(obj)
        obj = setUpSCOPT_level2(obj)
        
        %% Solve Methods
        solveSCOPTProblem(obj); 
        problemTranscription(obj,timeInitial_n,timeFinal_n,X_n,U_n,radius,lambda)
        solverOutput = solveSOCPProblem(obj)
        
        %% Globalisation Methods
        Z = computeNewPoint(~,searchDirection,stepLength,Z_p)
        M = computeMeritFunction(obj,searchDirection,stepLength,Z_p)
        directionalDerivative = computeDirectionalDerivative(obj,searchDirection,Z_p,M_p)
        sufficientDecrease = evaluateSufficientDecrease(obj,searchDirection,stepLength,directionalDerivative,Z_p,M_p)
        [radius,lambda,rejected] = adaptTrustRegion(obj,radius,lambda,rho,numericalIssues)
        
        %% Objective Function Related Methods
        [J,L,Jreal,Lreal] = computeObjectiveFunctionsNormal(obj,time_n,X_n,U_n,time_n_p,X_n_p,U_n_p)
        [Jreal,Lreal] = computeObjectiveFunctions(obj,time,X,U,time_p,X_p,U_p)
        [Jaug,Laug] = computeAugmentedObjectiveFunctions(obj,time_n,X_n,U_n,time_n_p,X_n_p,U_n_p,linearToo)
        [Reg] = computeRegularisation(obj,~,X_n,U_n)
        [penaltyDefectsNonLinear,penaltyBuffersNonLinear,penaltyTrustRegion,penaltyTime,penaltyDefectsLinear,penaltyBuffersLinear] = computeAugmentedObjectiveFunctionTerms(obj,time_n,X_n,U_n,time_n_p,X_n_p,U_n_p,linearToo)
        %% Termination Auxiliary Methods
        varTol = terminationVariableConvergence(obj,timeInitial_n_opt,timeFinal_n_opt,X_n_opt,U_n_opt,timeInitial_n_prev,timeFinal_n_prev,X_n_prev,U_n_prev)
        %% Variables Scaling Methods 
        Xs_n = scaleStatesToNormalized(obj,ph,Xs,index)
        Xs   = scaleStatesToReal(obj,ph,Xs_n,index)
        Us_n = scaleControlsToNormalized(obj,ph,Us,index)
        Us   = scaleControlsToReal(obj,ph,Us_n,index)
        t_n  = scaleTimeToNormalized(obj,ph,t)
        t    = scaleTimeToReal(obj,ph,t_n) 
        %% Postprocessing Methods 
        generateFigures(obj) 
        [hfig1,hfig2,hfig3] =  plotVerificationDefectConstraints(obj,fig1,fig2,fig3)
        [trustRegionErrorP1_ALL,trustRegionErrorP2_ALL,hfig1,hfig2] =  plotTrustRegionErrorHistory(obj,fig1,fig2)
        [trustRegionP1_ALL,trustRegionP2_ALL,hfig1,hfig2] =  plotTrustRegionHistory(obj,fig1,fig2)
        [virtualBuffersPathErrorPs,virtualBuffersPathErrorP,virtualBuffersEventsErrorP,hfig1,hfig2] =  plotVirtualBuffersErrorHistory(obj,fig1,fig2)
        [virtualBuffersPathP1_ALL,virtualBuffersPathP2_ALL,virtualBuffersEventsP1_ALL,hfig1,hfig2] =  plotVirtualBuffersHistory(obj,fig1,fig2)
        [virtualControlsErrorP1_ALL,virtualControlsErrorP2_ALL,hfig1,hfig2] =  plotVirtualControlErrorHistory(obj,fig1,fig2)
        [virtualControlsP1,virtualControlsP2] =  computePenaltiesVirtualControls(obj,virtualControls)
        [virtualControlsP1_ALL,virtualControlsP2_ALL,hfig1,hfig2] =  plotVirtualControlHistory(obj,fig1,fig2)
        variablesNorm = computeConvergenceHistory(obj)
        hfig =  plotConvergenceHistory(obj,fig,variablesNorm)
        hfig =  plotRealObjective(obj,fig)
        hfig =  plotTime(obj,fig)
        hfig =  plotAdaptiveTrustRegionMonitor(obj,fig)
        hfig =  plotSolutionMatrices(obj) 
        scoptSave2pdf(obj,pdfFileName)
        
        
    end
    %% Static Methods
    methods (Static)
        %% Create Cone
        function [s] = initCone(dimensions,n,Nodes)
            % Initializes a cone with given deimensions, for the states,
            % controls and parameters and for multidimensional number of
            % nodes if necessary
            if ~exist('Nodes','var')
                Nodes = 1;
            end
            s.nodes          = Nodes;
            s.dimensions     = dimensions;
            s.norm.cons      = zeros(dimensions-1,Nodes);
            s.right.cons     = zeros(1,Nodes);
            
            s.norm.states      = zeros(n.states,dimensions-1,Nodes);
            s.norm.controls    = zeros(n.controls,dimensions-1,Nodes);
            s.norm.parameters  = zeros(n.parameters,dimensions-1,Nodes);
            s.right.states     = zeros(n.states,Nodes);
            s.right.controls   = zeros(n.controls,Nodes);
            s.right.parameters = zeros(n.parameters,Nodes);
        end
        %% Postprocessing Static Methods
        plotNonZeroEntries(MAT,figNum,index)
        plotNonZeroEntriesAB(A,b,figNum,index)
    end
end

