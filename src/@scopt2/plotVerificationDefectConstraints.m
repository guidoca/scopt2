function [hfig1,hfig2,hfig3] =  plotVerificationDefectConstraints(obj,fig1,fig2,fig3)
% verifies the implementartion of the dynamic equations,
% figurwes should show values with small errors
states_IG                  = obj.problem.phases.initialGuess.states;
controls_IG                = obj.problem.phases.initialGuess.controls;
time_IG                    = obj.problem.phases.initialGuess.time;
timeFinal_IG               = obj.problem.phases.initialGuess.timeFinal;
timeInitial_IG             = obj.problem.phases.initialGuess.timeInitial;
states_ALL                 = cat(3,states_IG,obj.debugging.states_ALL);
controls_ALL               = cat(3,controls_IG,obj.debugging.controls_ALL);
time_ALL                   = cat(3,time_IG,obj.debugging.time_ALL);
timeFinal_ALL              = cat(2,timeFinal_IG,obj.debugging.timeFinal_ALL);
timeInitial_ALL            = cat(2,timeInitial_IG,obj.debugging.timeInitial_ALL);


iterations = obj.solution.iterations;

%             lineOpacity_vector = linspace(1,0,iterations);


Nodes           = size(states_ALL,2);
nDefects        = size(states_ALL,1);
defectsjColor   = jet(nDefects);
iterationsColor = autumn(iterations);
lineOpacity  = linspace(0,1,iterations+1);lineOpacity(1) = [];

statesDiff  = diff(states_ALL,1,2);
dtau        = 2/(Nodes-1);
dt          = dtau*(timeFinal_ALL-timeInitial_ALL)/2;
if obj.algorithm.virtualControl.phases(1).include
    virtualControls_ALL        = obj.debugging.virtualControls;
    valEvirtualControls        = obj.algorithm.virtualControl.phases(1).valE;
    rowEvirtualControls        = obj.algorithm.virtualControl.phases(1).rowE;
    colEvirtualControls        = obj.algorithm.virtualControl.phases(1).colE;
    scaleE                     = obj.algorithm.virtualControl.phases(1).scale;
    Evirtual                   = full(sparse(rowEvirtualControls,colEvirtualControls,valEvirtualControls,nDefects,max(colEvirtualControls)));
else
    virtualControls_ALL        = zeros(0,Nodes-1,iterations);
    scaleE                     = 0;
    Evirtual                   = zeros(nDefects,0);
end
stateDerivativeFunction = obj.problem.phases.dynamics.stateDerivativeFunction;
stateMatrixFunction     = obj.problem.phases.dynamics.stateMatrixFunction;
controlMatrixFunction   = obj.problem.phases.dynamics.controlMatrixFunction;
timeDerivativeFunction  = obj.problem.phases.dynamics.timeDerivativeFunction;

obj.problem.phases(1).BodyMap.evaluated = 0;

defectConstraints = zeros(nDefects,Nodes-1,iterations);
for kk = 1:iterations
    if kk>1
        if obj.debugging.rejected(kk-1)>0
            kklin = find(obj.debugging.rejected(1:kk-1)==0,1,'last')+1;
        else
            kklin = kk;
        end
    else
        kklin = kk;
    end
    f    = stateDerivativeFunction(time_ALL(:,:,kklin),states_ALL(:,:,kklin),controls_ALL(:,:,kklin),obj.problem.phases.BodyMap,1:Nodes);
    dfdt = timeDerivativeFunction(time_ALL(:,:,kklin),states_ALL(:,:,kklin),controls_ALL(:,:,kklin),obj.problem.phases.BodyMap,1:Nodes);
    switch obj.algorithm.collocationMethod
        case {'euler','trapezoidal'}
            if obj.problem.phases(1).dynamics.typeStateMatrix
                [rowA,colA,valA] = stateMatrixFunction(time_ALL(:,:,kklin),states_ALL(:,:,kklin),controls_ALL(:,:,kklin),obj.problem.phases.BodyMap,1:Nodes);
            else
                A = stateMatrixFunction(time_ALL(:,:,kklin),states_ALL(:,:,kklin),controls_ALL(:,:,kklin),obj.problem.phases.BodyMap,1:Nodes);
            end
            if obj.problem.phases(1).dynamics.typeControlMatrix
                [rowB,colB,valB] = controlMatrixFunction(time_ALL(:,:,kklin),states_ALL(:,:,kklin),controls_ALL(:,:,kklin),obj.problem.phases.BodyMap,1:Nodes);
            else
                B = controlMatrixFunction(time_ALL(:,:,kklin),states_ALL(:,:,kklin),controls_ALL(:,:,kklin),obj.problem.phases.BodyMap,1:Nodes);
            end
        case 'exact'
            if obj.problem.phases(1).dynamics.typeStateMatrix
                [rowA,colA,valA,~,rowA1,colA1,valA1] = stateMatrixFunction(time_ALL(:,:,kklin),states_ALL(:,:,kklin),controls_ALL(:,:,kklin),obj.problem.phases.BodyMap,1:Nodes);
            else
                [A,A1] = stateMatrixFunction(time_ALL(:,:,kklin),states_ALL(:,:,kklin),controls_ALL(:,:,kklin),obj.problem.phases.BodyMap,1:Nodes);
            end
            if obj.problem.phases(1).dynamics.typeControlMatrix
                [rowB,colB,valB,~,rowB1,colB1,valB1] = controlMatrixFunction(time_ALL(:,:,kklin),states_ALL(:,:,kklin),controls_ALL(:,:,kklin),obj.problem.phases.BodyMap,1:Nodes);
            else
                [B,B1] = controlMatrixFunction(time_ALL(:,:,kklin),states_ALL(:,:,kklin),controls_ALL(:,:,kklin),obj.problem.phases.BodyMap,1:Nodes);
            end
    end
    for ii = 1:Nodes-1 % Depends on discretisation option used
        switch obj.algorithm.collocationMethod
            case 'euler'
                fii = f(:,ii);
                if obj.problem.phases(1).dynamics.typeStateMatrix
                    Aii = full(sparse(rowA(:,ii),colA(:,ii),valA(:,ii),nDefects,nDefects));
                else
                    Aii = A(:,:,ii);
                end
                if obj.problem.phases(1).dynamics.typeControlMatrix
                    Bii = full(sparse(rowB(:,ii),colB(:,ii),valB(:,ii),nDefects,size(controls_IG,1)));
                else
                    Bii = B(:,:,ii);
                end
                dfdtii = dfdt(:,ii);
                
                defectConstraints(:,ii,kk) = obj.problem.phases.scale.states.^(-1).*statesDiff(:,ii,kk+1) - obj.problem.phases.scale.states.^(-1).*( ...
                    dt(kklin)*(...
                    fii...
                    + Aii*( states_ALL(:,ii,kk+1)   - states_ALL(:,ii,kklin)   )...
                    + Bii*( controls_ALL(:,ii,kk+1) - controls_ALL(:,ii,kklin) ) )...
                    + (dtau/2*fii+dt(kklin)*dfdtii)*(timeFinal_ALL(kk+1)-timeFinal_ALL(kklin)))...
                    - scaleE*Evirtual*virtualControls_ALL(:,ii,kk);
            case 'trapezoidal'
                fii  = f(:,ii);
                fii1 = f(:,ii+1);
                if obj.problem.phases(1).dynamics.typeStateMatrix
                    Aii  = full(sparse(rowA(:,ii),colA(:,ii),valA(:,ii),nDefects,nDefects));
                    Aii1 = full(sparse(rowA(:,ii+1),colA(:,ii+1),valA(:,ii+1),nDefects,nDefects));
                else
                    Aii  = A(:,:,ii);
                    Aii1 = A(:,:,ii+1);
                end
                if obj.problem.phases(1).dynamics.typeControlMatrix
                    Bii  = full(sparse(rowB(:,ii),colB(:,ii),valB(:,ii),nDefects,size(controls_IG,1)));
                    Bii1 = full(sparse(rowB(:,ii+1),colB(:,ii+1),valB(:,ii+1),nDefects,size(controls_IG,1)));
                else
                    Bii  = B(:,:,ii);
                    Bii1 = B(:,:,ii+1);
                end
                dfdtii  = dfdt(:,ii);
                dfdtii1 = dfdt(:,ii+1);
                
                defectConstraints(:,ii,kk) = obj.problem.phases.scale.states.^(-1).*statesDiff(:,ii,kk+1) - obj.problem.phases.scale.states.^(-1).*(...
                    dt(kklin)*(...
                    (fii+fii1)/2 ...
                    + ( Aii*( states_ALL(:,ii,kk+1)   - states_ALL(:,ii,kklin) )   + Aii1*( states_ALL(:,ii+1,kk+1)   - states_ALL(:,ii+1,kklin) ) )/2 ...
                    + ( Bii*( controls_ALL(:,ii,kk+1) - controls_ALL(:,ii,kklin) ) + Bii1*( controls_ALL(:,ii+1,kk+1) - controls_ALL(:,ii+1,kklin) ) )/2 ) ...
                    + ( dtau/2*(fii+fii1)/2+dt(kklin)*(dfdtii+dfdtii1)/2)*(timeFinal_ALL(kk+1)-timeFinal_ALL(kklin)))...
                    - scaleE*Evirtual*virtualControls_ALL(:,ii,kk);
            case 'exact'
                fii     = f(:,ii);
                if obj.problem.phases(1).dynamics.typeStateMatrix
                    Aii  = full(sparse(rowA(:,ii),colA(:,ii),valA(:,ii),nDefects,nDefects));
                    Aii1 = full(sparse(rowA1(:,ii),colA1(:,ii),valA1(:,ii),nDefects,nDefects));
                else
                    Aii  = A(:,:,ii);
                    Aii1 = A1(:,:,ii);
                end
                if obj.problem.phases(1).dynamics.typeControlMatrix
                    Bii     = full(sparse(rowB(:,ii),colB(:,ii),valB(:,ii),nDefects,size(controls_IG,1)));
                    Bii1    = full(sparse(rowB1(:,ii),colB1(:,ii),valB1(:,ii),nDefects,size(controls_IG,1)));
                else
                    Bii  = B(:,:,ii);
                    Bii1 = B1(:,:,ii);
                end
                dfdtii  = dfdt(:,ii);
                
                defectConstraints(:,ii,kk) = obj.problem.phases.scale.states.^(-1).*statesDiff(:,ii,kk+1) - obj.problem.phases.scale.states.^(-1).*( fii ...
                    + Aii*( states_ALL(:,ii,kk+1)   - states_ALL(:,ii,kklin) )    + Aii1*( states_ALL(:,ii+1,kk+1)   - states_ALL(:,ii+1,kklin) ) ...
                    + Bii*( controls_ALL(:,ii,kk+1) - controls_ALL(:,ii,kklin) )  + Bii1*( controls_ALL(:,ii+1,kk+1) - controls_ALL(:,ii+1,kklin) ) ...
                    + dfdtii*dtau/2*( timeFinal_ALL(kk+1) - timeFinal_ALL(kklin)  + timeInitial_ALL(kk+1) + timeInitial_ALL(kklin) ) )...
                    - scaleE*Evirtual*virtualControls_ALL(:,ii,kk);
            otherwise
                error('SCOPT Error: Collocation method does not exist')
        end
    end
end
%             defectConstraints = defectConstraints/scaleE;

switch obj.algorithm.sequential.trustRegion.defectsCompErrorPenalty
    case 0
        defectConstraintsP1 = sum(defectConstraints.*defectConstraints,1);
    case {1,2,Inf}
        defectConstraintsP1 = vecnorm(defectConstraints,obj.algorithm.sequential.trustRegion.defectsCompErrorPenalty,1);
end

switch obj.algorithm.sequential.trustRegion.defectsNodeErrorPenalty
    case {1,2,Inf}
        defectConstraintsP2 = vecnorm(defectConstraintsP1,obj.algorithm.sequential.trustRegion.defectsNodeErrorPenalty,2);
end

% Figure one is error at each state component at each node at
% each iteration
%   - (colours contour for each component)
%   - (opacity for iteration)
%   - (x axis is each node)

% Figure two is norm2 error at each node at each iteration
%   - (colours contour for each iteration)
%   - (x axis is each node)
% Figure three is norm2 error for concatenate full defect
% constraints vector at each iteration
%   - (x axis is each iteration)


if nargin >= 2
    if fig1~=0
        figure(fig1)
        dontPlot = 0;
    else
        dontPlot = 1;
    end
else
    figure
    dontPlot = 0;
end
if ~dontPlot
    labelLegend = {};
    for jj = 1:nDefects
        for kk = 1:iterations
            if kk ==iterations
                markerEdgeColor = [0 0 0];
            else
                markerEdgeColor = defectsjColor(jj,:);
            end
            hplotAux   = semilogy(1:Nodes-1,abs(defectConstraints(jj,:,kk)),'o-','linewidth',1,'color',[1,1,1]*(1-lineOpacity(kk)) + lineOpacity(kk)*defectsjColor(jj,:),'MarkerFaceColor',[1,1,1]*(1-lineOpacity(kk)) + lineOpacity(kk)*defectsjColor(jj,:),'MarkerEdgeColor',markerEdgeColor); hold on
        end
        hplot(jj)   = hplotAux;
        labelLegend = [labelLegend , ['Defect $j$ = ' num2str(jj)]];
    end
    ax = gca;
    ax.YGrid      = 'on';
    ax.YMinorGrid = 'on';
    xlabel('Defect $i$ [-]','interpreter','latex')
    ylabel('Error [-]','interpreter','latex')
    legend(hplot,labelLegend,'location','best','interpreter','latex')
    scoptSave2pdf(obj,'verificationDefectsj')
    hfig1 = gcf;
else
    hfig1 = nan;
end
if nargin >= 3
    if fig2~=0
        figure(fig2)
        dontPlot = 0;
    else
        dontPlot = 1;
    end
else
    figure
    dontPlot = 0;
end
if ~dontPlot
    labelLegend = {};
    for kk = 1:iterations
        hplot(kk)   = semilogy(1:Nodes-1,defectConstraintsP1(1,:,kk),'-*','linewidth',2,'color',iterationsColor(kk,:)); hold on
        labelLegend = [labelLegend , ['$k$ = ' num2str(kk)]];
    end
    ax = gca;
    ax.YGrid = 'on';
    ax.YMinorGrid = 'on';
    xlabel('Defect $i$ [-]','interpreter','latex')
    ylabel('Error [-]','interpreter','latex')
    legend(hplot,labelLegend,'location','best','interpreter','latex')
    scoptSave2pdf(obj,'verificationDefectsPC')
    hfig2 = gcf;
else
    hfig2 = nan;
end
if nargin >= 4
    if fig3~=0
        figure(fig3)
        dontPlot = 0;
    else
        dontPlot = 1;
    end
else
    figure
    dontPlot = 0;
end
if ~dontPlot
    labelLegend = {};
    hplotP2(1) = semilogy(1:iterations,squeeze(defectConstraintsP2),'-*','linewidth',2);  hold on
    ax = gca;
    ax.YGrid = 'on';
    ax.YMinorGrid = 'on';
    xlabel('Iteration $k$ [-]','interpreter','latex')
    ylabel('Error [-]','interpreter','latex')
    scoptSave2pdf(obj,'verificationDefectsPN')
    hfig3 = gcf;
else
    hfig3 = nan;
end
end


