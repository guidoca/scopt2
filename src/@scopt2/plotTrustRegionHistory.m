function [trustRegionP1_ALL,trustRegionP2_ALL,hfig1,hfig2] =  plotTrustRegionHistory(obj,fig1,fig2)

states_n_IG   = obj.problem.phases.initialGuess.states_n;
controls_n_IG   = obj.problem.phases.initialGuess.controls_n;
states_n_ALL   = cat(3,states_n_IG,obj.debugging.states_n_ALL);
controls_n_ALL = cat(3,controls_n_IG,obj.debugging.controls_n_ALL);
radius          = obj.debugging.radius;

iterations = obj.solution.iterations;

%             lineOpacity_vector = linspace(1,0,iterations);
defaultColor = autumn(iterations);

Nodes      = size(states_n_ALL,2);
variables_ALL   = zeros(0,Nodes,iterations+1);


if obj.algorithm.sequential.trustRegion.include.states
    variables_ALL = cat(1,variables_ALL,states_n_ALL);
end

if obj.algorithm.sequential.trustRegion.include.controls
    variables_ALL = cat(1,variables_ALL,controls_n_ALL);
end

trustRegionAux  = diff(variables_ALL,1,3);

%             IndexRejected = find(obj.debugging.rejected>0);
for kk = 2:iterations+0
    if obj.debugging.rejected(kk-1)
        kkValid = find(obj.debugging.rejected(1:kk-1)==0,1,'last');
        trustRegionAux(:,:,kk) = variables_ALL(:,:,kk+1)-variables_ALL(:,:,kkValid+1);
    end
    %                 trustRegionAux(:,:,IndexRejected+1) = variables_ALL(:,:,IndexRejected+2)-variables_ALL(:,:,IndexRejected-0);
end

switch obj.algorithm.sequential.trustRegion.variablesPenalty
    case 0
        trustRegionP1_ALL = sum(trustRegionAux.*trustRegionAux,1);
    case {1,2,inf}
        trustRegionP1_ALL = vecnorm(trustRegionAux,obj.algorithm.sequential.trustRegion.variablesPenalty,1);
    otherwise
        error('SCOPT ERROR \n Trust Region with variable penalty type %i not implemented \n',obj.algorithm.sequential.trustRegion.variablesPenalty )
end

switch obj.algorithm.sequential.trustRegion.nodePenalty
    case 0
        error('SCOPT ERROR \n Trust Region with node penalty type %i not implemented \n',obj.algorithm.sequential.trustRegion.nodePenalty )
    case {1,2,inf}
        trustRegionP2_ALL = reshape(vecnorm(trustRegionP1_ALL,obj.algorithm.sequential.trustRegion.nodePenalty,2),1,[]);
end

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
    for kk = 1:iterations
        hplot(kk)   = semilogy(1:Nodes,trustRegionP1_ALL(1,:,kk),'-*','linewidth',2,'color',defaultColor(kk,:)); hold on
        labelLegend = [labelLegend , ['$k$ = ' num2str(kk)]];
    end
    ax = gca;
    ax.YGrid = 'on';
    ax.YMinorGrid = 'on';
    xlabel('Node $i$ [-]','interpreter','latex')
    ylabel('$P^c\left(y^k_i-y^{k-1}_i\right)$ [-]','interpreter','latex')
    legend(hplot,labelLegend,'location','best','interpreter','latex')
    scoptSave2pdf(obj,'trustRegionPenaltyC')
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
    hplotP2(1) = semilogy(1:iterations,squeeze(trustRegionP2_ALL),'-*','linewidth',2);  hold on
    legendLabelP2{1}     = '$P^n\left(P^c\left(y^k-y^{k-1}\right)\right)$';
    if strcmp(obj.algorithm.sequential.trustRegion.type,'hard') && obj.algorithm.sequential.activate
        hplotP2(2) = semilogy(1:iterations,radius,'k--o','linewidth',2); hold on
        legendLabelP2{2} = '$\eta_{y}$';
    end
    ax = gca;
    ax.YGrid = 'on';
    ax.YMinorGrid = 'on';
    xlabel('Iteration $k$ [-]','interpreter','latex')
    ylabel('Trust Region Radius [-]','interpreter','latex')
    legend(hplotP2,legendLabelP2,'location','best','interpreter','latex')
    scoptSave2pdf(obj,'trustRegionPenaltyN')
    hfig2 = gcf;
else
    hfig2 = nan;
end
end


