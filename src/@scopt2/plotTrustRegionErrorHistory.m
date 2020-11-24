function [trustRegionErrorP1_ALL,trustRegionErrorP2_ALL,hfig1,hfig2] =  plotTrustRegionErrorHistory(obj,fig1,fig2)
% plots to verify trust regions
if strcmp(obj.algorithm.sequential.type,'trust-region')
    trustRegionP1   = obj.debugging.trustRegionP1;
    trustRegionP2   = obj.debugging.trustRegionP2;
    radius          = obj.debugging.radius;
    
    iterations = obj.solution.iterations;
    
    %             lineOpacity_vector = linspace(1,0,iterations);
    defaultColor = autumn(iterations);
    
    Nodes      = size(trustRegionP1,2);
    [trustRegionP1_ALL,trustRegionP2_ALL] =  plotTrustRegionHistory(obj,0,0);
    
    trustRegionErrorP1_ALL = abs(trustRegionP1 - trustRegionP1_ALL);
    trustRegionErrorP2_ALL = abs(trustRegionP2 - trustRegionP2_ALL);
    
    if obj.algorithm.sequential.activate
        radiusErrorP2_ALL      = abs(radius - trustRegionP2_ALL);
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
        %                 labelLegend = strsplit(num2str(1:iterations));
        labelLegend = {};
        for kk = 1:iterations
            hplot(kk)   = semilogy(1:Nodes,trustRegionErrorP1_ALL(1,:,kk),'-*','linewidth',2,'color',defaultColor(kk,:)); hold on
            labelLegend = [labelLegend , ['$k$ = ' num2str(kk)]];
        end
        ax = gca;
        ax.YGrid = 'on';
        ax.YMinorGrid = 'on';
        xlabel('Node $i$ [-]','interpreter','latex')
        ylabel('$|P^c\left(y^k_i-y^{k-1}_i\right)-\eta^{P1}|$ [-]','interpreter','latex')
        legend(hplot,labelLegend,'location','best','interpreter','latex')
        scoptSave2pdf(obj,'trustRegionPenaltyCErrors')
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
        hplotP2(1) = semilogy(1:iterations,squeeze(trustRegionErrorP2_ALL),'-*','linewidth',2); hold on
        legendLabelP2{1}     = '$|P^n\left(P^c\left(y^k-y^{k-1}\right)\right)-\eta^{P2}|$';
        if strcmp(obj.algorithm.sequential.trustRegion.type,'hard') && obj.algorithm.sequential.activate
            hplotP2(2) = semilogy(1:iterations,radiusErrorP2_ALL,'k--','linewidth',2); hold on
            legendLabelP2{2} = '$|P^n\left(P^c\left(y^k-y^{k-1}\right)\right)-\eta_{y}|$';
        end
        ax = gca;
        ax.YGrid = 'on';
        ax.YMinorGrid = 'on';
        xlabel('Iteration $k$ [-]','interpreter','latex')
        ylabel('Trust Region Error [-]','interpreter','latex')
        legend(hplotP2,legendLabelP2,'location','best','interpreter','latex')
        scoptSave2pdf(obj,'trustRegionPenaltyNErrors')
        hfig2 = gcf;
    else
        hfig2 = nan;
    end
else
    hfig1 = nan;
    hfig2 = nan;
end
end


