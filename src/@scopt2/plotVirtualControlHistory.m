function [virtualControlsP1_ALL,virtualControlsP2_ALL,hfig1,hfig2] =  plotVirtualControlHistory(obj,fig1,fig2)
% Plots virtual controls internal slack variable. User has to
% verify their values are close to zero or nonsignificant
% during the last iterations
if obj.algorithm.virtualControl.phases(1).include
    [virtualControlsP1_ALL,virtualControlsP2_ALL] =  computePenaltiesVirtualControls(obj,obj.debugging.virtualControls);
    
    defaultColor = autumn(obj.solution.iterations);
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
        for kk = 1:obj.solution.iterations
            semilogy(1:(obj.solution.Nodes-1),virtualControlsP1_ALL(1,:,kk),'-*','linewidth',2,'color',defaultColor(kk,:)); hold on
            labelLegend = [labelLegend , ['$k$ = ' num2str(kk)]];
        end
        ax = gca;
        ax.YGrid = 'on';
        ax.YMinorGrid = 'on';
        xlabel('Node $i$ [-]','interpreter','latex')
        ylabel('$P^c\left(E\nu_i\right)$ [-]','interpreter','latex')
        legend(labelLegend,'location','best','interpreter','latex')
        scoptSave2pdf(obj,'virtualControlsPenaltyC')
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
        semilogy(1:obj.solution.iterations,squeeze(virtualControlsP2_ALL),'-*','linewidth',2);
        ax = gca;
        ax.YGrid = 'on';
        ax.YMinorGrid = 'on';
        xlabel('Iteration $k$ [-]','interpreter','latex')
        ylabel('$P^n\left(P^c\left(E\nu\right)\right)$ [-]','interpreter','latex')
        legend(labelLegend,'location','best','interpreter','latex')
        scoptSave2pdf(obj,'virtualControlsPenaltyN')
        hfig2 = gcf;
    else
        hfig2 = nan;
    end
else
    hfig1 = nan;
    hfig2 = nan;
end
end


