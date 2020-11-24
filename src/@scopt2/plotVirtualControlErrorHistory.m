function [virtualControlsErrorP1_ALL,virtualControlsErrorP2_ALL,hfig1,hfig2] =  plotVirtualControlErrorHistory(obj,fig1,fig2)
% generates figures to verify virtual control errors, values
% should have small numbers
if obj.algorithm.virtualControl.phases.include
    virtualControlsP1   = obj.debugging.virtualControlsP1;
    virtualControlsP2   = obj.debugging.virtualControlsP2;
    iterations = obj.solution.iterations;
    defaultColor = autumn(iterations);
    
    Nodes      = size(virtualControlsP1,2);
    [virtualControlsP1_ALL,virtualControlsP2_ALL] =  computePenaltiesVirtualControls(obj,obj.debugging.virtualControls);
    
    virtualControlsErrorP1_ALL = abs(virtualControlsP1 - virtualControlsP1_ALL);
    virtualControlsErrorP2_ALL = abs(virtualControlsP2 - virtualControlsP2_ALL);
    
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
            hplot(kk)   = semilogy(1:Nodes,virtualControlsErrorP1_ALL(1,:,kk),'-*','linewidth',2,'color',defaultColor(kk,:)); hold on
            labelLegend = [labelLegend , ['$k$ = ' num2str(kk)]];
        end
        ax = gca;
        ax.YGrid = 'on';
        ax.YMinorGrid = 'on';
        xlabel('Node $i$ [-]','interpreter','latex')
        ylabel('$|P^c\left(E\nu_i\right)-\nu^{P1}|$ [-]','interpreter','latex')
        legend(hplot,labelLegend,'location','best','interpreter','latex')
        scoptSave2pdf(obj,'virtualControlsPenaltyCErrors')
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
        semilogy(1:iterations,squeeze(virtualControlsErrorP2_ALL),'-*','linewidth',2); hold on
        ax = gca;
        ax.YGrid = 'on';
        ax.YMinorGrid = 'on';
        xlabel('Iteration $k$ [-]','interpreter','latex')
        ylabel('$|P^n\left(P^c\left(E\nu\right)\right)-\nu^{P2}|$ [-]','interpreter','latex')
        scoptSave2pdf(obj,'virtualControlsPenaltyNErrors')
        hfig2 = gcf;
    else
        hfig2 = nan;
    end
end
end


