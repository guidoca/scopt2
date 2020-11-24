function [virtualBuffersPathP1_ALL,virtualBuffersPathP2_ALL,virtualBuffersEventsP1_ALL,hfig1,hfig2] =  plotVirtualBuffersHistory(obj,fig1,fig2)
%  plots virtual buffers slack variable values
if obj.algorithm.virtualBuffer.phases(1).include
    
    virtualBuffersPath     = obj.debugging.virtualBuffers.path;
    virtualBuffersEvent    = obj.debugging.virtualBuffers.events;
    
    iterations = obj.solution.iterations;
    Nodes      = size(virtualBuffersPath,2);
    NvbPath    = size(virtualBuffersPath,1);
    NvbEvent   = size(virtualBuffersEvent,1);
    
    defaultColorPath   = autumn(NvbPath);
    defaultColorEvent  = winter(NvbEvent);
    lineOpacity  = linspace(0,1,iterations+1);lineOpacity(1) = [];
    
    virtualBuffersPathP1_ALL   = max(0,virtualBuffersPath);
    virtualBuffersEventsP1_ALL = max(0,virtualBuffersEvent);
    virtualBuffersPathP2_ALL   = reshape(sum(virtualBuffersPathP1_ALL,2),NvbPath,iterations);
    
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
    if ~dontPlot && NvbPath>0
        labelLegend = {};
        for jj = 1:NvbPath
            for kk = 1:iterations
                hplotAux   = semilogy(1:Nodes,virtualBuffersPathP1_ALL(jj,:,kk),'o-','linewidth',1,'color',[1,1,1]*(1-lineOpacity(kk)) + lineOpacity(kk)*defaultColorPath(jj,:),'MarkerFaceColor',[1,1,1]*(1-lineOpacity(kk)) + lineOpacity(kk)*defaultColorPath(jj,:),'MarkerEdgeColor',defaultColorPath(jj,:)); hold on
            end
            hplot(jj)   = hplotAux;
            labelLegend = [labelLegend , ['Path Buffer Zone $j$ = ' num2str(jj)]];
        end
        
        ax = gca;
        ax.YGrid      = 'on';
        ax.YMinorGrid = 'on';
        xlabel('Node $i$ [-]','interpreter','latex')
        ylabel('$\max{\left(0,\xi_{g,i,j}\right)}$ [-]','interpreter','latex')
        legend(hplot,labelLegend,'location','best','interpreter','latex')
        scoptSave2pdf(obj,'virtualBuffersPathPenaltyMax')
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
    if ~dontPlot && (NvbPath+NvbEvent)>0
        labelLegend = {};
        for jj = 1:NvbPath
            hplot2(jj)   = semilogy(1:iterations,squeeze(virtualBuffersPathP2_ALL(jj,:)),'o-','linewidth',1,'color',defaultColorPath(jj,:)); hold on
            labelLegend = [labelLegend , ['Path Buffer Zone $j$ = ' num2str(jj)]];
        end
        
        for jj = 1:NvbEvent
            hplot2(jj+NvbPath)   = semilogy(1:iterations,squeeze(virtualBuffersEventsP1_ALL(jj,:)),'*-','linewidth',1,'color',defaultColorEvent(jj,:)); hold on
            labelLegend = [labelLegend , ['Event Buffer Zone $j$ = ' num2str(jj)]];
        end
        ax = gca;
        ax.YGrid = 'on';
        ax.YMinorGrid = 'on';
        xlabel('Iteration $k$ [-]','interpreter','latex')
        ylabel('$P^n\left(\max{\left(0,\xi_{g,j}\right)}\right),\max{\left(0,\xi_{e,j}\right)}$ [-]','interpreter','latex')
        legend(hplot2,labelLegend,'location','best','interpreter','latex')
        scoptSave2pdf(obj,'virtualBufferZonesFinalPenalty')
        hfig2 = gcf;
    else
        hfig2 = nan;
    end
else
    hfig1 = figure(fig1);set(hfig1, 'Visible','off')
    hfig2 = figure(fig2);set(hfig2, 'Visible','off')
    virtualBuffersPathP1_ALL = nan;
    virtualBuffersPathP2_ALL = nan;
    virtualBuffersEventsP1_ALL = nan;
end
end


