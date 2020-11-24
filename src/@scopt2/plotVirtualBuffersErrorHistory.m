function [virtualBuffersPathErrorPs,virtualBuffersPathErrorP,virtualBuffersEventsErrorP,hfig1,hfig2] =  plotVirtualBuffersErrorHistory(obj,fig1,fig2)
% generates plots to verify virtual buffers
if obj.algorithm.virtualBuffer.phases(1).include
    
    virtualBuffersPathPs  = obj.debugging.virtualBuffers.pathPs;
    virtualBuffersPathP   = obj.debugging.virtualBuffers.pathP;
    virtualBuffersEventsP  = obj.debugging.virtualBuffers.eventsP;
    
    [virtualBuffersPathPs_ALL,virtualBuffersPathP_ALL,virtualBuffersEventsP_ALL] =  plotVirtualBuffersHistory(obj,0,0);
    
    iterations = obj.solution.iterations;
    Nodes      = size(virtualBuffersPathPs,2);
    NvbPath    = size(virtualBuffersPathPs,1);
    NvbEvent   = size(virtualBuffersEventsP,1);
    
    
    defaultColorPath   = autumn(NvbPath);
    defaultColorEvent  = winter(NvbEvent);
    lineOpacity  = linspace(0,1,iterations+1);lineOpacity(1) = [];
    
    virtualBuffersPathErrorPs = abs(virtualBuffersPathPs - virtualBuffersPathPs_ALL);
    virtualBuffersPathErrorP  = abs(virtualBuffersPathP  - virtualBuffersPathP_ALL);
    virtualBuffersEventsErrorP = abs(virtualBuffersEventsP  - virtualBuffersEventsP_ALL);
    
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
    if ~dontPlot && (NvbPath)>0
        labelLegend = {};
        for jj = 1:NvbPath
            for kk = 1:iterations
                hplotAux   = semilogy(1:Nodes,virtualBuffersPathErrorPs(jj,:,kk),'o-','linewidth',1,'Color',[1,1,1]*(1-lineOpacity(kk)) + lineOpacity(kk)*defaultColorPath(jj,:),'MarkerFaceColor',[1,1,1]*(1-lineOpacity(kk)) + lineOpacity(kk)*defaultColorPath(jj,:),'MarkerEdgeColor',defaultColorPath(jj,:)); hold on
            end
            hplot(jj)   = hplotAux;
            labelLegend = [labelLegend , ['Path Buffer Zone $j$ = ' num2str(jj)]];
        end
        
        ax = gca;
        ax.YGrid      = 'on';
        ax.YMinorGrid = 'on';
        xlabel('Node $i$ [-]','interpreter','latex')
        ylabel('$|\max{\left(0,\xi_{g,i,j}\right)} - \xi^{P1}_{g,i,j}|$ [-]','interpreter','latex')
        legend(hplot,labelLegend,'location','best','interpreter','latex')
        scoptSave2pdf(obj,'virtualBuffersPathPenaltyMaxErrors')
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
            hplot(jj)   = semilogy(1:iterations,squeeze(virtualBuffersPathErrorP(jj,:)),'o-','linewidth',1,'color',defaultColorPath(jj,:)); hold on
            labelLegend = [labelLegend , ['Path Buffer Zone $j$ = ' num2str(jj)]];
        end
        
        for jj = 1:NvbEvent
            hplot(jj+NvbPath)   = semilogy(1:iterations,squeeze(virtualBuffersEventsErrorP(jj,:)),'*-','linewidth',1,'color',defaultColorEvent(jj,:)); hold on
            labelLegend = [labelLegend , ['Event Buffer Zone $j$ = ' num2str(jj)]];
        end
        ax = gca;
        ax.YGrid = 'on';
        ax.YMinorGrid = 'on';
        xlabel('Iteration $k$ [-]','interpreter','latex')
        ylabel('$|P^n\left(\max{\left(0,\xi_{g,j}\right)}\right)-\xi^{P2}_{e,j}|,\max{\left(0,\xi_{e,j}\right)} - \xi^{P1}_{e,j}$ [-]','interpreter','latex')
        legend(hplot,labelLegend,'location','best','interpreter','latex')
        scoptSave2pdf(obj,'virtualBufferZonesFinalPenalty')
        hfig2 = gcf;
    else
        hfig2 = nan;
    end
else
    hfig1 = figure(fig1);set(hfig1, 'Visible','off')
    hfig2 = figure(fig2);set(hfig2, 'Visible','off')
    virtualBuffersPathErrorPs = nan;
    virtualBuffersPathErrorP = nan;
    virtualBuffersEventsErrorP = nan;
end
end


