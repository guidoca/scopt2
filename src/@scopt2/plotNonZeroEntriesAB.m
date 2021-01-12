
function plotNonZeroEntriesAB(A,b,figNum,index)
% Function to plot nonzero entries in matrix A and column
% vector b side by side
figure(figNum)
n = 7;
subplot(1,n,1:n-1)
scopt2.plotNonZeroEntries(A);
if nargin ==5
    hold on
    %     plot(index.states(1)*ones(1,2),ylim,'r','linewidth',2)
    plot(index.controls(1)*ones(1,2)-0.5,ylim+[0.5,-0.5],'r','linewidth',1)
    if ~isempty(index.virtualControls)
        plot(index.virtualControls(1)*ones(1,2)-0.5,ylim+[0.5,-0.5],'r','linewidth',1)
    end
    if ~isempty(index.virtualControlsP1s)
        plot(index.virtualControlsP1s(1)*ones(1,2)-0.5,ylim+[0.5,-0.5],'r','linewidth',1)
    end
    if ~isempty(index.virtualControlsP1)
        plot(index.virtualControlsP1(1)*ones(1,2)-0.5,ylim+[0.5,-0.5],'r','linewidth',1)
    end
    if ~isempty(index.virtualControlsP2)
        plot(index.virtualControlsP2*ones(1,2)-0.5,ylim+[0.5,-0.5],'r','linewidth',1)
    end
    if ~isempty(index.virtualBuffers.events)
        plot(index.virtualBuffers.events(1)*ones(1,2)-0.5,ylim+[0.5,-0.5],'r','linewidth',1)
    end
    if ~isempty(index.virtualBuffers.eventsP)
        plot(index.virtualBuffers.eventsP(1)*ones(1,2)-0.5,ylim+[0.5,-0.5],'r','linewidth',1)
    end
    if ~isempty(index.virtualBuffers.path)
        plot(index.virtualBuffers.path(1)*ones(1,2)-0.5,ylim+[0.5,-0.5],'r','linewidth',1)
    end
    if ~isempty(index.virtualBuffers.pathPs)
        plot(index.virtualBuffers.pathPs(1)*ones(1,2)-0.5,ylim+[0.5,-0.5],'r','linewidth',1)
    end
    if ~isempty(index.virtualBuffers.pathP)
        plot(index.virtualBuffers.pathP(1)*ones(1,2)-0.5,ylim+[0.5,-0.5],'r','linewidth',1)
    end
    if ~isempty(index.objectiveLagrange)
        plot(index.objectiveLagrange(1)*ones(1,2)-0.5,ylim+[0.5,-0.5],'r','linewidth',1)
    end
    if ~isempty(index.objectiveSlack)
        plot(index.objectiveSlack*ones(1,2)-0.5,ylim+[0.5,-0.5],'r','linewidth',1)
    end
    
    if ~isempty(index.trustRegionP1)
        plot(index.trustRegionP1(1)*ones(1,2)-0.5,ylim+[0.5,-0.5],'r','linewidth',1)
    end
    
    if ~isempty(index.trustRegionP2)
        plot(index.trustRegionP2*ones(1,2)-0.5,ylim+[0.5,-0.5],'r','linewidth',1)
    end
    
    if ~isempty(index.timeFinalSlack)
        plot(index.timeFinalSlack*ones(1,2)-0.5,ylim+[0.5,-0.5],'r','linewidth',1)
    end
    if ~isempty(index.timeFinal)
        plot(index.timeFinal*ones(1,2)-0.5,ylim+[0.5,-0.5],'r','linewidth',1)
    end
end

subplot(1,n,n)
scopt2.plotNonZeroEntries(b);
end

