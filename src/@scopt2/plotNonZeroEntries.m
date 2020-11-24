function plotNonZeroEntries(MAT,figNum,index)
% Function to plot nonzero entries in a matrix
[m,n]=size(MAT);
[row,col]  = find(abs(MAT)>0.00000001);

if nargin >=3
    figure(figNum)
end

plot(col,row,'o') ;
hold on
plot([0.5 n+0.5],[0.5 0.5],'k');
plot([0.5 n+0.5],[m+0.5 m+0.5],'k');
plot([0.5 0.5],[0.5 m+0.5],'k');
plot([n+0.5 n+0.5],[0.5 m+0.5],'k');


axis([0 n+1 0 m+1]);

if nargin==4
    plot(xlim+[0.5,-0.5],index.controls(1)*ones(1,2)-0.5,'r','linewidth',1)
    if ~isempty(index.virtualControls)
        plot(xlim+[0.5,-0.5],index.virtualControls(1)*ones(1,2)-0.5,'r','linewidth',1)
    end
    if ~isempty(index.virtualControlsP1s)
        plot(xlim+[0.5,-0.5],index.virtualControlsP1s(1)*ones(1,2)-0.5,'r','linewidth',1)
    end
    if ~isempty(index.virtualControlsP1)
        plot(xlim+[0.5,-0.5],index.virtualControlsP1(1)*ones(1,2)-0.5,'r','linewidth',1)
    end
    if ~isempty(index.virtualControlsP2)
        plot(xlim+[0.5,-0.5],index.virtualControlsP2*ones(1,2)-0.5,'r','linewidth',1)
    end
    if ~isempty(index.virtualBuffers.events)
        plot(xlim+[0.5,-0.5],index.virtualBuffers.events(1)*ones(1,2)-0.5,'r','linewidth',1)
    end
    if ~isempty(index.virtualBuffers.eventsP)
        plot(xlim+[0.5,-0.5],index.virtualBuffers.eventsP(1)*ones(1,2)-0.5,'r','linewidth',1)
    end
    if ~isempty(index.virtualBuffers.path)
        plot(xlim+[0.5,-0.5],index.virtualBuffers.path(1)*ones(1,2)-0.5,'r','linewidth',1)
    end
    if ~isempty(index.virtualBuffers.pathPs)
        plot(xlim+[0.5,-0.5],index.virtualBuffers.pathPs(1)*ones(1,2)-0.5,'r','linewidth',1)
    end
    if ~isempty(index.virtualBuffers.pathP)
        plot(xlim+[0.5,-0.5],index.virtualBuffers.pathP(1)*ones(1,2)-0.5,'r','linewidth',1)
    end
    if ~isempty(index.objectiveLagrange)
        plot(xlim+[0.5,-0.5],index.objectiveLagrange(1)*ones(1,2)-0.5,'r','linewidth',1)
    end
    if ~isempty(index.objectiveSlack)
        plot(xlim+[0.5,-0.5],index.objectiveSlack*ones(1,2)-0.5,'r','linewidth',1)
    end
    
    if ~isempty(index.trustRegionP1)
        plot(xlim+[0.5,-0.5],index.trustRegionP1(1)*ones(1,2)-0.5,'r','linewidth',1)
    end
    if ~isempty(index.trustRegionP2)
        plot(xlim+[0.5,-0.5],index.trustRegionP2*ones(1,2)-0.5,'r','linewidth',1)
    end
    
    if ~isempty(index.timeFinalSlack)
        plot(xlim+[0.5,-0.5],index.timeFinalSlack*ones(1,2)-0.5,'r','linewidth',1)
    end
    if ~isempty(index.timeFinal)
        plot(xlim+[0.5,-0.5],index.timeFinal*ones(1,2)-0.5,'r','linewidth',1)
    end
end

grid on
set (gca,'Ydir','reverse')
end


