function hfig =  plotTime(obj,fig)
% Plots the free final time values
if obj.problem.phases(1).freeTimeFinal
    vectorIterations = 0:obj.solution.iterations;
    
    if nargin == 2
        figure(fig)
    else
        figure
    end
    
    plot(vectorIterations,[obj.problem.phases.initialGuess.timeFinal,obj.debugging.timeFinal_ALL],'-*','linewidth',2);hold on
    plot(vectorIterations,[obj.problem.phases.initialGuess.timeInitial,obj.debugging.timeInitial_ALL],'-o','linewidth',2);hold on
    
    xlabel('Iteration $k$ [-]','interpreter','latex')
    ylabel('Time [-]','interpreter','latex')
    legend({'Final Time','Initial Time'},'interpreter','latex','location','best')
    % legend({'Non-linear Cost Function','Cost Function','Non-linear Penalized Cost Function','Penalized Cost Function'},'interpreter','latex')
    grid on
    
    scoptSave2pdf(obj,'monitorTime')
    hfig = gcf;
    
    drawnow
else
    hfig = nan;
end
end

