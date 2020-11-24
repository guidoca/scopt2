function hfig =  plotRealObjective(obj,fig)
% Plots the real objective values
vectorIterations = 1:obj.solution.iterations;

if nargin == 2
    figure(fig)
else
    figure
end
plot([0,vectorIterations],[obj.problem.initialGuess.Jreal,obj.debugging.Jreal],'-*','linewidth',2);hold on
plot(vectorIterations,obj.debugging.Lreal,'-o','linewidth',2);hold on

xlabel('Iteration $k$ [-]','interpreter','latex')
ylabel('Objective [-]','interpreter','latex')
legend({'Non-convex','Convex'},'interpreter','latex','location','best')
% legend({'Non-linear Cost Function','Cost Function','Non-linear Penalized Cost Function','Penalized Cost Function'},'interpreter','latex')
grid on

scoptSave2pdf(obj,'realObjective')
hfig = gcf;

drawnow
end


