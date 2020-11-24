function hfig =  plotConvergenceHistory(obj,fig,variablesNorm)
% plots convergence history of the variables. user should
% verify it decreases during successive iterations
if nargin <3
    variablesNorm = computeConvergenceHistory(obj);
end
if nargin >= 2
    if ~ishandle(fig)
        figure(fig)
    else
        hold on
    end
else
    figure
end
semilogy(1:obj.solution.iterations,variablesNorm,'-*','linewidth',2);
ax = gca;
ax.YGrid = 'on';
ax.YMinorGrid = 'on';
xlabel('Iteration $k$','interpreter','latex')
ylabel('${||y^k-y^{k-1}||}_2$','interpreter','latex')
scoptSave2pdf(obj,'trustRegionConvergenceMonitor')

hfig = gcf;
end


