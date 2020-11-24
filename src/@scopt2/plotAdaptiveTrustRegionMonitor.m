function hfig =  plotAdaptiveTrustRegionMonitor(obj,fig)
% Plots a monitor of the adaptive trust region and objective
% values
vectorIterations = 1:obj.solution.iterations;

if nargin == 2
    figure(fig)
else
    figure
end

set(gcf,'units','normalized','position',[0 0 0.7 1])
subplot 311
plot([0 vectorIterations],[obj.problem.initialGuess.J obj.debugging.J],'-*','linewidth',2);hold on
plot(vectorIterations,obj.debugging.L,'-.*','linewidth',2);
plot(vectorIterations,obj.debugging.Jaug,'--*','linewidth',2);
plot(vectorIterations,obj.debugging.Laug,':*','linewidth',2);

xlabel('Iteration $k$ [-]','interpreter','latex')
ylabel('Objective [-]','interpreter','latex')
legend({'Non-linear Cost Function','Cost Function','Non-linear Penalised Cost Function','Penalised Cost Function'},'interpreter','latex')

if sum(obj.debugging.J<=0)<=0
    set(gca, 'YScale', 'log')
end
grid on

subplot 312
indexRejected = find(obj.debugging.rejected>0);
plot(vectorIterations,obj.debugging.rho,'r-*','linewidth',2);hold on
plot(vectorIterations,obj.algorithm.sequential.trustRegion.rho0*ones(1,obj.solution.iterations),'k--','linewidth',1);
plot(vectorIterations,obj.algorithm.sequential.trustRegion.rho1*ones(1,obj.solution.iterations),'k-.','linewidth',1);
plot(vectorIterations,obj.algorithm.sequential.trustRegion.rho2*ones(1,obj.solution.iterations),'k:','linewidth',1);
labels = {'$\rho$','$\rho_0$','$\rho_1$','$\rho_2$'};
if ~isempty(indexRejected)
    plot(vectorIterations(indexRejected),obj.debugging.rejected(indexRejected).*obj.debugging.rho(indexRejected),'bo');
    labels = [labels,'Rejected Step'];
end
xlabel('Iteration [-]','interpreter','latex')
ylabel('Ratio $\Delta J$/$\Delta L$ [-]','interpreter','latex')
legend(labels,'interpreter','latex','location','best');
grid on

subplot 313
indexRejected = find(obj.debugging.rejected>0);

if strcmp(obj.algorithm.sequential.trustRegion.type,'hard')
    plot(vectorIterations,obj.debugging.radius,'r-*','linewidth',2);hold on
    plot(vectorIterations,obj.debugging.radius(1)*ones(1,obj.solution.iterations),'k--','linewidth',1);
    labels = {'$\eta_y$','$\eta_{y1}$'};
    if ~isempty(indexRejected)
        plot(vectorIterations(indexRejected),obj.debugging.rejected(indexRejected).*obj.debugging.radius(indexRejected),'bo');
        labels = [labels,'Rejected Step'];
    end
    if find(obj.algorithm.sequential.trustRegion.radius_l==obj.debugging.radius)
        plot(vectorIterations,obj.algorithm.sequential.trustRegion.radius_l *ones(1,obj.solution.iterations),'b-.','linewidth',1);
        labels = [labels,'$\eta_l$'];
    end
    if find(obj.algorithm.sequential.trustRegion.radius_u==obj.debugging.radius)
        plot(vectorIterations,obj.algorithm.sequential.trustRegion.radius_u *ones(1,obj.solution.iterations),'g-.','linewidth',1);
        labels = [labels,'$\eta_u$'];
    end
    ylabel('$\eta_y$ [-]','interpreter','latex')
elseif strcmp(obj.algorithm.sequential.trustRegion.type,'soft')
    plot(vectorIterations,obj.debugging.lambda,'r-*','linewidth',2);hold on
    plot(vectorIterations,obj.debugging.lambda(1)*ones(1,obj.solution.iterations),'k--','linewidth',1);
    labels = {'$\lambda_y$','$\lambda_{y1}$'};
    if ~isempty(indexRejected)
        plot(vectorIterations(indexRejected),obj.debugging.rejected(indexRejected).*obj.debugging.lambda(indexRejected),'bo');
        labels = [labels,'Rejected Step'];
    end
    if find(obj.algorithm.sequential.trustRegion.lambda_l==obj.debugging.lambda)
        plot(vectorIterations,obj.algorithm.sequential.trustRegion.lambda_l *ones(1,obj.solution.iterations),'b-.','linewidth',1);
        labels = [labels,'$\lambda_l$'];
    end
    if find(obj.algorithm.sequential.trustRegion.lambda_u==obj.debugging.radius)
        plot(vectorIterations,obj.algorithm.sequential.trustRegion.lambda_u *ones(1,obj.solution.iterations),'g-.','linewidth',1);
        labels = [labels,'$\lambda_u$'];
    end
    ylabel('Weight [-]','interpreter','latex')
end

xlabel('Iteration [-]','interpreter','latex')
legend(labels,'interpreter','latex','location','best');
set(gca, 'YScale', 'log')
grid on

scoptSave2pdf(obj,'adaptiveTrustRegionMonitor')
hfig = gcf;

drawnow

end


