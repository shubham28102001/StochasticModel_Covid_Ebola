%A script to load the Stochastic data for Ebola and Covid and plot them and
%compare the same

clear,clc;
load EVD_Simulation_data
load COVID_Simulation_data
U=EVD_Simulation_data;
V=COVID_Simulation_data;

figure(1)
set(gca,'FontSize',20)
e=plot(U(:,1),U(:,4),'g')
hold on
c=plot(V(:,1),V(:,4),'r')
hold on
axis([0 1000 0 100])
legend([e,c],'Ebola','Covid');
xlabel('Optimal Path Points')
ylabel('Number of Infected Individuals')
title 'Extinction Time'
saveas(gcf,'ExtinctionTime_EVD_COVID.pdf')