%A script to load the Stochastic Ebola Data and plot it

clear
load COVID_Simulation_data
U=COVID_Simulation_data;

figure(30)
set(gca,'FontSize',20)
plot(U(:,1),U(:,4),'k')
axis([0 1000 0 100])
title 'Extinction Time for Covid'
xlabel 'Optimal Path Points'
ylabel 'Number of Infected Individuals'
print -f30 -dpdf InfectiousOverTime_COVID.pdf   