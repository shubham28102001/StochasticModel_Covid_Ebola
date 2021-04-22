%A script to load the Stochastic Ebola Data and plot it

clear
load EVD_Simulation_data
U=EVD_Simulation_data;

figure(30)
set(gca,'FontSize',20)
plot(U(:,1),U(:,4),'k')
axis([0 1000 0 100])
title 'Extinction Time of Ebola'
ylabel 'Number of Infected Individuals'
xlabel 'Optimal Path Points'
print -f30 -dpdf InfectiousOverTime_EVD.pdf   