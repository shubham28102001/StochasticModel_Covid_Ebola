clear,clc;
load EVD_Simulation_data
U=EVD_Simulation_data;

mu=0.00005;
Bi=0.5;
Bd=0.6;
Bh=0.00016;
N=500000;
sigma=0.10;
gir=0.07;
mue=0.12;
tau=0.20;
delta=0.33;
ghr=0.10;
k=0;

Rzero = sigma*(Bi+Bd*mue/(delta+mu)+Bh*tau/(ghr+mue+mu))/((gir+mue+tau+mu)*(mu+sigma)); 
format long
disp('Basic Reproduction Number for Ebola =')
Rzero

[numRows,numCols] = size(U);

Sa = zeros(numRows,1);
Ea = zeros(numRows,1);
Ia = zeros(numRows,1);
Da = zeros(numRows,1);
Ha = zeros(numRows,1);
Ra = zeros(numRows,1);
Ta = zeros(numRows,1);

for i = 1:numRows
    X=[U(i,2) U(i,3) U(i,4) U(i,5) U(i,6) U(i,7)];
    Sa(i,1)=(mu*N-mu*X(1)-Bi*X(3)*X(1)/N-Bd*X(4)*X(1)/N-Bh*X(5)*X(1)/N)/8;
    Ea(i,1)=(Bi*X(3)*X(1)/N+Bd*X(4)*X(1)/N+Bh*X(5)*X(1)/N-(sigma+mu)*X(2))/10000;
    Ia(i,1)=(sigma*X(2)-(gir+mue+tau)*X(3)-mu*X(3))/20000;
    Da(i,1)=(mue*X(3)-delta*X(4)-mu*X(4))/50000;
    Ha(i,1)=(tau*X(3)-(ghr+mue)*X(5)-mu*X(5))/20000;
    Ra(i,1)=(ghr*X(5)+gir*X(3)-mu*X(6))/5;
    Ta(i,1)=i;
end

figure(1)
plot(Ta(1:1500), Sa(1:1500));
title('Variation in Susceptible population with resivoiur transmission 0')
xlabel('Optimal path points')
ylabel('Rate of change in Susceptible Popultion')
hold on
figure(2)
plot(Ta(1:1500), Ea(1:1500));
title('Variation in Exposed population with resivoiur transmission 0')
xlabel('Optimal path points')
ylabel('Rate of change in Exposed Popultion')
hold on
figure(3)
plot(Ta(1:1500), Ia(1:1500));
title('Variation in Infected population with resivoiur transmission 0')
xlabel('Optimal path points')
ylabel('Rate of change in Infected Popultion')
hold on
figure(4)
plot(Ta(1:1500), Da(1:1500));
title('Variation in Deceased population with resivoiur transmission 0')
xlabel('Optimal path points')
ylabel('Rate of change in Deceased Popultion')
hold on
figure(5)
plot(Ta(1:1500), Ha(1:1500));
title('Variation in Hospitalized population with resivoiur transmission 0')
xlabel('Optimal path points')
ylabel('Rate of change in Hospitalized Popultion')
hold on
figure(6)
plot(Ta(1:1500), Ra(1:1500));
title('Variation in Recovered population with resivoiur transmission 0')
xlabel('Optimal path points')
ylabel('Rate of change in Recovered Popultion')
hold on