%A script to compare the change in S, E, I, D, H, R classes during Ebola
%and Corona Virus

clear,clc;
load EVD_Simulation_data
load COVID_Simulation_data
U=EVD_Simulation_data;
V=COVID_Simulation_data;

N=500000;
k=0;

%Paramaters for EVD
e_mu=0.00005;
e_Bi=0.5;
e_Bd=0.6;
e_Bh=0.00016;
e_sigma=0.10;
e_gir=0.07;
e_mue=0.12;
e_tau=0.20;
e_delta=0.33;
e_ghr=0.10;

%Paramaters for COVID
c_mu=0.00005;
c_Bi=0.65322;
c_Bd=0.64845;
c_Bh=0.00029;
c_sigma=0.10;
c_gir=0.092;
c_mue=0.14548;
c_tau=0.19013;
c_delta=0.33;
c_ghr=0.1103;

[numRows,numCols] = size(U);

Ta = zeros(numRows,1);

e_Sa = zeros(numRows,1);
e_Ea = zeros(numRows,1);
e_Ia = zeros(numRows,1);
e_Da = zeros(numRows,1);
e_Ha = zeros(numRows,1);
e_Ra = zeros(numRows,1);

c_Sa = zeros(numRows,1);
c_Ea = zeros(numRows,1);
c_Ia = zeros(numRows,1);
c_Da = zeros(numRows,1);
c_Ha = zeros(numRows,1);
c_Ra = zeros(numRows,1);
c_Ta = zeros(numRows,1);

for i = 1:numRows
    X=[U(i,2) U(i,3) U(i,4) U(i,5) U(i,6) U(i,7)];
    e_Sa(i,1)=(e_mu*N-e_mu*X(1)-e_Bi*X(3)*X(1)/N-e_Bd*X(4)*X(1)/N-e_Bh*X(5)*X(1)/N)/8;
    e_Ea(i,1)=(e_Bi*X(3)*X(1)/N+e_Bd*X(4)*X(1)/N+e_Bh*X(5)*X(1)/N-(e_sigma+e_mu)*X(2))/10000;
    e_Ia(i,1)=(e_sigma*X(2)-(e_gir+e_mue+e_tau)*X(3)-e_mu*X(3))/20000;
    e_Da(i,1)=(e_mue*X(3)-e_delta*X(4)-e_mu*X(4))/50000;
    e_Ha(i,1)=(e_tau*X(3)-(e_ghr+e_mue)*X(5)-e_mu*X(5))/20000;
    e_Ra(i,1)=(e_ghr*X(5)+e_gir*X(3)-e_mu*X(6))/5;
    Ta(i,1)=i;
end

for i = 1:numRows
    Y=[V(i,2) V(i,3) V(i,4) V(i,5) V(i,6) V(i,7)];
    c_Sa(i,1)=(c_mu*N-c_mu*Y(1)-c_Bi*Y(3)*Y(1)/N-c_Bd*Y(4)*Y(1)/N-c_Bh*Y(5)*Y(1)/N)/8;
    c_Ea(i,1)=(c_Bi*Y(3)*Y(1)/N+c_Bd*Y(4)*Y(1)/N+c_Bh*Y(5)*Y(1)/N-(c_sigma+c_mu)*Y(2))/10000;
    c_Ia(i,1)=(c_sigma*Y(2)-(c_gir+c_mue+c_tau)*Y(3)-c_mu*Y(3))/20000;
    c_Da(i,1)=(c_mue*Y(3)-c_delta*Y(4)-c_mu*Y(4))/50000;
    c_Ha(i,1)=(c_tau*Y(3)-(c_ghr+c_mue)*Y(5)-c_mu*Y(5))/20000;
    c_Ra(i,1)=(c_ghr*Y(5)+c_gir*Y(3)-c_mu*Y(6))/5;
end

figure(1)
e_s = plot(Ta(1:1500), e_Sa(1:1500), 'g');
hold on
c_s = plot(Ta(1:1500), c_Sa(1:1500), 'r');
hold on
legend([e_s,c_s],'Ebola','Covid');
title('Variation in Susceptible population with 0 reserviour transmission')
xlabel('Optimal path points')
ylabel('Rate of change in Susceptible Population')
hold on
figure(2)
e_e = plot(Ta(1:1500), e_Ea(1:1500), 'g');
hold on
c_e = plot(Ta(1:1500), c_Ea(1:1500), 'r');
hold on
legend([e_e,c_e],'Ebola','Covid');
title('Variation in Exposed population with 0 reserviour transmission')
xlabel('Optimal path points')
ylabel('Rate of change in Exposed Population')
hold on
figure(3)
e_i = plot(Ta(1:1500), e_Ia(1:1500), 'g');
hold on
c_i = plot(Ta(1:1500), c_Ia(1:1500), 'r');
hold on
legend([e_i,c_i],'Ebola','Covid');
title('Variation in Infected population with 0 reserviour transmission')
xlabel('Optimal path points')
ylabel('Rate of change in Infected Population')
hold on
figure(4)
e_d = plot(Ta(1:1500), e_Da(1:1500), 'g');
hold on
c_d = plot(Ta(1:1500), c_Da(1:1500), 'r');
hold on
legend([e_d,c_d],'Ebola','Covid');
title('Variation in Deceased population with 0 reserviour transmission')
xlabel('Optimal path points')
ylabel('Rate of change in Deceased Population')
hold on
figure(5)
e_h = plot(Ta(1:1500), e_Ha(1:1500), 'g');
hold on
c_h = plot(Ta(1:1500), c_Ha(1:1500), 'r');
hold on
legend([e_h,c_h],'Ebola','Covid');
title('Variation in Hospitalized population with 0 reserviour transmission')
xlabel('Optimal path points')
ylabel('Rate of change in Hospitalized Population')
hold on
figure(6)
e_r = plot(Ta(1:1500), e_Ra(1:1500), 'g');
hold on
c_r = plot(Ta(1:1500), c_Ra(1:1500), 'r');
hold on
legend([e_r,c_r],'Ebola','Covid');
title('Variation in Recovered population with 0 reserviour transmission')
xlabel('Optimal path points')
ylabel('Rate of change in Recovered Population')
hold on