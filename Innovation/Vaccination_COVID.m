%A script to predict compare change in Susceptible population due to different vaccination rate during Covid 

function Vaccination_plot

params.mu = 0.00005; 
params.bi = 0.65322; 
params.bd = 0.64845; 
params.bh = 0.00029; 
params.sigma = 0.1; 
params.mue = 0.14548; 
params.gir = 0.092; 
params.tau = 0.19013; 
params.delta = 0.33; 
params.ghr = 0.1103; 
params.v = 0.0; 
params.v2 = 0.05; 
params.v3 = 0.1; 


initial.S = 100; 
initial.E = 9;
initial.I = 6; 
initial.D = 2;
initial.H = 3;
initial.R = 1;

end_time = 3; % end of simulation time span starting a 0

[t, y] = ode45(@(t, x) derivative1(t, x, params), ...
[0 end_time], ...
[initial.S;initial.E; initial.I;initial.D;initial.H;initial.R], ...
[]);

outS = y(:,1);
outE = y(:,2);
outI = y(:,3);
outD = y(:,4);
outH = y(:,5);
outR = y(:,6);
xx = plot(t, outS,'r','Linewidth',0.5);
hold on

[t, y] = ode45(@(t, x) derivative2(t, x, params), ...
[0 end_time], ...
[initial.S;initial.E; initial.I;initial.D;initial.H;initial.R], ...
[]);

outS2 = y(:,1);
outE2 = y(:,2);
outI2 = y(:,3);
outD2 = y(:,4);
outH2 = y(:,5);
outR2 = y(:,6);
yy = plot(t, outS2,'b','Linewidth',0.5);
title('Variation in Susceptible population by COVID vaccination')
xlabel('Time in days')
ylabel('Number of Individuals')
hold on

[t, y] = ode45(@(t, x) derivative3(t, x, params), ...
[0 end_time], ...
[initial.S;initial.E; initial.I;initial.D;initial.H;initial.R], ...
[]);

outS3 = y(:,1);
outE3 = y(:,2);
outI3 = y(:,3);
outD3 = y(:,4);
outH3 = y(:,5);
outR3 = y(:,6);

zz = plot(t, outS3,'g','Linewidth',0.5);

legend([xx,yy,zz],'vaccination rate = 0','vaccination rate = 0.05','vaccination rate = 0.1');
saveas(gcf,'Vaccination_COVID.pdf')

function f = derivative1 (~, x, params)

S = x(1);
E = x(2);
I = x(3);
D = x(4);
H = x(5);
R = x(6);

ds = params.mu - params.bi*I*S - params.bd*D*S - params.bh*H*S - params.mu*S - params.v*S ;
de = params.bi*I*S + params.bd*D*S + params.bh*H*S - params.mu*E - params.sigma*E - params.v*E;
di = params.sigma*E - params.gir*I - params.mue*I - params.tau*I - params.mu*I;
dd = params.mue*I - params.delta*D - params.mu*D;
dh = params.tau*I - params.ghr*H - params.mue*H - params.mu*H;
dr = params.gir * I - params.mu * R + params.ghr*H + params.v*S + params.v*E;

f = [ds; de; di; dd; dh; dr];


function f = derivative2 (~, x, params)


S = x(1);
E = x(2);
I = x(3);
D = x(4);
H = x(5);
R = x(6);

ds = params.mu - params.bi*I*S - params.bd*D*S - params.bh*H*S - params.mu*S - params.v2*S*5 ;
de = params.bi*I*S + params.bd*D*S + params.bh*H*S - params.mu*E - params.sigma*E - params.v2*E;
di = params.sigma*E - params.gir*I - params.mue*I - params.tau*I - params.mu*I;
dd = params.mue*I - params.delta*D - params.mu*D;
dh = params.tau*I - params.ghr*H - params.mue*H - params.mu*H;
dr = params.gir * I - params.mu * R + params.ghr*H + params.v2*S + params.v2*E;

f = [ds; de; di; dd; dh; dr];

function f = derivative3 (~, x, params)

S = x(1);
E = x(2);
I = x(3);
D = x(4);
H = x(5);
R = x(6);

ds = params.mu - params.bi*I*S - params.bd*D*S - params.bh*H*S - params.mu*S - params.v3*S*10 ;
de = params.bi*I*S + params.bd*D*S + params.bh*H*S - params.mu*E - params.sigma*E - params.v3*E;
di = params.sigma*E - params.gir*I - params.mue*I - params.tau*I - params.mu*I;
dd = params.mue*I - params.delta*D - params.mu*D;
dh = params.tau*I - params.ghr*H - params.mue*H - params.mu*H;
dr = params.gir * I - params.mu * R + params.ghr*H + params.v3*S + params.v3*E;

f = [ds; de; di; dd; dh; dr];