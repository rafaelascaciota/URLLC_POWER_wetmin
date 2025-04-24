%% Simulation Parameters
EH = 1;             % Number of EH devices
N = 56;             % packet length
pc = 1e-6;          % control outage
B = 150e3;          % bandwidth
alpha = 3;        % path-loss exponent
fwit = 2.45e6;      % carrier frequency of the Source %%%-----> use the same freq. for both processes
fwet = 915e6;       % carrier frequency of PB %%%-----> use the same freq. for both processes
c = 3e8;            % speed of light
eta = 0.33;         % transmiter circuit efficiency
Pcirc = 1.33e-3;    % transceiver power transmission constant
Pmax_dBm = 3.3;     % transceiver max transmission power
Pmax = 10^((Pmax_dBm-30)/10);
Pb_dBm = 45;
Pb = 10.^((Pb_dBm-30)./10);

% Non-linear Energy Harvesting
% From: Massive Wireless Energy Transfer with Multiple Power Beacons 
% for very large Internet of Things
c0 = 0.2308;             % EH unitless constants
c1 = 5.365;
w = 10.73;              % energy harvesting saturation level

%% Baseline framework equations 
xx_dBm = 0:70;
%xx = 0:5e-4:0.1;
xx = 10.^((xx_dBm-30)/10);
gf = 1e-3.*w.*(1 - exp(-c0*xx*1e3))./(1 + exp(-c0.*(xx.*1e3 - c1)));
n = 100*(gf./xx);

gfl = (10*log10(gf))+30;
xxl = (10*log10(xx))+30;
figure(1)
subplot(2,1,1);
plot(xx_dBm,gf,'-*','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w')
grid on;
ax = gca;
ax.YAxis.FontSize = 12 %for y-axis 
ay = gca;
ay.XAxis.FontSize = 12 %for y-axis
xlabel('Incident RF Power [mW], $\mathrm{P}_\mathrm{RF}$','FontSize',  16,'Interpreter','latex');  
ylabel('$\mathcal{G}(x)$ [mW]', 'FontSize',  12,'Interpreter','latex'); 

subplot(2,1,2);
plot(xx_dBm,n,'-*','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w')
grid on;
ax = gca;
ax.YAxis.FontSize = 12 %for y-axis 
ay = gca;
ay.XAxis.FontSize = 12 %for y-axis
xlabel('Incident RF Power [mW], $\mathrm{P}_\mathrm{RF}$','FontSize',  16,'Interpreter','latex');  
ylabel('$\eta$ [$\%$]', 'FontSize',  16,'Interpreter','latex'); 

figure(2)
subplot(2,1,1);
plot(xxl,gfl,'-*','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w')
grid on;
ax = gca;
ax.YAxis.FontSize = 12 %for y-axis 
ay = gca;
ay.XAxis.FontSize = 12 %for y-axis
xlabel('Incident RF Power [dBm], $\mathrm{P}_\mathrm{RF}$','FontSize',  16,'Interpreter','latex');  
ylabel('$\mathcal{G}(x)$ [dBm]', 'FontSize',  12,'Interpreter','latex'); 

subplot(2,1,2);
plot(xxl,n,'-*','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w')
grid on;
ax = gca;
ax.YAxis.FontSize = 12 %for y-axis 
ay = gca;
ay.XAxis.FontSize = 12 %for y-axis
xlabel('Incident RF Power [dBm], $\mathrm{P}_\mathrm{RF}$','FontSize',  16,'Interpreter','latex');  
ylabel('$\eta$ [$\%$]', 'FontSize',  16,'Interpreter','latex'); 


