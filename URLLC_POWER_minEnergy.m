close all; clear; clc

rand('state',123456789)
randn('state',123456789)
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
Pcirc = 1.33e-3;     % transceiver power transmission constant
Pmax_dBm = 3.3;     % transceiver max transmission power
Pmax = 10^((Pmax_dBm-30)/10);
Pb_dBm = 45;
Pb = 10.^((Pb_dBm-30)./10);

Rmax = 8;           % Eficiência espectral máxima

% Non-linear Energy Harvesting
% From: Massive Wireless Energy Transfer with Multiple Power Beacons 
% for very large Internet of Things
c0 = 0.2308;             % EH unitless constants
c1 = 5.365;
w = 10.73;              % energy harvesting saturation level

% Rician fading LOS WET phase
kwet_dB = 4;
kwet = 10.^(kwet_dB./10);
% Rician fading LOS WIT phase
kwit_dB = 2;
kwit = 10.^(kwit_dB./10);

% Noise PSD
N0_dB = -204;
N0 = 10.^(N0_dB./10);
% Noise figure
Nf_dB = 10;
Nf = 10.^(Nf_dB./10);

% Packet replication
K = 1:10;
% Total time for transmission
T = 0.1e-3:0.01e-3:1e-3;           

dwit = 100;          % distance Source-Destination
dwet = 3;           % distance PB-Source

%% Baseline framework equations 
% average power gain in the WET link
beta_wet = (c^2)/((4*pi*fwet)^2*(dwet^alpha));

% [Eq. 9] average power gain in the WIT link
beta_wit = (c^2)/(((4*pi*fwit)^2)*(dwit^alpha)*Nf*N0*B);

%% Inline functions to help with the simulations
% [Eq. 32] Saturation non-linear EH function
syms xx
gf(xx) = 1e-3.*w.*(1 - exp(-c0*xx*1e3))./(1 + exp(-c0.*(xx.*1e3 - c1)));
ginv = finverse(gf);

% limites de integração
L1 = double(ginv(Pcirc));
L2 = double(ginv(Pmax/eta + Pcirc));
Linf = double(ginv(0.99999999*w*1e-3));

g = @(x4) 1e-3.*w.*(1 - exp(-c0*x4*1e3))./(1 + exp(-c0.*(x4.*1e3 - c1)));

ang = 2*pi*rand(1, EH);
rand('seed',1)

%% Start simulation - Minization of the Energy Consumption
M = [4 8];              % Number of PB antennas

for m=1:length(M)
    disp([num2str(M(m)), ' antenas']);

    % Channels' generation
    % Random deployment
    hlos = sqrt(kwet/(1+kwet))*exp(1i*(-pi)*([0:M(m)-1]')*sin(repmat(ang',1,EH)));
    RR = eye(M(m))*1/(1+kwet);
    hRh = real(hlos'*RR*hlos);                                %LOS channel component
    hnlos = sqrt(1/(1+kwet))*(randn(M(m),EH)+1i*randn(M(m),EH));     %Instantaneous nlos %%---dimension changes

    %channel realizations
    h = hlos + hnlos;
    
    Pfcsi = max(beta_wet.*Pb.*(norm(h)^2), L2);

    if g(Pfcsi) - Pcirc < 0
        error('Pb menor que o mínimo')
    end

    outFCSI = ones(length(K), length(T), length(M));
    outACSI = ones(length(K), length(T), length(M));
    EF = inf*ones(length(K), length(T), length(M));
    EA = inf*ones(length(K), length(T), length(M));

    for i = 1:length(K)
        for j = 1:length(T)
            % [Eq.21] Parameter of CDF WET
            a = sqrt(2/hRh)*(norm(hlos))^2;
            % [Eq.22] Parameter of CDF WET
            b = sqrt(2/(hRh*beta_wet*Pb))*norm(hlos);

            % [Eq.23] PDF of the energy at the EH terminal
            fe = @(x3) (1/2)*(b^2).*exp(-(a^2 + x3.*b^2)/2).*besseli(0, a*b*sqrt(x3));
            % [Eq.19] CDF of the energy at the EH terminal
            Fe = @(x3) 1 - marcumq(a,b*sqrt(x3)); %%---->this is needed

            % Transmit rate
            R(i,j) = (N*K(i))/(B*T(j));
            % Variable for outage
            q(i,j) = (2^R(i,j) - 1)./beta_wit;

            %% Average CSI
            % [Eq. 30] CDF WIT - maximum ratio combining
            Fz_mrc = @(x2) (1 - marcumq(sqrt(2*K(i)*kwit), sqrt(2*(1+kwit)*x2), K(i)));
            % maximum ratio combining
            out_mrc = @(x1) fe(x1).*Fz_mrc(q(i,j)./(eta*(g(x1) - Pcirc)));
            out_mrc2 = @(x1) fe(x1).*Fz_mrc(q(i,j)./(Pmax));

            % [Eq. 32] outage of A-CSI
            outACSI(i,j,m) = integral(out_mrc,L1,L2) + integral(out_mrc2,L2,Linf);

            % Find the values that guarantee constraint
            if((outACSI(i,j,m) <= pc) && (R(i,j) <= Rmax))
                EA(i,j,m) = Pb*T(j);
            end
            
            %% FULL CSI
            % [Eq. 25] power harvested F-CSI
            Xfcsi(i,j) = q(i,j)/(eta*(g(Pfcsi) - Pcirc));

            % [Eq. 33] maximum ratio combining
            outFCSI(i,j,m) = (1 - marcumq(sqrt(2*K(i)*kwit), sqrt(2*(1+kwit)*Xfcsi(i,j)), K(i)));

            if((outFCSI(i,j,m) <= pc) && (R(i,j) <= Rmax))
                EF(i,j,m) = Pb*T(j);
            end
        end
    end

    % Min Energy vs K
    [minEF_K(:,m), pos] = min(EF(:,:,m),[],2);
    TminEF(:,m) = T(pos);
    [minEA_K(:,m), pos] = min(EA(:,:,m),[],2);
    TminEA(:,m) = T(pos);

    % Min Energy vs T
    [minEF_T(:,m), pos] = min(EF(:,:,m),[],1);
    KminEF(:,m) = K(pos);
    [minEA_T(:,m), pos] = min(EA(:,:,m),[],1);
    KminEA(:,m) = K(pos);
end
TminEF(minEF_K==inf) = inf;
TminEA(minEA_K==inf) = inf;
KminEF(minEF_T==inf) = inf;
KminEA(minEA_T==inf) = inf;

%% ---- Figure

figure(1)
plot(K,TminEF(:,1)*1000,'-s','Color',[0 0.4470 0.7410],'LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w')
hold on
plot(K,TminEF(:,2)*1000,'-X','Color',[0 0.4470 0.7410],'LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w')
plot(K,TminEA(:,1)*1000,'--s','Color',[0.8500 0.3250 0.0980],'LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w')
plot(K,TminEA(:,2)*1000,'--X','Color',[0.8500 0.3250 0.0980],'LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w')
hold off
grid on;
ax = gca;
ax.YAxis.FontSize = 12; %for y-axis 
ay = gca;
ay.XAxis.FontSize = 12; %for y-axis
legend('F-CSI (M=4)','F-CSI (M=8)','A-CSI (M=4)','A-CSI (M=8)','FontSize', 10); 
xlabel('Number of Replicated Packets, $K$','FontSize',  16,'Interpreter','latex');  
ylabel('Optimal Time [ms], $T_i^{\star}$', 'FontSize',  16,'Interpreter','latex');
xlim([1 10])
ylim([0.05 0.6])
