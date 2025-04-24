%% Paper: Finite Blocklength Information Theory: What is the
%% Practical Impact on Wireless Communications? 
close all; clear; clc

%% Simulation Parameters
N = 56;             % packet length
K = [2 15];
R = [0.2 0.8];
T = N/2;
SNR_dB = -10:10;
SNR = 10.^(SNR_dB./10);

for k=1:length(K)
    for r=1:length(R)
        for s=1:length(SNR)
            C = @(x1) log2(1+x1);
            V = @(x2) ((log2(exp(1))).^2).*x2.*(x2+2)./(2.*(x2 + 1).^2);
            SNR0 = 2^R(r) - 1;
            Pt1 = @(x6) qfunc(sqrt(N./V(x6)).*(C(x6)-R(r)));
            Pt2 = @(x4) (1/2)*(sqrt(2*(1+K(k))./SNR(s))^2).*exp(-(sqrt(2*K(k))^2 + x4.*sqrt(2*(1+K(k))./SNR(s))^2)/2).*besseli(0, sqrt(2*K(k))*sqrt(2*(1+K(k))./SNR(s))*sqrt(x4));
            PDFav = @(x5) Pt1(x5).*Pt2(x5);
            AvgError(k,r,s) = integral(PDFav,0,50);
            Out(k,r,s) = 1 - marcumq(sqrt(2*K(k)),sqrt(2*SNR0*(1+K(k))/SNR(s)));
        end
    end
end


AvgError_K2R2(1,:) = AvgError(1,1,:);
AvgError_K2R8(1,:) = AvgError(1,2,:);
AvgError_K15R2(1,:) = AvgError(2,1,:);
AvgError_K15R8(1,:) = AvgError(2,2,:);

Out_K2R2(1,:) = Out(1,1,:);
Out_K2R8(1,:) = Out(1,2,:);
Out_K15R2(1,:) = Out(2,1,:);
Out_K15R8(1,:) = Out(2,2,:);

figure(1)
semilogy(SNR_dB,AvgError_K2R2,'-s','Color',[0 0.4470 0.7410],'LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w')
hold on
semilogy(SNR_dB,AvgError_K2R8,'-X','Color',[0 0.4470 0.7410],'LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w')
semilogy(SNR_dB,AvgError_K15R2,'-*','Color',[0 0.4470 0.7410],'LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w')
hold on
semilogy(SNR_dB,AvgError_K15R8,'-o','Color',[0 0.4470 0.7410],'LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w')
semilogy(SNR_dB,Out_K2R2,'--s','Color',[0.8500 0.3250 0.0980],'LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w')
semilogy(SNR_dB,Out_K2R8,'--X','Color',[0.8500 0.3250 0.0980],'LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w')
semilogy(SNR_dB,Out_K15R2,'--*','Color',[0.8500 0.3250 0.0980],'LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w')
semilogy(SNR_dB,Out_K15R8,'--o','Color',[0.8500 0.3250 0.0980],'LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w')
hold off
grid on;
ax = gca;
ax.YAxis.FontSize = 12; %for y-axis 
ay = gca;
ay.XAxis.FontSize = 12; %for y-axis
legend('NA Probability (R=0.2, K=2)','NA Probability (R=0.8, K=2)','NA Probability (R=0.2, K=15)','NA Probability (R=0.8, K=15)','Outage Probability (R=0.2, K=2)','Outage Probability (R=0.8, K=2)','Outage Probability (R=0.2, K=15)','Outage Probability (R=0.8, K=15)','FontSize', 10); 
% legend('NA Probability (R=0.2, K=15)','NA Probability (R=0.8, K=15)','Outage Probability (R=0.2, K=15)','Outage Probability (R=0.8, K=15)','FontSize', 10); 
xlabel('SNR [dB]','FontSize',  16,'Interpreter','latex');  
ylabel('Outage and Non Asymptotic Error Probabilities', 'FontSize',  16,'Interpreter','latex');
% xlim([1 10])
ylim([10^-8 5])