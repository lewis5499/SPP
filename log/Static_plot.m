%created by hz Liu, Dec.18th, 2023
figure;
set(gcf,'Position',[50 50 900 600])
box on
grid on 
hold on;
% %LS
%plot(Static_LS.EAST,Static_LS.mNORTH,LineWidth=1);
%EKF
plot(Static_256_36_01mm_16.EASTm,Static_256_36_01mm_16.NORTHm,LineWidth=1)
plot(StaticKF.EASTm,StaticKF.NORTHm,LineWidth=1)
%scatter(StaticKF.EASTm,StaticKF.NORTHm,'filled')
legend('$\bf{EKF-2.56-36.0-0.1mm-16.0}$','$\bf{EKF-6Sats}$')

set(legend,'LineWidth',1,'Interpreter','latex','FontSize',10);
xlabel('$\bf{East(m)}$','interpreter','latex','FontSize', 16)
ylabel('$\bf{North(m)}$','interpreter','latex','FontSize', 16) 
set(gca,'linewidth',1.2,'fontsize',12,'fontname','Times','FontWeight','bold')
title({'$\bf{NovAtelReceiver-Trajectory-Static(2D)}$'}, 'interpreter','latex','FontSize', 14);

% subPlot3(Static_LS.SOW,Static_LS.EAST,Static_LS.mNORTH,Static_LS.mUPm,Static_256_36_01mm_16.SOW, ...
%     Static_256_36_01mm_16.EASTm,Static_256_36_01mm_16.NORTHm,Static_256_36_01mm_16.UPm)

function []=subPlot3(t,value1,value2,value3,tkf,value1kf,value2kf,value3kf)
figure; 
set(gcf,'Position',[0 0 1000 1500])

subplot(3,1,1)
plot(t,value1,'Color',[0.62 0.49 0.31],LineWidth=1);hold on
plot(tkf,value1kf,'Color',[0.31 0.62 0.49],LineWidth=1);
set(gca,'linewidth',1.2,'fontsize',14,'fontname','Times','FontWeight','bold')
legend1=legend('$\bf{LS}$','$\bf{EKF}$','interpreter','latex','FontSize',10.5); 
set(legend1,'LineWidth',1,'Interpreter','latex','FontSize',10.5);
ylabel('$\bf{East(m)}$','interpreter','latex','FontSize', 17)  
title({'$\bf{Residual(E-N-U):Least-Squares-VS-Extended-Kalman-Filter}$'}, 'interpreter','latex','FontSize', 18);%
%subtitle('$\bf{Config:initVar=9.0,\sigma_x =1.5,D_\delta=81.0,allSats}$', 'interpreter','latex','FontSize', 14)
box on
grid on 

subplot(3,1,2)
plot(t,value2,'Color',[0.4660 0.6740 0.1880],LineWidth=1);hold on
plot(tkf,value2kf,'Color',[0.6740 0.1880 0.1880],LineWidth=1);
set(gca,'linewidth',1.2,'fontsize',14,'fontname','Times','FontWeight','bold')
legend1=legend('$\bf{LS}$','$\bf{EKF}$','interpreter','latex','FontSize',10.5); 
set(legend1,'LineWidth',1,'Interpreter','latex','FontSize',10.5);
ylabel('$\bf{North(m)}$','interpreter','latex','FontSize', 17)  
box on
grid on 

subplot(3,1,3)
plot(t,value3,'Color',[0.31 0.49 0.62],LineWidth=1);hold on
plot(tkf,value3kf,'Color',[0.8500 0.4250 0.2980],LineWidth=1);
set(gca,'linewidth',1.2,'fontsize',14,'fontname','Times','FontWeight','bold')
legend1=legend('$\bf{LS}$','$\bf{EKF}$','interpreter','latex','FontSize',10.5); 
set(legend1,'LineWidth',1,'Interpreter','latex','FontSize',10.5);
ylabel('$\bf{Up(m)}$','interpreter','latex','FontSize', 17)  
xlabel('$\bf{SOW(sec)}$','interpreter','latex','FontSize', 17)
box on
grid on 
hold off
end