%created by hz Liu, Dec.18th, 2023
%Plot2();


t=kinematic.SOW;
value1=kinematic.EASTm;
value2=kinematic.NORTHm;
value3=kinematic.UPm;

tkf=kinematicKF.SOW;
value1kf=kinematicKF.EASTm;
value2kf=kinematicKF.NORTHm;
value3kf=kinematicKF.UPm;

subPlot3(t,value1,value2,value3,tkf,value1kf,value2kf,value3kf);

function []=Plot2()
figure;
set(gcf,'Position',[50 50 900 600])
box on
grid on 
hold on;
% %LS
% plot(kinematic_LS_256.Bdeg,kinematic_LS_256.Ldeg,LineWidth=1);
%EKF
plot(kinematic_256_9_15_81.Bdeg,kinematic_256_9_15_81.Ldeg,LineWidth=1)
plot(kinematic.Bdeg,kinematic.Ldeg,LineWidth=1)
%scatter(kinematic.Bdeg,kinematic.Ldeg,'filled')
%plot(kinematic_225_9_15_81.Bdeg,kinematic_225_9_15_81.Ldeg,LineWidth=1)
%plot(kinematic_225_9_6_16.Bdeg,kinematic_225_9_6_16.Ldeg,'Color',[0.4660 0.6740 0.1880],LineWidth=1);
% plot(kinematic_225_9_6_81.Bdeg,kinematic_225_9_6_81.Ldeg,LineWidth=1);

% legend('$\bf{LS-2.56}$', ...
%     '$\bf{EKF-2.25-9.0-1.5-81.0}$', ...%'$\bf{EKF-2.56-9.0-1.5-81.0}$', ...
%     '$\bf{EKF-2.25-9.0-6.0-16.0}$')% ...'$\bf{EKF-2.25-9.0-6.0-81.0}$');
% 
legend('$\bf{EKF-2.56-9.0-1.5-81.0}$', ...
    '$\bf{EKF-6-Sats}$')
    

set(legend,'LineWidth',1,'Interpreter','latex','FontSize',10);
xlabel('$\bf{lat(deg)}$','interpreter','latex','FontSize', 16)
ylabel('$\bf{lon(deg)}$','interpreter','latex','FontSize', 16) 
set(gca,'linewidth',1.2,'fontsize',12,'fontname','Times','FontWeight','bold')
title({'$\bf{Nav-Trajectory-Kinematic(2D)}$'}, 'interpreter','latex','FontSize', 14);
end

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