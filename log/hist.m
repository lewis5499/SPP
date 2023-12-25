% Y=[std(LS.ECEFXm),std(EKF.ECEFXm);
%    std(LS.ECEFYm),std(EKF.ECEFYm);
%    std(LS.ECEFZm),std(EKF.ECEFZm);
%    rms(LS.EASTm), rms(EKF.EASTm);
%    rms(LS.NORTHm),rms(EKF.NORTHm);
%    rms(LS.UPm),   rms(EKF.UPm)];

Y=[std(LS.ECEFXm(301:end)),std(EKF.ECEFXm(301:end));
   std(LS.ECEFYm(301:end)),std(EKF.ECEFYm(301:end));
   std(LS.ECEFZm(301:end)),std(EKF.ECEFZm(301:end));
   rms(LS.EASTm(301:end)), rms(EKF.EASTm(301:end));
   rms(LS.NORTHm(301:end)),rms(EKF.NORTHm(301:end));
   rms(LS.UPm(301:end)),   rms(EKF.UPm(301:end))];

histogram6(Y);

function [] = histogram6(Y)
figure;
set(gcf,'Position',[50 50 1000 660])
box on
grid on 
hold on;

X = categorical({'stdX','stdY','stdZ','rmsE','rmsN','rmsU'});
X = reordercats(X,{'stdX','stdY','stdZ','rmsE','rmsN','rmsU'});

b=bar(X,Y);
b(1).FaceColor=[0.28 0.57 0.54];
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom','fontsize',10,'fontname','Times','FontWeight','bold');

b(2).FaceColor=[0.73 0.47 0.58];
xtips2 = b(2).XEndPoints;
ytips2 = b(2).YEndPoints;
labels2 = string(b(2).YData);
text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom','fontsize',10,'fontname','Times','FontWeight','bold');

ylim([0,1.5])
legend('$\bf{Single Epoch Least Squares}$','$\bf{Extended Kalman Filter}$','Orientation','vertical','Location','northwest')
set(legend,'LineWidth',1,'Interpreter','latex','FontSize',13);
xlabel('$\bf{Categories}$','interpreter','latex','FontSize', 19)
ylabel('$\bf{Values}$','interpreter','latex','FontSize', 19) 
set(gca,'linewidth',1.2,'fontsize',14,'fontname','Times','FontWeight','bold')
title({'$\bf{SPP-Accuracy-Assessment:LS-vs-EKF(after300s)}$'}, 'interpreter','latex','FontSize', 18);%

end

