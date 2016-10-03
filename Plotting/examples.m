%examples

%% Hit Rate from one subject
ss=11;
H=out.nHitsByPhase(ss,:);
M=out.nMissByPhase(ss,:);

xa = linspace(0,2*pi,1000);
x = cos(xa);

% Hits
figure(4); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[600 400],'paperposition',[0 0 600 400])
set(gcf,'position',[100,100,600,400])


a2 = axes('units','points','position',[0.12*600 0.15*400 0.8*600 0.15*400]);
axes(a2)
plot(xa./pi*180,x,'k','linewidth',4)
axis tight;
set(gca,'ytick',[],'ycolor','w','fontsize',24,'box','off','lineWidth',2)
set(gca,'xtick',[36:72:360])
xlabel(' Stimulation Phase (deg)')
xlim([-1 361])
grid on

a1 = axes('units','points','position',[0.12*600 0.32*400 0.8*600  0.65*400]);
axes(a1); hold on;
x1=21:72:360;
x2=51:72:360;
b1=bar(x1,H,0.4);
b1.FaceColor = [255 180 150]/255;
b1.EdgeColor = 'none';
b2=bar(x2,M,0.4);
b2.FaceColor = [150 220 220]/255;
b2.EdgeColor = 'none';
set(gca,'fontsize',24,'box','off','lineWidth',2)
set(gca,'xtick',[36:72:360],'xTickLabel','')
ylim([0 60])
xlim([-1 361])
set(gca,'ytick',0:20:60)
ylabel(' Trial Counts ' )
% plot([36:72:360],X','-','color',[180 180 180]/255)
% plot([36:72:360],mean(X), 'color',[100 100 100]/255,'linewidth',5)
% set(gca,'fontsize',16,'box','off','lineWidth',2)
% set(gca,'xtick',[36:72:360],'xTickLabel','')
% xlim([-1 361])
% ylim([0 1])
% ylabel(' Hit Rate ' )
 print(gcf, '-dpdf', ['../plots/tacs_enc_xdiva/Predictions/ExampleSubjsHMbyPhase']);
%% HC hit rate
ss=11;
H=out.nHitsByPhase(ss,:);
M=out.nMissByPhase(ss,:);

xa = linspace(0,2*pi,1000);
x = cos(xa);

% Hits
figure(4); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[600 400],'paperposition',[0 0 600 400])
set(gcf,'position',[100,100,600,400])


a2 = axes('units','points','position',[0.12*600 0.15*400 0.8*600 0.15*400]);
axes(a2)
plot(xa./pi*180,x,'k','linewidth',4)
axis tight;
set(gca,'ytick',[],'ycolor','w','fontsize',24,'box','off','lineWidth',2)
set(gca,'xtick',[36:72:360])
xlabel(' Stimulation Phase (deg)')
xlim([-1 361])
grid on

a1 = axes('units','points','position',[0.12*600 0.32*400 0.8*600  0.65*400]);
axes(a1); hold on;
x1=21:72:360;
x2=51:72:360;
b1=bar(x1,H,0.4);
b1.FaceColor = [255 180 150]/255;
b1.EdgeColor = 'none';
b2=bar(x2,M,0.4);
b2.FaceColor = [150 220 220]/255;
b2.EdgeColor = 'none';
set(gca,'fontsize',24,'box','off','lineWidth',2)
set(gca,'xtick',[36:72:360],'xTickLabel','')
ylim([0 60])
xlim([-1 361])
set(gca,'ytick',0:20:60)
ylabel(' Trial Counts ' )
% plot([36:72:360],X','-','color',[180 180 180]/255)
% plot([36:72:360],mean(X), 'color',[100 100 100]/255,'linewidth',5)
% set(gca,'fontsize',16,'box','off','lineWidth',2)
% set(gca,'xtick',[36:72:360],'xTickLabel','')
% xlim([-1 361])
% ylim([0 1])
% ylabel(' Hit Rate ' )
 print(gcf, '-dpdf', ['../plots/tacs_enc_xdiva/Predictions/ExampleSubjsHMbyPhase']);


%% Legend
colors = [255 180 150;150 220 220]/255;
txtStr = {'Rem','Forg'};
makeLegend(colors,txtStr);
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/Predictions/RemForgLegend']);

%% HitRate
X = H./(H+M);

xa = linspace(0,2*pi,1000);
x = cos(xa);

% Hits
figure(4); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[600 400],'paperposition',[0 0 600 400])
set(gcf,'position',[100,100,600,400])

a2 = axes('units','points','position',[0.15*600 0.15*400 0.8*600 0.15*400]);
axes(a2)
plot(xa./pi*180,x,'k','linewidth',4)
axis tight;
set(gca,'ytick',[],'ycolor','w','fontsize',24,'box','off','lineWidth',2)
set(gca,'xtick',[36:72:360])
xlabel(' Stimulation Phase (deg)')
xlim([-1 361])
grid on

a1 = axes('units','points','position',[0.15*600 0.32*400 0.8*600  0.65*400]);
axes(a1); hold on;
x1=36:72:360;
b1=bar(x1,X,0.8);
b1.FaceColor = [150 150 150]/255;
b1.EdgeColor = 'none';
set(gca,'fontsize',24,'box','off','lineWidth',2)
set(gca,'xtick',[36:72:360],'xTickLabel','')
ylim([0.4 0.8])
xlim([-1 361])
set(gca,'ytick',0:0.2:1)
ylabel(' Trial Counts ' )
ylabel(' Recognition Rate ' )
% 
 print(gcf, '-dpdf', ['../plots/tacs_enc_xdiva/Predictions/ExampleSubjsHitRatebyPhase']);

 %% Polar Plot 1
exNum=1
if exNum==1
dx = pi/2; 
 
th1 = pi/3;
th2 = th1+dx;

rho1 = 1;
rho2 = 0.5;
elseif exNum==2
dx = 0.3; 
 
th1 = -pi/5;
th2 = th1+dx;

rho1 = 0.75;
rho2 = 0.75; 
end

opts = [];
opts.meanVecs   = 0;
opts.markerSize = [200 200];
opts.colors     = [255 180 150; 150 220 220]/255;
opts.maxR       = 4/3;
opts.alpha      = 0.8;

opts.meanVecs   = 1;
opts.alpha      = 1;
opts.magText     = 0;
%opts.scatterDots = 0;
opts.connect      = 1;
opts.connectLinWidth=4;
% raw vectors
%nn = numel(th1_r);
han = PolarPlot([[th1;th1],[th2;th2]],[[rho1;rho1],[rho2;rho2]],opts);
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[600 600],'paperposition',[0 0 600 600])
set(gcf,'position',[100,100,400,400])
set(gca,'units','points','position',[50 50 300 300]); 
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Predictions/ExamplePolarDist' num2str(exNum)]) ....

 
%% single HitRate Vector
%X = H./(H+M);
%Z = mean(X.*exp(1j*out.StimPhases));
X=squeeze(out.HM_Conf_MeVects(11,1,3,:));
th = X(1);
rho = X(2);

opts = [];
opts.meanVecs   = 0;
opts.markerSize = [250];
opts.colors     = [255 180 150]/255;
%opts.maxR       = 4/3;
opts.magText     =0;
opts.alpha      = 0.8;

opts.meanVecs   = 1;
opts.alpha      = 1;
han = PolarPlot(th,rho,opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Predictions/SingleSubjsHR_Polar']) ....
%     num2str(round(dx/pi*180)) 'kp1_' num2str(kp1) 'kp2_' num2str(kp2) ]);

