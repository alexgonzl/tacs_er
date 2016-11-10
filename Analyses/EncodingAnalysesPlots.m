% Analyses for Encoding Data.

close all
clearvars

dataPath    = '~/Google Drive/Research/tACS/tACS_ER_task/data/tacs_enc_xdiva/';
load([dataPath 'Summary/BehavSummary.mat'])
load([dataPath 'Summary/PhaseDependentAnalyses.mat'])

uniqSeqSelection = 2;
if uniqSeqSelection==1
    subjs  = find(out.SubjWithUniqueSeq);
    SubjSelectStr = 'SSuniqueSeq';
elseif uniqSeqSelection==2
    subjs1  = find(out.SubjWithUniqueSeq);
    subjs2  = find(behav_out.encSummary.goodSubj);
    subjs   = intersect(subjs1,subjs2);
    SubjSelectStr = 'SSuniqueSeq_ValidEnc';
else
    subjs = 1:numel(behav_out.retSubj);
    SubjSelectStr = 'all';
end
nSubjs  = numel(subjs);
PhasesDeg = 36:72:359;
PhasesRad = PhasesDeg./180*pi;

%% Encoding Results: Categorization Task
clearvars -except out behav_out subjs nSubjs dataPath SubjSelectStr PhasesDeg
rng(1); % for location reproducibility
HR = [ behav_out.encSummary.meanAcc(subjs) behav_out.encSummary.FaceHR(subjs) behav_out.encSummary.SceneHR(subjs)]*100;
Dstrs = {'ACC','Face','Scn'};
disp(array2table(mean(HR),'variablenames',Dstrs))
[~,p,~,t]=ttest(HR(:,2),HR(:,3));
disp(table(t.tstat,p,'rownames',{'Face vs Scn ACC'},'variablenames',{'T','P'}))

% Figure : (a) ACC
figure(1); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[600 400],'paperposition',[0 0 600 400])
set(gcf,'position',[200,200,600,400])

a = axes('units','points','position',[100 100 80 200]); hold on;
x = randn(nSubjs,1)*0.05+0.5;
y = HR(:,1);
plot([0.2 0.8], ones(1,2)*mean(y),'linewidth',4,'color','k')
s=scatter(x,y);
s.MarkerFaceAlpha=0.5;
s.MarkerEdgeAlpha=0;
s.SizeData          = 120;
s.MarkerEdgeColor = [119,136,153]/255;
s.MarkerFaceColor = [119,136,153]/255;

set(gca,'fontsize',20,'xtick',[])
ylim([90 100])
xlim([0 1])
set(gca,'ytick',[90 100])
ylabel(' ACC (%) ')
set(gca,'LineWidth',2)

% Figure 1: (b) acc by category
a = axes('units','points','position',[300 100 180 200]); hold on;
y = HR(:,2:3);
%
for ii =1:nSubjs
    plot([x(ii) x(ii)+1], y(ii,:),'-','color',[0.6 0.6 0.6])
end

plot([0.2 0.8], ones(1,2)*mean(y(:,1)),'linewidth',4,'color','k')
plot([0.2 0.8]+1, ones(1,2)*mean(y(:,2)),'linewidth',4,'color','k')

% Faces
s=scatter(x,y(:,1));
s.MarkerFaceAlpha=0.5;
s.MarkerEdgeAlpha=0;
s.SizeData          = 120;
s.MarkerEdgeColor = [100 200 100]/255;
s.MarkerFaceColor = [100 200 100]/255;

% Scenes
s=scatter(x+1,y(:,2));
s.MarkerFaceAlpha   = 0.5;
s.MarkerEdgeAlpha   = 0;
s.SizeData          = 120;
s.MarkerEdgeColor   = [200 100 200]/255;
s.MarkerFaceColor   = [200 100 200]/255;

xlim([-0.1 2.1])
ylim([90 100])
set(gca,'fontsize',20,'xtick',[0.5 1.5],'xticklabel',{'F','S'})
set(gca,'ytick',[90 100])
set(gca,'LineWidth',2)
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/Encoding/' 'CategorizationPerf_' SubjSelectStr])

%%  Categorization by Phase
clearvars -except out behav_out subjs nSubjs dataPath SubjSelectStr PhasesDeg
rng(1); % for location reproducibility
HR = behav_out.encSummary.HRByPhase(subjs,:,:)*100;

% Figure : ACC by Phase
figure(1); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[600 400],'paperposition',[0 0 600 400])
set(gcf,'position',[200,200,600,400])
a = axes('units','points','position',[100 150 400 200]); hold on;
x = randn(nSubjs,1)*0.05+0.5;
y = squeeze(HR(:,1,:));
nBars = 5;
for ii =1:nBars
    xx = x+(ii-1);
    yy = y(:,ii);
    plot([0.2 0.8]+(ii-1), ones(1,2)*mean(yy),'linewidth',4,'color','k')
    s = scatter(xx,yy);
    s.MarkerFaceAlpha=0.5;
    s.MarkerEdgeAlpha=0.4;
    s.SizeData          = 120;
    s.MarkerEdgeColor = [119,136,153]/255;
    s.MarkerFaceColor = [119,136,153]/255;
end
set(gca,'fontsize',20,'xtick',[0.5:5],'xticklabel',[])%,'xtick',[0.5:5],'xticklabel',PhasesDeg)
ylim([85 100])
xlim([0 5])
set(gca,'ytick',[90 100])
ylabel(' ACC (%) ')
set(gca,'LineWidth',2)


a2 = axes('units','points','position',[100 80 400 50]); hold on;
xa = linspace(0,2*pi,1000); x = cos(xa);
axes(a2)
plot(xa./pi*180,x,'k','linewidth',4)
axis tight
set(gca,'ytick',[],'ycolor','w','fontsize',20,'box','off','lineWidth',2)
set(gca,'xtick',PhasesDeg)
xlabel(' Encoding Phase (deg)')
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/Encoding/' 'CategorizationPerfByPhase_' SubjSelectStr])

thRad = (PhasesDeg-36)/180*pi;
Z = repmat(exp(1j*thRad),nSubjs,1);
yZ = mean(y.*Z,2);
th = angle(yZ);

opts = [];
opts.colors = [119,136,153]/255;
opts.maxR = 4/3;
opts.alpha = 0.8;
opts.markerSize=300;
opts.magText    = '';
PolarPlot(th,[],opts)
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/Encoding/' 'CategorizationPerfByPhasePolar_' SubjSelectStr])

[p,r] = circ_rtest(th);
fprintf('Rayleight Test for modulation across subjects:\n')
fprintf('p=%g ; r = %g \n',p,r)

%% Categorization By Phase / Stim Type
clearvars -except out behav_out subjs nSubjs dataPath SubjSelectStr PhasesDeg
rng(2); % for location reproducibility
HR = behav_out.encSummary.HRByPhase(subjs,:,:)*100;

% Figure : ACC by Phase
figure(1); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[600 400],'paperposition',[0 0 600 400])
set(gcf,'position',[200,200,600,400])
a = axes('units','points','position',[100 150 400 200]); hold on;
x1 = randn(nSubjs,1)*0.02+0.3;
x2 = randn(nSubjs,1)*0.02+0.7;
y1 = squeeze(HR(:,2,:));
y2 = squeeze(HR(:,3,:));
nBars = 5;
for ii =1:nBars
    xx1 = x1+(ii-1);
    xx2 = x2+(ii-1);
    yy1 = y1(:,ii);
    yy2 = y2(:,ii);
    
    for ss =1:nSubjs
        plot([xx1(ss) xx2(ss)], [yy1(ss) yy2(ss)],'-','color',[0.6 0.6 0.6])
    end
    
    % Faces
    plot([0.15 0.45]+(ii-1), ones(1,2)*mean(yy1),'linewidth',4,'color','k')
    s = scatter(xx1,yy1);
    s.MarkerFaceAlpha=0.5;
    s.MarkerEdgeAlpha=0.4;
    s.SizeData          = 120;
    s.MarkerEdgeColor = [100 200 100]/255;
    s.MarkerFaceColor = [100 200 100]/255;
    
    % Scenes
    plot([0.55 0.85]+(ii-1), ones(1,2)*mean(yy2),'linewidth',4,'color','k')
    s = scatter(xx2,yy2);
    s.MarkerFaceAlpha=0.5;
    s.MarkerEdgeAlpha=0;
    s.SizeData          = 120;
    s.MarkerEdgeColor = [200 100 200]/255;
    s.MarkerFaceColor = [200 100 200]/255;
end
set(gca,'fontsize',20,'xtick',[0.5:5],'xticklabel',[])%,'xtick',[0.5:5],'xticklabel',PhasesDeg)
ylim([85 100])
xlim([0 5])
set(gca,'ytick',[90 100])
ylabel(' ACC (%) ')
set(gca,'LineWidth',2)

a2 = axes('units','points','position',[100 80 400 50]); hold on;
xa = linspace(0,2*pi,1000); x = cos(xa);
axes(a2)
plot(xa./pi*180,x,'k','linewidth',4)
axis tight
set(gca,'ytick',[],'ycolor','w','fontsize',20,'box','off','lineWidth',2)
set(gca,'xtick',PhasesDeg)
xlabel(' Encoding Phase (deg)')

print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/Encoding/' 'CategorizationPerfByPhaseStimType_' SubjSelectStr])

% Polar Plots
thRad = (PhasesDeg-36)/180*pi;
Z = repmat(exp(1j*thRad),nSubjs,1);
y1Z = mean(y1.*Z,2);
th1 = angle(y1Z);
y2Z = mean(y2.*Z,2);
th2 = angle(y2Z);

opts = [];
opts.colors = [100 200 100; 200 100 200]/255;
opts.maxR = 4/3;
opts.alpha = 0.8;
opts.markerSize=250;
opts.magText    = '';
PolarPlot([th1 th2],[],opts)
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/Encoding/' 'CategorizationPerfByPhaseByCatPolar_' SubjSelectStr])

[p,r] = circ_rtest(th1);
fprintf('Rayleight Test for modulation across subjects for Faces:\n')
fprintf('p=%g ; r = %g \n',p,r)
[p,r] = circ_rtest(th2);
fprintf('Rayleight Test for modulation across subjects for Scenes:\n')
fprintf('p=%g ; r = %g \n',p,r)
%% make legend
colors = [100 200 100;200 100 200]/255;
txtStr = {'Faces','Scenes'};
makeLegend(colors,[200 200],txtStr);
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/Encoding/FaceScnLegend']);

%% Categorization RTs:
clearvars -except out behav_out subjs nSubjs dataPath SubjSelectStr PhasesDeg
rng(1); % for location reproducibility
RTs = [behav_out.encSummary.medianRTs(subjs)  ...
    behav_out.encSummary.medianFaceRTs(subjs) behav_out.encSummary.medianSceneRTs(subjs)];
Dstrs = {'RTs','Face','Scn'};
disp(array2table(mean(RTs),'variablenames',Dstrs))
[~,p,~,t]=ttest(RTs(:,2),RTs(:,3));
disp(table(t.tstat,p,'rownames',{'Face vs Scn RTs'},'variablenames',{'T','P'}))
fprintf('\n Sign Test \n')
fprintf('\n p=%g\n',signtest(RTs(:,2),RTs(:,3)))
% Figure : (a) ACC
figure(1); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[600 400],'paperposition',[0 0 600 400])
set(gcf,'position',[200,200,600,400])

a = axes('units','points','position',[100 100 80 200]); hold on;
x = randn(nSubjs,1)*0.05+0.5;
y = RTs(:,1);
plot([0.2 0.8], ones(1,2)*mean(y),'linewidth',4,'color','k')
s=scatter(x,y);
s.MarkerFaceAlpha=0.5;
s.MarkerEdgeAlpha=0;
s.SizeData          = 120;
s.MarkerEdgeColor = [119,136,153]/255;
s.MarkerFaceColor = [119,136,153]/255;

set(gca,'fontsize',20,'xtick',[])
ylim([0.5 1])
xlim([0 1])
set(gca,'ytick',[0.5:0.2:1])
ylabel(' RTs (s) ')
set(gca,'LineWidth',2)

%
% Figure 1: (b) acc by category
a = axes('units','points','position',[300 100 180 200]); hold on;
y = RTs(:,2:3);

for ii =1:nSubjs
    if y(ii,1)>y(ii,2)
        plot([x(ii) x(ii)+1], y(ii,:),'-','color',[0.8 0.2 0.2])
    elseif y(ii,1)<=y(ii,2)
        plot([x(ii) x(ii)+1], y(ii,:),'-','color',[0.6 0.6 0.6])
    end
end

plot([0.2 0.8], ones(1,2)*mean(y(:,1)),'linewidth',4,'color','k')
plot([0.2 0.8]+1, ones(1,2)*mean(y(:,2)),'linewidth',4,'color','k')

% Faces
s=scatter(x,y(:,1));
s.MarkerFaceAlpha=0.5;
s.MarkerEdgeAlpha=0.4;
s.SizeData          = 120;
s.MarkerEdgeColor = [100 200 100]/255;
s.MarkerFaceColor = [100 200 100]/255;

% Scenes
s=scatter(x+1,y(:,2));
s.MarkerFaceAlpha   = 0.5;
s.MarkerEdgeAlpha   = 0;
s.SizeData          = 120;
s.MarkerEdgeColor   = [200 100 200]/255;
s.MarkerFaceColor   = [200 100 200]/255;

xlim([-0.1 2.1])
ylim([0.5 1])
set(gca,'fontsize',20,'xtick',[0.5 1.5],'xticklabel',{'F','S'})
set(gca,'ytick',[0.5:0.2:1])
set(gca,'LineWidth',2)
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/Encoding/' 'CategorizationRTs_' SubjSelectStr])

%% Categorization RTs by Phase
clearvars -except out behav_out subjs nSubjs dataPath SubjSelectStr PhasesDeg
rng(1); % for location reproducibility
RTs = squeeze(behav_out.encSummary.medianRTsByPhase(subjs,1,:));

% Figure : ACC by Phase
figure(1); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[600 400],'paperposition',[0 0 600 400])
set(gcf,'position',[200,200,600,400])
a = axes('units','points','position',[100 150 400 200]); hold on;
x = randn(nSubjs,1)*0.05+0.5;
y = RTs;
nBars = 5;
for ii =1:nBars
    xx = x+(ii-1);
    yy = y(:,ii);
    plot([0.2 0.8]+(ii-1), ones(1,2)*mean(yy),'linewidth',4,'color','k')
    s = scatter(xx,yy);
    s.MarkerFaceAlpha=0.5;
    s.MarkerEdgeAlpha=0;
    s.SizeData          = 120;
    s.MarkerEdgeColor = [119,136,153]/255;
    s.MarkerFaceColor = [119,136,153]/255;
end
set(gca,'fontsize',20,'xtick',[0.5:5],'xticklabel',[])%,'xtick',[0.5:5],'xticklabel',PhasesDeg)
ylim([0.6 1])
xlim([0 5])
set(gca,'ytick',[0.5:0.2:1])
ylabel(' RTs (s) ')
set(gca,'LineWidth',2)

a2 = axes('units','points','position',[100 80 400 50]); hold on;
xa = linspace(0,2*pi,1000); x = cos(xa);
axes(a2)
plot(xa./pi*180,x,'k','linewidth',4)
axis tight
set(gca,'ytick',[],'ycolor','w','fontsize',20,'box','off','lineWidth',2)
set(gca,'xtick',PhasesDeg)
xlabel(' Encoding Phase (deg)')

print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/Encoding/' 'CategorizationRTsByPhase_' SubjSelectStr])

% compute modulation by phase
thRad = (PhasesDeg-36)/180*pi;
Z = repmat(exp(1j*thRad),nSubjs,1);
y = -log(RTs);
yZ = mean(y.*Z,2);
th = angle(yZ);

opts = [];
opts.colors = [119,136,153]/255;
opts.maxR = 4/3;
opts.alpha = 0.8;
opts.markerSize=300;
opts.magText    = '';
PolarPlot(th,[],opts)
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/Encoding/' 'CategorizationRTsByPhasePolar_' SubjSelectStr])

[p,r] = circ_rtest(th);
th_hat = mod(angle(mean(exp(1j*th))),2*pi)/pi*180;
fprintf('Rayleight Test for Categorization RT modulation across subjects:\n')
fprintf('p=%g ; r = %g , angle = %g \n',p,r,th_hat)

%% Categorization RTs by Phase and category
clearvars -except out behav_out subjs nSubjs dataPath SubjSelectStr PhasesDeg
rng(2); % for location reproducibility
RTs = behav_out.encSummary.medianRTsByPhase(subjs,:,:);

% Figure : ACC by Phase
figure(1); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[600 400],'paperposition',[0 0 600 400])
set(gcf,'position',[200,200,600,400])
a = axes('units','points','position',[100 150 400 200]); hold on;
x1 = randn(nSubjs,1)*0.02+0.3;
x2 = randn(nSubjs,1)*0.02+0.7;
y1 = squeeze(RTs(:,2,:));
y2 = squeeze(RTs(:,3,:));
nBars = 5;
for ii =1:nBars
    xx1 = x1+(ii-1);
    xx2 = x2+(ii-1);
    yy1 = y1(:,ii);
    yy2 = y2(:,ii);
    
    for jj =1:nSubjs
        if yy1(jj)>yy2(jj)
            plot([xx1(ii) xx2(ii)], [yy1(jj) yy2(jj)],'-','color',[0.8 0.2 0.2])
        else
            plot([xx1(ii) xx2(ii)], [yy1(jj) yy2(jj)],'-','color',[0.6 0.6 0.6])
        end
    end
    
    % Faces
    plot([0.15 0.45]+(ii-1), ones(1,2)*mean(yy1),'linewidth',4,'color','k')
    s = scatter(xx1,yy1);
    s.MarkerFaceAlpha=0.5;
    s.MarkerEdgeAlpha=0;
    s.SizeData          = 120;
    s.MarkerEdgeColor = [100 200 100]/255;
    s.MarkerFaceColor = [100 200 100]/255;
    
    % Scenes
    plot([0.55 0.85]+(ii-1), ones(1,2)*mean(yy2),'linewidth',4,'color','k')
    s = scatter(xx2,yy2);
    s.MarkerFaceAlpha=0.5;
    s.MarkerEdgeAlpha=0.4;
    s.SizeData          = 120;
    s.MarkerEdgeColor = [200 100 200]/255;
    s.MarkerFaceColor = [200 100 200]/255;
end
set(gca,'fontsize',20,'xtick',[0.5:5],'xticklabel',[])%,'xtick',[0.5:5],'xticklabel',PhasesDeg)
ylim([0.6 1])
xlim([0 5])
set(gca,'ytick',[0.5:0.2:1])
ylabel(' RTs (s) ')
set(gca,'LineWidth',2)

a2 = axes('units','points','position',[100 80 400 50]); hold on;
xa = linspace(0,2*pi,1000); x = cos(xa);
axes(a2)
plot(xa./pi*180,x,'k','linewidth',4)
axis tight
set(gca,'ytick',[],'ycolor','w','fontsize',20,'box','off','lineWidth',2)
set(gca,'xtick',PhasesDeg)
xlabel(' Encoding Phase (deg)')

print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/Encoding/' 'CategorizationRTsByPhaseStimType_' SubjSelectStr])

% compute modulation by phase
thRad = (PhasesDeg-36)/180*pi;
Z = repmat(exp(1j*thRad),nSubjs,1);
y1 = -log(y1); y2 = -log(y2);
y1Z = mean(y1.*Z,2); y2Z = mean(y2.*Z,2);
th1 = angle(y1Z); th2 = angle(y2Z);

opts = [];
opts.colors = [100 200 100;200 100 200 ]/255;
opts.maxR = 4/3;
opts.alpha = 0.8;
opts.markerSize= 250;
opts.magText    = '';
PolarPlot([th1 th2],[],opts)
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/Encoding/' 'CategorizationRTsByPhaseByCatPolar_' SubjSelectStr])

[p,r] = circ_rtest(th1);
th1_hat = mod(angle(mean(exp(1j*th1))),2*pi)/pi*180;
fprintf('Rayleight Test for RT modulation across subjects for Faces:\n')
fprintf('p=%g ; r = %g , angle = %g \n',p,r,th1_hat)
[p,r] = circ_rtest(th2);
th2_hat = mod(angle(mean(exp(1j*th2))),2*pi)/pi*180;
fprintf('Rayleight Test for RT modulation across subjects for Scenes:\n')
fprintf('p=%g ; r = %g , angle = %g \n',p,r,th2_hat)
% print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/Encoding/' 'CategorizationRTsByPhasePolar_' SubjSelectStr])
%
% [p,r] = circ_rtest(th);
% th_hat = mod(angle(mean(exp(1j*th))),2*pi)/pi*180;
% fprintf('Rayleight Test for Categorization RT modulation across subjects:\n')
% fprintf('p=%g ; r = %g , angle = %g \n',p,r,th_hat)
