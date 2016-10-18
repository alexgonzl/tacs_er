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
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/Behavior/' 'CategorizationPerf_' SubjSelectStr])

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
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/Behavior/' 'CategorizationPerfByPhase_' SubjSelectStr])

thRad = (PhasesDeg-36)/180*pi;
Z = repmat(exp(1j*thRad),nSubjs,1);
yZ = mean(y.*Z,2);
th = angle(yZ);

opts = [];
opts.colors = [119,136,153]/255;
opts.maxR = 4/3;
opts.alpha = 0.8;
opts.markerSize=300;
PolarPlot(th,[],opts)
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/Behavior/' 'CategorizationPerfByPhasePolar_' SubjSelectStr])

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

print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/Behavior/' 'CategorizationPerfByPhaseStimType_' SubjSelectStr])

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
PolarPlot([th1 th2],[],opts)
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/Behavior/' 'CategorizationPerfByPhaseByCatPolar_' SubjSelectStr])

[p,r] = circ_rtest(th1);
fprintf('Rayleight Test for modulation across subjects for Faces:\n')
fprintf('p=%g ; r = %g \n',p,r)
[p,r] = circ_rtest(th2);
fprintf('Rayleight Test for modulation across subjects for Scenes:\n')
fprintf('p=%g ; r = %g \n',p,r)
%% make legend
colors = [100 200 100;200 100 200]/255;
txtStr = {'Faces','Scenes'};
makeLegend(colors,txtStr);
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/Behavior/FaceScnLegend']);

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
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/Behavior/' 'CategorizationRTs_' SubjSelectStr])

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

print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/Behavior/' 'CategorizationRTsByPhase_' SubjSelectStr])

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

print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/Behavior/' 'CategorizationRTsByPhaseStimType_' SubjSelectStr])

%% dPrimes
clearvars -except out behav_out subjs nSubjs dataPath SubjSelectStr
rng(1); % for location reproducibility

Dstrs   = {'dPrime','Face_dPrime','Scene_dPrime',...
    'dPrime_C','Face_dPrime_C','Scene_dPrime_C'};

Dstrs2  = {'dP', 'Face dP', 'Scn dP', 'C', 'Face C', 'Scn C'};
Dstrs3  = {'dP', 'FaceDP', 'ScnDP', 'C', 'FaceC', 'ScnC'};
nD = numel(Dstrs);
D       = zeros(nSubjs,nD);
for ii = 1:nD
    D(:,ii) = behav_out.retSummary.(Dstrs{ii})(subjs);
end
disp(array2table(mean(D),'variablenames',Dstrs3))
[~,p,~,t]=ttest(D(:,2),D(:,3));
disp(table(t.tstat,p,'rownames',{'Face vs Scn DP'}))

yLims = [-0.1 1.5; 0 2; 0 2; -1.2 1.2; -1.2 1.2; -1.2 1.2];
yTicks1 = [0 1 2];

% Figure 1: (a) d-prime
figure(1); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[600 400],'paperposition',[0 0 600 400])
set(gcf,'position',[200,200,600,400])

a = axes('units','points','position',[100 100 80 200]); hold on;
x = randn(nSubjs,1)*0.05+0.5;
y = D(:,1);
plot([0.2 0.8], ones(1,2)*mean(y),'linewidth',4,'color','k')
s=scatter(x,y);

s.MarkerFaceAlpha=0.5;
s.MarkerEdgeAlpha=0.4;
s.SizeData          = 120;
s.MarkerEdgeColor = [119,136,153]/255;
s.MarkerFaceColor = [119,136,153]/255;

set(gca,'fontsize',20,'xtick',[])
ylim(yLims(1,:))
xlim([0 1])
set(gca,'ytick',yTicks1)
ylabel('dP')
set(gca,'LineWidth',2)

% Figure 1: (b) d-prime by category
a = axes('units','points','position',[300 100 180 200]); hold on;
y = D(:,2:3);

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
s.MarkerEdgeAlpha   = 0.4;
s.SizeData          = 120;
s.MarkerEdgeColor   = [200 100 200]/255;
s.MarkerFaceColor   = [200 100 200]/255;

xlim([-0.1 2.1])
ylim([-0.1 2.1])
set(gca,'fontsize',20,'xtick',[0.5 1.5],'xticklabel',{'F','S'})
set(gca,'ytick',[0 1 2])
ylabel(' dP ')
set(gca,'LineWidth',2)
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/Behavior/' 'dPrimes' SubjSelectStr])

%% Confidence
clearvars -except out behav_out subjs nSubjs dataPath SubjSelectStr

% dPrime by confidence:
DPC = behav_out.retSummary.dPrimeConf(subjs,:);
figure(2); clf;
set(gcf,'paperpositionmode','auto','color','white')
AR = [500 300];
set(gcf,'paperUnits','points','papersize',AR,'paperposition',[0 0 AR])
set(gcf,'position',[50,500,AR])

hold on;
for ss = 1:nSubjs
    p=plot(1:3,DPC(ss,:),'linewidth',1,'color',[0.8 0.8 0.8]);
end
plot(1:3,mean(DPC),'linewidth',4,'color',[0.1 0.1 0.1]);
set(gca,'fontsize',30,'xTick',1:3,'xticklabel',{'Low','Med','High'})
xlim([0.8 3.2])
ylim([-0.5 3.5])
set(gca,'LineWidth',3,'ytick',[0:3])
ylabel(' dP ')

print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/Behavior/'  'Confidence_dPrime' SubjSelectStr])


disp(array2table(mean(DPC),'variablenames',{'Low','Med','High'}))
[~,p1,~,t1]=ttest(DPC(:,2),DPC(:,1));
[~,p2,~,t2]=ttest(DPC(:,3),DPC(:,2));
[~,p3,~,t3]=ttest(DPC(:,3),DPC(:,1));
p = [p1;p2;p3];
t = [t1.tstat;t2.tstat;t3.tstat];
disp(table(t,p,'rownames',{'Mid>Low','Hi>Mid','Hi>Lo'}))

anova1(DPC(:),[ones(nSubjs,1);2*ones(nSubjs,1);3*ones(nSubjs,1)])

%% Confidence Face/Scenes
clearvars -except out behav_out subjs nSubjs dataPath SubjSelectStr
close all
% Faces
DPCF = behav_out.retSummary.Face_dPrimeConf(subjs,:);
figure(3); clf;
AR = [500 300;]
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',AR,'paperposition',[0 0 AR])
set(gcf,'position',[50,500,AR])
hold on;
for ss = 1:nSubjs
    %p=plot(1:3,DPCF(ss,:),'linewidth',1,'color',[0.8 0.8 0.8]);
    p=plot(1:3,DPCF(ss,:),'linewidth',1,'color',[100 200 100]/255);
end
%plot(1:3,nanmean(DPCF),'linewidth',3,'color',[0.1 0.1 0.1]);
plot(1:3,nanmean(DPCF),'linewidth',4,'color',[0.1 0.1 0.1]);
set(gca,'fontsize',30,'xTick',1:3,'xticklabel',{'Low','Med','High'})
xlim([0.8 3.2])
ylim([-0.5 3.5])
set(gca,'LineWidth',3,'ytick',[0:3])
ylabel(' dP ')

print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/Behavior/' 'Confidence_dPrimeFaces' SubjSelectStr])
disp('Face dPrime by Confidence ')
disp(array2table(nanmean(DPCF),'variablenames',{'Low','Med','High'}))

% Scenes
DPCS = behav_out.retSummary.Scene_dPrimeConf(subjs,:);
figure(4); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[AR],'paperposition',[0 0 AR])
set(gcf,'position',[50,500,AR])
hold on;
for ss = 1:nSubjs
    p=plot(1:3,DPCS(ss,:),'linewidth',1,'color',[200 100 200]/255);
end
plot(1:3,nanmean(DPCS),'linewidth',4,'color',[0.1 0.1 0.1]);
set(gca,'fontsize',30,'xTick',1:3,'xticklabel',{'Low','Med','High'})
xlim([0.8 3.2])
ylim([-0.5 3.5])
set(gca,'LineWidth',3,'ytick',[0:3])
ylabel(' dP ')
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/Behavior/' 'Confidence_dPrimeScn' SubjSelectStr])

disp('Scene dPrime by Confidence ')
disp(array2table(nanmean(DPCS),'variablenames',{'Low','Med','High'}))

%% RTs
clearvars -except out behav_out subjs nSubjs dataPath SubjSelectStr
rng(1)
strs    = {'medianHit_RTs','medianMiss_RTs','medianFA_RTs','medianCRs_RTs'};
strs2   = {'Hits','Misses','FA','CRs'};
nRTConds = numel(strs);

RTs       = zeros(nSubjs,nRTConds);
for ii = 1:nRTConds
    RTs(:,ii) = behav_out.retSummary.(strs{ii})(subjs);
end
figure(1); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[600 400],'paperposition',[0 0 600 400])
set(gcf,'position',[100,100,600,400])

dx = 80;
xPosCore = [150:dx:500];
xPos = [xPosCore ];
a     = zeros(nRTConds,1);
for ii =1:nRTConds
    a(ii)=axes('units','points','position',[xPos(ii) 100 60 200]);
end
yLims = [0.8 2];
yTicks1 = [0:0.5:2.5];
x = randn(nSubjs,1)*0.1+0.5;
for ii=1:nRTConds
    axes(a(ii))
    %set(gca,'Position',[xPos(ii) 100 80 250])
    y = RTs(:,ii);
    s = scatter(x,y); hold on;
    s.MarkerFaceAlpha=0.5;
    s.MarkerEdgeAlpha=0.4;
    s.SizeData          = 120;
    s.MarkerEdgeColor = [119,136,153]/255;
    s.MarkerFaceColor = [119,136,153]/255;
    
    
    set(gca,'fontsize',20,'xTick','')
    plot([0.2 0.8], ones(1,2)*mean(y),'linewidth',4,'color','k')
    xlabel(strs2{ii})
    ylim(yLims)
    xlim([0 1])
    if ii==1
        set(gca,'ytick',yTicks1)
        ylabel(' RTs (s) ')
    end
    if ii>1
        set(gca,'ycolor','none')
    end
    set(gca,'LineWidth',2)
end

print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/Behavior/' 'RTs' SubjSelectStr])

disp(array2table(mean(RTs),'variablenames',strs2))
[~,p1,~,t1]=ttest(RTs(:,1),RTs(:,2));
[~,p2,~,t2]=ttest(RTs(:,1),RTs(:,3));
[~,p3,~,t3]=ttest(RTs(:,1),RTs(:,4));
p = [p1;p2;p3];
t = [t1.tstat;t2.tstat;t3.tstat];
disp(table(t,p,'rownames',{'HvsM','HvFA','HvCRs'}))

%% Retrieval RTs by Stimulys Category 
clearvars -except out behav_out subjs nSubjs dataPath SubjSelectStr
rng(1)
strs    = {'medianFaceScnHit_RTs','medianFaceScnMiss_RTs','medianFaceScnFA_RTs','medianFaceScnCRs_RTs'};
strs2   = {'Hits','Misses','FA','CRs'};
nRTConds = numel(strs);

RTsF       = zeros(nSubjs,nRTConds);
RTsS       = zeros(nSubjs,nRTConds);
for ii = 1:nRTConds
    RTsF(:,ii) = behav_out.retSummary.(strs{ii})(subjs,1);
    RTsS(:,ii) = behav_out.retSummary.(strs{ii})(subjs,2);
end
figure(1); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[600 400],'paperposition',[0 0 600 400])
set(gcf,'position',[100,100,600,400])

dx = 110;
xPosCore = [100:dx:500];
xPos = [xPosCore ];
a     = zeros(nRTConds,1);
for ii =1:nRTConds
    a(ii)=axes('units','points','position',[xPos(ii) 100 100 200]);
end
yLims = [0.8 2];
yTicks1 = [0:0.5:2.5];
x1 = randn(nSubjs,1)*0.03+0.3;
x2 = randn(nSubjs,1)*0.03+0.7;
for ii=1:nRTConds
    axes(a(ii)); hold on;
    y1  = RTsF(:,ii);
    y2  = RTsS(:,ii);
    
    for jj =1:nSubjs
    if y1(jj)>y2(jj)
        plot([x1(jj) x2(jj)], [y1(jj) y2(jj)],'-','color',[0.8 0.2 0.2])
    else
        plot([x1(jj) x2(jj)], [y1(jj) y2(jj)],'-','color',[0.6 0.6 0.6])
    end
    end

    % Faces
    s   = scatter(x1,y1); 
    plot([0.2 0.4], ones(1,2)*mean(y1),'linewidth',4,'color','k');    
    s.MarkerFaceAlpha   = 0.5;
    s.MarkerEdgeAlpha   = 0.4;
    s.SizeData          = 120;
    s.MarkerEdgeColor = [100 200 100]/255;
    s.MarkerFaceColor = [100 200 100]/255;
    
    set(gca,'fontsize',20,'xTick','')    
    xlabel(strs2{ii})
    ylim(yLims)
    xlim([0 1])
    
    % Scenes    
    s   = scatter(x2,y2);
    plot([0.6 0.8], ones(1,2)*mean(y2),'linewidth',4,'color','k')    
    s.MarkerFaceAlpha   = 0.5;
    s.MarkerEdgeAlpha   = 0.4;
    s.SizeData          = 120;
    s.MarkerEdgeColor = [200 100 200]/255;
    s.MarkerFaceColor = [200 100 200]/255;
    
    set(gca,'fontsize',20,'xTick','')
    
    xlabel(strs2{ii})
    ylim(yLims)
    xlim([0 1])
    
    if ii==1
        set(gca,'ytick',yTicks1)
        ylabel(' RTs (s) ')
    end
    if ii>1
        set(gca,'ycolor','none')
    end
    set(gca,'LineWidth',2)
end

print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/Behavior/' 'RetRTsCategory' SubjSelectStr])

%% RTs by Confidence
clearvars -except out behav_out subjs nSubjs dataPath SubjSelectStr
strs    = {'medianHit_RTsConf','medianMiss_RTsConf','medianFA_RTsConf','medianCRs_RTsConf'};
strs2   = {'Hits','Misses','FA','CRs'};
nRTConds = numel(strs);

RTsConf = zeros(nSubjs,3,nRTConds);
for ii = 1:nRTConds
    RTsConf(:,:,ii) = behav_out.retSummary.(strs{ii})(subjs,:);
end
figure(1); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[600 400],'paperposition',[0 0 600 400])
set(gcf,'position',[100,150,600,400])

dx = 120;
xPosCore = [100:dx:600];
xPos = [xPosCore ];
a     = zeros(nRTConds,1);
for ii =1:nRTConds
    a(ii)=axes('units','points','position',[xPos(ii) 100 100 200]);
end
yLims = [0.5 2.5];
yTicks1 = [0:1:2.5];
%
for ii=1:nRTConds
    axes(a(ii));
    hold on;
    %set(gca,'Position',[xPos(ii) 100 80 250])
    Y = RTsConf(:,:,ii);
    
    for ss=1:nSubjs
        p = plot(1:3,Y(ss,:));
        p.Color = [0.7 0.7 0.7];
        p.LineWidth = 1;
    end
    set(gca,'fontsize',20,'xTick',1:3, 'xticklabel',{'L','M','H'} )
    plot(1:3,nanmean(Y),'linewidth',4,'color','k')
    xlabel(strs2{ii})
    ylim(yLims)
    xlim([0.8 3.2])
    if ii==1
        set(gca,'ytick',yTicks1)
        ylabel(' RTs (s) ')
    end
    if ii>1
        set(gca,'ycolor','none')
    end
    set(gca,'LineWidth',2)
end

print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/Behavior/' 'RTsByConf' SubjSelectStr])

%% PropHits PropMisses and Difference by Phase
clearvars -except out behav_out subjs nSubjs dataPath SubjSelectStr

xa = linspace(0,2*pi,1000);
x = cos(xa);

% Hits
figure(4); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[1000 600],'paperposition',[0 0 1000 600])
set(gcf,'position',[100,100,600,400])

a2 = axes('units','points','position',[0.12*600 0.1*400 0.8*600 0.15*400]);
axes(a2)
plot(xa./pi*180,x,'k','linewidth',4)
axis tight;
set(gca,'ytick',[],'ycolor','w','fontsize',16,'box','off','lineWidth',2)
set(gca,'xtick',[0:72:360])
xlabel(' Encoding Phase (deg)')
grid on

a1 = axes('position', [0.12 0.3 0.8 0.6]);
axes(a1); hold on;
X = out.propHitsByPhase(subjs,:);
plot([36:72:360],X','-','color',[255 180 150]/255)
plot([36:72:360],mean(X), 'color',[240 80 40]/255,'linewidth',5)
set(gca,'fontsize',16,'box','off','lineWidth',2)
set(gca,'xtick',[36:72:360],'xTickLabel','')
xlim([0 360])
ylim([0.1 0.3])
ylabel(' proportion (hits) ' )

print(gcf, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_HitMiss/' 'HitsPropByPhase'  SubjSelectStr]);

% Misses
figure(5); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[1000 600],'paperposition',[0 0 1000 600])
set(gcf,'position',[100,100,600,400])

a2 = axes('units','points','position',[0.12*600 0.1*400 0.8*600 0.15*400]);
axes(a2)
plot(xa./pi*180,x,'k','linewidth',4)
axis tight;
set(gca,'ytick',[],'ycolor','w','fontsize',16,'box','off','lineWidth',2)
set(gca,'xtick',[0:72:360])
xlabel(' Encoding Phase (deg)')
grid on

a1 = axes('position', [0.12 0.3 0.8 0.6]);
axes(a1); hold on;
X = out.propMissByPhase(subjs,:);
plot([36:72:360],X','-','color',[150 220 220]/255)
plot([36:72:360],mean(X), 'color',[120 200 200]/255,'linewidth',5)
set(gca,'fontsize',16,'box','off','lineWidth',2)
set(gca,'xtick',[36:72:360],'xTickLabel','')
xlim([0 360])
ylim([0.1 0.3])
ylabel(' proportion (misses) ' )

print(gcf, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_HitMiss/' 'MissPropByPhase'  SubjSelectStr]);

% Difference
figure(6); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[1000 600],'paperposition',[0 0 1000 600])
set(gcf,'position',[100,100,600,400])

a2 = axes('units','points','position',[0.12*600 0.1*400 0.8*600 0.15*400]);
axes(a2)
plot(xa./pi*180,x,'k','linewidth',4)
axis tight;
set(gca,'ytick',[],'ycolor','w','fontsize',16,'box','off','lineWidth',2)
set(gca,'xtick',[0:72:360])
xlabel(' Encoding Phase (deg)')
grid on

a1 = axes('position', [0.12 0.3 0.8 0.6]);
axes(a1); hold on;
X = out.propHitsByPhase-out.propMissByPhase;
X = X(subjs,:);
plot([36:72:360],X','-','color',[180 180 180]/255)
plot([36:72:360],mean(X), 'color',[100 100 100]/255,'linewidth',5)
set(gca,'fontsize',16,'box','off','lineWidth',2)
set(gca,'xtick',[36:72:360],'xTickLabel','')
xlim([0 360])
ylim([-0.15 0.15])
ylabel(' p(h)-p(m) ' )

print(gcf, '-dpdf',['../plots/tacs_enc_xdiva/Phase_HitMiss/' 'H-MPropByPhase' SubjSelectStr]);

%% Hit Rate by phase
clearvars -except out behav_out subjs nSubjs dataPath SubjSelectStr

X = out.HRByPhase(subjs,:);
xa = linspace(0,2*pi,1000);
x = cos(xa);

figure(4); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[1000 600],'paperposition',[0 0 1000 600])
set(gcf,'position',[100,100,600,400])

a2 = axes('units','points','position',[0.12*600 0.1*400 0.8*600 0.15*400]);
axes(a2)
plot(xa./pi*180,x,'k','linewidth',4)
axis tight;
set(gca,'ytick',[],'ycolor','w','fontsize',16,'box','off','lineWidth',2)
set(gca,'xtick',[36:72:360])
xlabel(' Encoding Phase (deg)')
xlim([-1 361])
grid on

a1 = axes('position', [0.12 0.3 0.8 0.6]);
axes(a1); hold on;
plot([36:72:360],X','-','color',[180 180 180]/255)
plot([36:72:360],mean(X), 'color',[100 100 100]/255,'linewidth',5)
set(gca,'fontsize',16,'box','off','lineWidth',2)
set(gca,'xtick',[36:72:360],'xTickLabel','')
xlim([-1 361])
ylim([0 1])
ylabel(' Hit Rate ' )
print(gcf, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_HitMiss/' 'HitsRatepByPhase' SubjSelectStr]);

% mean Vectors
Z= X.*exp(1j*repmat([36:72:360]./180*pi,[nSubjs,1]));
mZ = mean(Z,2);
th = mod(angle(mZ),2*pi); rho = abs(mZ);
opts = [];
opts.colors = [180 180 180]/255;
opts.markerSize=300;
%opts.markerSize=behav_out.retSummary.dPrime(subjs)*150;
han = PolarPlot(th,rho,opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_HitMiss/' 'HitsRateMeanVec' SubjSelectStr]);
disp(table(circ_rtest(th),'variablenames',{'rhao_test_pval'}))

% 
opts.maxR = 4/3;
opts.markerSize=500;
opts.colors = [255 180 150]/255;
han = PolarPlot(th,[],opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_HitMiss/' 'HitsRateMeanVecAng' SubjSelectStr]);


% scatter of radii to dPrime
x = rho;
y = behav_out.retSummary.dPrime(subjs);
opts =[];
opts.colors = [150 150 150]/255;
opts.ylabel = ' dPrime ';
opts.xlabel = ' \rho ';
opts.polyfitN = 1;
opts.text       =['R = ' num2str(round(corr(x,y,'type','spearman')*100)/100)];
opts.xytext     = [0.04 1.1];
han = xyScatter(x,y,opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_HitMiss/' 'HitRateScatterRho-Dprime' SubjSelectStr]);

%% MemScore by Phase
clearvars -except out behav_out subjs nSubjs dataPath SubjSelectStr

xa = linspace(0,2*pi,1000);
x = cos(xa);

X = out.MemScore2ByPhase(subjs,:);
figure(7); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[1000 600],'paperposition',[0 0 1000 600])
set(gcf,'position',[100,100,600,400])

a2 = axes('units','points','position',[0.12*600 0.1*400 0.8*600 0.15*400]);
axes(a2)
plot(xa/pi*180,x,'k','linewidth',4)
axis tight;
set(gca,'ytick',[],'ycolor','w','fontsize',16,'box','off','lineWidth',2)
set(gca,'xtick',[36:72:360])
xlabel(' Encoding Phase (deg)')
xlim([-1 361])
grid on

a1 = axes('position', [0.12 0.3 0.8 0.6]);
axes(a1); hold on;
Xm=(X-repmat(mean(X,2),[1,5]));
plot([36:72:360],Xm','-','color',[180 180 180]/255)
plot([36:72:360],mean(Xm), 'color',[100 100 100]/255,'linewidth',5)
set(gca,'fontsize',16,'box','off','lineWidth',2)
set(gca,'xtick',[36:72:360],'xTickLabel','')
xlim([-1 361])
ylim([-1 1])
ylabel(' MemScore ' )
print(gcf, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_HitMiss/' 'MemScoreByPhase' SubjSelectStr] );

% polar plot
Z=X.*exp(1j*repmat([36:72:360]./180*pi,[nSubjs,1]));
mZ = mean(Z,2);
th = mod(angle(mZ),2*pi); rho = abs(mZ);
opts = [];
opts.colors = [180 180 180]/255;
opts.markerSize=300;%behav_out.retSummary.dPrime(subjs)*150;
han = PolarPlot(th,rho,opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_HitMiss/'  'MemScoreMeanVec' SubjSelectStr]);
disp(table(circ_rtest(th),'variablenames',{'rhao_test_pval'}))

opts.maxR = 4/3;
opts.markerSize=300;
han = PolarPlot(th,[],opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_HitMiss/' 'MemScoreMeanVecAng' SubjSelectStr]);

% scatter of radii to dPrime
x = rho;
y = behav_out.retSummary.dPrime(subjs);
opts =[];
opts.colors = [150 150 150]/255;
opts.ylabel = ' dPrime ';
opts.xlabel = ' \rho ';
opts.polyfitN = 1;
opts.text       =['R = ' num2str(round(corr(x,y,'type','spearman')*100)/100)];
opts.xytext     = [0.04 1.1];
han = xyScatter(x,y,opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_HitMiss/' 'MemScoreRho-Dprime' SubjSelectStr]);

%% Hit Miss Phase distributions
clearvars -except out behav_out subjs nSubjs dataPath SubjSelectStr

opts = [];
opts.colors = [255 180 150; 150 220 220]/255;
opts.markerSize=squeeze(sum(behav_out.retSummary.nH_nMiss_nFA_nCRs(subjs,:,1:2),2));

% raw vectors
th  = mod(out.HM_MeVects(subjs,:,1),2*pi);
rho = out.HM_MeVects(subjs,:,2);
opts.markerSize=300;
han = PolarPlot(th,rho,opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_HitMiss/' 'MeanPhaseVectsByCond' SubjSelectStr]);

%test and disply circular uniformity
disp('Rayleigh Test:')
p = zeros(2,1); u = p;
for ii = 1:2
    [p(ii),u(ii)]=circ_rtest(th(:,ii));
end
disp(table(p,u,'VariableNames',{'P_Val','Z'},'rownames',{'Hits','Misses'}))

% raw only theta (ignores individuals SS strength)
opts.maxR = 4/3;
opts.markerSize=300;
han = PolarPlot(th,ones(nSubjs,2),opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_HitMiss/''MeanPhaseVectsByCond-Theta' SubjSelectStr]);

%% line between points
clearvars -except out behav_out subjs nSubjs dataPath SubjSelectStr

th  = mod(out.HM_MeVects(subjs,:,1),2*pi);
rho = out.HM_MeVects(subjs,:,2);
z   = rho.*exp(1i*th);

opts = [];
%opts.markerSize =squeeze(sum(behav_out.retSummary.nH_nMiss_nFA_nCRs(subjs,:,1:2),2));
opts.markerSize = 300;
opts.colors = [255 180 150;150 220 220]/255;
opts.connect = 1; opts.polarGrid =0; opts.meanVecs = 0;
opts.magText =0;
han = PolarPlot(th,rho,opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_HitMiss/' 'MeanPhaseVectsByCond-Connected'  SubjSelectStr]);

% line between hit-miss angles without magnitude
opts.connect = 1; opts.meanVecs = 0; opts.polarGrid=0;
han = PolarPlot(th,ones(nSubjs,2),opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_HitMiss/' 'MeanPhaseVectsByCond-ThetaConnected'  SubjSelectStr]);

% figure with angles re-centered
thc  = th-repmat(th(:,1),[1 2]);
opts.polarGrid = 1;
han = PolarPlot(thc,rho,opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_HitMiss/' 'MeanPhaseVectsByCond-ConnectedReCenter'  SubjSelectStr]);

% mean Vecs re-centered
opts.meanVecs = 1;
opts.connect  = 0;
han = PolarPlot(thc,rho,opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_HitMiss/' 'MeanPhaseVectsByCond-ReCenter'  SubjSelectStr]);

%% difference in hit phase and miss phase
clearvars -except out behav_out subjs nSubjs dataPath SubjSelectStr
% mean differences
th = mod(out.HM_MeVects(subjs,:,1),2*pi);
rho = out.HM_MeVects(subjs,:,2);
dTh= th(:,1)-th(:,2);

z = rho.*exp(1j*th);
R = abs(z(:,1)-z(:,2));

opts = [];
opts.colors = [150 150 150]/255;
opts.markerSize = 300;
han = PolarPlot(dTh,R,opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_HitMiss/' 'MeanPhaseVectDiff' SubjSelectStr]);

% difference in hit phase and miss phase Rstat
rhoR = out.HM_R(subjs,:);
z = rhoR.*exp(1j*th); RR = abs(z(:,1)-z(:,2));
han = PolarPlot(dTh,RR,opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_HitMiss/' 'MeanPhaseVectDiff-Rstat' SubjSelectStr]);

% only Delta in Theta
opts.maxR = 4/3;
han = PolarPlot(dTh,ones(nSubjs,1),opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_HitMiss/' 'MeanPhaseVectDiffTheta' SubjSelectStr]);

% variance decreases as mean vector length increases
%d = 0.5*(1-cos(dTh));
x = R;
d = abs(pi-abs(dTh));
opts =[];
opts.colors = [150 150 150]/255;
opts.xlabel = ' \rho ';
opts.ylabel = ' |\pi-|\Delta \theta| |';
opts.polyfitN = 1;
opts.text       =['R = ' num2str(round(corr(x,d,'type','spearman')*100)/100)];
opts.xytext     = [0.1 0.5];
han = xyScatter(x,d,opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_HitMiss/' 'DiffAbAllHitMissesScatter' SubjSelectStr]);

% variance decreases as mean vector length increases R-stat
%d = 0.5*(1-cos(th(:,1)-th(:,2)));
x = RR;
opts.text       =['R = ' num2str(round(corr(x,d,'type','spearman')*100)/100)];
opts.xytext     = [1 0.5];
han = xyScatter(x,d,opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_HitMiss/' 'Diff-RstatAllHitMissesScatter' SubjSelectStr]);

[p,u]=circ_rtest(th(:,1)-th(:,2));
disp( 'Rayleigh test for Hit-Miss Angle: ')
disp(table(p,u,'variablenames',{'P_Val','Z'},'rownames',{'Hit-Miss'}))

%% HM vec Length By Confidence
Ns = squeeze(sum(out.HM_Conf_N(subjs,:,:),2));
lowNsubjs = find(sum(squeeze(sum(Ns<10,2)),2));
subjs2 = setdiff(subjs,lowNsubjs);
Ns(lowNsubjs,:)      = [];
nSubjs2 = numel(subjs2);
X=out.HM_Conf_Z(subjs2,:);
Y=Ns;
%
figure(); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[500 300],'paperposition',[0 0 500 300])
set(gcf,'position',[50,500,500,300])

hold on;
for ss = 1:nSubjs2
    p=plot(1:3,X(ss,:),'linewidth',1,'color',[0.8 0.8 0.8]);
end
plot(1:3,nanmedian(X),'linewidth',3,'color',[0.1 0.1 0.1]);
set(gca,'fontsize',20,'xTick',1:3,'xticklabel',{'Low','Med','High'})
xlim([0.8 3.2])
ylim([0 0.55])
set(gca,'LineWidth',2,'ytick',[0 0.2 0.4])
ylabel(' P ')

print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/Phase_HitMiss/VecLengthByConf' SubjSelectStr])
%% High Confidence Hits.
HitsConf    = squeeze(out.HM_Conf_MeVects(subjs,1,:,:));
HitsConfR   = squeeze(out.HM_Conf_R(subjs,1,:));
HitsConfN   = squeeze(out.HM_Conf_N(subjs,1,:));
Conf        = squeeze(out.Conf_MeVects(subjs,:,:));
ConfR       = squeeze(out.Conf_R(subjs,:));
ConfN       = squeeze(out.ConfPhases_N(subjs,1,:,:));

opts = [];
opts.colors = [255 180 120; 180 200 120]/255;
magMult = 25;

th  = [HitsConf(:,3,1) HitsConf(:,1,1)];
%rho = [HitsConfR(:,3) HitsConfR(:,1)];
opts.maxR = 4/3;
opts.alpha = 0.9;
opts.markerSize = [HitsConfR(:,3) HitsConfR(:,1)]*magMult;
han = PolarPlot(th,[],opts); 
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_HitMiss/High_LowConfHits' SubjSelectStr]);

txtStr = {'Hits HC','Hits LC'};
makeLegend(opts.colors,[200;200],txtStr);
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/HitsHC_HitsLC']);


th = [HitsConf(:,3,1)-HitsConf(:,1,1)];
opts.markerSize = 300;
opts.colors = [0.8 0.85 0.9]/1.2;
han = PolarPlot(th,[],opts); 
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_HitMiss/High-LowConfHits' SubjSelectStr]);

% HC hits vs low conf resp
th  = [HitsConf(:,3,1) Conf(:,1,1)];
opts.maxR = 4/3;

% only the angle
opts.alpha = 0.9;
opts.markerSize=300;
opts.colors = [255 180 120; 100 120 150]/255;
opts.markerSize = [HitsConfR(:,3) ConfR(:,1)]*magMult;
han = PolarPlot(th,[],opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_HitMiss/HitsHighConf_AllLowConf' SubjSelectStr]);

txtStr = {'Hits HC','All LC'};
makeLegend(opts.colors,[200;200],txtStr);
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/HitsHC_AllLC']);

colors = repmat([0 0 0],[4,1]);
qq = [3:3:12];
%qq = quantile(opts.markerSize(:),[0.1:0.25:1]);
%txtStr = cellfun(@num2str,num2cell(round(qq/magMult)),'UniformOutput',0);
txtStr = cellfun(@num2str,num2cell(qq),'UniformOutput',0);
makeLegend(colors,qq*magMult,txtStr);
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/MagnitudeLeg']);

%% Difference of HC hits - all low confidence
th  = [HitsConf(:,3,1)-Conf(:,1,1)];
opts.markerSize = 300;
opts.colors = [0.8 0.85 0.9]/1.2;
han = PolarPlot(th,[],opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_HitMiss/HitsHighConf-AllLowConf' SubjSelectStr]);

%% all confidence hits and misses

% HC Hits only
th = HitsConf(:,3,1);
opts.maxR = 4/3;
opts.alpha = 0.8;
opts.markerSize = HitsConfR(:,3)*magMult;
opts.colors = [255 180 150]/255;
han = PolarPlot(th,[],opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_HitMiss/HitsHighConf' SubjSelectStr]);

th = HitsConf(:,:,1);
opts.maxR = 4/3;
opts.alpha = 0.8;
opts.markerSize = HitsConfR*magMult;
opts.colors = brewermap(3,'set2');
han = PolarPlot(th,[],opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_HitMiss/HitsConf' SubjSelectStr]);


MissConf    = squeeze(out.HM_Conf_MeVects(subjs,2,:,1));
MissConfR   = squeeze(out.HM_Conf_R(subjs,2,:));
th = MissConf(:,:);
opts.maxR = 4/3;
opts.alpha = 0.8;
opts.colors = brewermap(3,'set2');
opts.markerSize = MissConfR*magMult;
han = PolarPlot(th,[],opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_HitMiss/MissConf' SubjSelectStr]);

%% Confidence independent of memory
th  = Conf(:,:,1);
opts.colors = brewermap(3,'set2');
opts.maxR = 4/3;
opts.alpha = 0.8;
%opts.markerSize=300;
opts.markerSize = ConfR*magMult;
han = PolarPlot(th,[],opts);
print(han,'-dpdf',['../plots/tacs_enc_xdiva/Phase_Confidence/ConfidencePolar']);
txtStr = {'Low','Med','High'};
makeLegend(opts.colors,[300;300;300],txtStr);
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/Phase_Confidence/ConfidenceLeg']);
%% Hits High Confidence by Face Scene
FH_HC=squeeze(out.HM_Conf_FaScn_MeVects(subjs,1,3,1,:));
FH_HC_R = out.HM_Conf_FaScn_R(subjs,1,3,1);
SH_HC=squeeze(out.HM_Conf_FaScn_MeVects(subjs,1,3,2,:));
SH_LC=squeeze(out.HM_Conf_FaScn_MeVects(subjs,1,1,2,:));
SH_HC_R = out.HM_Conf_FaScn_R(subjs,1,3,2);
SH_LC_R = out.HM_Conf_FaScn_R(subjs,1,1,2);

magMult = 25;
opts = [];
opts.colors = [255 180 150]/255;
opts.maxR = 4/3;
opts.alpha = 0.8;
opts.markerSize=300;
opts.markerSize = [FH_HC_R]*magMult;


% FAces High Confidence Hits 
th  = [FH_HC(:,1)];
han = PolarPlot(th,[],opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_FaceScene/FaceHighConf' SubjSelectStr]);

% Scenes High Confidence Hits 
th  = SH_HC(:,1);
opts.markerSize = [SH_HC_R]*magMult;
han = PolarPlot(th,[],opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_FaceScene/ScnHighConf' SubjSelectStr]);

% both on the same plot
th = [FH_HC(:,1),SH_HC(:,1)];
opts.markerSize = [FH_HC_R SH_HC_R]*magMult;
opts.colors = [100 200 100; 200 100 200]/255;
han = PolarPlot(th,[],opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_FaceScene/FaScnHighConf' SubjSelectStr]);
%% magnitude of HC effects
FH_HC=squeeze(out.HM_Conf_FaScn_MeVects(subjs,1,3,1,:));
FH_HC_R = out.HM_Conf_FaScn_R(subjs,1,3,1);
FH_LC_R = out.HM_Conf_FaScn_R(subjs,1,1,1);
SH_HC=squeeze(out.HM_Conf_FaScn_MeVects(subjs,1,3,2,:));
SH_LC=squeeze(out.HM_Conf_FaScn_MeVects(subjs,1,1,2,:));
SH_HC_R = out.HM_Conf_FaScn_R(subjs,1,3,2);
SH_LC_R = out.HM_Conf_FaScn_R(subjs,1,1,2);
ConfR       = squeeze(out.Conf_R(subjs,:));

x = [FH_HC_R./out.HM_FaScn_N(subjs,1,1) SH_HC_R./out.HM_FaScn_N(subjs,1,2)];
x = [FH_HC_R SH_HC_R];

figure(); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[300 200],'paperposition',[0 0 300 200])
set(gcf,'position',[50,500,300,200])

hold on;
for ss = 1:nSubjs
    
    if x(ss,1)<x(ss,2)
        p=plot(1:2,x(ss,:),'linewidth',1,'color',[0.8 0.8 0.8]);
    else
        p=plot(1:2,x(ss,:),'linewidth',1,'color',[0.8 0.2 0.2]);
    end
end
plot(1:2,nanmean(x),'linewidth',3,'color',[0.1 0.1 0.1]);
set(gca,'fontsize',20,'xTick',1:2,'xticklabel',{'Face','Scene'})
xlim([0.8 2.2])
%ylim([0 0.15])
set(gca,'LineWidth',2,'ytick',[0:3:15])
ylabel(' nTrials ')

print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/Phase_FaceScene/MagNormHighConf_FaceScn' SubjSelectStr])

%
x = randn(nSubjs,1)*0.05+0.5;
y = [FH_HC_R SH_HC_R];

figure(); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[300 200],'paperposition',[0 0 300 200])
set(gcf,'position',[50,500,300,200]); hold on;

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
s.MarkerEdgeAlpha   = 0.4;
s.SizeData          = 120;
s.MarkerEdgeColor   = [200 100 200]/255;
s.MarkerFaceColor   = [200 100 200]/255;

xlim([-0.1 2.1])
ylim([-0.1 15])
set(gca,'fontsize',20,'xtick',[0.5 1.5],'xticklabel',{'F','S'})
set(gca,'LineWidth',2,'ytick',[0:3:15])
ylabel(' nTrials ')
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/Phase_FaceScene/' 'MagNormHighConf_FaceScn2' SubjSelectStr])

%% Face / Scene differences independent of memory
clearvars -except out behav_out subjs nSubjs dataPath SubjSelectStr

%th      = mod(out.FaScnMePhases(subjs,:),2*pi);
%rho     = out.FaScnAbPhases(subjs,:);
th      = out.FaScn_MeVects(subjs,:,1);
rho      = out.FaScn_MeVects(subjs,:,2);
opts            = [];
opts.colors     = [100 200 100; 200 100 200]/255;
%opts.markerSize = out.FaScnHitN(subjs,:)+out.FaScnMissN(subjs,:);
opts.maxR       = 1/100;
han = PolarPlot(th,rho,opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_FaceScene/' 'MeanPhaseVectsFaceScn' SubjSelectStr]);

%tests of circular uniformity
disp('Rayleigh Test:')
p = zeros(2,1); u = p;
for ii = 1:2
    [p(ii),u(ii)]=circ_rtest(th(:,ii));
end
disp(table(p,u,'VariableNames',{'P_Val','Z'},'rownames',{'Faces','Scenes'}))

% difference
dTh  = th(:,1)-th(:,2);
z = rho.*exp(1j*th); R = abs(z(:,1)-z(:,2));

opts.colors = [150 150 150]/255;
han = PolarPlot(dTh,R,opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_FaceScene/' 'MeanPhaseVectsFaceScnDiff' SubjSelectStr]);

[p,u]=circ_rtest(dTh);
disp( 'Rayleigh test for Face-Scene Angle: ')
disp(table(p,u,'variablenames',{'P_Val','Z'},'rownames',{'Face-Scene'}))

%% Faces difference by memory
clearvars -except out behav_out subjs nSubjs dataPath SubjSelectStr

% Hit/Miss for Faces
%thF = mod([out.FaScnHitMePhases(subjs,1) , out.FaScnMissMePhases(subjs,1)],2*pi);
%rhoF = [out.FaScnHitAbPhases(subjs,1), out.FaScnMissAbPhases(subjs,1)];
thF  = [out.HM_FaScn_MeVects(subjs,1,1,1) out.HM_FaScn_MeVects(subjs,2,1,1)];
rhoF = [out.HM_FaScn_MeVects(subjs,1,1,2) out.HM_FaScn_MeVects(subjs,2,1,2)];
dThF = mod(thF(:,1)-thF(:,2),2*pi);
zF   = rhoF.*exp(1j*thF); RF = abs(zF(:,1)-zF(:,2));

opts        = [];
opts.colors = [50 100 50;100 200 100]/255;
%opts.markerSize = [out.FaScnHitN(subjs,1) out.FaScnMissN(subjs,1)];
opts.markerSize= out.HM_FaScn_N(subjs,:,1);
han = PolarPlot(thF,rhoF,opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_FaceScene/' 'MeanPhaseVectsFacesHitMiss' SubjSelectStr]);

% hit/miss for faces projected to the unit circle
opts.maxR = 4/3;
han = PolarPlot(thF,ones(nSubjs,2),opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_FaceScene/' 'MeanPhaseVectsFacesHitMissTheta' SubjSelectStr]);

% hit/miss for faces lines
opts.maxR = 4/3;
opts.connect = 1; opts.meanVecs = 0; opts.polarGrid =0;
han = PolarPlot(thF,ones(nSubjs,2),opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_FaceScene/' 'MeanPhaseVectsFacesHitMissThetaConnected' SubjSelectStr]);

%tests of circular uniformity
disp('Rayleigh Test Faces:')
p = zeros(2,1); u = p;
for ii = 1:2
    [p(ii),u(ii)]=circ_rtest(thF(:,ii));
end
disp(table(p,u,'VariableNames',{'P_Val','Z'},'rownames',{'Hits','Misses'}))

% Hit Miss Difference for Faces
opts        = []; opts.colors = [150 150 150]/255;
han = PolarPlot(dThF,RF,opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_FaceScene/' 'MeanPhaseVectsFacesHitMissDiff' SubjSelectStr]);

[p,u]=circ_rtest(thF(:,1)-thF(:,2));
disp( 'Rayleigh test for Hit-Miss Angle Faces: ')
disp(table(p,u,'variablenames',{'P_Val','Z'},'rownames',{'Hit-Miss'}))

%% Scenes difference by memory
clearvars -except out behav_out subjs nSubjs dataPath SubjSelectStr

% Hit/Miss for Scenes
thS = mod([out.FaScnHitMePhases(subjs,2) , out.FaScnMissMePhases(subjs,2)],2*pi);
rhoS= [out.FaScnHitAbPhases(subjs,2), out.FaScnMissAbPhases(subjs,2)];
dThS = mod(thS(:,1)-thS(:,2),2*pi);
zS = rhoS.*exp(1j*thS); RS = abs(zS(:,1)-zS(:,2));

opts        = [];
opts.colors = [100 50 100;200 100 200]/255;
opts.markerSize = [out.FaScnHitN(subjs,2) out.FaScnMissN(subjs,2)];
han = PolarPlot(thS,rhoS,opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_FaceScene/' 'MeanPhaseVectsScnHitMiss' SubjSelectStr]);

% hit/miss for scenes projected to the unit circle
opts.maxR = 4/3;
han = PolarPlot(thS,ones(nSubjs,2),opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_FaceScene/' 'MeanPhaseVectsScnHitMissTheta' SubjSelectStr]);

% hit/miss for scenes lines
opts.maxR = 4/3;
opts.connect = 1; opts.meanVecs = 0; opts.polarGrid =0;
han = PolarPlot(thS,ones(nSubjs,2),opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_FaceScene/' 'MeanPhaseVectsFacesHitMissThetaConnected' SubjSelectStr]);

%tests of circular uniformity
disp('Rayleigh Test Scenes:')
p = zeros(2,1); u = p;
for ii = 1:2
    [p(ii),u(ii)]=circ_rtest(thS(:,ii));
end
disp(table(p,u,'VariableNames',{'P_Val','Z'},'rownames',{'Hits','Misses'}))

% Hit Miss Difference for Scenes
opts        = []; opts.colors = [150 150 150]/255;
han = PolarPlot(dThS,RS,opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_FaceScene/' 'MeanPhaseVectsScnHitMissDiff' SubjSelectStr]);

[p,u]=circ_rtest(dThS);
disp( 'Rayleigh test for Hit-Miss Scenes: ')
disp(table(p,u,'variablenames',{'P_Val','Z'},'rownames',{'Hit-Miss'}))

%% Face / Scene Hits
close all
clearvars -except out behav_out subjs nSubjs dataPath SubjSelectStr

th  = mod(out.FaScnHitMePhases(subjs,:),2*pi);
rho = out.FaScnHitAbPhases(subjs,:);
dTh = mod(th(:,1)-th(:,2),2*pi);
z   = rho.*exp(1j*th); R = abs(z(:,1)-z(:,2));

opts = [];
opts.colors = [50 100 50; 100 50 100]/255;
opts.markerSize = out.FaScnHitN(subjs,:);

han = PolarPlot(th,rho,opts);
print(han, '-dpdf',['../plots/tacs_enc_xdiva/Phase_FaceScene/' 'MeanPhaseVectsFaScnHits' SubjSelectStr]);

% hits projected to the unit circle
opts.maxR = 4/3;
han = PolarPlot(th,ones(nSubjs,2),opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_FaceScene/' 'MeanPhaseVectsFaScnHitsTheta' SubjSelectStr]);

% hits for connected lines
opts.maxR = 4/3;
opts.connect = 1; opts.meanVecs = 0; opts.polarGrid =0;
han = PolarPlot(th,ones(nSubjs,2),opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_FaceScene/' 'MeanPhaseVectsFaScnHitsThetaConnected' SubjSelectStr]);

% Hits Difference for Faces and Scenes
opts        = []; opts.colors = [150 150 150]/255;
han = PolarPlot(dTh,R,opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_FaceScene/' 'MeanPhaseVectsFaScnHitsDiff' SubjSelectStr]);

[p,u]=circ_rtest(dTh);
disp( 'Rayleigh test for Hits Angle Faces-Scn: ')
disp(table(p,u,'variablenames',{'P_Val','Z'},'rownames',{'Fa-Scn'}))

%% Face / Scene Misses
close all
clearvars -except out behav_out subjs nSubjs dataPath SubjSelectStr

th  = mod(out.FaScnMissMePhases(subjs,:),2*pi);
rho = out.FaScnMissAbPhases(subjs,:);
dTh = mod(th(:,1)-th(:,2),2*pi);
z   = rho.*exp(1j*th); R = abs(z(:,1)-z(:,2));

opts = [];
opts.colors = [100 200 100; 200 100 200]/255;
opts.markerSize = out.FaScnMissN(subjs,:);

han = PolarPlot(th,rho,opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_FaceScene/' 'MeanPhaseVectsFaScnMiss' SubjSelectStr]);

% hits projected to the unit circle
opts.maxR = 4/3;
han = PolarPlot(th,ones(nSubjs,2),opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_FaceScene/' 'MeanPhaseVectsFaScnMissTheta' SubjSelectStr]);

% hits for connected lines
opts.maxR = 4/3;
opts.connect = 1; opts.meanVecs = 0; opts.polarGrid =0;
han = PolarPlot(th,ones(nSubjs,2),opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_FaceScene/' 'MeanPhaseVectsFaScnMissThetaConnected' SubjSelectStr]);

% Hits Difference for Faces and Scenes
opts        = []; opts.colors = [150 150 150]/255;
han = PolarPlot(dTh,R,opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_FaceScene/' 'MeanPhaseVectsFaScnMissDiff' SubjSelectStr]);

[p,u]=circ_rtest(dTh);
disp( 'Rayleigh test for Miss Angle Faces-Scn: ')
disp(table(p,u,'variablenames',{'P_Val','Z'},'rownames',{'Fa-Scn'}))



