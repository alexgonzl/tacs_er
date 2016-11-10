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
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/BehavRetrieval/' 'dPrimes' SubjSelectStr])

m = mean(D(:,1));
[~,p,~,t] = ttest(D(:,1));
fprintf('Across subject overall dprime = %0.2f, t=%0.2f, p=%0.2g\n',m,t.tstat,p)
%% dP by Confidence
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

print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/BehavRetrieval/'  'Confidence_dPrime' SubjSelectStr])


disp(array2table(mean(DPC),'variablenames',{'Low','Med','High'}))
[~,p1,~,t1]=ttest(DPC(:,2),DPC(:,1));
[~,p2,~,t2]=ttest(DPC(:,3),DPC(:,2));
[~,p3,~,t3]=ttest(DPC(:,3),DPC(:,1));
p = [p1;p2;p3];
t = [t1.tstat;t2.tstat;t3.tstat];
disp(table(t,p,'rownames',{'Mid>Low','Hi>Mid','Hi>Lo'}))

%anova1(DPC(:),[ones(nSubjs,1);2*ones(nSubjs,1);3*ones(nSubjs,1)])
% Repeated measures Anova table
t=array2table(DPC);
conf = dataset([1:3]','varnames',{'confidence'});
rm = fitrm(t,'DPC1-DPC3~1','withindesign',conf);
ranova(rm,'withinmodel','confidence-1')

%% dP Confidence Face/Scenes
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

print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/BehavRetrieval/' 'Confidence_dPrimeFaces' SubjSelectStr])
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
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/BehavRetrieval/' 'Confidence_dPrimeScn' SubjSelectStr])

disp('Scene dPrime by Confidence ')
disp(array2table(nanmean(DPCS),'variablenames',{'Low','Med','High'}))

% Statistics
fprintf('Model for dP Confidence Faces: \n')
t = array2table(DPCF);
within = dataset([1:3]','varnames',{'confidence'});
rm = fitrm(t,'DPCF1-DPCF3~1','withindesign',within);
ranova(rm,'withinmodel','confidence-1')

fprintf('Model for dP Confidence Scenes: \n')
t = array2table(DPCS);
within = dataset([1:3]','varnames',{'confidence'});
rm = fitrm(t,'DPCS1-DPCS3~1','withindesign',within);
ranova(rm,'withinmodel','confidence-1')

fprintf('Model for dP Confidence X Face/Scenes: \n')
t = array2table([DPCF DPCS]);
within = dataset([1:3,1:3]',['F','F','F','S','S','S']','varnames',{'confidence','category'});
rm = fitrm(t,'Var1-Var6~1','withindesign',within);
ranova(rm,'withinmodel','confidence*category-1')
% t=array2table(DPC);
% conf = dataset([1:3]','varnames',{'confidence'});
% rm = fitrm(t,'DPC1-DPC3~1','withindesign',conf);
% ranova(rm)

%

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

print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/BehavRetrieval/' 'RTs' SubjSelectStr])

% statistics
disp(array2table(mean(RTs),'variablenames',strs2))
y = -log(RTs);
[~,p1,~,t1]=ttest(y(:,1),y(:,2));
[~,p2,~,t2]=ttest(y(:,1),y(:,3));
[~,p3,~,t3]=ttest(y(:,1),y(:,4));
p = [p1;p2;p3];
t = [t1.tstat;t2.tstat;t3.tstat];
disp(table(t,p,'rownames',{'HvsM','HvFA','HvCRs'}))

t=array2table(y);
conds = dataset([1:4]','varnames',{'MemoryCondition'});
rm = fitrm(t,'y1-y4~1','withindesign',conds);
ranova(rm,'withinmodel','MemoryCondition-1')

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

print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/BehavRetrieval/' 'RetRTsCategory' SubjSelectStr])

% Statistics
fprintf('Model for RTs Category X Memory Conidtion: \n')
y1 = -log(RTsF);
y2 = -log(RTsS);
t = array2table([y1 y2]);
within = dataset([1:4,1:4]',['F','F','F','F','S','S','S','S']','varnames',{'memory','category'});
rm = fitrm(t,'Var1-Var8~1','withindesign',within);
ranova(rm,'withinmodel','memory*category-1')

% statistics
disp(array2table(mean([y1 y2]),'variablenames',[strcat('F_',strs2) strcat('S_',strs2)]))
for ii = 1:4
    [~,p,~,t]=ttest(y1(:,ii),y2(:,ii));
    fprintf('RTs %s Face vs Scene t=%0.2f, p=%0.2g \n',strs2{ii},t.tstat,p)
end

%% RTs by Confidence and memory
clearvars -except out behav_out subjs nSubjs dataPath SubjSelectStr
strs    = {'medianHit_RTsConf','medianMiss_RTsConf','medianFA_RTsConf','medianCRs_RTsConf'};
strs2   = {'Hits','Misses','FA','CRs'};
nRTConds = numel(strs);

RTsConf = zeros(nSubjs,3,nRTConds);
for ii = 1:nRTConds
    RTsConf(:,:,ii) = behav_out.retSummary.(strs{ii})(subjs,:);
end


figure(1); clf; AR = [550 300];
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',AR,'paperposition',[0 0 AR])
set(gcf,'position',[100,100,AR])

dx = 120;
xPosCore = [80:dx:550];
xPos = [xPosCore ];
a     = zeros(nRTConds,1);
for ii =1:nRTConds
    a(ii)=axes('units','points','position',[xPos(ii) 100 100 150]);
end
yLims = [0.5 2.5];
yTicks1 = [0:0.5:2.5];
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

print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/BehavRetrieval/' 'RTsByConf' SubjSelectStr])

fprintf('Model for RTs Confidence X Memory Conidtion: \n')
% Statistics
y = -log10(RTsConf(:,:));
t = array2table(y);
within = dataset(repmat((1:3)',[4,1]),['H','M','F','C','H','M','F','C'...
    ,'H','M','F','C',]','varnames',{'confidence','memory'});
rm = fitrm(t,'y1-y12~1','withindesign',within);
ranova(rm,'withinmodel','memory*confidence-1')


%% proportion of HC hits to to other categories:

X=behav_out.retSummary.nH_nMiss_nFA_nCRs(subjs,:,:);
HC = squeeze(X(:,3,:));
MC = squeeze(X(:,2,:));
LC = squeeze(X(:,1,:));

Ns=squeeze(sum(X,2));
propHC = HC./Ns;
fprintf('Proportion of HC responses by memory:\n')
disp(array2table(mean(propHC),'variablenames',{'H','M','FA','CRs'}))
[~,p,~,t]=ttest(propHC(:,1),propHC(:,2));
fprintf('Prop HC Hits > Prop HC M: t=%0.2f,p=%0.2g\n',t.tstat,p)
[~,p,~,t]=ttest(propHC(:,1),propHC(:,3));
fprintf('Prop HC Hits > Prop HC FA: t=%0.2f,p=%0.2g\n',t.tstat,p)
[~,p,~,t]=ttest(propHC(:,1),propHC(:,4));
fprintf('Prop HC Hits > Prop HC CRs: t=%0.2f,p=%0.2g\n',t.tstat,p)

% Proportion of HC Hits to MC to LC Hits
propLC = LC./Ns;
propMC = MC./Ns;

% Plot # responses
figure(1); clf;
AR = [550 300];
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[AR],'paperposition',[0 0 AR])
set(gcf,'position',[0,0,AR])

dx = 120;
xPosCore = [80:dx:600];
xPos = [xPosCore ];
a     = zeros(nRTConds,1);
for ii =1:4
    a(ii)=axes('units','points','position',[xPos(ii) 100 100 150]);
end
Y=X;
for ii=1:4
    yLims = [0 120];
    yTicks1 = [0:25:100];
    axes(a(ii)); cla;
    hold on;
    Yii = Y(:,:,ii);
    
    for ss=1:nSubjs
        p = plot(1:3,Yii(ss,:));
        p.Color = [0.7 0.7 0.7];
        p.LineWidth = 1;
    end
    set(gca,'fontsize',20,'xTick',1:3, 'xticklabel',{'L','M','H'} )
    plot(1:3,nanmean(Yii),'linewidth',4,'color','k')
    xlabel(strs2{ii})
    ylim(yLims)
    xlim([0.8 3.2])
    if ii==1
        set(gca,'ytick',yTicks1)
        ylabel(' # Trials ')
    end
    if ii>1
        set(gca,'ycolor','none')
    end
    set(gca,'LineWidth',2)
end
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/BehavRetrieval/' 'nRespConf' SubjSelectStr])
%statistcs for Hits%

% Plot proportions
Y = cat(3,propLC,propMC,propHC)*100;
for ii=1:4
    yLims = [0 90];
    yTicks1 = [0:25:75];
    axes(a(ii)); cla;
    hold on;
    Yii = squeeze(Y(:,ii,:));
    
    for ss=1:nSubjs
        p = plot(1:3,Yii(ss,:));
        p.Color = [0.7 0.7 0.7];
        p.LineWidth = 1;
    end
    set(gca,'fontsize',20,'xTick',1:3, 'xticklabel',{'L','M','H'} )
    plot(1:3,nanmean(Yii),'linewidth',4,'color','k')
    xlabel(strs2{ii})
    ylim(yLims)
    xlim([0.8 3.2])
    if ii==1
        set(gca,'ytick',yTicks1)
        ylabel(' % Responses ')
    end
    if ii>1
        set(gca,'ycolor','none')
    end
    set(gca,'LineWidth',2)
end
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/BehavRetrieval/' 'propRespConf' SubjSelectStr])

%% RTs by confidence/category
RTs_Faces_Conf      = zeros(nSubjs,3);
RTs_Scenes_Conf     = zeros(nSubjs,3);

for ss = 1:nSubjs
    s = subjs(ss);
    RTs = behav_out.retSubj{s}.RTs;
    
    % confidence
    face_trials     = behav_out.retSubj{s}.FaceTrials;
    scene_trials    = behav_out.retSubj{s}.SceneTrials;
    for cc = 1:3
        trials = behav_out.retSubj{s}.Confidence==cc;
        RTs_Faces_Conf(ss,cc)   = median(RTs(trials & face_trials));
        RTs_Scenes_Conf(ss,cc)  = median(RTs(trials & scene_trials));
    end
end
%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--
%Faces
figure(3); clf;
AR = [500 300];
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',AR,'paperposition',[0 0 AR])
set(gcf,'position',[50,500,AR])
hold on;
for ss = 1:nSubjs
    p=plot(1:3,RTs_Faces_Conf(ss,:),'linewidth',1,'color',[100 200 100]/255);
end
plot(1:3,nanmean(RTs_Faces_Conf),'linewidth',4,'color',[0.1 0.1 0.1]);
set(gca,'fontsize',30,'xTick',1:3,'xticklabel',{'Low','Med','High'})
xlim([0.8 3.2])
ylim([0.5 2.1])
set(gca,'LineWidth',3,'ytick',[0.5:0.5:2])
ylabel(' RTs (s) ')
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/BehavRetrieval/' 'RTs_Faces_Confidence_' SubjSelectStr])
% Stats
disp('Face dPrime by Confidence ')
disp(array2table(nanmean(RTs_Faces_Conf),'variablenames',{'Low','Med','High'}))
% Statistics
fprintf('Model for dP Confidence Faces: \n')
y=-log10(RTs_Faces_Conf);
t = array2table(y);
within = dataset([1:3]','varnames',{'confidence'});
rm = fitrm(t,'y1-y3~1','withindesign',within);
ranova(rm,'withinmodel','confidence-1')
%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--
% Scenes
figure(4); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[AR],'paperposition',[0 0 AR])
set(gcf,'position',[50,500,AR])
hold on;
for ss = 1:nSubjs
    p=plot(1:3,RTs_Scenes_Conf(ss,:),'linewidth',1,'color',[200 100 200]/255);
end
plot(1:3,nanmean(RTs_Scenes_Conf),'linewidth',4,'color',[0.1 0.1 0.1]);
set(gca,'fontsize',30,'xTick',1:3,'xticklabel',{'Low','Med','High'})
xlim([0.8 3.2])
ylim([0.5 2.1])
set(gca,'LineWidth',3,'ytick',[0:0.5:3])
ylabel(' RTs(s) ')
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/BehavRetrieval/' 'RTs_Scenes_Confidence_' SubjSelectStr])

% Stats
disp('Scene dPrime by Confidence ')
disp(array2table(nanmean(RTs_Scenes_Conf),'variablenames',{'Low','Med','High'}))
fprintf('Model for dP Confidence Faces: \n')
y=-log10(RTs_Scenes_Conf);
t = array2table(y);
within = dataset([1:3]','varnames',{'confidence'});
rm = fitrm(t,'y1-y3~1','withindesign',within);
ranova(rm,'withinmodel','confidence-1')

%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--%--
% Interaction Stats
fprintf('Model for dP Confidence X Face/Scenes: \n')
y1=-log10(RTs_Faces_Conf);
y2=-log10(RTs_Scenes_Conf);
t = array2table([y1 y2]);
within = dataset([1:3,1:3]',['F','F','F','S','S','S']','varnames',{'confidence','category'});
rm = fitrm(t,'Var1-Var6~1','withindesign',within);
ranova(rm,'withinmodel','confidence*category-1')

%% RTs by confidence/category and memory
RTs_Faces       = zeros(nSubjs,3,4);
nFace_trials    =zeros(nSubjs,3,4);
RTs_Scenes      = zeros(nSubjs,3,4);
nScene_trials    =zeros(nSubjs,3,4);
memConds        = {'Hits','Misses','FA','CRs'};
for ss = 1:nSubjs
    s = subjs(ss);
    % confidence
    RTs             = behav_out.retSubj{s}.RTs;
    face_trials     = behav_out.retSubj{s}.FaceTrials;
    scene_trials    = behav_out.retSubj{s}.SceneTrials;
    for cc = 1:3
        % memory
        conf_trials = behav_out.retSubj{s}.Confidence==cc;
        for mm =1:4
            memCond_trials          = behav_out.retSubj{s}.(memConds{mm});
            % Faces
            trials          = face_trials & conf_trials & memCond_trials;
            nFace_trials(ss,cc,mm)  = sum(trials);
            RTs_Faces(ss,cc,mm)     = median(RTs(trials));
            % Scenes
            trials          = scene_trials & conf_trials & memCond_trials;
            nScene_trials(ss,cc,mm) = sum(trials);
            RTs_Scenes(ss,cc,mm)    = median(RTs(trials));
        end
    end
end

% Figure:
figure(1); clf; AR = [550 300];
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',AR,'paperposition',[0 0 AR])
set(gcf,'position',[100,100,AR])

dx = 120;
xPos = [80:dx:550];
a     = zeros(4,1);
for ii =1:nRTConds
    a(ii)=axes('units','points','position',[xPos(ii) 100 100 150]);
end
yLims = [0.5 2.5];
yTicks1 = [0.5:0.5:2.5];

strs = {'Faces','Scenes'};
% RTs for scenes/faces:
Y = cat(4,RTs_Faces,RTs_Scenes);
Colors = [100 200 100; 200 100 200]/255;
for jj = 1:2
    for ii=1:4
        axes(a(ii)); cla;
        hold on;        
        Yjj = Y(:,:,ii,jj);
        for ss=1:nSubjs
            p = plot(1:3,Yjj(ss,:));
            p.Color = Colors(jj,:);
            p.LineWidth = 1;
        end
        set(gca,'fontsize',20,'xTick',1:3, 'xticklabel',{'L','M','H'} )
        plot(1:3,nanmean(Yjj),'linewidth',4,'color','k')
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
    print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/BehavRetrieval/' 'RTs_' strs{jj} '_Conf_Mem' SubjSelectStr])
end

% Statistics for Hits: Interaction Stats
fprintf('Model for RTs Confidence X Face/Scenes: \n')
y =  -log([RTs_Faces(:,:,1) RTs_Scenes(:,:,1)]);
t = array2table(y);
within = dataset([1:3,1:3]',['F','F','F','S','S','S']','varnames',{'confidence','category'});
rm = fitrm(t,'y1-y6~1','withindesign',within);
ranova(rm,'withinmodel','confidence*category-1')

y1 = -log(RTs_Faces); y2 = -log(RTs_Scenes);
% Paired test between RTs of hits for scenes and faces by confidence
[~,p,~,t]=ttest(y1(:,:,1) ,y2(:,:,1));
for ii=1:3    
    fprintf(['RTs Confidence Level %i; Faces %0.2f, Scenes %0.2f\n' ...
    ' Face vs Scene: t=%0.2f, p = %0.2g \n\n'],ii, mean(RTs_Faces(:,ii,1)), ...
    mean(RTs_Scenes(:,ii,1)),t.tstat(ii), p(ii) )
end
% Paired test between LC/HC hits for scenes and for faces
[~,p,~,t]=ttest(y1(:,3,1) ,y1(:,1,1));
fprintf(' RTs Faces HC vs LC:  t=%0.2f, p = %0.2g \n\n',t.tstat,p)
[~,p,~,t]=ttest(y2(:,3,1) ,y2(:,1,1));
fprintf(' RTs Scenes HC vs LC:  t=%0.2f, p = %0.2g \n\n',t.tstat,p)

%% nTrials for Scenes and Faces
Y = cat(4,nFace_trials,nScene_trials);
yLims = [0 95];
yTicks1 = [0: 25: 100]
for jj = 1:2
    for ii=1:4
        axes(a(ii)); cla;
        hold on;        
        Yjj = Y(:,:,ii,jj);
        for ss=1:nSubjs
            p = plot(1:3,Yjj(ss,:));
            p.Color = Colors(jj,:);
            p.LineWidth = 1;
        end
        set(gca,'fontsize',20,'xTick',1:3, 'xticklabel',{'L','M','H'} )
        plot(1:3,nanmean(Yjj),'linewidth',4,'color','k')
        xlabel(strs2{ii})
        ylim(yLims)
        xlim([0.8 3.2])
        if ii==1
            set(gca,'ytick',yTicks1)
            ylabel(' # Trials ')
        end
        if ii>1
            set(gca,'ycolor','none')
        end
        set(gca,'LineWidth',2)
    end
    print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/BehavRetrieval/' 'nTrials_' strs{jj} '_Conf_Mem' SubjSelectStr])
end

% Hit Statistics in number of trials.
y = [nFace_trials(:,:,1) nScene_trials(:,:,1)];
t = array2table(y);
within = dataset([1:3,1:3]',['F','F','F','S','S','S']','varnames',{'confidence','category'});
rm = fitrm(t,'y1-y6~1','withindesign',within);
ranova(rm,'withinmodel','confidence*category-1')

% Paired test between # of hits for scenes and faces
[~,p,~,t]=ttest(sum(nFace_trials(:,:,1),2) ,sum(nScene_trials(:,:,1),2));
fprintf('#Trials Face vs Scene t=%0.2f, p=%0.2g \n',t.tstat,p)

[~,p,~,t]=ttest(nScene_trials(:,3,1), nScene_trials(:,1,1));
fprintf('#Trials HCH Scene vs LCH Scene t=%0.2f, p=%0.2g \n',t.tstat,p)

[~,p,~,t]=ttest(nFace_trials(:,3,1), nFace_trials(:,1,1));
fprintf('#Trials HCH FAce vs LCH FAce t=%0.2f, p=%0.2g \n',t.tstat,p)
