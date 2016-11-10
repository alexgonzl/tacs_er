% Analyses of Retrieval behavior

close all
clearvars

dataPath    = '~/Google Drive/Research/tACS/tACS_ER_task/data/tacs_er_objstim/';
plotPath    = '~/Google Drive/Research/tACS/tACS_ER_task/plots/tacs_er_objstim/BehavRetrieval/';
load([dataPath 'Summary/BehavSummary.mat'])
load([dataPath 'Summary/DataMatrix.mat'])

% Exclusion of subjects based on poor behavioral performance
badSubjs =  [3 6];

nTotalSubjs      = 22;
subjs           = setdiff(1:nTotalSubjs,badSubjs);
nSubjs          = numel(subjs);

f = @(th)(mean(exp(1j*th)));
f2 = @(th,mag)(mean(mag.*exp(1j*th)));
temp= brewermap(6,'accent');

% data colors
SmallColor  = temp(5,:);
BigColor    = temp(6,:);
greyColor   = [119,136,153]/255;

hitColor  = [255 180 150]/255;
missColor = [150 220 220]/255;

temp=brewermap(2,'set1');
HCcolor = temp(1,:);
LCcolor = temp(2,:);

%% dPrimes
rng(1); % for location reproducibility

Dstrs   = {'dPrime','Small_dPrime','Big_dPrime',...
    'dPrime_C','Small_dPrime_C','Big_dPrime_C'};

Dstrs2  = {'dP', 'Small dP', 'Big dP', 'C', 'Small C', 'Big C'};
Dstrs3  = {'dP', 'SmallDP', 'BigDP', 'C', 'SmallC', 'BigC'};
nD      = numel(Dstrs);
D       = zeros(nSubjs,nD);
for ii = 1:nD
    D(:,ii) = behav_out.retSummary.(Dstrs{ii})(subjs);
end
disp(array2table(mean(D),'variablenames',Dstrs3))
[~,p,~,t]=ttest(D(:,2),D(:,3));
disp(table(t.tstat,p,'rownames',{'Small vs Big DP'}))

yLims = [-0.1 2.5];
yTicks1 = [0 1 2];

% Figure 1: (a) d-prime
figure(1); clf;
AR = [450 300];
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',AR,'paperposition',[0 0 AR])
set(gcf,'position',[200,200,AR])

a = axes('units','points','position',[50 50 100 200]); hold on;
x = randn(nSubjs,1)*0.05+0.5;
y = D(:,1);
plot([0.2 0.8], ones(1,2)*mean(y),'linewidth',4,'color','k')
s=scatter(x,y);

s.MarkerFaceAlpha   = 0.5;
s.MarkerEdgeAlpha   = 0.4;
s.SizeData          = 120;
s.MarkerEdgeColor   = greyColor;
s.MarkerFaceColor   = greyColor;

set(gca,'fontsize',20,'xtick',[])
ylim(yLims)
xlim([0 1])
set(gca,'ytick',yTicks1)
ylabel(' d'' ')
set(gca,'LineWidth',2)

% Figure 1: (b) d-prime by category
a = axes('units','points','position',[220 50 200 200]); hold on;
y = D(:,2:3);

for ii =1:nSubjs    
    plot([x(ii) x(ii)+1], y(ii,:),'-','color',[0.6 0.6 0.6])    
end

plot([0.2 0.8], ones(1,2)*mean(y(:,1)),'linewidth',4,'color','k')
plot([0.2 0.8]+1, ones(1,2)*mean(y(:,2)),'linewidth',4,'color','k')

% Small
s=scatter(x,y(:,1));
s.MarkerFaceAlpha=0.5;
s.MarkerEdgeAlpha=0.4;
s.SizeData          = 120;
s.MarkerEdgeColor = SmallColor;
s.MarkerFaceColor = SmallColor;

% Big
s=scatter(x+1,y(:,2));
s.MarkerFaceAlpha   = 0.5;
s.MarkerEdgeAlpha   = 0.4;
s.SizeData          = 120;
s.MarkerEdgeColor   = BigColor;
s.MarkerFaceColor   = BigColor;

xlim([-0.1 2.1])
ylim(ylim)
set(gca,'fontsize',20,'xtick',[0.5 1.5],'xticklabel',{'Small','Big'})
set(gca,'ytick',[0 1 2])
ylabel(' d'' ')
set(gca,'LineWidth',2)
print(gcf,'-dpdf',[plotPath 'dPrimes'])

m = mean(D(:,1));
[~,p,~,t] = ttest(D(:,1));
fprintf('Across subject overall dprime = %0.2f, t=%0.2f, p=%0.2g\n',m,t.tstat,p)
%% dP by Confidence

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
ylim([-1 4])
set(gca,'LineWidth',3,'ytick',[0:3])
ylabel(' d'' ')

print(gcf,'-dpdf',[plotPath  'Confidence_dPrime'])


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
close all

ylims = [-1.5 4.1];
% Small
DPC_Sm = behav_out.retSummary.Small_dPrimeConf(subjs,:);
figure(3); clf;
AR = [500 300;];
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',AR,'paperposition',[0 0 AR])
set(gcf,'position',[50,500,AR])
hold on;
for ss = 1:nSubjs    
    p=plot(1:3,DPC_Sm(ss,:),'linewidth',1,'color',SmallColor);
end
plot(1:3,nanmean(DPC_Sm),'linewidth',4,'color',[0.1 0.1 0.1]);
set(gca,'fontsize',30,'xTick',1:3,'xticklabel',{'Low','Med','High'})
xlim([0.8 3.2])
ylim(ylims)
set(gca,'LineWidth',3,'ytick',[-2:2:6])
ylabel(' d'' ')

print(gcf,'-dpdf',[plotPath 'Confidence_dPrimeSmall'])
disp('Small dPrime by Confidence ')
disp(array2table(nanmean(DPC_Sm),'variablenames',{'Low','Med','High'}))

% Big
DPC_Big = behav_out.retSummary.Big_dPrimeConf(subjs,:);
figure(4); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[AR],'paperposition',[0 0 AR])
set(gcf,'position',[50,500,AR])
hold on;
for ss = 1:nSubjs
    p=plot(1:3,DPC_Big(ss,:),'linewidth',1,'color',BigColor);
end
plot(1:3,nanmean(DPC_Big),'linewidth',4,'color',[0.1 0.1 0.1]);
set(gca,'fontsize',30,'xTick',1:3,'xticklabel',{'Low','Med','High'})
xlim([0.8 3.2])
ylim(ylims)
set(gca,'LineWidth',3,'ytick',[-2:2:6])
ylabel(' d'' ')
print(gcf,'-dpdf',[plotPath 'Confidence_dPrimeBig'])

disp('Big dPrime by Confidence ')
disp(array2table(nanmean(DPC_Big),'variablenames',{'Low','Med','High'}))

% Statistics
fprintf('Model for dP Confidence Small: \n')
t = array2table(DPC_Sm);
within = dataset([1:3]','varnames',{'confidence'});
rm = fitrm(t,'DPC_Sm1-DPC_Sm3~1','withindesign',within);
ranova(rm,'withinmodel','confidence-1')

fprintf('Model for dP Confidence Big: \n')
t = array2table(DPC_Big);
within = dataset([1:3]','varnames',{'confidence'});
rm = fitrm(t,'DPC_Big1-DPC_Big3~1','withindesign',within);
ranova(rm,'withinmodel','confidence-1')

fprintf('Model for dP Confidence X Small/Big: \n')
t = array2table([DPC_Sm DPC_Big]);
within = dataset([1:3,1:3]',['S','S','S','B','B','B']','varnames',{'confidence','category'});
rm = fitrm(t,'Var1-Var6~1','withindesign',within);
ranova(rm,'withinmodel','confidence*category-1')
% t=array2table(DPC);
% conf = dataset([1:3]','varnames',{'confidence'});
% rm = fitrm(t,'DPC1-DPC3~1','withindesign',conf);
% ranova(rm)

%

%% RTs
rng(1)
strs    = {'medianHit_RTs','medianMiss_RTs','medianFA_RTs','medianCRs_RTs'};
strs2   = {'Hits','Misses','FA','CRs'};
nRTConds = numel(strs);
RTs       = zeros(nSubjs,nRTConds);
for ii = 1:nRTConds
    RTs(:,ii) = behav_out.retSummary.(strs{ii})(subjs);
end
figure(1); clf;
AR = [450 300];
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',AR,'paperposition',[0 0 AR])
set(gcf,'position',[200,200,AR])

dx = 80;
xPosCore = [100:dx:500];
xPos = [xPosCore ];
a     = zeros(nRTConds,1);
for ii =1:nRTConds
    a(ii)=axes('units','points','position',[xPos(ii) 80 60 160]);
end
yLims = [0.5 2.1];
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

print(gcf,'-dpdf',[plotPath 'RetRTs'])

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

%% Retrieval RTs by Stimulus Category
rng(1)
strs    = {'medianSmallBigHit_RTs','medianSmallBigMiss_RTs',...
    'medianSmallBigFA_RTs','medianSmallBigCRs_RTs'};
strs2   = {'Hits','Misses','FA','CRs'};
nRTConds = numel(strs);

RTs_Small       = zeros(nSubjs,nRTConds);
RTs_Big       = zeros(nSubjs,nRTConds);
for ii = 1:nRTConds
    RTs_Small(:,ii) = behav_out.retSummary.(strs{ii})(subjs,1);
    RTs_Big(:,ii) = behav_out.retSummary.(strs{ii})(subjs,2);
end
figure(1); clf;
AR = [450 300];
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',AR,'paperposition',[0 0 AR])
set(gcf,'position',[200,200,AR])

dx = 90;
xPosCore = [80:dx:500];
xPos = [xPosCore ];
a     = zeros(nRTConds,1);
for ii =1:nRTConds
    a(ii)=axes('units','points','position',[xPos(ii) 75 75 160]);
end
yLims = [0.5 2.5];
yTicks1 = [0:0.5:2.5];
x1 = randn(nSubjs,1)*0.03+0.3;
x2 = randn(nSubjs,1)*0.03+0.7;
for ii=1:nRTConds
    axes(a(ii)); hold on;
    y1  = RTs_Small(:,ii);
    y2  = RTs_Big(:,ii);
    
    for jj =1:nSubjs
        plot([x1(jj) x2(jj)], [y1(jj) y2(jj)],'-','color',[0.6 0.6 0.6])
    end
    
    % Small
    s   = scatter(x1,y1);
    plot([0.2 0.4], ones(1,2)*mean(y1),'linewidth',4,'color','k');
    s.MarkerFaceAlpha   = 0.5;
    s.MarkerEdgeAlpha   = 0.4;
    s.SizeData          = 120;
    s.MarkerEdgeColor = SmallColor;
    s.MarkerFaceColor = SmallColor;
    
    set(gca,'fontsize',20,'xTick','')
    xlabel(strs2{ii})
    ylim(yLims)
    xlim([0 1])
    
    % Big
    s   = scatter(x2,y2);
    plot([0.6 0.8], ones(1,2)*mean(y2),'linewidth',4,'color','k')
    s.MarkerFaceAlpha   = 0.5;
    s.MarkerEdgeAlpha   = 0.4;
    s.SizeData          = 120;
    s.MarkerEdgeColor   = BigColor;
    s.MarkerFaceColor   = BigColor;
    
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

print(gcf,'-dpdf',[plotPath 'RetRTsCategory'])

% Statistics
fprintf('Model for RTs Category X Memory Conidtion: \n')
y1 = -log(RTs_Small);
y2 = -log(RTs_Big);
t = array2table([y1 y2]);
within = dataset([1:4,1:4]',['S','S','S','S','B','B','B','B',]','varnames',{'memory','category'});
rm = fitrm(t,'Var1-Var8~1','withindesign',within);
ranova(rm,'withinmodel','memory*category-1')

% statistics
disp(array2table(mean([y1 y2]),'variablenames',[strcat('S_',strs2) strcat('B_',strs2)]))
for ii = 1:4
    [~,p,~,t]=ttest(y1(:,ii),y2(:,ii));
    fprintf('RTs %s Small vs Big t=%0.2f, p=%0.2g \n',strs2{ii},t.tstat,p)
end

%% RTs by Confidence and memory
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

print(gcf,'-dpdf',[plotPath 'RTsByConf'])

fprintf('Model for RTs Confidence X Memory Conidtion: \n')
% Statistics
y = -log(RTsConf(:,:));
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
    yLims = [0 250];
    yTicks1 = [0:50:250];
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
print(gcf,'-dpdf',[plotPath 'nRespConf'])
%statistcs for Hits%


figure(2); clf;
AR = [550 300];
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[AR],'paperposition',[0 0 AR])
set(gcf,'position',[0,0,AR])
for ii =1:4
    a(ii)=axes('units','points','position',[xPos(ii) 100 100 150]);
end
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
print(gcf,'-dpdf',[plotPath 'propRespConf'])

%% RTs by confidence/category and memory
RTs_Small       = zeros(nSubjs,3,4);
nSmall_trials    =zeros(nSubjs,3,4);
RTs_Big      = zeros(nSubjs,3,4);
nBig_trials    =zeros(nSubjs,3,4);
memConds        = {'Hits','Misses','FA','CRs'};
for ss = 1:nSubjs
    s = subjs(ss);
    % confidence
    RTs             = behav_out.retSubj{s}.RTs;
    small_trials    = behav_out.retSubj{s}.SmallTrials ;
    big_trials      = behav_out.retSubj{s}.BigTrials;
    for cc = 1:3
        % memory
        conf_trials = behav_out.retSubj{s}.Confidence==cc;
        for mm =1:4
            memCond_trials          = behav_out.retSubj{s}.(memConds{mm});
            % Small
            trials                  = small_trials & conf_trials & memCond_trials;
            nSmall_trials(ss,cc,mm)  = sum(trials);
            RTs_Small(ss,cc,mm)     = nanmedian(RTs(trials));
            % Big
            trials                  = big_trials & conf_trials & memCond_trials;
            nBig_trials(ss,cc,mm) = sum(trials);
            RTs_Big(ss,cc,mm)       = nanmedian(RTs(trials));
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
Y = cat(4,RTs_Small,RTs_Big);
Colors = [SmallColor; BigColor];
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
    print(gcf,'-dpdf',[plotPath 'RTs_' strs{jj} '_Conf_Mem'])
end

% Statistics for Hits: Interaction Stats
fprintf('Model for RTs Confidence X Small/Big: \n')
y =  -log([RTs_Small(:,:,1) RTs_Big(:,:,1)]);
t = array2table(y);
within = dataset([1:3,1:3]',['S','S','S','B','B','B',]','varnames',{'confidence','category'});
rm = fitrm(t,'y1-y6~1','withindesign',within);
ranova(rm,'withinmodel','confidence*category-1')

y1 = -log(RTs_Small); y2 = -log(RTs_Big);
% Paired test between RTs of hits for scenes and faces by confidence
[~,p,~,t]=ttest(y1(:,:,1) ,y2(:,:,1));
for ii=1:3    
    fprintf(['RTs Confidence Level %i; Small %0.2f, Big %0.2f\n' ...
    ' Small vs Big: t=%0.2f, p = %0.2g \n\n'],ii, nanmean(RTs_Small(:,ii,1)), ...
    nanmean(RTs_Big(:,ii,1)),t.tstat(ii), p(ii) )
end
% Paired test between LC/HC hits for scenes and for faces
[~,p,~,t]=ttest(y1(:,3,1) ,y1(:,1,1));
fprintf(' RTs Small HC vs LC:  t=%0.2f, p = %0.2g \n\n',t.tstat,p)
[~,p,~,t]=ttest(y2(:,3,1) ,y2(:,1,1));
fprintf(' RTs Big HC vs LC:  t=%0.2f, p = %0.2g \n\n',t.tstat,p)

%% nTrials for Scenes and Faces
Y = cat(4,nSmall_trials,nBig_trials);
yLims = [0 150];
yTicks1 = [0: 50: 200];
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
    print(gcf,'-dpdf',[plotPath 'nTrials_' strs{jj} '_Conf_Mem'])
end

% Hit Statistics in number of trials.
y = [nSmall_trials(:,:,1) nBig_trials(:,:,1)];
t = array2table(y);
within = dataset([1:3,1:3]',['S','S','S','B','B','B']','varnames',{'confidence','category'});
rm = fitrm(t,'y1-y6~1','withindesign',within);
ranova(rm,'withinmodel','confidence*category-1')

% Paired test between # of hits for scenes and faces
[~,p,~,t]=ttest(sum(nSmall_trials(:,:,1),2) ,sum(nBig_trials(:,:,1),2));
fprintf('#Trials Small vs Big t=%0.2f, p=%0.2g \n',t.tstat,p)

[~,p,~,t]=ttest(nBig_trials(:,3,1), nBig_trials(:,1,1));
fprintf('#Trials HCH Big vs LCH Big t=%0.2f, p=%0.2g \n',t.tstat,p)

[~,p,~,t]=ttest(nSmall_trials(:,3,1), nSmall_trials(:,1,1));
fprintf('#Trials HCH Small vs LCH Small t=%0.2f, p=%0.2g \n',t.tstat,p)
