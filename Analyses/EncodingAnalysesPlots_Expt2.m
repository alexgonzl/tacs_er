% Analyses for Encoding Data.

close all
clearvars

dataPath    = '~/Google Drive/Research/tACS/tACS_ER_task/data/tacs_er_objstim/';
plotPath    = '~/Google Drive/Research/tACS/tACS_ER_task/plots/tacs_er_objstim/Encoding/';
load([dataPath 'Summary/BehavSummary.mat'])
load([dataPath 'Summary/DataMatrix.mat'])

% Exclusion of subjects based on poor behavioral performance or no
% entrainment
badSubjs =  [3 6 8 27 29 32 34];

nTotalSubjs      = 37;
subjs           = setdiff(1:nTotalSubjs,badSubjs);
nSubjs          = numel(subjs);

f = @(th)(nansum(exp(1j*th)));
f2 = @(th,mag)(nanmean(mag.*exp(1j*th)));
temp= brewermap(6,'accent');
SmallColor  = temp(5,:);
BigColor    = temp(6,:);
greyColor = [119,136,153]/255;

PhasesCol       = 4;
CorrectRespCol  = 2;
StimTypeCol     = 1; 
EncRTsCol       = 3;
%% Encoding Results: Categorization Task
rng(1); % for location reproducibility
HR = [ behav_out.encSummary.meanAcc(subjs) behav_out.encSummary.SmallHR(subjs) behav_out.encSummary.BigHR(subjs)]*100;
Dstrs = {'ACC','Small','Big'};
disp(array2table(mean(HR),'variablenames',Dstrs))
[~,p,~,t]=ttest(HR(:,2),HR(:,3));
disp(table(t.tstat,p,'rownames',{'Small vs Big ACC'},'variablenames',{'T','P'}))

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
ylim([60 100])
xlim([0 1])
set(gca,'ytick',[60:20:100])
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

% Small
s=scatter(x,y(:,1));
s.MarkerFaceAlpha=0.5;
s.MarkerEdgeAlpha=0;
s.SizeData          = 120;
s.MarkerEdgeColor = SmallColor;
s.MarkerFaceColor = SmallColor;

% Big
s=scatter(x+1,y(:,2));
s.MarkerFaceAlpha   = 0.5;
s.MarkerEdgeAlpha   = 0;
s.SizeData          = 120;
s.MarkerEdgeColor   = BigColor;
s.MarkerFaceColor   = BigColor;

xlim([-0.1 2.1])
ylim([60 100])
set(gca,'fontsize',20,'xtick',[0.5 1.5],'xticklabel',{'Small','Big'})
set(gca,'ytick',[60:20:100])
set(gca,'LineWidth',2)
print(gcf,'-dpdf',[plotPath 'CategorizationPerf'])

%% Categorization by Phase
rng(1); % for location reproducibility

% Modulation by Correct Responses
Z =zeros(nSubjs,1);
for ss = 1:nSubjs
    s = subjs(ss);
    trials = out.datMat(s,:,CorrectRespCol)==1;
    phases = out.datMat(s,:,PhasesCol);
    Z(ss) = f(phases(trials));
    [~,R(ss)] = circ_rtest(phases(trials));
end

th  = angle(Z);
rho = abs(Z);

opts = [];
opts.colors = greyColor;
%opts.maxR = 4/3;
opts.alpha = 0.8;
opts.markerSize=300;
%opts.magText    = '';
PolarPlot(th,rho,opts)
print(gcf,'-dpdf',[plotPath 'CategorizationPerfByPhasePolarRho'])

opts = [];
opts.colors = greyColor;
opts.maxR = 4/3;
opts.alpha = 0.8;
opts.markerSize=300;
opts.magText    = '';
PolarPlot(th,[],opts)
print(gcf,'-dpdf',[plotPath 'CategorizationPerfByPhasePolar'])

[p,r] = circ_rtest(th);
fprintf('Rayleight Test for modulation across subjects:\n')
fprintf('p=%g ; r = %g \n',p,r)
circ_mean(th)*180/pi
%% Categorization By Phase / Stim Type
close all; rng(1); % for location reproducibility

% Modulation by Correct Responses
Z1 =zeros(nSubjs,1);
Z2   =zeros(nSubjs,1);
for ss = 1:nSubjs
    s       = subjs(ss);
    phases  = out.datMat(s,:,PhasesCol);
    
    trials  = (out.datMat(s,:,CorrectRespCol)==1) & (out.datMat(s,:,StimTypeCol)==1);    
    Z1(ss)   = f(phases(trials));
    
    trials  = (out.datMat(s,:,CorrectRespCol)==1) & (out.datMat(s,:,StimTypeCol)==2);    
    Z2(ss)   = f(phases(trials));
end

% Polar Plots
th1 = angle(Z1); rho1 = abs(Z1);
th2 = angle(Z2); rho2 = abs(Z2);

% Polar w magnitude
opts = [];
opts.colors = [SmallColor; BigColor];
opts.alpha = 0.4;
opts.markerSize=300;
PolarPlot([th1 th2],[rho1 rho2],opts)
print(gcf,'-dpdf',[plotPath 'EncPerfByCategoryPolarRho'])

% only angle
opts = [];
opts.colors = [SmallColor; BigColor];
opts.maxR = 4/3;
opts.alpha = 0.8;
opts.markerSize=300;
opts.magText    = '';
PolarPlot([th1 th2],[],opts)
print(gcf,'-dpdf',[plotPath 'EncPerfByCategoryPolar'])

% difference
Z3 = Z2-Z1;
th3 = angle(Z3); rho3 = abs(Z3);
% Polar w magnitude
opts = [];
opts.colors = greyColor;
opts.alpha = 0.8;
opts.markerSize=300;
PolarPlot(th3,rho3,opts)
print(gcf,'-dpdf',[plotPath 'EncPerfCatDiffPolarRho'])

% only angle
opts = [];
opts.colors = greyColor;
opts.maxR = 4/3;
opts.alpha = 0.8;
opts.markerSize=300;
opts.magText    = '';
PolarPlot(th3,[],opts)
print(gcf,'-dpdf',[plotPath 'EncPerfCatDiffPolar'])

[p,r] = circ_rtest(th1);
fprintf('Rayleight Test for modulation across subjects for Small:\n')
th_hat = circ_mean(th1)/pi*180;
fprintf('angle = %g; r = %g; p=%g  \n',th_hat,r,p)

[p,r] = circ_rtest(th2);
fprintf('Rayleight Test for modulation across subjects for Big:\n')
th_hat = circ_mean(th2)/pi*180;
fprintf('angle = %g; r = %g; p=%g  \n',th_hat,r,p)

[p,r] = circ_rtest(th3);
fprintf('Rayleight Test for modulation across subjects for Big-Sma;;:\n')
th_hat = circ_mean(th3)/pi*180;
fprintf('angle = %g; r = %g; p=%g  \n',th_hat,r,p)

%% make legend
colors = [SmallColor; BigColor];
txtStr = {'Small','Big'};
makeLegend(colors,[300 300],txtStr);
print(gcf,'-dpdf',[plotPath 'SmallBigLegend']);

%% Categorization RTs:
rng(1); % for location reproducibility
RTs = [behav_out.encSummary.medianRTs(subjs)  ...
    behav_out.encSummary.medianSmallRTs(subjs) behav_out.encSummary.medianBigRTs(subjs)];

Dstrs = {'RTs','Small','Big'};
disp(array2table(mean(RTs),'variablenames',Dstrs))
[~,p,~,t]=ttest(RTs(:,2),RTs(:,3));
disp(table(t.tstat,p,'rownames',{'Small vs Big RTs'},'variablenames',{'T','P'}))
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
ylim([0.4 1.21])
xlim([0 1])
set(gca,'ytick',[0.4:0.2:1.2])
ylabel(' RTs (s) ')
set(gca,'LineWidth',2)

% Figure 1: (b) acc by category
a = axes('units','points','position',[300 100 180 200]); hold on;
y = RTs(:,2:3);

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
s.MarkerEdgeAlpha   = 0;
s.SizeData          = 120;
s.MarkerEdgeColor   = BigColor;
s.MarkerFaceColor   = BigColor;

xlim([-0.1 2.1])
ylim([0.4 1.21])
set(gca,'fontsize',20,'xtick',[0.5 1.5],'xticklabel',{'Small','Big'})
set(gca,'ytick',[0.4:0.2:1.2])
set(gca,'LineWidth',2)
print(gcf,'-dpdf',[plotPath 'EncRTs'])

%% Categorization RTs by Phase
close all;
rng(1); % for location reproducibility

% Modulation by Correct Responses
Z  = zeros(nSubjs,1);
for ss = 1:nSubjs
    s       = subjs(ss);
    phases  = out.datMat(s,:,PhasesCol);
    RTs     = -log(out.datMat(s,:,EncRTsCol));
    
    trials  = out.datMat(s,:,CorrectRespCol)==1 ;   
    Z(ss)   = f2(phases(trials),RTs(trials));
end
th  = angle(Z);
rho = abs(Z);

opts = [];
opts.colors = greyColor;
opts.alpha = 0.8;
opts.markerSize=300;
PolarPlot(th,rho,opts)
print(gcf,'-dpdf',[plotPath 'EncRTsPolarRho'])

opts = [];
opts.colors = greyColor;
opts.maxR = 4/3;
opts.alpha = 0.8;
opts.markerSize=300;
opts.magText    = '';
PolarPlot(th,[],opts)
print(gcf,'-dpdf',[plotPath 'EncRTsPolar'])

[p,r] = circ_rtest(th);
th_hat = circ_mean(th);
fprintf('Rayleight Test for Categorization RT modulation across subjects:\n')
fprintf('angle = %g;  r = %g; p=%g \n',th_hat,r,p)

%% Categorization RTs by Phase and category
close all; rng(1); % for location reproducibility

StimTypeCol     = 1;  
CorrectRespCol  = 2;
EncRTsCol       = 3;
PhasesCol       = 4;

% Modulation by Correct Responses
Z1  = zeros(nSubjs,1);
Z2   =zeros(nSubjs,1);
for ss = 1:nSubjs
    s       = subjs(ss);
    phases  = out.datMat(s,:,PhasesCol);
    RTs     = -log(out.datMat(s,:,EncRTsCol));
    
    trials  = (out.datMat(s,:,CorrectRespCol)==1) & (out.datMat(s,:,StimTypeCol)==1);    
    Z1(ss)   = f2(phases(trials),RTs(trials));
    
    trials  = (out.datMat(s,:,CorrectRespCol)==1) & (out.datMat(s,:,StimTypeCol)==2);    
    Z2(ss)   = f2(phases(trials),RTs(trials));
end

% Polar Plots
th1 = angle(Z1); rho1 = abs(Z1);
th2 = angle(Z2); rho2 = abs(Z2);

% Polar w magnitude
opts = [];
opts.colors = [SmallColor; BigColor];
opts.alpha = 0.4;
opts.markerSize=300;
PolarPlot([th1 th2],[rho1 rho2],opts)
print(gcf,'-dpdf',[plotPath 'EncRTsByCategoryPolarRho'])

% only angle
opts = [];
opts.colors = [SmallColor; BigColor];
opts.maxR = 4/3;
opts.alpha = 0.8;
opts.markerSize=300;
opts.magText    = '';
PolarPlot([th1 th2],[],opts)
print(gcf,'-dpdf',[plotPath 'EncRTsByCategoryPolar'])

% difference
Z3 = Z2-Z1;
th3 = angle(Z3); rho3 = abs(Z3);
% Polar w magnitude
opts = [];
opts.colors = greyColor;
opts.alpha = 0.8;
opts.markerSize=300;
PolarPlot(th3,rho3,opts)
print(gcf,'-dpdf',[plotPath 'EncRTsCatDiffPolarRho'])

% only angle
opts = [];
opts.colors = greyColor;
opts.maxR = 4/3;
opts.alpha = 0.8;
opts.markerSize=300;
opts.magText    = '';
PolarPlot(th3,[],opts)
print(gcf,'-dpdf',[plotPath 'EncRTsCatDiffPolar'])

[p,r] = circ_rtest(th1);
fprintf('Rayleight Test for modulation across subjects for Small:\n')
th_hat = circ_mean(th1)/pi*180;
fprintf('angle = %g; r = %g; p=%g  \n',th_hat,r,p)

[p,r] = circ_rtest(th2);
fprintf('Rayleight Test for modulation across subjects for Big:\n')
th_hat = circ_mean(th2)/pi*180;
fprintf('angle = %g; r = %g; p=%g  \n',th_hat,r,p)

[p,r] = circ_rtest(th3);
fprintf('Rayleight Test for modulation across subjects for Big-Sma;;:\n')
th_hat = circ_mean(th3)/pi*180;
fprintf('angle = %g; r = %g; p=%g  \n',th_hat,r,p)

