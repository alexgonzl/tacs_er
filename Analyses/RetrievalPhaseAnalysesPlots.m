close all
clearvars
addpath Plotting/
addpath Analyses/
addpath Library/

dataPath    = '~/Google Drive/Research/tACS/tACS_ER_task/data/tacs_enc_xdiva/';
load([dataPath 'Summary/BehavSummary.mat'])
load([dataPath 'Summary/PhaseDependentAnalyses.mat'])

uniqSeqSelection = 2;
if uniqSeqSelection==1 % unique sequences
    subjs  = find(out.SubjWithUniqueSeq);
    SubjSelectStr = 'SSuniqueSeq';
elseif uniqSeqSelection==2 % unique sequences with good encoding data*****
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
thRad     = (PhasesDeg-36)/180*pi;
% nHits, nMisses by Phase
% Raw # by phase
nHitP = out.nHitsByPhase(subjs,:);
nMissP = out.nMissByPhase(subjs,:);

% Number of hits and misses by confidence
nHitP_C = out.nHitsConfByPhase(subjs,:,:);
nHitP_HC = squeeze(nHitP_C(:,3,:));
nHitP_MC = squeeze(nHitP_C(:,2,:));
nHitP_LC = squeeze(nHitP_C(:,1,:));

nMissP_C = out.nMissConfByPhase(subjs,:,:);
nMissP_HC = squeeze(nMissP_C(:,3,:));
nMissP_MC = squeeze(nMissP_C(:,2,:));
nMissP_LC = squeeze(nMissP_C(:,1,:));

nConfP_C    = nHitP_C+nMissP_C;
nConfP_HC   = squeeze(nConfP_C(:,3,:));
nConfP_MC   = squeeze(nConfP_C(:,2,:));
nConfP_LC   = squeeze(nConfP_C(:,1,:));


% Phase Matrix and conversion function
PhaseMat = repmat(exp(1j*thRad),nSubjs,1);
f = @(X)(mean(X.*PhaseMat,2));

% functions
m_th     = @(th)(angle(mean(exp(1j*th))));
m_rho    = @(th)(abs(mean(exp(1j*th))));
deltaTh = @(th1,th2)(mod(abs(m_th(th1)-m_th(th2)),2*pi));
% Colors
hitColor  = [255 180 150]/255;
missColor = [150 220 220]/255;
faceColor = [100 200 100]/255;
sceneColor =[200 100 200]/255;
greyColor = [119,136,153]/255;
greyColor2 = [0.25 0.25 0.25];
temp=brewermap(2,'set1');
HCcolor = temp(1,:);
LCcolor = temp(2,:);
%% all Hits by Phase
close all
hitColor  = [255 180 150]/255;
X = nHitP;
% Hits scatter by phase
figure(2); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[600 400],'paperposition',[0 0 600 400])
set(gcf,'position',[100,100,600,400])

axes('units','points','position',[100 150 400 200]); hold on;
x = randn(nSubjs,1)*0.05+0.5;
y = X;
nBars = 5;
for ii =1:nBars
    xx = x+(ii-1);
    yy = y(:,ii);
    plot([0.2 0.8]+(ii-1), ones(1,2)*median(yy),'linewidth',4,'color','k')
    s = scatter(xx,yy);
    s.MarkerFaceAlpha=0.8;
    s.MarkerEdgeAlpha=0.7;
    s.SizeData          = 120;
    s.MarkerEdgeColor = hitColor;
    s.MarkerFaceColor = hitColor;
end
set(gca,'fontsize',20,'xtick',[0.5:5],'xticklabel',[])
ylim([10 50])
xlim([0 5])
set(gca,'ytick',[0:10:100])
ylabel(' # Trials ')
set(gca,'LineWidth',2)

% sinwave (x-axis)
axes('units','points','position',[100 80 400 50]); hold on;
xa = linspace(0,2*pi,1000);
x  = cos(xa);
plot(xa./pi*180,x,'k','linewidth',4)
axis tight;
set(gca,'ytick',[],'ycolor','w','fontsize',20,'box','off','lineWidth',2)
set(gca,'xtick',[36:72:360])
xlabel(' Encoding Phase (deg)')
xlim([-1 361])
print(gcf, '-dpdf', ['../plots/tacs_enc_xdiva/PhaseRetrieval/' 'nHitsByPhase' SubjSelectStr]);

% Polar Plots for nHits
z   = f(X);
th  = angle(z);
rho = abs(z);

% with magnitude of modulation
opts = [];
opts.colors = hitColor;
opts.maxR   = 4;
opts.markerSize=300;
PolarPlot(th,rho,opts)
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/PhaseRetrieval/' 'nHitsPolarRho' SubjSelectStr])

opts = [];
opts.colors = hitColor;
opts.maxR = 4/3;
opts.alpha = 0.8;
opts.markerSize=300;
opts.magText    = '';
PolarPlot(th,[],opts)
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/PhaseRetrieval/' 'nHitsPolar' SubjSelectStr])

[p,r] = circ_rtest(th);
th_hat = mod(angle(mean(exp(1j*th))),2*pi)/pi*180;
fprintf('Rayleight Test for modulation on Hits across subjects:\n')
fprintf(' angle = %g, r = %g ; p = %g \n',th_hat,r,p)

%% Misses by Phase
close all
Y = nMissP;
missColor = [150 220 220]/255;

% Misses scatter by phase
figure(); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[600 400],'paperposition',[0 0 600 400])
set(gcf,'position',[100,100,600,400])

axes('units','points','position',[100 150 400 200]); hold on;
x = randn(nSubjs,1)*0.05+0.5;
y = Y;
nBars = 5;
for ii =1:nBars
    xx = x+(ii-1);
    yy = y(:,ii);
    plot([0.2 0.8]+(ii-1), ones(1,2)*median(yy),'linewidth',4,'color','k')
    s = scatter(xx,yy);
    s.MarkerFaceAlpha=0.8;
    s.MarkerEdgeAlpha=0.7;
    s.SizeData          = 120;
    s.MarkerEdgeColor = missColor;
    s.MarkerFaceColor = missColor;
end
set(gca,'fontsize',20,'xtick',[0.5:5],'xticklabel',[])
xlim([0 5])
ylim([10 50])
set(gca,'ytick',[0:10:100])
ylabel(' Hit Rate ')
set(gca,'LineWidth',2)
% sinwave (x-axis)
axes('units','points','position',[100 80 400 50]); hold on;
xa = linspace(0,2*pi,1000);
x  = cos(xa);
plot(xa./pi*180,x,'k','linewidth',4)
axis tight;
set(gca,'ytick',[],'ycolor','w','fontsize',20,'box','off','lineWidth',2)
set(gca,'xtick',[36:72:360])
xlabel(' Encoding Phase (deg)')
xlim([-1 361])
print(gcf, '-dpdf', ['../plots/tacs_enc_xdiva/PhaseRetrieval/' 'MissesByPhase' SubjSelectStr]);

% Polar Plots for Misses
z   = f(Y);
th  = angle(z);
rho = abs(z);

% with magnitude of modulation
opts = [];
opts.colors = missColor;
opts.maxR   = 4;
opts.markerSize=300;
PolarPlot(th,rho,opts)
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/PhaseRetrieval/' 'MissesPolarRho' SubjSelectStr])

opts = [];
opts.colors = missColor;
opts.maxR = 4/3;
opts.alpha = 0.8;
opts.markerSize=300;
opts.magText    = '';
PolarPlot(th,[],opts)
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/PhaseRetrieval/' 'MissesPolar' SubjSelectStr])

[p,r] = circ_rtest(th);
fprintf('Rayleight Test for Misses modulation across subjects:\n')
th_hat = mod(angle(mean(exp(1j*th))),2*pi)/pi*180;
fprintf(' angle = %g, r = %g ; p = %g \n',th_hat,r,p)

%% HitRate by Phase
close all
Z = nHitP./(nHitP+nMissP);
greyColor = [119,136,153]/255;

% Hit Rate scatter by phase
figure(3); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[600 400],'paperposition',[0 0 600 400])
set(gcf,'position',[100,100,600,400])

axes('units','points','position',[100 150 400 200]); hold on;
x = randn(nSubjs,1)*0.05+0.5;
y = Z;
nBars = 5;
for ii =1:nBars
    xx = x+(ii-1);
    yy = y(:,ii);
    plot([0.2 0.8]+(ii-1), ones(1,2)*median(yy),'linewidth',4,'color','k')
    s = scatter(xx,yy);
    s.MarkerFaceAlpha=0.8;
    s.MarkerEdgeAlpha=0.7;
    s.SizeData          = 120;
    s.MarkerEdgeColor = greyColor;
    s.MarkerFaceColor = greyColor;
end
set(gca,'fontsize',20,'xtick',[0.5:5],'xticklabel',[])
xlim([0 5])
ylim([0.2 0.9])
set(gca,'ytick',[0:0.2:1])
ylabel(' Hit Rate ')
set(gca,'LineWidth',2)
% sinwave (x-axis)
axes('units','points','position',[100 80 400 50]); hold on;
xa = linspace(0,2*pi,1000);
x  = cos(xa);
plot(xa./pi*180,x,'k','linewidth',4)
axis tight;
set(gca,'ytick',[],'ycolor','w','fontsize',20,'box','off','lineWidth',2)
set(gca,'xtick',[36:72:360])
xlabel(' Encoding Phase (deg)')
xlim([-1 361])
print(gcf, '-dpdf', ['../plots/tacs_enc_xdiva/PhaseRetrieval/' 'HitRateByPhase' SubjSelectStr]);

% Polar Plots for Hit Rate
z   = f(Z);
th  = angle(z);
rho = abs(z);

% with magnitude of modulation
opts = [];
opts.colors = greyColor;
opts.markerSize=300;
PolarPlot(th,rho,opts)
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/PhaseRetrieval/' 'HitRatePolarRho' SubjSelectStr])

opts = [];
opts.colors = greyColor;
opts.maxR = 4/3;
opts.alpha = 0.8;
opts.markerSize=300;
opts.magText    = '';
PolarPlot(th,[],opts)
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/PhaseRetrieval/' 'HitRatePolar' SubjSelectStr])

[p,r] = circ_rtest(th);
th_hat = mod(angle(mean(exp(1j*th))),2*pi)/pi*180;
fprintf(' angle = %g, r = %g ; p = %g \n',th_hat,r,p)

%% Hits Misses on Same Polar
close all
X = nHitP;
Y = nMissP;
missColor = [150 220 220]/255;
hitColor  = [255 180 150]/255;

% Polar Plots
z1   = f(X);
z2   = f(Y);
th1  = angle(z1);
rho1 = abs(z1);
th2  = angle(z2);
rho2  = abs(z2);

% with magnitude of modulation
opts = [];
opts.colors = [hitColor;missColor];
opts.maxR   = 4;
opts.markerSize=300;
PolarPlot([th1 th2],[rho1 rho2],opts)
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/PhaseRetrieval/' 'HM_PolarRho' SubjSelectStr])

% v2
z3 = f(X)-f(Y);
th3 = angle(z3); rho3 = abs(z3);
opts = [];
opts.colors = [greyColor];
%opts.maxR   = 4;
opts.markerSize=300;
PolarPlot(th3,rho3,opts)
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/PhaseRetrieval/' 'H-M_PolarRho' SubjSelectStr])


opts = [];
opts.colors = [hitColor;missColor];
opts.maxR = 4/3;
opts.alpha = 0.8;
opts.markerSize=300;
opts.magText    = '';
PolarPlot([th1 th2],[],opts)
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/PhaseRetrieval/' 'HM_Polar' SubjSelectStr])

[p,r] = circ_rtest(th1-th2);
fprintf('Rayleight Test for H-M modulation across subjects:\n')
fprintf('r = %g ; p = %g \n',r,p)

makeLegend(opts.colors,[300 300],{'Hits','Misses'});
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/PhaseRetrieval/HitMissLegend']);

%% High Confidence vs Low Condidence Hits.
close all

temp=brewermap(2,'set1');
HCcolor = temp(1,:);
LCcolor = temp(2,:);

% Polar Plots
z1   = f(nHitP_LC);
z2   = f(nHitP_HC);
th1  = angle(z1);
rho1 = abs(z1);
th2  = angle(z2);
rho2  = abs(z2);

% with magnitude of modulation
opts = [];
opts.colors = [LCcolor;HCcolor;];
opts.maxR   = 4;
opts.markerSize=300;
opts.alpha  = 0.5;
PolarPlot([th1 th2],[rho1 rho2],opts);
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/PhaseConfRetrieval/' 'HCLCHits_PolarRho' SubjSelectStr])

opts = [];
opts.colors = [LCcolor;HCcolor;];
opts.maxR   = 4/3;
opts.markerSize=300;
opts.connect = 1;
opts.meanVecs = 0;
opts.alpha  = 0.5;
opts.polarGrid = 0;
opts.magText   ='';
PolarPlot([th1 th2],[],opts);
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/PhaseConfRetrieval/' 'HCLCHits_PolarRho_v2' SubjSelectStr])


% Projected to the circle
opts = [];
opts.colors = [LCcolor;HCcolor;];
opts.markerSize=300;
opts.maxR   = 4/3;
opts.magText ='';
opts.alpha  = 0.8;
PolarPlot([th1 th2],[],opts);
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/PhaseConfRetrieval/' 'HCLCHits_Polar' SubjSelectStr])

% Difference
z3 = z2-z1;
th3 = angle(z3);
rho3 = abs(z3);

% with magnitude
opts = [];
opts.colors = hitColor;
opts.maxR   = 4;
opts.markerSize=300;
opts.alpha  = 0.8;
PolarPlot(th3,rho3,opts);
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/PhaseConfRetrieval/' 'HC-LCHits_PolarRho' SubjSelectStr])

% Projected to the circle
opts.maxR   = 4/3;
opts.magText ='';
PolarPlot(th3,[],opts);
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/PhaseConfRetrieval/' 'HC-LCHits_Polar' SubjSelectStr])

% HC Hits Stats
[p,r] = circ_rtest(th1);
[a,b,c]=circ_mean(th1); a = mod(a*180/pi,360); b = mod(b*180/pi,360); c= mod(c*180/pi,360);
th_hat = mod(angle(mean(exp(1j*th1))),2*pi)/pi*180;
fprintf('Rayleight Test for LC Hits modulation across subjects:\n')
fprintf(' angle = %g, r = %g ; p = %g \n',th_hat,r,p)
fprintf('lower bound = %g, upper bound = %g ; mean = %g \n',c,b,a)

% LC Hits Stats
[p,r] = circ_rtest(th2);
[a,b,c]=circ_mean(th2); a = mod(a*180/pi,360); b = mod(b*180/pi,360); c= mod(c*180/pi,360);
th_hat = mod(angle(mean(exp(1j*th2))),2*pi)/pi*180;
fprintf('Rayleight Test for HC Hits modulation across subjects:\n')
fprintf(' angle = %g, r = %g ; p = %g \n',th_hat,r,p)
fprintf('lower bound = %g, upper bound = %g ; mean = %g \n',c,b,a)
% HC-LC Hits Stats
[p,r] = circ_rtest(th3);
[a,b,c]=circ_mean(th3); a = mod(a*180/pi,360); b = mod(b*180/pi,360); c= mod(c*180/pi,360);

th_hat = mod(angle(mean(exp(1j*th3))),2*pi)/pi*180;
fprintf('Rayleight Test for HC-LC Hits modulation across subjects:\n')
fprintf('angle = %g, r = %g ; p = %g \n',th_hat,r,p)
fprintf('lower bound = %g, upper bound = %g ; mean = %g \n',c,b,a)

makeLegend([LCcolor;HCcolor;hitColor],[300 300 300],{'LC Hits','HC Hits','HC-LC Hits'})
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/PhaseConfRetrieval/' 'HCLC_Legend'])

%% HC Hits to all Hits proportion
z = f(nHitP_HC./nHitP);
th = angle(z);
rho = abs(z);

% with magnitude
opts = [];
opts.colors = greyColor;
opts.markerSize=300;
opts.alpha  = 0.8;
PolarPlot(th,rho,opts)
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/PhaseConfRetrieval/' 'HCH_Prop_PolarRho' SubjSelectStr])

% Projected to the circle
opts.maxR   = 4/3;
opts.magText ='';
PolarPlot(th,[],opts)
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/PhaseConfRetrieval/' 'HCH_Prop_Polar' SubjSelectStr])

% HC Hits Prop Stats
[p,r] = circ_rtest(th);
th_hat = mod(angle(mean(exp(1j*th))),2*pi)/pi*180;
fprintf('Rayleight Test for Prp[ Hits modulation across subjects:\n')
fprintf(' angle = %g, r = %g ; p = %g \n',th_hat,r,p)

% Relationship to dprime
d=behav_out.retSummary.dPrime(subjs);
opts = [];
opts.alpha = 1;
opts.colors = greyColor;
opts.markerSize = 200;
opts.xlabel = ' \rho ';
opts.ylabel = ' d'' ' ;
opts.polyfitN = 1;

opts.text       =['R = ' num2str(round(corr(rho,d,'type','spearman')*100)/100)];
opts.xytext    = [0.05 0.4];
han = xyScatter(rho,d,opts);
fprintf('Correlation of magnitude of modulation with d'' :\n')
[c,p]=corr(rho,d,'type','spearman');
fprintf('r = %g ; p = %g \n',c,p)

print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/PhaseConfRetrieval/' 'HCH_PropMagDprimeScatter_' SubjSelectStr])

%% Misses HCvsLC 
close all; clc;
temp=brewermap(8,'dark2');
HCcolor = temp(3,:);
LCcolor = temp(4,:);
greyColor2 = [0.25 0.25 0.25];

% LC Misses
z1 = f(nMissP_LC);
th1 = angle(z1);
rho1 = abs(z1);

% HC Misses
z2 = f(nMissP_HC);
th2 = angle(z2);
rho2 = abs(z2);

% with magnitude of modulation
opts = [];
opts.colors = [LCcolor;HCcolor;];
opts.maxR   = 4;
opts.markerSize=300;
opts.alpha  = 0.5;
PolarPlot([th1 th2],[rho1 rho2],opts);
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/PhaseConfRetrieval/' 'HCLCMiss_PolarRho' SubjSelectStr])

% along the unit circle
opts = [];
opts.colors = [LCcolor;HCcolor;];
opts.maxR   = 4/3;
opts.markerSize=300;
opts.alpha  = 0.8;
opts.magText = '';
PolarPlot([th1 th2],[],opts);
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/PhaseConfRetrieval/' 'HCLCMiss_Polar' SubjSelectStr])

% Difference
z3 = z2-z1;
th3 = angle(z3);
rho3 = abs(z3);
opts        =[];
opts.maxR   = 4;
opts.markerSize=300;
opts.alpha  = 0.8;
opts.colors = missColor;
PolarPlot(th3,rho3,opts);
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/PhaseConfRetrieval/' 'HCLCMissRho_Polar' SubjSelectStr])

% Projected to the circle
opts.maxR   = 4/3;
opts.magText ='';
opts.alpha = 0.8;
PolarPlot(th3,[],opts);
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/PhaseConfRetrieval/' 'HC-LCMiss_Polar' SubjSelectStr])

% Proportion of HC Misses
z4 = f(nMissP_HC./nMissP);
th4 = angle(z4);
rho4 = abs(z4);

% Projected to the circle
opts.maxR   = 4/3;
opts.magText ='';
opts.alpha = 0.8;
opts.colors = greyColor;
PolarPlot(th4,[],opts);
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/PhaseConfRetrieval/' 'propHCMisses_Polar' SubjSelectStr])

% HC Misses Stats
[p,r] = circ_rtest(th1);
th_hat = mod(angle(mean(exp(1j*th1))),2*pi)/pi*180;
fprintf('Rayleight Test for LC Misses modulation across subjects:\n')
fprintf(' angle = %g, r = %g ; p = %g \n',th_hat,r,p)

% LC Misses Stats
[p,r] = circ_rtest(th2);
th_hat = mod(angle(mean(exp(1j*th2))),2*pi)/pi*180;
fprintf('Rayleight Test for HC Misses modulation across subjects:\n')
fprintf(' angle = %g, r = %g ; p = %g \n',th_hat,r,p)

% HC-LC Misses Stats
[p,r] = circ_rtest(th3);
th_hat = mod(angle(mean(exp(1j*th3))),2*pi)/pi*180;
fprintf('Rayleight Test for HC-LC Misses modulation across subjects:\n')
fprintf(' angle = %g, r = %g ; p = %g \n',th_hat,r,p)

% HC/all Misses Stats
[p,r] = circ_rtest(th4);
th_hat = mod(angle(mean(exp(1j*th4))),2*pi)/pi*180;
fprintf('Rayleight Test for HC/all hits Misses modulation across subjects:\n')
fprintf(' angle = %g, r = %g ; p = %g \n',th_hat,r,p)

makeLegend([LCcolor; HCcolor ;missColor],[300 300 300],{'LC Misses','HC Misses','HC-LC Misses'});
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/PhaseConfRetrieval/' 'HCLCMiss_Legend'])

%% All Confidence
close all; clc;
% confidence independent of memory:
z1 = f(nConfP_LC); % low
z2 = f(nConfP_MC); % mid
z3 = f(nConfP_HC); % mid
th1 = angle(z1); th2 = angle(z2); th3=angle( z3);

confColors = brewermap(3,'set3');

% all confidence
opts = [];
opts.colors = [confColors];
opts.maxR   = 4/3;
opts.markerSize=300;
opts.alpha  = 0.7;
opts.magText = '';
PolarPlot([th1 th2 th3],[],opts)
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/PhaseConfRetrieval/' 'allConf_Polar' SubjSelectStr])

% proportion of high confidence 
th = angle(f(nConfP_HC./squeeze(sum(nConfP_C,2))));
opts = [];
opts.colors = greyColor;
opts.maxR   = 4/3;
opts.markerSize=300;
opts.alpha  = 0.7;
opts.magText = '';
PolarPlot(th,[],opts)
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/PhaseConfRetrieval/' 'propHC_Polar' SubjSelectStr])

% Low Conf Stats
[p,r] = circ_rtest(th1);
th_hat = mod(angle(mean(exp(1j*th1))),2*pi)/pi*180;
fprintf('Rayleight Test for LC modulation across subjects:\n')
fprintf(' angle = %g, r = %g ; p = %g \n',th_hat,r,p)

% Mid Conf Stats
[p,r] = circ_rtest(th2);
th_hat = mod(angle(mean(exp(1j*th2))),2*pi)/pi*180;
fprintf('Rayleight Test for MC modulation across subjects:\n')
fprintf(' angle = %g, r = %g ; p = %g \n',th_hat,r,p)

% High Conf Stats
[p,r] = circ_rtest(th3);
th_hat = mod(angle(mean(exp(1j*th3))),2*pi)/pi*180;
fprintf('Rayleight Test for HC modulation across subjects:\n')
fprintf(' angle = %g, r = %g ; p = %g \n',th_hat,r,p)

% Prop of High Confidence
% High Conf Stats
[p,r] = circ_rtest(th);
th_hat = mod(angle(mean(exp(1j*th))),2*pi)/pi*180;
fprintf('Rayleight Test for proportion HC modulation across subjects:\n')
fprintf(' angle = %g, r = %g ; p = %g \n',th_hat,r,p)

makeLegend([confColors],[300 300 300],{'Low','Mid','High'})
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/PhaseConfRetrieval/' 'Conf_Legend'])

%% Splits for HC Hits for Faces and Scenes

% Get data splits
nFaceHitsConf = zeros(nSubjs,5,3);
nFaceHits     = zeros(nSubjs,5);
nSceneHitsConf = zeros(nSubjs,5,3);
nSceneHits     = zeros(nSubjs,5);
for ss = 1:nSubjs
    s       = subjs(ss);
    Faces   = squeeze(out.datMat(s,:,1)==1);    
    Hits    = squeeze(out.datMat(s,:,7)==1);
    Conf    = squeeze(out.datMat(s,:,9));
    Phase   = squeeze(out.datMat(s,:,5));   
    for pp = 1:5
        trials              = (Faces & Hits) & (Phase==pp);
        nFaceHits(ss,pp)    = sum(trials);
        trials              = (~Faces & Hits) & (Phase==pp);
        nSceneHits(ss,pp)   = sum(trials);
        for cc = 1:3        
            trials = (Faces & Hits) & (Conf==cc) & (Phase==pp);
            nFaceHitsConf(ss,pp,cc) = sum(trials);
            trials = (~Faces & Hits) & (Conf==cc) & (Phase==pp);
            nSceneHitsConf(ss,pp,cc) = sum(trials);
        end        
    end
end
%% Plots for HC Scene Hits 
close all; clc;
% Polar w Mag Plot of HC vs LC Hit Scenes
z1  = f(nSceneHitsConf(:,:,1)); % LC
z2 = f(nSceneHitsConf(:,:,3));  % HC
th1 = angle(z1);    th2 = angle(z2);
rho1= abs(z1);      rho2 = abs(z2);

temp=brewermap(2,'set1');
HCcolor = temp(1,:);
LCcolor = temp(2,:);

opts = [];
opts.colors = [LCcolor;HCcolor;];
opts.maxR   = 4;
opts.markerSize=300;
opts.alpha  = 0.5;
PolarPlot([th1 th2],[rho1 rho2],opts)
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/PhaseConfRetrieval/' 'HCLCSceneHits_PolarRho' SubjSelectStr])

% Projected to the circle
opts = [];
opts.colors = [LCcolor;HCcolor;];
opts.markerSize=300;
opts.maxR   = 4/3;
opts.magText ='';
opts.alpha  = 0.8;
PolarPlot([th1 th2],[],opts)
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/PhaseConfRetrieval/' 'HCLCSceneHits_Polar' SubjSelectStr])

% Polar Plot of HC - LC Hit Scenes
z3   = z2-z1;
th3  = angle(z3);
rho3 = abs(z3);
opts        =[];
opts.maxR   = 4;
opts.markerSize=300;
opts.alpha  = 0.8;
opts.colors = sceneColor;
PolarPlot(th3,rho3,opts);
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/PhaseConfRetrieval/' 'HC-LCSceneHits_PolarRho' SubjSelectStr])

opts        =[];
opts.maxR   = 4/3;
opts.markerSize=300;
opts.magText ='';
opts.colors = sceneColor;
PolarPlot(th3,[],opts)
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/PhaseConfRetrieval/' 'HC-LCSceneHits_Polar' SubjSelectStr])

% HC Hit Scene Stats
[p,r] = circ_rtest(th1);
th_hat = mod(angle(mean(exp(1j*th1))),2*pi)/pi*180;
fprintf('Rayleight Test for LC Hit Scene modulation across subjects:\n')
fprintf(' angle = %g, r = %g ; p = %g \n',th_hat,r,p)
[a,b,c]=circ_mean(th1); a = mod(a*180/pi,360); b = mod(b*180/pi,360); c= mod(c*180/pi,360);
fprintf('lower bound = %g, upper bound = %g ; mean = %g \n',floor(c),ceil(b),round(a))

% LC Hit Scene Stats
[p,r] = circ_rtest(th2);
th_hat = mod(angle(mean(exp(1j*th2))),2*pi)/pi*180;
fprintf('Rayleight Test for HC Hit Scene modulation across subjects:\n')
fprintf(' angle = %g, r = %g ; p = %g \n',th_hat,r,p)
[a,b,c]=circ_mean(th2); a = mod(a*180/pi,360); b = mod(b*180/pi,360); c= mod(c*180/pi,360);
fprintf('lower bound = %g, upper bound = %g ; mean = %g \n',floor(c),ceil(b),round(a))

% HC-LC Hit Scene Stats
[p,r] = circ_rtest(th3);
th_hat = mod(angle(mean(exp(1j*th3))),2*pi)/pi*180;
fprintf('Rayleight Test for HC-LC Hit Scene modulation across subjects:\n')
fprintf(' angle = %g, r = %g ; p = %g \n',th_hat,r,p)
[a,b,c]=circ_mean(th3); a = mod(a*180/pi,360); b = mod(b*180/pi,360); c= mod(c*180/pi,360);
fprintf('lower bound = %g, upper bound = %g ; mean = %g \n',floor(c),ceil(b),round(a))

makeLegend([LCcolor; HCcolor; sceneColor],[300 300 300],{'LC Scn Hits','HC Scn Hits','HC-LC Scn Hits'})
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/PhaseConfRetrieval/' 'HCLCScnHits_Legend'])

%% Plots for HC Face Hits 
close all; clc;
% Polar w Mag Plot of HC vs LC Hit Scenes
z1  = f(nFaceHitsConf(:,:,1)); % LC
z2 = f(nFaceHitsConf(:,:,3));  % HC
th1 = angle(z1);    th2 = angle(z2);
rho1= abs(z1);      rho2 = abs(z2);

temp=brewermap(2,'set1');
HCcolor = temp(1,:);
LCcolor = temp(2,:);
FaceColor = [100 200 100]/255;
colors = [100 200 100;200 100 200]/255;
txtStr = {'Faces','Scenes'};

opts = [];
opts.colors = [LCcolor;HCcolor;];
opts.maxR   = 4;
opts.markerSize=300;
opts.alpha  = 0.5;
PolarPlot([th1 th2],[rho1 rho2],opts);
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/PhaseConfRetrieval/' 'HCLCFaceHits_PolarRho' SubjSelectStr])

% Projected to the circle
opts = [];
opts.colors = [LCcolor;HCcolor;];
opts.markerSize=300;
opts.maxR   = 4/3;
opts.magText ='';
opts.alpha  = 0.8;
PolarPlot([th1 th2],[],opts);
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/PhaseConfRetrieval/' 'HCLCFaceHits_Polar' SubjSelectStr])

% Polar Plot of HC - LC Hit Faces
z3   = z2-z1;
th3  = angle(z3);
rho3 = abs(z3);
opts        =[];
opts.maxR   = 4;
opts.markerSize=300;
opts.alpha  = 0.8;
opts.colors = faceColor;
PolarPlot(th3,rho3,opts);
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/PhaseConfRetrieval/' 'HC-LCFaceHits_PolarRho' SubjSelectStr])

opts        =[];
opts.maxR   = 4/3;
opts.markerSize=300;
opts.magText ='';
opts.colors = [faceColor];
PolarPlot(th3,[],opts);
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/PhaseConfRetrieval/' 'HC-LCFaceHits_Polar' SubjSelectStr])

% HC Hit Face Stats
[p,r] = circ_rtest(th1);
th_hat = mod(angle(mean(exp(1j*th1))),2*pi)/pi*180;
fprintf('Rayleight Test for LC Hit Face modulation across subjects:\n')
fprintf(' angle = %g, r = %g ; p = %g \n',th_hat,r,p)
[a,b,c]=circ_mean(th1); a = mod(a*180/pi,360); b = mod(b*180/pi,360); c= mod(c*180/pi,360);
fprintf('lower bound = %g, upper bound = %g ; mean = %g \n',floor(c),ceil(b),round(a))

% LC Hit Face Stats
[p,r] = circ_rtest(th2);
th_hat = mod(angle(mean(exp(1j*th2))),2*pi)/pi*180;
fprintf('Rayleight Test for HC Hit Face modulation across subjects:\n')
fprintf(' angle = %g, r = %g ; p = %g \n',th_hat,r,p)
[a,b,c]=circ_mean(th2); a = mod(a*180/pi,360); b = mod(b*180/pi,360); c= mod(c*180/pi,360);
fprintf('lower bound = %g, upper bound = %g ; mean = %g \n',floor(c),ceil(b),round(a))

% HC-LC Hit Face Stats
[p,r] = circ_rtest(th3);
th_hat = mod(angle(mean(exp(1j*th3))),2*pi)/pi*180;
fprintf('Rayleight Test for HC-LC Hit Face modulation across subjects:\n')
fprintf(' angle = %g, r = %g ; p = %g \n',th_hat,r,p)
[a,b,c]=circ_mean(th3); a = mod(a*180/pi,360); b = mod(b*180/pi,360); c= mod(c*180/pi,360);
fprintf('lower bound = %g, upper bound = %g ; mean = %g \n',floor(c),ceil(b),round(a))

makeLegend([LCcolor; HCcolor; faceColor],[300 300 300],{'LC Face Hits','HC Face Hits','HC-LC Face Hits'})
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/PhaseConfRetrieval/' 'HCLCFaceHits_Legend'])

%% %% Bootstrap difference between HC and LC Hits;
close all; clc; rng(1)

nBoot = 1000;
lowerB = round(nBoot*0.024);
upperB = round(nBoot*0.976);

% both categories
z1  = f(nHitP_LC); % LC
z2 = f(nHitP_HC);  % HC
th1 = mod(angle(z1),2*pi);    th2 = mod(angle(z2),2*pi);
B = bootstrp(nBoot,deltaTh,th1,th2);
medianB = circ_median(B);
[sB,i]=sort(B-medianB);
bounds = [sB(lowerB) sB(upperB)]+medianB;

opts        =[];
opts.maxR   = 4/3;
opts.markerSize=50;
opts.alpha = 0.3;
opts.magText ='';
opts.meanVecs = 1;
opts.colors = hitColor;
PolarPlot(B,1+0.05*randn(nBoot,1),opts)   
for jj = 1:2
    x = 1.2*cos(bounds(jj));
    y = 1.2*sin(bounds(jj));
    plot([0 x],[0 y],'linewidth',5,'color',greyColor)    
end

print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/PhaseConfRetrieval/' 'HC-LCHitsBoot_Polar' SubjSelectStr])
%

% Faces
z1  = f(nFaceHitsConf(:,:,1)); % LC
z2 = f(nFaceHitsConf(:,:,3));  % HC
th1 = mod(angle(z1),2*pi);    th2 = mod(angle(z2),2*pi);
B = bootstrp(nBoot,deltaTh,th1,th2);
medianB = circ_median(B);
[sB,i]=sort(B-medianB);
bounds = [sB(lowerB) sB(upperB)]+medianB;

opts        =[];
opts.maxR   = 4/3;
opts.markerSize=50;
opts.alpha = 0.3;
opts.magText ='';
opts.meanVecs = 1;
opts.colors = faceColor;
PolarPlot(B,1+0.05*randn(nBoot,1),opts)   
for jj = 1:2
    x = 1.2*cos(bounds(jj));
    y = 1.2*sin(bounds(jj));
    plot([0 x],[0 y],'linewidth',5,'color',greyColor)    
end
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/PhaseConfRetrieval/' 'HC-LCFaceHitsBoot_Polar' SubjSelectStr])

% Scenes
% Polar w Mag Plot of HC vs LC Hit Scenes
z1  = f(nSceneHitsConf(:,:,1)); % LC
z2 = f(nSceneHitsConf(:,:,3));  % HC
th1 = mod(angle(z1),2*pi);    th2 = mod(angle(z2),2*pi);
B = bootstrp(nBoot,deltaTh,th1,th2);
medianB = circ_median(B);
[sB,i]=sort(B-medianB);
bounds = [sB(lowerB) sB(upperB)]+medianB;

opts        =[];
opts.maxR   = 4/3;
opts.markerSize=50;
opts.alpha = 0.3;
opts.magText ='';
opts.meanVecs = 1;
opts.colors = sceneColor;
PolarPlot(B,1+0.05*randn(nBoot,1),opts)   
for jj = 1:2
    x = 1.2*cos(bounds(jj));
    y = 1.2*sin(bounds(jj));
    plot([0 x],[0 y],'linewidth',5,'color',greyColor)    
end
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/PhaseConfRetrieval/' 'HC-LCSceneHitsBoot_Polar' SubjSelectStr])

makeLegend([hitColor; sceneColor; faceColor],[300 300 300],{'All Hits','Scene Hits','Face Hits'})
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/PhaseConfRetrieval/' 'allHits_Legend'])

%% Retrieval RT Modulation
% Get data splits
meRTsHits           = zeros(nSubjs,5);
meRTsHitsConf       = zeros(nSubjs,5,3);
meRTsFaceHitsConf   = zeros(nSubjs,5,3);
meRTsFaceHits       = zeros(nSubjs,5);
meRTsSceneHitsConf   = zeros(nSubjs,5,3);
meRTsSceneHits     = zeros(nSubjs,5);
for ss = 1:nSubjs
    s       = subjs(ss);
    RTs     = -log(squeeze(out.datMat(subjs,:,11)));
    Faces   = squeeze(out.datMat(s,:,1)==1);    
    Hits    = squeeze(out.datMat(s,:,7)==1);
    Conf    = squeeze(out.datMat(s,:,9));
    Phase   = squeeze(out.datMat(s,:,5));   
    for pp = 1:5
        % all hits
        trials              = (Hits) & (Phase==pp);
        meRTsHits(ss,pp)    = nanmedian(RTs(trials));
        % face hits
        trials              = (Faces & Hits) & (Phase==pp);        
        meRTsFaceHits(ss,pp)    = nanmedian(RTs(trials));
        % scene hits
        trials              = (~Faces) & Hits & (Phase==pp);        
        meRTsSceneHits(ss,pp)    = nanmedian(RTs(trials));
        for cc = 1:3
            % all hits
            trials = Hits & (Conf==cc) & (Phase==pp);
            meRTsHitsConf(ss,pp,cc)    = nanmedian(RTs(trials));
            % face hits
            trials = (Faces & Hits) & (Conf==cc) & (Phase==pp);
            meRTsFaceHitsConf(ss,pp,cc) = nanmedian(RTs(trials));
            % scene hits
            trials = (~Faces & Hits) & (Conf==cc) & (Phase==pp);
            meRTsSceneHitsConf(ss,pp,cc) = nanmedian(RTs(trials));
        end        
    end
end

% rows-> hits, face hits, scene hits
% columns -> all hits, HC and LC projected, HC-LC

% create data structure for easy iteration by row
X = cell(3,2);
X{1,1} = f(meRTsHits);      X{1,2} = [f(meRTsHitsConf(:,:,1))     f(meRTsHitsConf(:,:,3))];
X{2,1} = f(meRTsFaceHits);  X{2,2} = [f(meRTsFaceHitsConf(:,:,1)) f(meRTsFaceHitsConf(:,:,3))];
X{3,1} = f(meRTsSceneHits); X{3,2} = [f(meRTsSceneHitsConf(:,:,1)) f(meRTsSceneHitsConf(:,:,3))];
Y = cell(3,3);
for ii = 1:3    
    Y{ii,1} = angle(X{ii,1});
    Y{ii,2} = [angle(X{ii,2}(:,1)) angle(X{ii,2}(:,2))];
    Y{ii,3} = angle( X{ii,2}(:,2)-X{ii,2}(:,1) );     
end

% compute statistics for hits, and hc-lc hits
P = zeros(3,2); R = zeros(3,2);
strs1 = {'hits', 'face hits', 'scene hits'};
strs2 = {'', 'HC-LC'};
for ii = 1:3
    cnt = 1;
    for jj = [1 3]
        [P(ii,cnt),R(ii,cnt)]=circ_rtest(Y{ii,jj});                
        fprintf(['Rayleight Test for %s %s modulation across subjects:\n'],strs2{cnt},strs1{ii})
        fprintf('r = %g ; p = %g \n',R(ii,cnt),P(ii,cnt))
        cnt = cnt +1;        
    end
end

% create 3x3 figure: 
close all
han=figure(1); clf;
Dims = [550 550];
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',Dims,'paperposition',[0 0 Dims])
set(gcf,'position',[300,100,Dims])

dxy = 50; margins = 0; axDims = 150;
iiPixCnt = margins; jjPixCnt = margins;
axHan    = cell(3,3);
for ii = 3:-1:1 
    for jj = 1:3        
        axHan{ii,jj}=axes('units','points','position',[iiPixCnt,jjPixCnt,axDims,axDims]);
        iiPixCnt = iiPixCnt+dxy+axDims;
    end
    iiPixCnt = margins;
    jjPixCnt = jjPixCnt+dxy+axDims;
end

colors = cell(3,3);
colors{1,1} = hitColor;     colors{1,3} = greyColor;
colors{2,1} = faceColor;    colors{2,3} = greyColor;
colors{3,1} = sceneColor;   colors{3,3} = greyColor;
colors{1,2} = [LCcolor; HCcolor]; colors{2,2} = [LCcolor; HCcolor]; colors{3,2} = [LCcolor; HCcolor];

opts            =[];
opts.maxR       = 4/3;
opts.markerSize = 150;
opts.alpha      = 0.7;
opts.magText    ='';
for ii = 1:3
    for jj = 1:3
        opts.handle = axHan{ii,jj};
        opts.colors = colors{ii,jj};
        PolarPlot(Y{ii,jj},[],opts);
    end
end
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/PhaseConfRetrieval/' 'RTsHits_Polar' SubjSelectStr])
makeLegend([greyColor],[300 300 300],{'HC-LC Hits'})
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/PhaseConfRetrieval/' 'HC-LCHitsGrey_Legend'])

