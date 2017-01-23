close all
clearvars
addpath Plotting/
addpath Analyses/
addpath Library/
close all
clearvars

dataPath        = '~/Google Drive/Research/tACS/tACS_ER_task/data/tacs_er_objstim/';
plotPath1       = '~/Google Drive/Research/tACS/tACS_ER_task/plots/tacs_er_objstim/PhaseRetrieval/';
plotPath2       = '~/Google Drive/Research/tACS/tACS_ER_task/plots/tacs_er_objstim/PhaseConfRetrieval/';
load([dataPath 'Summary/BehavSummary.mat'])
load([dataPath 'Summary/DataMatrix.mat'])

% Exclusion of subjects based on poor behavioral performance
badSubjs =  [3 6 8 27 29 32 34];

nTotalSubjs      = 37;
subjs           = setdiff(1:nTotalSubjs,badSubjs);
nSubjs          = numel(subjs);

%% data colors
temp        = brewermap(6,'accent');
SmallColor  = temp(5,:);
BigColor    = temp(6,:);
temp        = brewermap(2,'set2');
HCcolor     = temp(1,:);
LCcolor     = temp(2,:);

greyColor   = [119,136,153]/255;
hitColor    = [255 180 150]/255;
missColor   = [150 220 220]/255;

% Trial Types
P = mod(out.datMat(subjs,:,4),2*pi); % Phase
H = out.datMat(subjs,:,5); % Hits
M = out.datMat(subjs,:,6); % Misses
C = out.datMat(subjs,:,7); % Confidence
RR = -log(out.datMat(subjs,:,8)); % Retrieval RTs
ER = out.datMat(subjs,:,3); % Encoding RTs
ST = out.datMat(subjs,:,1); % stimulus names
RRr = (RR*6)*2*pi; % retrieval rts in radians
ERr = (ER*6)*2*pi; % retrieval rts in radians

% Phase Selection
% Phase Bins
PhaseDelta  = 2*pi/5;
[~,PhaseEdges]=histcounts(P(1,:),'BinWidth',PhaseDelta);
PhaseCenters = PhaseEdges+PhaseDelta/2;
PhaseCenters(end) =[];
nPhaseBins = numel(PhaseCenters);

NsPhase  = zeros(nSubjs, nPhaseBins);
NHsPhase = zeros(nSubjs, nPhaseBins);
NMsPhase = zeros(nSubjs, nPhaseBins);
for ss = 1:nSubjs
    NsPhase(ss,:)       = histcounts(P(ss,:),PhaseEdges);
    NHsPhase(ss,:)      = histcounts(P(ss,H(ss,:)==1),PhaseEdges);
    NHsHCPhase(ss,:)    = histcounts(P(ss,H(ss,:)==1 & C(ss,:)==3),PhaseEdges);
    NHsMCPhase(ss,:)    = histcounts(P(ss,H(ss,:)==1 & C(ss,:)==2),PhaseEdges);
    NHsLCPhase(ss,:)    = histcounts(P(ss,H(ss,:)==1 & C(ss,:)==1),PhaseEdges);
    NMsPhase(ss,:)      = histcounts(P(ss,M(ss,:)==1),PhaseEdges);
end

% Functions
PhaseMat = repmat(exp(1j*PhaseCenters),nSubjs,1);
fQuant = @(X)(mean(X.*PhaseMat./NsPhase,2));

f1 = @(th)(nansum(exp(1j*th)));
f1b = @(th)(nanmean(exp(1j*th)));
f2 = @(th,mag)(nanmean(mag.*exp(1j*th)));
m_th     = @(th)(angle(mean(exp(1j*th))));
m_rho    = @(th)(abs(mean(exp(1j*th))));
deltaTh = @(th1,th2)(mod(abs(m_th(th1)-m_th(th2)),2*pi));

%% all Hits by Phase (phase bins)
close all
Z_qH = fQuant(NHsPhase);
Z_qM = fQuant(NMsPhase);

% Polar Plots for nHits
th1  = angle(Z_qH); th2 = angle(Z_qM);
rho1 = abs(Z_qH);   rho2= abs(Z_qM);

% with magnitude of modulation
opts = [];
opts.colors = hitColor ;
opts.markerSize=300;
PolarPlot(th1,rho1,opts)
print(gcf,'-dpdf',[plotPath1 'HitProbPolarRhoQuant'])

opts = [];
opts.colors = hitColor;
opts.maxR   = 4/3;
opts.alpha  = 0.8;
opts.markerSize=300;
opts.magText    = '';
PolarPlot(th1,[],opts)
print(gcf,'-dpdf',[plotPath1 'HitProbPolarQuant'])

[p,r] = circ_rtest(th1);
th_hat = circ_mean(th1);
fprintf('Rayleight Test for modulation on Hits across subjects:\n')
fprintf(' angle = %g, r = %g ; p = %g \n',th_hat,r,p)

%% Misses by Phase (phase bins)
close all

% with magnitude of modulation
opts = [];
opts.colors = missColor;
opts.markerSize=300;
PolarPlot(th2,rho2,opts)
print(gcf,'-dpdf',[plotPath1 'MissProbPolarRho'])

opts = [];
opts.colors = missColor;
opts.maxR = 4/3;
opts.alpha = 0.8;
opts.markerSize=300;
opts.magText    = '';
PolarPlot(th2,[],opts)
print(gcf,'-dpdf',[plotPath1 'MissProbPolar'])

[p,r] = circ_rtest(th2);
fprintf('Rayleight Test for Misses modulation across subjects:\n')
th_hat = circ_mean(th2);
fprintf(' angle = %g, r = %g ; p = %g \n',th_hat,r,p)

%% hits by phase
clc; close all;
Z1 = zeros(nSubjs,1);
Z2 = zeros(nSubjs,1);
for ss = 1:nSubjs
    Z1(ss) = f1(P(ss,H(ss,:)==1));
    
    Z2(ss) = f1(P(ss,M(ss,:)==1));
    
    Z3(ss) = f1([P(ss,H(ss,:)==1) -P(ss,M(ss,:)==1)]) ;
    [P3r(ss),Z3r(ss)] = circ_rtest([P(ss,H(ss,:)==1) -P(ss,M(ss,:)==1)]);
end
% Polar Plots for nHits
th1  = angle(Z1); th2 = angle(Z2);
rho1 = abs(Z1);   rho2= abs(Z2);

% with magnitude of modulation
opts = [];
opts.colors = hitColor ;
opts.markerSize=300;
PolarPlot(th1,rho1,opts)
print(gcf,'-dpdf',[plotPath1 'HitProbPolarRho'])

opts = [];
opts.colors = hitColor;
opts.maxR   = 4/3;
opts.alpha  = 0.8;
opts.markerSize=300;
opts.magText    = '';
PolarPlot(th1,[],opts)
print(gcf,'-dpdf',[plotPath1 'HitProbPolar'])

[p,r] = circ_rtest(th1);
th_hat = circ_mean(th1);
fprintf('Rayleight Test for modulation on Hits across subjects:\n')
fprintf(' angle = %g, r = %g ; p = %g \n',th_hat,r,p)

%% misses by phase
close all

% with magnitude of modulation
opts = [];
opts.colors = missColor;
opts.markerSize=300;
PolarPlot(th2,rho2,opts)
print(gcf,'-dpdf',[plotPath1 'MissProbPolarRho'])

opts = [];
opts.colors = missColor;
opts.maxR = 4/3;
opts.alpha = 0.8;
opts.markerSize=300;
opts.magText    = '';
PolarPlot(th2,[],opts)
print(gcf,'-dpdf',[plotPath1 'MissProbPolar'])

[p,r] = circ_rtest(th2);
fprintf('Rayleight Test for Misses modulation across subjects:\n')
th_hat = circ_mean(th2);
fprintf(' angle = %g, r = %g ; p = %g \n',th_hat,r,p)

%% Hits Misses on Same Polar
close all
clc; close all;
Z1 = zeros(nSubjs,1);
Z2 = zeros(nSubjs,1);
for ss = 1:nSubjs
    Z1(ss) = f1(P(ss,H(ss,:)==1));
    Z2(ss) = f1(P(ss,M(ss,:)==1));
    Z1b(ss) = f1b(P(ss,H(ss,:)==1));
    Z2b(ss) = f1b(P(ss,M(ss,:)==1));
     [~,R1(ss)] = circ_rtest(P(ss,H(ss,:)==1));
     [~,R2(ss)] = circ_rtest(P(ss,M(ss,:)==1));
end
th1     = angle(Z1);
rho1    = abs(Z1);
th2     = angle(Z2);
rho2    = abs(Z2);

% with magnitude of modulation
opts = [];
opts.colors = [hitColor;missColor];
%opts.maxR   = 4;
opts.markerSize=300;
PolarPlot([th1 th2],[rho1 rho2],opts);
print(gcf,'-dpdf',[plotPath1 'HM_PolarRho'])

% v2
z3 = Z1-Z2;
th3 = angle(z3); rho3 = abs(z3);
opts = [];
opts.colors = [greyColor];
%opts.maxR   = 4;
opts.markerSize=300;
PolarPlot(th3,rho3,opts);
print(gcf,'-dpdf',[plotPath1 'H-M_PolarRho'])

opts = [];
opts.colors = [hitColor;missColor];
opts.maxR = 4/3;
opts.alpha = 0.8;
opts.markerSize=300;
opts.magText    = '';
PolarPlot([th1 th2],[],opts);
print(gcf,'-dpdf',[plotPath1 'HM_Polar'])

opts = [];
opts.colors = [greyColor];
opts.maxR = 4/3;
opts.alpha = 0.8;
opts.markerSize=300;
opts.magText    = '';
PolarPlot([th3],[],opts);
print(gcf,'-dpdf',[plotPath1 'H-M_Polar'])

[p,r] = circ_rtest(th3);
fprintf('Rayleight Test for H-M modulation across subjects:\n')
fprintf('r = %g ; p = %g \n',r,p)

fprintf('RankSum Test for greater modulation across subjects:\n')
[p,~,k]=signrank(abs(Z1),abs(Z2));
fprintf('Hits Median = %g; Miss Median = %g; Z = %g ; p = %g \n',median(rho1),median(rho2),k.zval,p)

opts.colors = [hitColor;missColor;greyColor];
makeLegend(opts.colors,[300 300 300],{'Hits','Misses','H-M'});
print(gcf,'-dpdf',[plotPath1 'HM_Legend'])

%% Hits-Misses by subject score
Z3 = zeros(nSubjs,1);
for ss = 1:nSubjs
    Z3(ss) = f1([P(ss,H(ss,:)==1) -P(ss,M(ss,:)==1)]) ;
    [P3r(ss),Z3r(ss)] = circ_rtest([P(ss,H(ss,:)==1) -P(ss,M(ss,:)==1)]);        
end
th3 = angle(Z3); rho3 = abs(Z3);
opts = [];
opts.colors = [greyColor];
%opts.maxR   = 4;
opts.markerSize=300;
opts.alpha = 0.6;
PolarPlot(th3,rho3,opts);
print(gcf,'-dpdf',[plotPath1 'H-M_PolarRho2'])

opts = [];
opts.colors = [greyColor];
opts.maxR = 4/3;
opts.alpha = 0.8;
opts.markerSize=300;
opts.magText    = '';
PolarPlot([th3],[],opts);
print(gcf,'-dpdf',[plotPath1 'HM_Polar2'])

figure(); clf;
AR = [450 300];
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',AR,'paperposition',[0 0 AR])
set(gcf,'position',[200,200,AR])
a = axes('units','points','position',[50 50 100 200]); hold on;

%% High Confidence vs Low Condidence Hits.
close all

temp=brewermap(3,'set1');
HCcolor = temp(1,:);
LCcolor = temp(2,:);
MCcolor = temp(3,:);

% Polar Plots
clc; close all;
Z_H_HC = zeros(nSubjs,1);
Z_H_MC = zeros(nSubjs,1);
Z_H_LC = zeros(nSubjs,1);
Z_H_MLC= zeros(nSubjs,1);
for ss = 1:nSubjs
    Z_H_HC(ss)  = f1(P(ss,H(ss,:)==1 & C(ss,:)==3));
    Z_H_MC(ss)  = f1(P(ss,H(ss,:)==1 & C(ss,:)==2));
    Z_H_LC(ss)  = f1(P(ss,H(ss,:)==1 & C(ss,:)==1));
    Z_H_MLC(ss) = f1(P(ss,H(ss,:)==1 & C(ss,:)>=2));
    Z_LC(ss)    = f1(P(ss,C(ss,:)<=2));
    
    [~,R1(ss)] = circ_rtest(P(ss,H(ss,:)==1 & C(ss,:)==3));
    [~,R2(ss)] = circ_rtest(P(ss,H(ss,:)==1 & C(ss,:)==1));
    
    pp = [P(ss,H(ss,:)==1 & C(ss,:)==3) -P(ss,H(ss,:)==1 & C(ss,:)==1)];
    Z3_H_HC_LC(ss) = f1(pp);
    pp2 = [P(ss,H(ss,:)==1 & C(ss,:)==3) -P(ss,M(ss,:)==1)];
    Z3_H_HC_M(ss) = f1(pp2);
    
    Z1b(ss) = f1b(P(ss,H(ss,:)==1 & C(ss,:)==3));
    Z2b(ss) = f1b(P(ss,H(ss,:)==1 & C(ss,:)==1));
end
th_H  = angle(Z_H_HC);
rho_H = abs(Z_H_HC);
th_M  = angle(Z_H_MC);
rho_M = abs(Z_H_MC);
th_L  = angle(Z_H_LC);
rho_L  = abs(Z_H_LC);

% with magnitude of modulation
opts = [];
opts.colors = [HCcolor; LCcolor;];
%opts.maxR   = 4;
opts.markerSize=300;
opts.alpha  = 0.5;
PolarPlot([th_H th_L],[rho_H rho_L],opts);
%PolarPlot([th_H th_M th_L],[rho_H rho_M rho_L],opts);
print(gcf,'-dpdf',[plotPath1 'HCLCHits_PolarRho'])

opts = [];
opts.colors = [HCcolor;LCcolor;];
opts.maxR   = 4/3;
opts.markerSize=300;
opts.connect = 1;
opts.meanVecs = 0;
opts.alpha  = 0.5;
opts.polarGrid = 0;
opts.magText   ='';
PolarPlot([th_H th_L],[],opts);
print(gcf,'-dpdf',[plotPath1 'HCLCHits_PolarRho_v2'])

% Projected to the circle
opts = [];
opts.colors = [HCcolor;LCcolor;];
opts.markerSize=300;
opts.maxR   = 4/3;
opts.magText ='';
opts.alpha  = 0.8;
PolarPlot([th_H th_L ],[],opts);
print(gcf,'-dpdf',[plotPath1 'HCLCHits_Polar'])

% all hit types
opts = [];
opts.colors = [HCcolor; MCcolor; LCcolor;];
opts.markerSize = 300;
opts.maxR       = 4/3;
opts.magText    = '';
opts.alpha      = 0.8;
PolarPlot([th_H th_M th_L ],[],opts);
print(gcf,'-dpdf',[plotPath1 'HCMCLCHits_Polar'])

% Difference
Z3 = Z_H_HC-Z_H_LC;
th3 = angle(Z3);
rho3 = abs(Z3);

% with magnitude
opts = [];
opts.colors = hitColor;
%opts.maxR   = 4;
opts.markerSize=300;
opts.alpha  = 0.8;
PolarPlot(th3,rho3,opts);
print(gcf,'-dpdf',[plotPath1 'HC-LCHits_PolarRho'])

% Projected to the circle
opts.maxR   = 4/3;
opts.magText ='';
PolarPlot(th3,[],opts);
print(gcf,'-dpdf',[plotPath1 'HC-LCHits_Polar'])

% HC Hits Stats
[p,r] = circ_rtest(th_H);
[a,b,c]=circ_mean(th_H); a = mod(a*180/pi,360); b = mod(b*180/pi,360); c= mod(c*180/pi,360);
th_hat = mod(angle(mean(exp(1j*th_H))),2*pi)/pi*180;
fprintf('Rayleight Test for HC Hits modulation across subjects:\n')
fprintf(' angle = %g, r = %g ; p = %g \n',th_hat,r,p)
fprintf('lower bound = %g, upper bound = %g ; mean = %g \n',c,b,a)

% LC Hits Stats
[p,r] = circ_rtest(th_L);
[a,b,c]=circ_mean(th_L); a = mod(a*180/pi,360); b = mod(b*180/pi,360); c= mod(c*180/pi,360);
th_hat = mod(angle(mean(exp(1j*th_L))),2*pi)/pi*180;
fprintf('Rayleight Test for LC Hits modulation across subjects:\n')
fprintf(' angle = %g, r = %g ; p = %g \n',th_hat,r,p)
fprintf('lower bound = %g, upper bound = %g ; mean = %g \n',c,b,a)

% HC-LC Hits Stats
[p,r] = circ_rtest(th3);
[a,b,c]=circ_mean(th3); a = mod(a*180/pi,360); b = mod(b*180/pi,360); c= mod(c*180/pi,360);

th_hat = mod(angle(mean(exp(1j*th3))),2*pi)/pi*180;
fprintf('Rayleight Test for HC-LC Hits modulation across subjects:\n')
fprintf('angle = %g, r = %g ; p = %g \n',th_hat,r,p)
fprintf('lower bound = %g, upper bound = %g ; mean = %g \n',c,b,a)

fprintf('RankSum Test for greater modulation across subjects:\n')
[p,~,k]=signrank(rho_H,rho_L);
fprintf('HCH Median = %g; LCH Median = %g; Z = %g ; p = %g \n',median(rho_H),median(rho_L),k.zval,p)

makeLegend([LCcolor;HCcolor;hitColor],[300 300 300],{'LC Hits','HC Hits','HC-LC Hits'})
print(gcf,'-dpdf',[plotPath1 'HCLC_Legend'])

%% Misses HCvsLC 
close all; clc;
temp=brewermap(8,'dark2');
HCcolor = temp(3,:);
LCcolor = temp(4,:);
MCcolor = temp(5,:);
greyColor2 = [0.25 0.25 0.25];

% Polar Plots
clc; close all;
Z_M_HC = zeros(nSubjs,1);
Z_M_MC = zeros(nSubjs,1);
Z_M_LC = zeros(nSubjs,1);
for ss = 1:nSubjs
    Z_M_HC(ss) = f1(P(ss,M(ss,:)==1 & C(ss,:)==3));
    Z_M_MC(ss) = f1(P(ss,M(ss,:)==1 & C(ss,:)==2));
    Z_M_LC(ss) = f1(P(ss,M(ss,:)==1 & C(ss,:)==1));
end

% HC Misses
th_M_HC  = angle(Z_M_HC);
rho_M_HC = abs(Z_M_HC);

% LC Misses
th_M_LC  = angle(Z_M_LC);
rho_M_LC = abs(Z_M_LC);


% with magnitude of modulation
opts = [];
opts.colors = [HCcolor;LCcolor;];
%opts.maxR   = 4;
opts.markerSize=300;
opts.alpha  = 0.5;
PolarPlot([th_M_HC th_M_LC],[rho_M_HC rho_M_LC],opts);
print(gcf,'-dpdf',[plotPath1 'HCLCMiss_PolarRho'])

% along the unit circle
opts = [];
opts.colors = [HCcolor; LCcolor;];
opts.maxR   = 4/3;
opts.markerSize=300;
opts.alpha  = 0.8;
opts.magText = '';
PolarPlot([th_M_HC th_M_LC],[],opts);
print(gcf,'-dpdf',[plotPath1 'HCLCMiss_PolarRho'])

% Difference
z3 = Z_M_HC-Z_M_LC;
th3 = angle(z3);
rho3 = abs(z3);
opts        =[];
opts.maxR   = 4;
opts.markerSize=300;
opts.alpha  = 0.8;
opts.colors = missColor;
PolarPlot(th3,rho3,opts);
print(gcf,'-dpdf',[plotPath1 'HCLCMissRho_Polar'])

% Projected to the circle
opts.maxR   = 4/3;
opts.magText ='';
opts.alpha = 0.8;
PolarPlot(th3,[],opts);
print(gcf,'-dpdf',[plotPath1 'HC-LCMiss_Polar'])

% HC Misses Stats
[p,r] = circ_rtest(th_M_HC);
th_hat = mod(angle(mean(exp(1j*th_M_HC))),2*pi)/pi*180;
fprintf('Rayleight Test for LC Misses modulation across subjects:\n')
fprintf(' angle = %g, r = %g ; p = %g \n',th_hat,r,p)

% LC Misses Stats
[p,r] = circ_rtest(th_M_LC);
th_hat = mod(angle(mean(exp(1j*th_M_LC))),2*pi)/pi*180;
fprintf('Rayleight Test for HC Misses modulation across subjects:\n')
fprintf(' angle = %g, r = %g ; p = %g \n',th_hat,r,p)

% HC-LC Misses Stats
[p,r] = circ_rtest(th3);
th_hat = mod(angle(mean(exp(1j*th3))),2*pi)/pi*180;
fprintf('Rayleight Test for HC-LC Misses modulation across subjects:\n')
fprintf(' angle = %g, r = %g ; p = %g \n',th_hat,r,p)

makeLegend([LCcolor; HCcolor ;missColor],[300 300 300],{'LC Misses','HC Misses','HC-LC Misses'});
print(gcf,'-dpdf',[plotPath1 'HCLCMiss_Legend'])

%% All Confidence
close all; clc;
% confidence independent of memory:
Z_LC = zeros(nSubjs,1);
Z_MC = zeros(nSubjs,1);
Z_HC = zeros(nSubjs,1);
for ss = 1:nSubjs
    Z_LC(ss) = f1(P(ss,C(ss,:)==1));
    Z_MC(ss) = f1(P(ss,C(ss,:)==2));
    Z_HC(ss) = f1(P(ss,C(ss,:)==3));
end
th1 = angle(Z_LC); th2 = angle(Z_MC); th3=angle(Z_HC);

confColors = brewermap(3,'set3');
% all confidence
opts = [];
opts.colors = [confColors];
opts.maxR   = 4/3;
opts.markerSize=300;
opts.alpha  = 0.7;
opts.magText = '';
PolarPlot([th1 th2 th3],[],opts)
print(gcf,'-dpdf',[plotPath1 'allConf_Polar'])

% all confidence
opts = [];
opts.colors = [confColors];
%opts.maxR   = 4/3;
opts.markerSize=300;
opts.alpha  = 0.7;
opts.magText = '';
PolarPlot([th1 th2 th3],[rho1 rho2 rho3],opts)
print(gcf,'-dpdf',[plotPath1 'allConf_PolarRho'])

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

makeLegend([confColors],[300 300 300],{'Low','Mid','High'})
print(gcf,'-dpdf',[plotPath1 'Conf_Legend'])

%% Splits for Big and Small Obj

% Get data splits
Z_H_SmallBigConf    = zeros(nSubjs,2,3);
Z_H_SmallBig        = zeros(nSubjs,2);
Z_M_SmallBigConf    = zeros(nSubjs,2,3);
Z_M_SmallBig        = zeros(nSubjs,2);
Ns                  = zeros(nSubjs,2,3,2); % stimtype, conf, hit/miss
for ss = 1:nSubjs % subjects
    for st = 1:2 % stimulus type
        trials1 = (ST(ss,:)==st) & (H(ss,:)==1);
        Z_H_SmallBig(ss,st) = f1(P(ss,trials1));
        trials2 = (ST(ss,:)==st) & (M(ss,:)==1);
        Z_M_SmallBig(ss,st) = f1(P(ss,trials2));
        for cc = 1:3 % confidence        
            H_trials = trials1 & (C(ss,:)==cc);
            M_trials = trials2 & (C(ss,:)==cc);
            
            Z_H_SmallBigConf(ss,st,cc) = f1(P(ss,H_trials));
            Z_M_SmallBigConf(ss,st,cc) = f1(P(ss,M_trials));
            Ns(ss,st,cc,1) = sum(H_trials);
            Ns(ss,st,cc,2) = sum(M_trials);
        end
    end 
end
%% Plots for HC Small Hits 
close all; clc;
% Polar w Mag Plot of HC vs LC Hit Small
z1  = squeeze(Z_H_SmallBigConf(:,1,1)); %LC 
z2  = squeeze(Z_H_SmallBigConf(:,1,3)); %HC
th1 = angle(z1);    th2 = angle(z2);
rho1= abs(z1);      rho2 = abs(z2);

opts = [];
opts.colors = [LCcolor;HCcolor;];
opts.markerSize=300;
opts.alpha  = 0.5;
PolarPlot([th1 th2],[rho1 rho2],opts)
print(gcf,'-dpdf',[plotPath1 'HCLCSmallHits_PolarRho'])

% Projected to the circle
opts = [];
opts.colors = [LCcolor;HCcolor;];
opts.markerSize=300;
opts.maxR   = 4/3;
opts.magText ='';
opts.alpha  = 0.8;
PolarPlot([th1 th2],[],opts)
print(gcf,'-dpdf',[plotPath1 'HCLCSmallHits_Polar'])

% Polar Plot of HC - LC Hit Small
z3   = z2-z1;
th3  = angle(z3);
rho3 = abs(z3);
opts        =[];
opts.markerSize=300;
opts.alpha  = 0.8;
opts.colors = SmallColor;
PolarPlot(th3,rho3,opts);
print(gcf,'-dpdf',[plotPath1 'HC-LCSmallHits_PolarRho'])


opts        =[];
opts.maxR   = 4/3;
opts.markerSize=300;
opts.magText ='';
opts.colors = SmallColor;
PolarPlot(th3,[],opts)
print(gcf,'-dpdf',[plotPath1 'HC-LCSmallHits_Polar'])

%
% LC Hit Small Stats
[p,r] = circ_rtest(th1);
th_hat = mod(angle(mean(exp(1j*th1))),2*pi)/pi*180;
fprintf('Rayleight Test for LC Hit Small modulation across subjects:\n')
fprintf(' angle = %g, r = %g ; p = %g \n',th_hat,r,p)
[a,b,c]=circ_mean(th1); a = mod(a*180/pi,360); b = mod(b*180/pi,360); c= mod(c*180/pi,360);
fprintf('lower bound = %g, upper bound = %g ; mean = %g \n',floor(c),ceil(b),round(a))

% HC Hit Small Stats
[p,r] = circ_rtest(th2);
th_hat = mod(angle(mean(exp(1j*th2))),2*pi)/pi*180;
fprintf('Rayleight Test for HC Hit Small modulation across subjects:\n')
fprintf(' angle = %g, r = %g ; p = %g \n',th_hat,r,p)
[a,b,c]=circ_mean(th2); a = mod(a*180/pi,360); b = mod(b*180/pi,360); c= mod(c*180/pi,360);
fprintf('lower bound = %g, upper bound = %g ; mean = %g \n',floor(c),ceil(b),round(a))

% HC-LC Hit Small Stats
[p,r] = circ_rtest(th3);
th_hat = mod(angle(mean(exp(1j*th3))),2*pi)/pi*180;
fprintf('Rayleight Test for HC-LC Hit Small modulation across subjects:\n')
fprintf(' angle = %g, r = %g ; p = %g \n',th_hat,r,p)
[a,b,c]=circ_mean(th3); a = mod(a*180/pi,360); b = mod(b*180/pi,360); c= mod(c*180/pi,360);
fprintf('lower bound = %g, upper bound = %g ; mean = %g \n',floor(c),ceil(b),round(a))

makeLegend([LCcolor; HCcolor; SmallColor],[300 300 300],{'LC Scn Hits','HC Scn Hits','HC-LC Scn Hits'})
print(gcf,'-dpdf',[plotPath1 'HCLCSmallHits_Legend'])

rho2s = rho2;
%% Plots for HC Big Hits 
close all; clc;
% Polar w Mag Plot of HC vs LC Hit Big
z1  = squeeze(Z_H_SmallBigConf(:,2,1)); %LC
z2  = squeeze(Z_H_SmallBigConf(:,2,3)); %HC
th1 = angle(z1);    th2 = angle(z2);
rho1= abs(z1);      rho2 = abs(z2);

temp=brewermap(2,'set2');
HCcolor = temp(1,:);
LCcolor = temp(2,:);

opts = [];
opts.colors = [LCcolor;HCcolor;];
opts.markerSize=300;
opts.alpha  = 0.5;
PolarPlot([th1 th2],[rho1 rho2],opts)
print(gcf,'-dpdf',[plotPath1 'HCLCBigHits_PolarRho'])

% Projected to the circle
opts = [];
opts.colors = [LCcolor;HCcolor;];
opts.markerSize=300;
opts.maxR   = 4/3;
opts.magText ='';
opts.alpha  = 0.8;
PolarPlot([th1 th2],[],opts)
print(gcf,'-dpdf',[plotPath1 'HCLCBigHits_Polar'])

% Polar Plot of HC - LC Hit Big
z3   = z2-z1;
th3  = angle(z3);
rho3 = abs(z3);
opts        =[];
opts.markerSize=300;
opts.alpha  = 0.8;
opts.colors = BigColor;
PolarPlot(th3,rho3,opts);
print(gcf,'-dpdf',[plotPath1 'HC-LCBigHits_PolarRho'])

opts        =[];
opts.maxR   = 4/3;
opts.markerSize=300;
opts.magText ='';
opts.colors = BigColor;
PolarPlot(th3,[],opts)
print(gcf,'-dpdf',[plotPath1 'HC-LCBigHits_Polar'])

% LC Hit Big Stats
[p,r] = circ_rtest(th1);
th_hat = mod(angle(mean(exp(1j*th1))),2*pi)/pi*180;
fprintf('Rayleight Test for LC Hit Big modulation across subjects:\n')
fprintf(' angle = %g, r = %g ; p = %g \n',th_hat,r,p)
[a,b,c]=circ_mean(th1); a = mod(a*180/pi,360); b = mod(b*180/pi,360); c= mod(c*180/pi,360);
fprintf('lower bound = %g, upper bound = %g ; mean = %g \n',floor(c),ceil(b),round(a))

% HC Hit Big Stats
[p,r] = circ_rtest(th2);
th_hat = mod(angle(mean(exp(1j*th2))),2*pi)/pi*180;
fprintf('Rayleight Test for HC Hit Big modulation across subjects:\n')
fprintf(' angle = %g, r = %g ; p = %g \n',th_hat,r,p)
[a,b,c]=circ_mean(th2); a = mod(a*180/pi,360); b = mod(b*180/pi,360); c= mod(c*180/pi,360);
fprintf('lower bound = %g, upper bound = %g ; mean = %g \n',floor(c),ceil(b),round(a))

% HC-LC Hit Big Stats
[p,r] = circ_rtest(th3);
th_hat = mod(angle(mean(exp(1j*th3))),2*pi)/pi*180;
fprintf('Rayleight Test for HC-LC Hit Big modulation across subjects:\n')
fprintf(' angle = %g, r = %g ; p = %g \n',th_hat,r,p)
[a,b,c]=circ_mean(th3); a = mod(a*180/pi,360); b = mod(b*180/pi,360); c= mod(c*180/pi,360);
fprintf('lower bound = %g, upper bound = %g ; mean = %g \n',floor(c),ceil(b),round(a))

makeLegend([LCcolor; HCcolor; BigColor],[300 300 300],{'LC Scn Hits','HC Scn Hits','HC-LC Scn Hits'})
print(gcf,'-dpdf',[plotPath1 'HCLCBigHits_Legend'])

makeLegend([HCcolor; LCcolor; BigColor; SmallColor],[300 300 300 300],{'HC Hits','LC Hits','HC-LC Big','HC-LC Small' })
print(gcf,'-dpdf',[plotPath1 'HCLC_Cat_Hits_Legend'])
%% Plots for Small Hits vs Misses
close all; clc;
% Polar w Mag Plot of Hits vs Misses
z1  = squeeze(Z_H_SmallBig(:,1)); %LC
z2  = squeeze(Z_M_SmallBig(:,1)); %HC
th1 = angle(z1);    th2 = angle(z2);
rho1= abs(z1);      rho2 = abs(z2);

opts = [];
opts.colors = [hitColor;missColor;];
opts.markerSize=300;
opts.alpha  = 0.8;
PolarPlot([th1 th2],[rho1 rho2],opts)
print(gcf,'-dpdf',[plotPath1 'HM_Small_PolarRho'])

% Projected to the circle
opts = [];
opts.colors = [hitColor;missColor;];
opts.markerSize=300;
opts.maxR   = 4/3;
opts.magText ='';
opts.alpha  = 0.8;
PolarPlot([th1 th2],[],opts)
print(gcf,'-dpdf',[plotPath1 'HM_Small_Polar'])

% Polar Plot of HC - LC Hit Small
z3   = z1-z2;
th3  = angle(z3);
rho3 = abs(z3);
opts        =[];
opts.markerSize=300;
opts.alpha  = 0.8;
opts.colors = greyColor;
PolarPlot(th3,rho3,opts);
print(gcf,'-dpdf',[plotPath1 'H-M_Small_PolarRho'])

opts        =[];
opts.maxR   = 4/3;
opts.markerSize=300;
opts.magText ='';
opts.colors = greyColor;
PolarPlot(th3,[],opts)
print(gcf,'-dpdf',[plotPath1 'H-M_Small_Polar'])

% Hit Small Stats
[p,r] = circ_rtest(th1);
th_hat = mod(angle(mean(exp(1j*th1))),2*pi)/pi*180;
fprintf('Rayleight Test for Hit Small modulation across subjects:\n')
fprintf(' angle = %g, r = %g ; p = %g \n',th_hat,r,p)
[a,b,c]=circ_mean(th1); a = mod(a*180/pi,360); b = mod(b*180/pi,360); c= mod(c*180/pi,360);
fprintf('lower bound = %g, upper bound = %g ; mean = %g \n',floor(c),ceil(b),round(a))

% Miss Small Stats
[p,r] = circ_rtest(th2);
th_hat = mod(angle(mean(exp(1j*th2))),2*pi)/pi*180;
fprintf('Rayleight Test for Miss Small modulation across subjects:\n')
fprintf(' angle = %g, r = %g ; p = %g \n',th_hat,r,p)
[a,b,c]=circ_mean(th2); a = mod(a*180/pi,360); b = mod(b*180/pi,360); c= mod(c*180/pi,360);
fprintf('lower bound = %g, upper bound = %g ; mean = %g \n',floor(c),ceil(b),round(a))

% Hit-Miss Small Stats
[p,r] = circ_rtest(th3);
th_hat = mod(angle(mean(exp(1j*th3))),2*pi)/pi*180;
fprintf('Rayleight Test for Hit-Miss Small modulation across subjects:\n')
fprintf(' angle = %g, r = %g ; p = %g \n',th_hat,r,p)
[a,b,c]=circ_mean(th3); a = mod(a*180/pi,360); b = mod(b*180/pi,360); c= mod(c*180/pi,360);
fprintf('lower bound = %g, upper bound = %g ; mean = %g \n',floor(c),ceil(b),round(a))

fprintf('(Small )RankSum Test for greater of Hits than Misses modulation across subjects:\n')
[p,~,k]=signrank(rho1,rho2);
fprintf('Hits Median = %g; Miss Median = %g; Z = %g ; p = %g \n',median(rho1),median(rho2),k.zval,p)
rho1s = rho1; rho2s = rho2;
%% Plots for Big Hits vs Misses
close all; clc;
% Polar w Mag Plot of Hits vs Misses
z1  = squeeze(Z_H_SmallBig(:,2)); %LC
z2  = squeeze(Z_M_SmallBig(:,2)); %HC
th1 = angle(z1);    th2 = angle(z2);
rho1= abs(z1);      rho2 = abs(z2);

opts = [];
opts.colors = [hitColor;missColor;];
opts.markerSize=300;
opts.alpha  = 0.8;
PolarPlot([th1 th2],[rho1 rho2],opts)
print(gcf,'-dpdf',[plotPath1 'HM_Big_PolarRho'])

% Projected to the circle
opts = [];
opts.colors = [hitColor;missColor;];
opts.markerSize=300;
opts.maxR   = 4/3;
opts.magText ='';
opts.alpha  = 0.8;
PolarPlot([th1 th2],[],opts)
print(gcf,'-dpdf',[plotPath1 'HM_Big_Polar'])

% Polar Plot of HC - LC Hit Big
z3   = z1-z2;
th3  = angle(z3);
rho3 = abs(z3);
opts        =[];
opts.markerSize=300;
opts.alpha  = 0.8;
opts.colors = greyColor;
PolarPlot(th3,rho3,opts);
print(gcf,'-dpdf',[plotPath1 'H-M_Big_PolarRho'])

opts        =[];
opts.maxR   = 4/3;
opts.markerSize=300;
opts.magText ='';
opts.colors = greyColor;
PolarPlot(th3,[],opts)
print(gcf,'-dpdf',[plotPath1 'H-M_Big_Polar'])

% Hit Big Stats
[p,r] = circ_rtest(th1);
th_hat = mod(angle(mean(exp(1j*th1))),2*pi)/pi*180;
fprintf('Rayleight Test for Hit Big modulation across subjects:\n')
fprintf(' angle = %g, r = %g ; p = %g \n',th_hat,r,p)
[a,b,c]=circ_mean(th1); a = mod(a*180/pi,360); b = mod(b*180/pi,360); c= mod(c*180/pi,360);
fprintf('lower bound = %g, upper bound = %g ; mean = %g \n',floor(c),ceil(b),round(a))

% Miss Big Stats
[p,r] = circ_rtest(th2);
th_hat = mod(angle(mean(exp(1j*th2))),2*pi)/pi*180;
fprintf('Rayleight Test for Miss Big modulation across subjects:\n')
fprintf(' angle = %g, r = %g ; p = %g \n',th_hat,r,p)
[a,b,c]=circ_mean(th2); a = mod(a*180/pi,360); b = mod(b*180/pi,360); c= mod(c*180/pi,360);
fprintf('lower bound = %g, upper bound = %g ; mean = %g \n',floor(c),ceil(b),round(a))

% Hit-Miss Big Stats
[p,r] = circ_rtest(th3);
th_hat = mod(angle(mean(exp(1j*th3))),2*pi)/pi*180;
fprintf('Rayleight Test for Hit-Miss Big modulation across subjects:\n')
fprintf(' angle = %g, r = %g ; p = %g \n',th_hat,r,p)
[a,b,c]=circ_mean(th3); a = mod(a*180/pi,360); b = mod(b*180/pi,360); c= mod(c*180/pi,360);
fprintf('lower bound = %g, upper bound = %g ; mean = %g \n',floor(c),ceil(b),round(a))

fprintf('(Big )RankSum Test for greater of Hits than Misses modulation across subjects:\n')
[p,~,k]=signrank(rho1,rho2);
fprintf('Hits Median = %g; Miss Median = %g; Z = %g ; p = %g \n',median(rho1),median(rho2),k.zval,p)

fprintf('(RankSum Test for greater Small Hits than Big hits modulation across subjects:\n')
[p,~,k]=signrank(rho1s,rho1);
fprintf('Small Hits = %g; Big Hits = %g; Z = %g ; p = %g \n',median(rho1s),median(rho1),k.zval,p)


makeLegend([hitColor; missColor; BigColor],[300 300 300],{'Hits','Miss','H-M'})
print(gcf,'-dpdf',[plotPath1 'HMBig_Legend'])

makeLegend([hitColor; missColor; greyColor],[300 300 300],{'Hits','Misses','H-M'})
print(gcf,'-dpdf',[plotPath1 'HM_Legend'])
%% Bootstrap difference between Hits and Misses;
close all; clc; rng(1)

nBoot = 1000;
lowerB = round(nBoot*0.024);
upperB = round(nBoot*0.976);

% both categories
for ss = 1:nSubjs
    Z1(ss) = f1(P(ss,H(ss,:)==1));
    Z2(ss) = f1(P(ss,M(ss,:)==1));
end
th1     = mod(angle(Z1),2*pi);
th2     = mod(angle(Z2),2*pi);

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

print(gcf,'-dpdf',[plotPath1 'Hits-MissBoot_Polar'])
%
% Small
z1  = Z_H_SmallBig(:,1); % Hits
z2 = Z_M_SmallBig(:,1);  % Miss
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
opts.colors = SmallColor;
PolarPlot(B,1+0.05*randn(nBoot,1),opts)   
for jj = 1:2
    x = 1.2*cos(bounds(jj));
    y = 1.2*sin(bounds(jj));
    plot([0 x],[0 y],'linewidth',5,'color',greyColor)    
end
print(gcf,'-dpdf',[plotPath1 'Hit-MissSmallBoot_Polar'])

% Big
% Polar w Mag Plot of HC vs LC Hit Scenes
z1  = Z_H_SmallBig(:,2); % Hits
z2 = Z_M_SmallBig(:,2);  % Miss

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
opts.colors = BigColor;
PolarPlot(B,1+0.05*randn(nBoot,1),opts)   
for jj = 1:2
    x = 1.2*cos(bounds(jj));
    y = 1.2*sin(bounds(jj));
    plot([0 x],[0 y],'linewidth',5,'color',greyColor)    
end
print(gcf,'-dpdf',[plotPath1 'Hit-MissBigBoot_Polar'])

makeLegend([hitColor; SmallColor; BigColor],[300 300 300],{'All','Small ','Big'})
print(gcf,'-dpdf',[plotPath1 'All_SmBig_Legend'])

 %% Retrieval RT Modulation
% Get data splits
ZR_H                = zeros(nSubjs,1);
ZR_M                = zeros(nSubjs,1);
ZR_H_Conf           = zeros(nSubjs,3);
ZR_M_Conf           =zeros(nSubjs,3);
ZR_H_SmallBigConf    = zeros(nSubjs,2,3);
ZR_H_SmallBig        = zeros(nSubjs,2);
ZR_M_SmallBigConf    = zeros(nSubjs,2,3);
ZR_M_SmallBig        = zeros(nSubjs,2);
Ns                  = zeros(nSubjs,2,3,2); % stimtype, conf, hit/miss
for ss = 1:nSubjs % subjects
    ZR_H(ss) = f2(P(ss,H(ss,:)==1),RR(ss,H(ss,:)==1));
    ZR_M(ss) = f2(P(ss,M(ss,:)==1),RR(ss,M(ss,:)==1));
    for st = 1:2 % stimulus type
        trials1 = (ST(ss,:)==st) & (H(ss,:)==1);
        ZR_H_SmallBig(ss,st) = f2(P(ss,trials1),RR(ss,trials1));
        trials2 = (ST(ss,:)==st) & (M(ss,:)==1);
        ZR_M_SmallBig(ss,st) = f2(P(ss,trials2),RR(ss,trials2));
        for cc = 1:3 % confidence        
            H_trials = trials1 & (C(ss,:)==cc);
            M_trials = trials2 & (C(ss,:)==cc);
            
            ZR_H_SmallBigConf(ss,st,cc) = f2(P(ss,H_trials),RR(ss,H_trials));
            ZR_M_SmallBigConf(ss,st,cc) = f2(P(ss,M_trials),RR(ss,M_trials));
        end
    end
    for cc = 1:3
       trials1 = (H(ss,:)==1) & (C(ss,:)==cc);
       trials2 = (M(ss,:)==1) & (C(ss,:)==cc);
       ZR_H_Conf(ss,cc) = f2(P(ss,trials1),RR(ss,trials1)); 
       ZR_M_Conf(ss,cc) = f2(P(ss,trials2),RR(ss,trials2)); 
    end
end
%% RT Hit Misses

% rows-> hits, big hits, small hits
% columns -> hits/misses, hits/misses projected, hits-misses projected

% create data structure for easy iteration by row
X = cell(3,1);
X{1} = [ZR_H ZR_M];               
X{2} = [ZR_H_SmallBig(:,2) ZR_M_SmallBig(:,2)];     
X{3} = [ZR_H_SmallBig(:,1) ZR_M_SmallBig(:,1)]; 
Y = cell(3,3);
for ii = 1:3    
    Y{ii,1} = [angle(X{ii}) abs(X{ii,1})];
    Y{ii,2} = angle([X{ii}]);
    Y{ii,3} = angle([X{ii}(:,1)-X{ii}(:,2)]);     
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

colors = cell(2,1);
colors{1} = [hitColor;missColor];  
colors{2} = greyColor;

for ii = 1:3
    opts            =[];
    opts.markerSize = 150;
    opts.alpha      = 0.7;    
    for jj = 1:3
        opts.handle = axHan{ii,jj};
        if jj==3
            opts.colors = colors{2};
        else
            opts.colors = colors{1};
        end
        
        if jj ==1;            
            PolarPlot(Y{ii,jj}(:,1:2),Y{ii,jj}(:,3:4),opts);            
        else
            opts.magText    ='';
            opts.maxR       = 4/3;
            PolarPlot(Y{ii,jj},[],opts);            
        end
    end
end

print(gcf,'-dpdf',[plotPath1 'HM_RT_Phase'])
% Hit Stats
th_H = angle(X{1}(:,1)); th_M=angle(X{1}(:,2));
[p,r] = circ_rtest(th_H);
th_hat = mod(angle(mean(exp(1j*th_H))),2*pi)/pi*180;
fprintf('Rayleight Test for Hit RT modulation across subjects:\n')
fprintf(' angle = %g, r = %g ; p = %g \n',th_hat,r,p)
[a,b,c]=circ_mean(th_H); a = mod(a*180/pi,360); b = mod(b*180/pi,360); c= mod(c*180/pi,360);
fprintf('lower bound = %g, upper bound = %g ; mean = %g \n',floor(c),ceil(b),round(a))

% Miss Stats
[p,r] = circ_rtest(th_M);
th_hat = mod(angle(mean(exp(1j*th_M))),2*pi)/pi*180;
fprintf('Rayleight Test for Miss RT modulation across subjects:\n')
fprintf(' angle = %g, r = %g ; p = %g \n',th_hat,r,p)
[a,b,c]=circ_mean(th_M); a = mod(a*180/pi,360); b = mod(b*180/pi,360); c= mod(c*180/pi,360);
fprintf('lower bound = %g, upper bound = %g ; mean = %g \n',floor(c),ceil(b),round(a))

% H-M
 th3= Y{1,3} ;
 [p,r] = circ_rtest(th3);
th_hat = circ_mean(th3);
fprintf('Rayleight Test for Hit-Miss RT modulation across subjects:\n')
fprintf(' angle = %g, r = %g ; p = %g \n',th_hat,r,p)
[a,b,c]=circ_mean(th3); a = mod(a*180/pi,360); b = mod(b*180/pi,360); c= mod(c*180/pi,360);
fprintf('lower bound = %g, upper bound = %g ; mean = %g \n',floor(c),ceil(b),round(a))

% H-M Big
 th3= Y{2,3} ;
 [p,r] = circ_rtest(th3);
th_hat = circ_mean(th3);
fprintf('Rayleight Test for Big Hit-Miss RT modulation across subjects:\n')
fprintf(' angle = %g, r = %g ; p = %g \n',th_hat,r,p)
[a,b,c]=circ_mean(th3); a = mod(a*180/pi,360); b = mod(b*180/pi,360); c= mod(c*180/pi,360);
fprintf('lower bound = %g, upper bound = %g ; mean = %g \n',floor(c),ceil(b),round(a))
% H-M Small
 th3= Y{3,3} ;
 [p,r] = circ_rtest(th3);
th_hat = circ_mean(th3);
fprintf('Rayleight Test for Small Hit-Miss RT modulation across subjects:\n')
fprintf(' angle = %g, r = %g ; p = %g \n',th_hat,r,p)
[a,b,c]=circ_mean(th3); a = mod(a*180/pi,360); b = mod(b*180/pi,360); c= mod(c*180/pi,360);
fprintf('lower bound = %g, upper bound = %g ; mean = %g \n',floor(c),ceil(b),round(a))

%
makeLegend([hitColor;missColor;greyColor],[300 300 300],{'Hits','Misses','H-M'})
print(gcf,'-dpdf',[plotPath1 'HM_Legend2'])
makeLegend([HCcolor;LCcolor;greyColor],[300 300 300],{'HC Hits','LC Hits','HC-LC'})
print(gcf,'-dpdf',[plotPath1 'HCLC_Legend2'])
%% RT HC -LC Hits
X = cell(3,1);
X{1} = [ZR_H_Conf(:,3) ZR_H_Conf(:,1) ];               
X{2} = [ZR_H_SmallBigConf(:,2,3) ZR_M_SmallBigConf(:,2,1)];     
X{3} = [ZR_H_SmallBigConf(:,1,3) ZR_M_SmallBigConf(:,1,1)]; 
Y = cell(3,3);
for ii = 1:3    
    Y{ii,1} = [angle(X{ii}) abs(X{ii,1})];
    Y{ii,2} = angle([X{ii}]);
    Y{ii,3} = angle([X{ii}(:,1)-X{ii}(:,2)]);     
end

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

colors = cell(2,1);
colors{1} = [HCcolor;LCcolor];  
colors{2} = greyColor;


for ii = 1:3
    opts            =[];
    opts.markerSize = 150;
    opts.alpha      = 0.7;    
    for jj = 1:3
        opts.handle = axHan{ii,jj};
        if jj==3
            opts.colors = colors{2};
        else
            opts.colors = colors{1};
        end
        
        if jj ==1;            
            PolarPlot(Y{ii,jj}(:,1:2),Y{ii,jj}(:,3:4),opts);            
        else
            opts.magText    ='';
            opts.maxR       = 4/3;
            PolarPlot(Y{ii,jj},[],opts);            
        end
    end
end
print(gcf,'-dpdf',[plotPath1 'HCLCHits_RT_Phase'])

% HC H Stats
th_H = angle(X{1}(:,1)); th_L=angle(X{1}(:,2)); 
[p,r] = circ_rtest(th_H);
th_hat = circ_mean(th_H);
fprintf('Rayleight Test for HC Hit RT modulation across subjects:\n')
fprintf(' angle = %g, r = %g ; p = %g \n',th_hat,r,p)
[a,b,c]=circ_mean(th_H); a = mod(a*180/pi,360); b = mod(b*180/pi,360); c= mod(c*180/pi,360);
fprintf('lower bound = %g, upper bound = %g ; mean = %g \n',floor(c),ceil(b),round(a))

% LC H Stats
[p,r] = circ_rtest(th_L);
th_hat = circ_mean(th_L);
fprintf('Rayleight Test for LC Hit RT modulation across subjects:\n')
fprintf(' angle = %g, r = %g ; p = %g \n',th_hat,r,p)
[a,b,c]=circ_mean(th_L); a = mod(a*180/pi,360); b = mod(b*180/pi,360); c= mod(c*180/pi,360);
fprintf('lower bound = %g, upper bound = %g ; mean = %g \n',floor(c),ceil(b),round(a))

% HC-LC H Stats
th3= Y{1,3};
[p,r] = circ_rtest(th3);
th_hat = circ_mean(th3);
fprintf('Rayleight Test for HC-LC Hit RT modulation across subjects:\n')
fprintf(' angle = %g, r = %g ; p = %g \n',th_hat,r,p)
[a,b,c]=circ_mean(th3); a = mod(a*180/pi,360); b = mod(b*180/pi,360); c= mod(c*180/pi,360);
fprintf('lower bound = %g, upper bound = %g ; mean = %g \n',floor(c),ceil(b),round(a))

% HC-LC H Big Stats
th3= Y{2,3}; 
th3(isnan(th3))=[];
[p,r] = circ_rtest(th3);
th_hat = circ_mean(th3);
fprintf('Rayleight Test for HC-LC Big Hit RT modulation across subjects:\n')
fprintf(' angle = %g, r = %g ; p = %g \n',th_hat,r,p)
[a,b,c]=circ_mean(th3); a = mod(a*180/pi,360); b = mod(b*180/pi,360); c= mod(c*180/pi,360);
fprintf('lower bound = %g, upper bound = %g ; mean = %g \n',floor(c),ceil(b),round(a))

% HC-LC H Big Stats
th3= Y{3,3};
[p,r] = circ_rtest(th3);
th_hat = circ_mean(th3);
fprintf('Rayleight Test for  HC-LC Small Hit RT modulation across subjects:\n')
fprintf(' angle = %g, r = %g ; p = %g \n',th_hat,r,p)
[a,b,c]=circ_mean(th3); a = mod(a*180/pi,360); b = mod(b*180/pi,360); c= mod(c*180/pi,360);
fprintf('lower bound = %g, upper bound = %g ; mean = %g \n',floor(c),ceil(b),round(a))


