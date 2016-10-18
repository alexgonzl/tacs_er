function vonMissesPlotExamples(mu,kappa)
%%
% plot examples of circular distributions

rng(3); close all;
% von misses distribution for continous distribution
vm = @(th,mu,k)(exp(k*cos(th-mu))/(2*pi*besseli(0,k)));
th = linspace(0,2*pi,1000);
%mu = pi/2;
mu = 5*pi/4;
kp  = 0.5;
rho = vm(th,mu,kp);

N       = 100;
th_hat  = randsample(th,N,true,rho)'; 


opts            = [];
opts.maxR       = 4/3;
opts.colors     = [255 180 150]/255;
opts.alpha      = 0.9;
han=PolarPlot(th_hat, ones(N,1)+randn(N,1)*0.05,opts);
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[500 500],'paperposition',[0 0 500 500])
set(gcf,'position',[50,50,400,400])
set(gca,'units','points','position',[0 0 400 400]); 
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Simulations/ProbabPlot_mu' ....
    num2str(round(mu/pi*180)) '_kp' strrep(num2str(kp),'.','p') '_v2']);


%plot for theta at stimulus onse presentation locations
th = 0:2*pi/5:(2*pi-0.1);
rho = vm(th,mu,kp);
th_hat  = randsample(th,N,true,rho)'; 
han=PolarPlot(th_hat+unifrnd(0,36/180*pi,N,1), ones(N,1)+randn(N,1)*0.05,opts);
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[500 500],'paperposition',[0 0 500 500])
set(gcf,'position',[50,50,400,400])
set(gca,'units','points','position',[0 0 400 400]); 
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Simulations/ProbabPlot_mu' ....
    num2str(round(mu/pi*180)) 'kp' strrep(num2str(kp),'.','p') '_v3']);

% plot the pdf
NN = 1000;
th = linspace(0,2*pi,NN);
rho = vm(th,mu,kp);

xlocks = 0:pi/2:2*pi;
xlockslabel = xlocks/pi*180;

figure(3);
plot(th,rho,'color',[255 180 150]/255,'linewidth',5)
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[500 300],'paperposition',[0 0 500 300])
set(gcf,'position',[0,0,500,300])
set(gca,'units','points','position',[80 80 400 200]); 
set(gca,'fontsize',20,'linewidth',2,'box','off')
set(gca,'xtick',xlocks,'XTickLabel',xlockslabel)
xlabel(' \theta (degrees) ')
ylabel([' f(\theta; \mu=' num2str(round(mu/pi*180)) ',\kappa=' num2str(kp) ')' ])
axis tight
grid on
ylim([-0.01 0.6])
print(gcf, '-dpdf', ['../plots/tacs_enc_xdiva/Simulations/ProbabPlot_mu' ....
    num2str(round(mu/pi*180)) 'kp' strrep(num2str(kp),'.','p') '_v1']);

%% Population Simulations
rng(3); close all;

nSimulations = 1000;
nSubjs = 20;
trueMu = 5*pi/8;
kappaRange = [0 0.8];









