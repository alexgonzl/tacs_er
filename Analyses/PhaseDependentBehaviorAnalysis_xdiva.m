
dataPath    = '~/Google Drive/Research/tACS/tACS_ER_task/data/tacs_enc_xdiva/';
nSubjs      = 20;

% load behavioral data
load([dataPath 'Summary/BehavSummary.mat']) 
%%
out                         = [];
out.nSubjs                  = nSubjs;
out.EventPhase              = behav_out.Trials.EncTrialPhase; % phase at encoding
out.EventPhaseCond          = behav_out.Trials.EncTrialPhaseCondition;
out.EventEncCondAtRet       = behav_out.Trials.EncCondAtRet; 
out.EncStimAtRet            = behav_out.Trials.EncStimIDAtRet;
out.RetStimAtEnc            = nan(nSubjs,300);
out.HitMissEncTrialIDs         = nan(nSubjs,2,300);
out.ValidEncTrialIDs           = nan(nSubjs,300);

out.HitMissPhases           = cell(nSubjs,2); % first column is hits, second misses
out.HitMissPhaseCond        = cell(nSubjs,2); % first column is hits, second misses
out.HitMisRTSubj            = cell(nSubjs,2); % first column is hits, second misses
out.ConfidenceByCond        = cell(nSubjs,2); % first column is hits, second misses
out.PhasesByConf            = cell(nSubjs,3); % high to low

% circular distribution stats
out.HCRDistFromUniform      = nan(nSubjs,2);
out.ConfDistFromUniform     = nan(nSubjs,2);
out.DiffInDist              = nan(nSubjs,2);

for ss = 1:out.nSubjs
        
    % Sort trials by their position at encoding.   
    EncIDsAtRet = out.EncStimAtRet(ss,:);
    [s,i] = sort(EncIDsAtRet);
    out.RetIDsAtEnc(ss,:) = i(s>0);
    
    % Trial IDs for hits and misses.
    hits                            = behav_out.retSubj{ss}.Hits;   % hits IDs @ retrieval
    misses                          = behav_out.retSubj{ss}.Misses; % miss IDs @ retrieval
    % get hits and misses IDs at encoding.
    out.HitMissEncTrialIDs(ss,1,:)  = hits(out.RetIDsAtEnc(ss,:));
    out.HitMissEncTrialIDs(ss,2,:)  = misses(out.RetIDsAtEnc(ss,:));
    out.ValidEncTrialIDs(ss,:)      = out.HitMissEncTrialIDs(ss,1,:) | out.HitMissEncTrialIDs(ss,2,:);

    % Get true phase for Hits and Misses
    out.HitMissPhases{ss,1} = out.EventPhase(ss,out.HitMissEncTrialIDs(ss,1,:));    
    out.HitMissPhases{ss,2} = out.EventPhase(ss,out.HitMissEncTrialIDs(ss,2,:));    
    
    % Get phase condition for Hits and Misses
    out.HitMissPhaseCond{ss,1} = out.EventPhaseCond(ss,out.HitMissEncTrialIDs(ss,1,:));    
    out.HitMissPhaseCond{ss,2} = out.EventPhaseCond(ss,out.HitMissEncTrialIDs(ss,2,:));    
    
    % Get Confidence and assign sign based on hit/miss
    Confidence                  = behav_out.retSubj{ss}.Confidence;
    out.ConfidenceByCond(ss,1)  = Confidence(hits);
    out.ConfidenceByCond(ss,2)  = Confidence(misses)*-1;
    
    %CondPerSubj{ss,1} = EventPhaseCond(ss,behav_out.EncRet.EncHitTrialIDs{ss});
    
    Confidence(~(hits|misses))  = [];            
    
    
   % EncStimAtRet(~(hits|misses))=[];
    %ConfBySubj{ss} = Confidence;
end
% %%
% % circularStatsByCond       = cell(nSubjs,2);
% % HitsFixedEffects        = []; 
% % CRsFixedEffects         = [];
% 
% % startEdge               = tacs_er.EncPhases;
% % binWidth                = 2*pi/10;
% % %phaseEdges              = (-pi:binWidth:pi);
% % nPhaseBins              = numel(startEdge);
% % conds                   = 1:5;
% % %phaseEdges
% % ConfBySubj              = cell(nSubjs,1);
% % PhasesBySubj            = cell(nSubjs,1);
% % PhasesBySubjCounts      = zeros(nSubjs,nPhaseBins);
% % PhaseBySubjCond1        = zeros(nSubjs,nPhaseBins);
% % PhaseBySubjCond2        = zeros(nSubjs,nPhaseBins);
% 
% % propHits        = nan(nSubjs,nPhaseBins);
% % HitRateByPhase  = nan(nSubjs,nPhaseBins);
% % P               = nan(nSubjs,nPhaseBins);
% % ConfByPhase     = cell(nSubjs,1);
% % ConfByPhaseNorm = cell(nSubjs,1);
% 
% for ss = 1:nSubjs
%     
%     EncStimAtRet    = behav_out.retSubj{ss}.EncStimIDAtRet;
%     hits            = behav_out.retSubj{ss}.Hits;
%     misses          = behav_out.retSubj{ss}.Misses;
%     
%     Confidence                  = behav_out.retSubj{ss}.Confidence;
%     Confidence(misses)          = Confidence(misses)*-1;
%     Confidence(~(hits|misses))  = [];            
%     
%     EncStimAtRet(~(hits|misses))=[];
%     ConfBySubj{ss} = Confidence;
%     
%     % Hits
%     PhasesPerSubj{ss,1} = EventPhase(ss,behav_out.EncRet.EncHitTrialIDs{ss});    
%     CondPerSubj{ss,1} = EventPhaseCond(ss,behav_out.EncRet.EncHitTrialIDs{ss});    
%     %temp = histc(PhasesPerCondSubj{ss,1},-pi:binWidth:pi);
%     %PhaseBySubjCond1(ss,:) = temp(1:end-1);    
%     
%     %Misses
%     %PhasesPerCondSubj{ss,2} = subjEventPhases{ss}.out.TrueAngleStims(behav_out.EncRet.EncMissTrialIDs{ss});
%     PhasesPerSubj{ss,2} = EventPhase(ss,behav_out.EncRet.EncMissTrialIDs{ss});
%     CondPerSubj{ss,2} = EventPhaseCond(ss,behav_out.EncRet.EncMissTrialIDs{ss});    
%     %temp = histc(PhasesPerCondSubj{ss,2},-pi:binWidth:pi);
%     %PhaseBySubjCond2(ss,:) = temp(1:end-1);
% 
%     %
%     PhasesBySubj{ss} = EventPhase(ss,EncStimAtRet);   
%     CondBySubj{ss}   = EventPhaseCond(ss,EncStimAtRet);    
%     %temp = histc(PhasesBySubj{ss},-pi:binWidth:pi);
%     %PhasesBySubjCounts(ss,:)   = temp(1:end-1);
%     
%     [~,HCRDistFromUniform(ss,1)]   = circ_rtest(PhasesPerSubj{ss,1});
%     [~,HCRDistFromUniform(ss,2)]   = circ_rtest(PhasesPerSubj{ss,2});
%     
%     circularStatsByCond{ss,1}      = circ_stats(PhasesPerSubj{ss,1});
%     circularStatsByCond{ss,2}      = circ_stats(PhasesPerSubj{ss,2});
%     
%     me(ss,1)= circularStatsByCond{ss,1}.mean;
%     me(ss,2)= circularStatsByCond{ss,2}.mean;    
%        
%     CondRTSubj{ss,1}    = behav_out.EncRet.EncHitRTs{ss};
%     CondRTSubj{ss,2}    = behav_out.EncRet.EncMissRTs{ss};
%     
%     nHitsByPhase(ss,:)        = histc(CondPerSubj{ss,1},conds);
%     nMissesByPhase(ss,:)      = histc(CondPerSubj{ss,2},conds);
%     temp = nHitsByPhase./(nHitsByPhase+nMissesByPhase);
%     HitRateByPhase(ss,:) = temp;
%     
%     temp = nHitsByPhase/numel(CondPerSubj{ss,1});
%     propHits(ss,:) = temp;
%     
% end
% 
% %% trial counts by phase
% 
% figure(1); clf;
% set(gcf,'position',[100 100 600 400],'paperpositionmode','auto','color','w',...
% 'paperposition',[0.1 0.2 0.6 0.4],'paperunits','normalized');
% 
% t = linspace(0,1,1000);
% x = cos(2*pi*t-pi);
% xa = angle(hilbert(x));
% 
% a2 = axes('position',[0.12 0.1 0.8 0.15]);
% 
% axes(a2)
% plot((xa+pi)./pi*180,x,'k','linewidth',4)
% axis tight; 
% set(gca,'ytick',[],'ycolor','w','fontsize',16,'box','off','lineWidth',2)
% set(gca,'xtick',[0:60:360])
% xlabel(' Encoding Phase (deg)')
% grid on
% 
% a1 = axes('position', [0.12 0.3 0.8 0.6]);
% axes(a1); hold on;
% plot([30:60:360],PhasesBySubjCounts','-','color',[0.8 0.8 0.8])
% plot([30:60:360],mean(PhasesBySubjCounts), 'k','linewidth',5)
% set(gca,'fontsize',16,'box','off','lineWidth',2)
% set(gca,'xtick',[0:30:360],'xTickLabel','')
% xlim([0 360])
% ylabel(' trial counts' )
% 
% print(gcf, '-dpdf', '../plots/trialCountsByPhase_xdiva');
% 
% %% hits and misses raw #s
% 
% figure(2); clf;
% set(gcf,'position',[100 100 600 400],'paperpositionmode','auto','color','w',...
% 'paperposition',[0.1 0.2 0.6 0.4],'paperunits','normalized');
% 
% t = linspace(0,1,1000);
% x = cos(2*pi*t-pi);
% xa = angle(hilbert(x));
% 
% a2 = axes('position',[0.12 0.1 0.8 0.15]);
% axes(a2);
% plot((xa+pi)./pi*180,x,'k','linewidth',4)
% axis tight; 
% set(gca,'ytick',[],'ycolor','w','fontsize',16,'box','off','lineWidth',2)
% set(gca,'xtick',[0:60:360])
% xlabel(' Encoding Phase (deg)')
% grid on
% 
% a1 = axes('position', [0.12 0.3 0.8 0.6]);
% axes(a1); hold on;
% plot([30:60:360],PhaseBySubjCond1','-','color',[255 180 150]/255)
% plot([30:60:360],mean(PhaseBySubjCond1), 'color',[240 80 40]/255,'linewidth',5)
% set(gca,'fontsize',16,'box','off','lineWidth',2)
% set(gca,'xtick',[0:30:360],'xTickLabel','')
% xlim([0 360])
% ylabel(' trial counts (hits) ' )
% 
% print(gcf, '-dpdf', '../plots/HitsCountsByPhase');
% 
% %
% 
% figure(3); clf;
% set(gcf,'position',[100 100 600 400],'paperpositionmode','auto','color','w',...
% 'paperposition',[0.1 0.2 0.6 0.4],'paperunits','normalized');
% 
% a2 = axes('position',[0.12 0.1 0.8 0.15]);
% axes(a2)
% plot((xa+pi)./pi*180,x,'k','linewidth',4)
% axis tight; 
% set(gca,'ytick',[],'ycolor','w','fontsize',16,'box','off','lineWidth',2)
% set(gca,'xtick',[0:60:360])
% xlabel(' Encoding Phase (deg)')
% grid on
% 
% a1 = axes('position', [0.12 0.3 0.8 0.6]);
% axes(a1); hold on;
% X = PhaseBySubjCond1./repmat(sum(PhaseBySubjCond1,2),[1,6]);
% plot([30:60:360],X','-','color',[255 180 150]/255)
% plot([30:60:360],mean(X), 'color',[240 80 40]/255,'linewidth',5)
% set(gca,'fontsize',16,'box','off','lineWidth',2)
% set(gca,'xtick',[0:30:360],'xTickLabel','')
% xlim([0 360])
% ylim([0 0.3])
% ylabel(' proportion (hits) ' )
% 
% print(gcf, '-dpdf', '../plots/HitsPropByPhase');
% 
% %% hits and misses, proportions
% figure(4); clf;
% set(gcf,'position',[100 100 600 400],'paperpositionmode','auto','color','w',...
% 'paperposition',[0.1 0.2 0.6 0.4],'paperunits','normalized');
% 
% t = linspace(0,1,1000);
% x = cos(2*pi*t-pi);
% xa = angle(hilbert(x));
% 
% a2 = axes('position',[0.12 0.1 0.8 0.15]);
% 
% axes(a2)
% plot((xa+pi)./pi*180,x,'k','linewidth',4)
% axis tight; 
% set(gca,'ytick',[],'ycolor','w','fontsize',16,'box','off','lineWidth',2)
% set(gca,'xtick',[0:60:360])
% xlabel(' Encoding Phase (deg)')
% grid on
% 
% a1 = axes('position', [0.12 0.3 0.8 0.6]);
% axes(a1); hold on;
% plot([30:60:360],PhaseBySubjCond2','-','color',[180 250 250]/255)
% plot([30:60:360],mean(PhaseBySubjCond2), 'color',[120 200 200]/255,'linewidth',5)
% set(gca,'fontsize',16,'box','off','lineWidth',2)
% set(gca,'xtick',[0:30:360],'xTickLabel','')
% xlim([0 360])
% ylabel(' trial counts (misses) ' )
% 
% print(gcf, '-dpdf', '../plots/MissesCountsByPhase');
% 
% figure(5); clf;
% set(gcf,'position',[100 100 600 400],'paperpositionmode','auto','color','w',...
% 'paperposition',[0.1 0.2 0.6 0.4],'paperunits','normalized');
% 
% a2 = axes('position',[0.12 0.1 0.8 0.15]);
% axes(a2)
% plot((xa+pi)./pi*180,x,'k','linewidth',4)
% axis tight; 
% set(gca,'ytick',[],'ycolor','w','fontsize',16,'box','off','lineWidth',2)
% set(gca,'xtick',[0:60:360])
% xlabel(' Encoding Phase (deg)')
% grid on
% 
% a1 = axes('position', [0.12 0.3 0.8 0.6]);
% axes(a1); hold on;
% X = PhaseBySubjCond2./repmat(sum(PhaseBySubjCond2,2),[1,6]);
% plot([30:60:360],X','-','color',[180 250 250]/255)
% plot([30:60:360],mean(X), 'color',[120 200 200]/255,'linewidth',5)
% set(gca,'fontsize',16,'box','off','lineWidth',2)
% set(gca,'xtick',[0:30:360],'xTickLabel','')
% xlim([0 360])
% ylim([0 0.3])
% ylabel(' proportion (misses) ' )
% 
% print(gcf, '-dpdf', '../plots/MissesPropByPhase');
% 
% figure(6); clf;
% 
% set(gcf,'position',[100 100 600 400],'paperpositionmode','auto','color','w',...
% 'paperposition',[0.1 0.2 0.6 0.4],'paperunits','normalized');
% 
% a2 = axes('position',[0.12 0.1 0.8 0.15]);
% axes(a2)
% plot((xa+pi)./pi*180,x,'k','linewidth',4)
% axis tight; 
% set(gca,'ytick',[],'ycolor','w','fontsize',16,'box','off','lineWidth',2)
% set(gca,'xtick',[0:60:360])
% xlabel(' Encoding Phase (deg)')
% grid on
% 
% a1 = axes('position', [0.12 0.3 0.8 0.6]);
% axes(a1); hold on;
% X = PhaseBySubjCond1./repmat(sum(PhaseBySubjCond1,2),[1,6]);
% Y = PhaseBySubjCond2./repmat(sum(PhaseBySubjCond2,2),[1,6]);
% plot([30:60:360],mean(X), 'color',[240 80 40]/255,'linewidth',5)
% plot([30:60:360],mean(Y), 'color',[120 200 200]/255,'linewidth',5)
% set(gca,'fontsize',16,'box','off','lineWidth',2)
% set(gca,'xtick',[0:30:360],'xTickLabel','')
% xlim([0 360])
% ylim([0.1 0.25])
% 
% ylabel(' proportion of trials ' )
% 
% print(gcf, '-dpdf', '../plots/HMPropByPhase');
% 
% %% distance from uniform
% figure(8); clf;
% set(gcf,'position',[100 100 600 400],'paperpositionmode','auto','color','w',...
% 'paperposition',[0.1 0.2 0.6 0.4],'paperunits','normalized');
% a2 = axes('position',[0.15 0.1 0.3 0.8]); hold on;
% for ss = 1:nSubjs
%     h1 = plot([1 2],[me(ss,1) me(ss,2)]/pi*180,'o-','linewidth',2,'color',[0.5 0.5 0.5],...
%         'MarkerFaceColor',[0.5 0.5 0.5],'markeredgecolor','k','MarkerSize',10);
% end
% set(a2,'fontSize',16,'lineWidth',3,'xtick',[1 2],'xticklabel',{'hits','misses'})
% set(a2,'ytick',[-pi:pi/3:pi]./pi*180)
% ylabel(' mean Phase (�) ')
% plot([0.8 1.2], [1 1]*mean(me(:,1))/pi*180,'linewidth',5,'color','k')
% plot([1.8 2.2], [1 1]*mean(me(:,2))/pi*180,'linewidth',5,'color','k')
% grid on
% xlim([0.6 2.4])
% 
% a3 = axes('position',[0.6 0.1 0.35 0.8]); hold on;
% for ss = 1:nSubjs
%     h1 = plot([1 2],[HCRDistFromUniform(ss,1) HCRDistFromUniform(ss,2)],'o-','linewidth',2,'color',[0.5 0.5 0.5],...
%         'MarkerFaceColor',[0.5 0.5 0.5],'markeredgecolor','k','MarkerSize',10);
% end
% set(a3,'fontSize',16,'lineWidth',3,'xtick',[1 2],'xticklabel',{'D(H)','D(M)'})
% set(a3,'ytick',[0:0.5:4])
% ylim([0 4])
% ylabel(' Z-Stat ')
% plot([0.8 1.2], [1 1]*mean(HCRDistFromUniform(:,1)),'linewidth',5,'color','k')
% plot([1.8 2.2], [1 1]*mean(HCRDistFromUniform(:,2)),'linewidth',5,'color','k')
% 
% xlim([0.6 2.4])
% grid on
% print(gcf, '-dpdf', '../plots/DistFromUniformByCond');
% 
% %% 
% 
% %% modulation hit rate and memory score
% 
% 
% %%
% stepSize   = pi/6;
% phaseEdges = 0:%(-pi:stepSize:pi);
% nPhaseBins = 5;%numel(phaseEdges);
% 
% propHits        = nan(nSubjs,nPhaseBins);
% propMisses      = nan(nSubjs,nPhaseBins);
% HitRateByPhase  = nan(nSubjs,nPhaseBins);
% H_Hits          = nan(nSubjs,1);
% P               = nan(nSubjs,nPhaseBins);
% ConfByPhase     = cell(nSubjs,1);
% ConfByPhaseNorm = cell(nSubjs,1);
% 
% confScores = [-3 -2 -1 1 2 3]; nConfScores = numel(confScores);
% phaseMemScore = nan(nSubjs-1,nPhaseBins-1);
% for ss=1:nSubjs
%     nHitsByPhase        = histc(PhasesPerCondSubj{ss,1},conds);
%     nMissesByPhase      = histc(PhasesPerCondSubj{ss,2},conds);
%     HitRateByPhase(ss-1,:) = nHitsByPhase./(nHitsByPhase+nMissesByPhase);
%     
%     propHits(ss-1,:) = nHitsByPhase/numel(PhasesPerSubj{ss,1});
%     P(ss-1,:) = HitRateByPhase(ss-1,:)-nanmean(HitRateByPhase(ss-1,:));    
%     
%     propMisses(ss-1,:) = nMissesByPhase/numel(PhasesPerSubj{ss,2});
%     
%     counter = 1;
%     for cc=confScores
%         ConfByPhase{ss-1}(counter,:) = histc(PhasesBySubj{ss}(ConfBySubj{ss}==cc),conds);                
%         counter = counter+1;
%     end    
%     [~, Idx ]=histc(PhasesBySubj{ss},phaseEdges);
%     for ph = 1:(nPhaseBins-1)
%         phaseMemScore(ss-1,ph)=sum(ConfBySubj{ss}(Idx==ph));
%     end
%     
%      ConfByPhase{ss-1}(:,nPhaseBins)=[];
%     ConfByPhaseNorm{ss-1}=ConfByPhase{ss-1}./repmat(sum(ConfByPhase{ss-1}),[nConfScores 1]);
%     ConfByPhaseNorm{ss-1} =  ConfByPhaseNorm{ss-1}./repmat(sum(ConfByPhaseNorm{ss-1},2),[1 nPhaseBins-1]);
% %     
% end
% 
% phaseMemScoreNorm=phaseMemScore-repmat(mean(phaseMemScore,2),[1 nPhaseBins-1]);
% 
% P(:,nPhaseBins)=[];
% HitRateByPhase(:,nPhaseBins)=[];
% 
% propHits(:,nPhaseBins)=[];
% propMisses(:,nPhaseBins)=[];
% 
% pu = 1/(nPhaseBins-1)*ones(nPhaseBins-1,1);
% sum(-pu.*log2(pu))
% %X = propHits-propMisses;
% 
% %plot(diff(cumsum(phaseEdges)),HitRateByPhase-repmat(mean(HitRateByPhase,2),[1,nPhaseBins-1]),'-*')
