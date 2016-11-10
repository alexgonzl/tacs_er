function behav_out = EncRetBehavAnalyses(expt,nSubjs)
% behavioral analysis function for tACS_ER task
% analyzes both encoding and retrieval data for tacs_enc
%
% tune for experiment tacs_er_objstim
% thePath input determines the experimental folder to load data from
%
%------------------------------------------------------------------------%
% Author:         Alex Gonzalez
% Created:        Dec 1, 2015
% LastUpdate:     Oct 26, 2016
%------------------------------------------------------------------------%
%
% encoding task:
% Perceptual Discrimination
% 1) Big hit rate
% 2) Small hit rate
% 3) Discrimination rate (accuracy) (stats across subjects)
% 4) RTs of for Bigs and Smalls; (stats across subjects)
%
% retrieval task:
% 1) hit rate, fa rate, cr rate, miss rate
% 2) hit rate, fa rate, cr rate, miss rate for Bigs and Smalls
% independently
% 3) RTs
% 4) Confidence
%

dataPath = ['~/Google Drive/Research/tACS/tACS_ER_task/data/' expt '/'];
behav_out     = [];
behav_out.encSubj = cell(nSubjs,1);
behav_out.retSubj = cell(nSubjs,1);

%% Pre-allocations
% encoding stats pre-allocation
behav_out.encSummary = [];
behav_out.encSummary.goodSubj           = nan(nSubjs,1);
behav_out.encSummary.meanAcc            = nan(nSubjs,1);
behav_out.encSummary.BigHR              = nan(nSubjs,1);
behav_out.encSummary.SmallHR            = nan(nSubjs,1);
behav_out.encSummary.medianRTs          = nan(nSubjs,1);
behav_out.encSummary.medianBigRTs       = nan(nSubjs,1);
behav_out.encSummary.medianSmallRTs     = nan(nSubjs,1);
behav_out.encSummary.BigSmall_RTs_TVals = nan(nSubjs,1);    
        
% retrieval stats pre-allocation
behav_out.retSummary =[];
behav_out.retSummary.dPrime             = nan(nSubjs,1);
behav_out.retSummary.dPrime_C           = nan(nSubjs,1);
behav_out.retSummary.Big_dPrime        = nan(nSubjs,1);
behav_out.retSummary.Big_dPrime_C      = nan(nSubjs,1);
behav_out.retSummary.Small_dPrime       = nan(nSubjs,1);
behav_out.retSummary.Small_dPrime_C     = nan(nSubjs,1);

% reaction time summaries
behav_out.retSummary.medianHit_RTs     = nan(nSubjs,1);
behav_out.retSummary.medianCRs_RTs     = nan(nSubjs,1);
behav_out.retSummary.medianFA_RTs      = nan(nSubjs,1);
behav_out.retSummary.medianMiss_RTs    = nan(nSubjs,1);
behav_out.retSummary.HitCRs_RTs_TVals   = nan(nSubjs,1);

% confidence summaries
behav_out.retSummary.dPrimeConf             = nan(nSubjs,3);
behav_out.retSummary.dPrimeConf_C           = nan(nSubjs,3);
behav_out.retSummary.meanAccuracyByConf     = nan(nSubjs,3);
behav_out.retSummary.nRespByConf            = nan(nSubjs,3);
behav_out.retSummary.nH_nMiss_nFA_nCRs      = nan(nSubjs,3,4);
behav_out.retSummary.medianHit_RTsConf     = nan(nSubjs,3);
behav_out.retSummary.medianCRs_RTsConf         = nan(nSubjs,3);
behav_out.retSummary.medianFA_RTsConf          = nan(nSubjs,3);
behav_out.retSummary.medianMiss_RTsConf        = nan(nSubjs,3);

% trials
behav_out.EncRet = [];
behav_out.EncRet.EncHitTrialIDs         = cell(nSubjs,1);
behav_out.EncRet.EncMissTrialIDs        = cell(nSubjs,1);
behav_out.EncRet.EncHitRTs              = cell(nSubjs,1);
behav_out.EncRet.EncMissRTs             = cell(nSubjs,1);

behav_out.Trials =[]; nEncTrials = 400; nRetTrials = 600;
behav_out.Trials.EncTrialPhaseCond      = nan(nSubjs,nEncTrials);
behav_out.Trials.EncTrialPhase          = nan(nSubjs,nEncTrials);
behav_out.Trials.RetStimIDAtEnc         = nan(nSubjs,nEncTrials);
behav_out.Trials.StimTypeEnc            = nan(nSubjs,nEncTrials); %
behav_out.Trials.EncCondAtRet           = nan(nSubjs,nRetTrials);
behav_out.Trials.EncStimIDAtRet         = nan(nSubjs,nRetTrials);

%%
for ss = 1:nSubjs
    %% Encoding
    try
        if strcmp(expt,'tacs_er_objstim')
            load([dataPath 's' num2str(ss) '/tacs_er_objstim.enc.mat'])
            load([dataPath 's' num2str(ss) '/tacs_er_objstim.task.mat'])
            behav_out.encSubj{ss}                  = EncAnalysisBySubj(enc_out);       
        else
            error('wrong expt')
        end        
        behav_out.encSummary.meanAcc(ss)            = behav_out.encSubj{ss}.meanAccuracy;
        behav_out.encSummary.BigHR(ss)              = behav_out.encSubj{ss}.BigHitRate;
        behav_out.encSummary.SmallHR(ss)            = behav_out.encSubj{ss}.SmallHitRate;
        behav_out.encSummary.medianRTs(ss)          = behav_out.encSubj{ss}.RTsStats.medianRT;
        behav_out.encSummary.medianBigRTs(ss)       = behav_out.encSubj{ss}.BigRTsStats.medianRT;
        behav_out.encSummary.medianSmallRTs(ss)     = behav_out.encSubj{ss}.SmallRTsStats.medianRT;
        behav_out.encSummary.BigSmall_RTs_TVals(ss) ...
            = behav_out.encSubj{ss}.RTsStats_Big_Small.tstat;
        
    catch msg   
        fprintf('issues processing subject %i \n',ss)
    end
    %% Retrieval
    try
        
        if strcmp(expt,'tacs_er_objstim')
            load([dataPath 's' num2str(ss) '/tacs_er_objstim.test.mat'])
            load([dataPath 's' num2str(ss) '/EventsPhase.mat']);
            behav_out.Trials.PhaseInfo          = out;
            behav_out.Trials.AllStimPhases(ss,:)    = out.FzFC2Phases;      
            [~,i]=intersect(out.AllStimIdx,out.SmallIdx);
            behav_out.Trials.SmallStimPhases(ss,:)  = out.FzFC2Phases(i);
            [~,i]=intersect(out.AllStimIdx,out.BigIdx);
            behav_out.Trials.BigStimPhases(ss,:)    = out.FzFC2Phases(i);
        end
         behav_out.retSubj{ss}                           = RetAnalysisBySubj(ret_out);
        
        % IDs of test items during encoding.
        behav_out.Trials.EncStimIDAtRet(ss,:)           = behav_out.retSubj{ss}.EncStimIDAtRet;
        [s,i]                                           = sort(behav_out.Trials.EncStimIDAtRet(ss,:));
        behav_out.Trials.RetStimIDAtEnc(ss,:)           = i(s>0);
        behav_out.Trials.StimTypeEnc(ss,:)              = ret_out.expInfo.EncStimType;
        
        % DPrimes and accuracy summary measures
        behav_out.retSummary.dPrime(ss)                 = behav_out.retSubj{ss}.dPrime;
        behav_out.retSummary.dPrime_C(ss)               = behav_out.retSubj{ss}.dPrime_C;
        behav_out.retSummary.Big_dPrime(ss)            = behav_out.retSubj{ss}.Big_dPrime;
        behav_out.retSummary.Big_dPrime_C(ss)          = behav_out.retSubj{ss}.Big_dPrime_C;
        behav_out.retSummary.Small_dPrime(ss)           = behav_out.retSubj{ss}.Small_dPrime;
        behav_out.retSummary.Small_dPrime_C(ss)         = behav_out.retSubj{ss}.Small_dPrime_C;
        behav_out.retSummary.meanAccuracyByConf(ss,:)   = behav_out.retSubj{ss}.meanAccuracyByConf;
        behav_out.retSummary.dPrimeConf(ss,:)           = behav_out.retSubj{ss}.dPrimeConf;
        behav_out.retSummary.dPrimeConf_C(ss,:)         = behav_out.retSubj{ss}.dPrimeConf_C;
        
        behav_out.retSummary.Big_dPrimeConf(ss,:)      = behav_out.retSubj{ss}.Big_dPrimeConf;
        behav_out.retSummary.Big_dPrimeConf_C(ss,:)    = behav_out.retSubj{ss}.Big_dPrimeConf_C;
        behav_out.retSummary.Small_dPrimeConf(ss,:)     = behav_out.retSubj{ss}.Small_dPrimeConf;
        behav_out.retSummary.Small_dPrimeConf_C(ss,:)   = behav_out.retSubj{ss}.Small_dPrimeConf_C;
        
        % numbers of trials per conditions
        behav_out.retSummary.nH_nMiss_nFA_nCRs(ss,:,:)  = behav_out.retSubj{ss}.nH_nMiss_nFA_nCRs;
        behav_out.retSummary.nRespByConf(ss,:)          = behav_out.retSubj{ss}.nRespByConf;
        
        % reaction times        
        behav_out.retSummary.medianHit_RTs(ss)          = behav_out.retSubj{ss}.Hit_RTsStats.medianRT;
        behav_out.retSummary.medianCRs_RTs(ss)          = behav_out.retSubj{ss}.CRs_RTsStats.medianRT;
        behav_out.retSummary.medianFA_RTs(ss)           = behav_out.retSubj{ss}.FA_RTsStats.medianRT;
        behav_out.retSummary.medianMiss_RTs(ss)         = behav_out.retSubj{ss}.Misses_RTsStats.medianRT;
        behav_out.retSummary.HitCRs_RTs_TVals(ss)       = behav_out.retSubj{ss}.HitCRs_RTs.tstat;
        
        behav_out.retSummary.medianSmallBigHit_RTs(ss,:) = [behav_out.retSubj{ss}.Small_Hit_RTsStats.medianRT,   behav_out.retSubj{ss}.Big_Hit_RTsStats.medianRT];
        behav_out.retSummary.medianSmallBigCRs_RTs(ss,:) = [behav_out.retSubj{ss}.Small_CRs_RTsStats.medianRT,   behav_out.retSubj{ss}.Big_CRs_RTsStats.medianRT];
        behav_out.retSummary.medianSmallBigFA_RTs(ss,:)  = [behav_out.retSubj{ss}.Small_FA_RTsStats.medianRT,    behav_out.retSubj{ss}.Big_FA_RTsStats.medianRT];
        behav_out.retSummary.medianSmallBigMiss_RTs(ss,:)= [behav_out.retSubj{ss}.Small_Misses_RTsStats.medianRT, behav_out.retSubj{ss}.Big_Misses_RTsStats.medianRT];
        
        behav_out.retSummary.medianHit_RTsConf(ss,:)    = [behav_out.retSubj{ss}.Hit_RTsStatsConf.medianRT];
        behav_out.retSummary.medianCRs_RTsConf(ss,:)    = [behav_out.retSubj{ss}.CRs_RTsStatsConf.medianRT];
        behav_out.retSummary.medianFA_RTsConf(ss,:)     = [behav_out.retSubj{ss}.FA_RTsStatsConf.medianRT];
        behav_out.retSummary.medianMiss_RTsConf(ss,:)   = [behav_out.retSubj{ss}.Misses_RTsStatsConf.medianRT];
      
        % trials
        behav_out.EncRet.EncHitTrialIDs{ss}         = behav_out.retSubj{ss}.EncHitTrialIDs;
        behav_out.EncRet.EncMissTrialIDs{ss}        = behav_out.retSubj{ss}.EncMissTrialIDs;
        behav_out.EncRet.EncHCTrials{ss}            = behav_out.retSubj{ss}.EncHCTrials;
        behav_out.EncRet.EncMCTrials{ss}            = behav_out.retSubj{ss}.EncMCTrials;
        behav_out.EncRet.EncLCTrials{ss}            = behav_out.retSubj{ss}.EncLCTrials;
        behav_out.EncRet.EncHitRTs{ss}              = behav_out.retSubj{ss}.EncHitRTs;
        behav_out.EncRet.EncMissRTs{ss}             = behav_out.retSubj{ss}.EncMissRTs;
        
    catch msg
        disp(msg)
    end
    
end
save([dataPath '/Summary/BehavSummary'],'behav_out')
end

%% Auxiliary functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Encoding Analyses by subject
function enc=EncAnalysisBySubj(enc_out)
enc =[];
enc.encCond = enc_out.exptInfo.EncStimType; % Big/Small trials
enc.RTs     = enc_out.TimingInfo.trialRT;
keyPress    = enc_out.TimingInfo.trialKeyPress;

cond1Key=enc_out.PresParams.RespToCond2;
cond2Key=enc_out.PresParams.RespToCond1;

% get simple responses information
binOut = tabulateBinaryResponses(keyPress,enc.encCond,cond1Key,cond2Key);

enc.CorrectResp     = binOut.Correct;
enc.InCorrectResp   = binOut.InCorrect;
enc.meanAccuracy    = binOut.mACC;
enc.BigHitRate     = binOut.Cond1HR;
enc.SmallHitRate    = binOut.Cond2HR;
enc.BigCorrectResp = binOut.Cond1Correct;
enc.SmallCorrectResp= binOut.Cond2Correct;
enc.otherEncInfo    = binOut;

% RT info (for correct)
enc.RTsStats        = analyzeRTs(enc.RTs(enc.CorrectResp));
enc.BigRTsStats    = analyzeRTs(enc.RTs(enc.BigCorrectResp));
enc.SmallRTsStats   = analyzeRTs(enc.RTs(enc.SmallCorrectResp));
[~,a1,~,a2] = ttest2(enc.RTs(enc.BigCorrectResp),enc.RTs(enc.SmallCorrectResp));
enc.RTsStats_Big_Small.p=a1;
enc.RTsStats_Big_Small.tstat=a2.tstat;

end

% Retrieval Analyses
function ret=RetAnalysisBySubj(ret_out)

ret = [];
ret.retCond     = ret_out.expInfo.RetCondTrialCode;
% 1 old Small, 2 old Big, 3 new Small, 4 new Big
ret.OldTrials   = ret.retCond==1 | ret.retCond==2;
ret.NewTrials   = ret.retCond==3 | ret.retCond==4;
ret.SmallTrials = ret.retCond==1 | ret.retCond==3;
ret.BigTrials   = ret.retCond==2 | ret.retCond==4;

[~,ret.EncStimIDAtRet]=ismember(ret_out.expInfo.RetStimNames,ret_out.expInfo.EncStimNames);


ret.RTs         = ret_out.TimingInfo.trialRT;

keyPress        = ret_out.TimingInfo.trialKeyPress;

cond1Key=ret_out.PresParams.RespButtons(1);
cond2Key=ret_out.PresParams.RespButtons(3);

% get overall response information
binOut = tabulateBinaryResponses(keyPress,1*ret.OldTrials+2*ret.NewTrials,cond1Key,cond2Key);

ret.CorrectResp     = binOut.Correct;
ret.InCorrectResp   = binOut.InCorrect;
ret.meanAccuracy    = binOut.mACC;

ret.HitRate         = binOut.Cond1HR;
ret.CorrRejRate     = binOut.Cond2HR;
ret.FARate          = binOut.Cond2FAR;
ret.MissRate        = binOut.Cond1FAR;

ret.Hits            = binOut.Cond1Correct;
ret.CRs             = binOut.Cond2Correct;
ret.Misses          = binOut.Cond1InCorrect;
ret.FA              = binOut.Cond2InCorrect;

[ret.dPrime, ret.dPrime_C] = calc_dPrime(sum(ret.Hits),sum(ret.Misses),sum(ret.FA),sum(ret.CRs));

% RT analyses
ret.Hit_RTsStats     = analyzeRTs(ret.RTs(ret.Hits));
ret.CRs_RTsStats     = analyzeRTs(ret.RTs(ret.CRs));
ret.Misses_RTsStats  = analyzeRTs(ret.RTs(ret.Misses));
ret.FA_RTsStats      = analyzeRTs(ret.RTs(ret.FA));

[~,a1,~,a2]          = ttest2(ret.RTs(ret.Hits),ret.RTs(ret.CRs));
ret.HitCRs_RTs.p     = a1;
ret.HitCRs_RTs.tstat =a2.tstat;

% separated by Bigs & Smalls
types = {'Big','Small'};
for tt = 1:numel(types)
    type = types{tt};
    trials = ret.([type 'Trials']);
    
    binOut = tabulateBinaryResponses(keyPress(trials),1*ret.OldTrials(trials)+2*ret.NewTrials(trials),cond1Key,cond2Key);
    
    ret.([type '_meanACC'])      = binOut.mACC;
    ret.([type '_HitRate'])      = binOut.Cond1HR;
    ret.([type '_CorrRejRate'])  = binOut.Cond2HR;
    ret.([type '_FARate'])       = binOut.Cond2FAR;
    ret.([type '_MissRate'])     = binOut.Cond1FAR;
    
    [ret.([type '_dPrime']), ret.([type '_dPrime_C'])]= calc_dPrime(sum(binOut.Cond1Correct),....
        sum(binOut.Cond1InCorrect),sum(binOut.Cond2InCorrect),sum(binOut.Cond2Correct));
    
    % RT analyses
    ret.([type '_Hit_RTsStats'])     = analyzeRTs(ret.RTs(ret.Hits & trials));
    ret.([type '_CRs_RTsStats'])     = analyzeRTs(ret.RTs(ret.CRs & trials));
    ret.([type '_Misses_RTsStats'])  = analyzeRTs(ret.RTs(ret.Misses & trials));
    ret.([type '_FA_RTsStats'])      = analyzeRTs(ret.RTs(ret.FA & trials));
end

% Confidence
ret.ConfidenceResp = ret_out.TimingInfo.ConfResp;
ret.Confidence     = 3*strcmp(ret.ConfidenceResp,'high')+ ...
    2*strcmp(ret.ConfidenceResp,'mid')+1*strcmp(ret.ConfidenceResp,'low');
ret.nRespByConf    = histc(ret.Confidence,1:3);

% Tabulate Hits, FA, CRs and Misses as a function of confidence
for cc = 1:3
    % get overall response information
    trials = ret.Confidence==cc;
    
    binOut = tabulateBinaryResponses(keyPress(trials),...
        1*ret.OldTrials(trials)+2*ret.NewTrials(trials),...
        cond1Key,cond2Key);
    
    ret.meanAccuracyByConf(cc)    = binOut.mACC;
    
    ret.HitRateByConf(cc)         = binOut.Cond1HR;
    ret.CorrRejRateByConf(cc)     = binOut.Cond2HR;
    ret.FARateByConf(cc)          = binOut.Cond2FAR;
    ret.MissRateByConf(cc)        = binOut.Cond1FAR;
    
    x1           = binOut.Cond1Correct; % hits
    x2           = binOut.Cond1InCorrect; % misses
    x3           = binOut.Cond2InCorrect; % # FA
    x4           = binOut.Cond2Correct; % # CRs

    % #s of hits/misses/fa/CRs
    ret.nH_nMiss_nFA_nCRs(cc,:)   = sum([x1 x2 x3 x4]);
    
    % compute dprime
    ret.dPrimeConf(cc)      = binOut.dPrime;
    ret.dPrimeConf_C(cc)    = binOut.dPrimeC;   
    
    % get RT stats per confidence.
    rts = ret.RTs(trials);
    ret.Hit_RTsStatsConf(cc)     = analyzeRTs(rts(x1));
    ret.CRs_RTsStatsConf(cc)     = analyzeRTs(rts(x2));
    ret.Misses_RTsStatsConf(cc)  = analyzeRTs(rts(x3));
    ret.FA_RTsStatsConf(cc)      = analyzeRTs(rts(x4));
    
    % compute dPrime for Bigs and Smalls
    types = {'Big','Small'};
    for tt = 1:numel(types)
        type = types{tt};
        trials2 = ret.([type 'Trials']) & trials;
        binOut = tabulateBinaryResponses(keyPress(trials2),...
        1*ret.OldTrials(trials2)+2*ret.NewTrials(trials2),...
        cond1Key,cond2Key);
        
        ret.([type '_dPrimeConf'])(cc)   = binOut.dPrime;
        ret.([type '_dPrimeConf_C'])(cc) = binOut.dPrimeC;
    end
end

ret.EncHitTrialIDs     = ret.EncStimIDAtRet(ret.Hits);
ret.EncHCHitTrialIDs   = ret.EncStimIDAtRet(ret.Hits & ret.Confidence==3);
ret.EncMCHitTrialIDs   = ret.EncStimIDAtRet(ret.Hits & ret.Confidence==2);
ret.EncLCHitTrialIDs   = ret.EncStimIDAtRet(ret.Hits & ret.Confidence==1);
ret.EncMissTrialIDs    = ret.EncStimIDAtRet(ret.Misses );
ret.EncHCMissTrialIDs  = ret.EncStimIDAtRet(ret.Misses & ret.Confidence==3);
ret.EncMCMissTrialIDs  = ret.EncStimIDAtRet(ret.Misses & ret.Confidence==2);
ret.EncLCMissTrialIDs  = ret.EncStimIDAtRet(ret.Misses & ret.Confidence==1);
ret.EncHCTrials        = ret.EncStimIDAtRet(ret.OldTrials & ret.Confidence==3);
ret.EncMCTrials        = ret.EncStimIDAtRet(ret.OldTrials & ret.Confidence==2);
ret.EncLCTrials        = ret.EncStimIDAtRet(ret.OldTrials & ret.Confidence==1);
ret.EncHitRTs          = ret.RTs(ret.Hits);
ret.EncMissRTs         = ret.RTs(ret.Misses);


end

% Other functions.
function out=tabulateBinaryResponses(responses,correctAns,key1,key2)

out = [];
fields = {'multipleResp','noResp','Cond1Correct','Cond2Correct',...
    'Cond1InCorrect','Cond2InCorrect','Cond1HR','Cond2HR','Cond1FAR',...
    'Cond2FAR','Correct','InCorrect','mACC','dPrime','dPrimeC'};
for ii = 1:numel(fields)
    out.(fields{ii}) = [];
end

% handle special cases
% for more than one answer, use the first button response
if any(cellfun(@iscell,responses))
    ids = find(cellfun(@iscell,responses));
    for ii = 1:numel(ids)
        responses(ids(ii))=responses{ids(ii)}(1);
    end
    out.multipleResp = ids;
end

% no response trials
if any(cellfun(@isempty,responses))
    out.noResp=find(cellfun(@isempty,responses));
end

try
    % tabulate correct/incorrect by cond
    out.Cond1Correct   = strcmp(responses,key1)&(correctAns==1); % hits
    out.Cond1InCorrect = strcmp(responses,key2)&(correctAns==1); % misses
    out.Cond2Correct   =  strcmp(responses,key2)&(correctAns==2); % FA
    out.Cond2InCorrect = strcmp(responses,key1)&(correctAns==2); % CRs
    
    out.Cond1HR  = sum(out.Cond1Correct)/sum(correctAns==1);  % correct response to cond 1
    out.Cond2HR  = sum(out.Cond2Correct)/sum(correctAns==2);   % correct response to cond 2
    out.Cond1FAR = sum(out.Cond1InCorrect)/(sum(out.Cond1InCorrect)+sum(out.Cond1Correct)); % responded as condition 1 (when it was condition 2)
    out.Cond2FAR = sum(out.Cond2InCorrect)/(sum(out.Cond2InCorrect)+sum(out.Cond2Correct)); % responded as condition 2 (when it was condition 1)
    
    out.Correct   = out.Cond1Correct | out.Cond2Correct;
    out.InCorrect = out.Cond1InCorrect | out.Cond2InCorrect;
    out.mACC = mean(out.Correct);
    
    x1           = out.Cond1Correct; % hits
    x1s          = sum(x1);
    x2           = out.Cond1InCorrect; % misses
    x2s          = sum(x2);
    x3           = out.Cond2InCorrect; % # FA
    x3s          = sum(x3); % # FA
    x4           = out.Cond2Correct; % # CRs
    x4s          = sum(x4);
    
    % compute dprime
    [out.dPrime, out.dPrimeC] = calc_dPrime(x1s,x2s,x3s,x4s);
    
catch msg
    keyboard
end

end

% version when resp and ans match.
function out=tabulateBinaryResponses2(resp,stim)

out = [];
fields = {'multipleResp','noResp','Cond1Correct','Cond2Correct',...
    'Cond1InCorrect','Cond2InCorrect','Cond1HR','Cond2HR','Cond1FAR',...
    'Cond2FAR','Correct','InCorrect','mACC','dPrime','dPrimeC'};
for ii = 1:numel(fields)
    out.(fields{ii}) = [];
end

% no response trials
out.noResp = (resp==0);

try
    % tabulate correct/incorrect by cond
    out.Cond1Correct   = (resp==1)&(stim==1); % hits
    out.Cond1InCorrect = (resp==2)&(stim==1); % misses
    out.Cond2Correct   = (resp==2)&(stim==2); % FA
    out.Cond2InCorrect = (resp==1)&(stim==2); % CRs
    
    out.Cond1HR  = sum(out.Cond1Correct)/sum(stim==1);  % correct response to cond 1
    out.Cond2HR  = sum(out.Cond2Correct)/sum(stim==2);   % correct response to cond 2
    out.Cond1FAR = sum(out.Cond1InCorrect)/(sum(out.Cond1InCorrect)+sum(out.Cond1Correct)); % responded as condition 1 (when it was condition 2)
    out.Cond2FAR = sum(out.Cond2InCorrect)/(sum(out.Cond2InCorrect)+sum(out.Cond2Correct)); % responded as condition 2 (when it was condition 1)
    
    out.Correct   = out.Cond1Correct | out.Cond2Correct;
    out.InCorrect = out.Cond1InCorrect | out.Cond2InCorrect;
    out.mACC = mean(out.Correct);
    
    x1           = out.Cond1Correct; % hits
    x1s          = sum(x1);
    x2           = out.Cond1InCorrect; % misses
    x2s          = sum(x2);
    x3           = out.Cond2InCorrect; % # FA
    x3s          = sum(x3); % # FA
    x4           = out.Cond2Correct; % # CRs
    x4s          = sum(x4);
    
    % compute dprime
    [out.dPrime, out.dPrimeC] = calc_dPrime(x1s,x2s,x3s,x4s);
    
catch msg
    keyboard
end

end
function out=analyzeRTs(RTs)
out.meanRT = nanmean(RTs);
out.sdRT = nanstd(RTs);
out.medianRT = nanmedian(RTs);
out.kurtosisRT = kurtosis(RTs);
end

