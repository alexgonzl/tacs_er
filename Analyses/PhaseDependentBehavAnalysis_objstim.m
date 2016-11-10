function out = DataMat_objstim()
expt = 'tacs_er_objstim';

% load behavioral data
dataPath    = ['~/Google Drive/Research/tACS/tACS_ER_task/data/' expt '/'];
load([dataPath 'Summary/BehavSummary.mat'])
addpath CircStats


%%
out                         = [];
out.nEncTrials              = 400;
out.nSubjs                  = 22;

% matrix with all conditions for all subjects
% 1st   column: small=1/big=2
% 2     column: correct small/big categorization (bool)
% 3     column: reaction times for perceptual decision
% 4     column: Phase at encoding(rads)
% 5     column: Hits [(binary) -> subsequently remember
% 6     column: Misses (binary -> subsequenty forgotten
% 7     column: Confidence-> subsequent confidence in hits /misses
% 8     column: reaction times on memory decision
% 9     column: position of encoding item at retrieval (for correctly
%               sorting retrieval data

out.datMatColumnNames      = ...
    {'StimType','EncCorrect','EncRTs','PhaseRad',...
    'Hit','Miss','Confidence','RetRTs','RetPos'};


out.EncEventPhases          = behav_out.Trials.AllStimPhases;
out.EncStimAtRet            = behav_out.Trials.EncStimIDAtRet;
out.RetStimAtEnc            = behav_out.Trials.RetStimIDAtEnc;
% Pre-allocation of data
out.datMat                  = nan(out.nSubjs,out.nEncTrials,numel(out.datMatColumnNames));
out.datMat(:,:,4)           = behav_out.Trials.AllStimPhases; %a all phases
out.datMat(:,:,12)          = behav_out.Trials.RetStimIDAtEnc;

%%
for ss = 1:out.nSubjs
    try
        %% Additional entries for data matrix
        
        % Correct encoding responses 
        out.datMat(ss,:,1)           = behav_out.encSubj{ss}.encCond;
        out.datMat(ss,:,3)           = behav_out.encSubj{ss}.CorrectResp;        
        out.datMat(ss,:,4)          = behav_out.encSubj{ss}.RTs; 
        %Get retrieval RTs
        out.datMat(ss,:,11)         = behav_out.retSubj{ss}.RTs(out.RetStimAtEnc(ss,:));
        
        % Trial IDs for hits and misses.        
        out.datMat(ss,:,7) = behav_out.retSubj{ss}.Hits(out.RetStimAtEnc(ss,:));   % hits IDs at encoding
        out.datMat(ss,:,8) = behav_out.retSubj{ss}.Misses(out.RetStimAtEnc(ss,:)); % miss IDs at encoding
        

        % Get Confidence and assign sign based on hit/miss
        Confidence                  = behav_out.retSubj{ss}.Confidence;
        Confidence(Confidence==0)   = nan;
        out.datMat(ss,:,9)          = Confidence(out.RetStimAtEnc(ss,:));
           
    catch msg
        warning(sprintf('on subject %i',ss))
        disp(msg)
    end
end

save([dataPath 'Summary/PhaseDependentAnalyses.mat'],'out')
end
