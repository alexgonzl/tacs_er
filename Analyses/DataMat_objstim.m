function out = DataMat_objstim(nsubjs)
expt = 'tacs_er_objstim';

% load behavioral data
dataPath    = ['~/Google Drive/Research/tACS/tACS_ER_task/data/' expt '/'];
load([dataPath 'Summary/BehavSummary.mat'])

%%
out                         = [];
out.nEncTrials              = 400;
out.nSubjs                  = nsubjs;

% matrix with all conditions for all subjects
% 1st   column: small=1/big=2
% 2     column: correct small/big categorization (bool)
% 3     column: reaction times for perceptual decision
% 4     column: Phase at encoding(rads) Pz
% 5     column: Hits [(binary) -> subsequently remember
% 6     column: Misses (binary -> subsequenty forgotten
% 7     column: Confidence-> subsequent confidence in hits /misses
% 8     column: reaction times on memory decision
% 9     column: position of encoding item at retrieval (for correctly
%               sorting retrieval data
% 10    column: Phase at encoding(rads) Fz-FC2 (this is based on the
% stimulation template file)

out.datMatColumnNames      = ...
    {'StimType','EncCorrect','EncRTs','PhaseRad',...
    'Hit','Miss','Confidence','RetRTs','RetPos','FzPhaseRad'};

% Pre-allocation of data
out.datMat                  = nan(out.nSubjs,out.nEncTrials,numel(out.datMatColumnNames));
out.datMat(:,:,4)           = mod(behav_out.Trials.AllStimPhases1,2*pi); % all phases
% position of encoding item at retrieval (for sorting hits/misses and RTs)
out.datMat(:,:,9)           = behav_out.Trials.RetStimIDAtEnc; 
out.datMat(:,:,10)          = mod(behav_out.Trials.AllStimPhases2,2*pi); 

%%
for ss = 1:out.nSubjs
    try
        %% Additional entries for data matrix        
        % Correct encoding responses 
        out.datMat(ss,:,1)           = behav_out.encSubj{ss}.encCond;
        out.datMat(ss,:,2)           = behav_out.encSubj{ss}.CorrectResp;        
        out.datMat(ss,:,3)           = behav_out.encSubj{ss}.RTs; 
        
        encPos = out.datMat(ss,:,9);
        % Trial IDs for hits and misses.        
        out.datMat(ss,:,5) = behav_out.retSubj{ss}.Hits(encPos);   % hits IDs at encoding
        out.datMat(ss,:,6) = behav_out.retSubj{ss}.Misses(encPos); % miss IDs at encoding
        
        % Get Confidence and assign sign based on hit/miss
        Confidence                  = behav_out.retSubj{ss}.Confidence;
        Confidence(Confidence==0)   = nan;
        out.datMat(ss,:,7)          = Confidence(encPos);
        
        %Get retrieval RTs
        out.datMat(ss,:,8)         = behav_out.retSubj{ss}.RTs(encPos);
           
    catch msg
        warning(sprintf('on subject %i',ss))
        disp(msg)
    end
end

save([dataPath 'Summary/DataMatrix.mat'],'out')
end
