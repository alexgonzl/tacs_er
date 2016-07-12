
function [tacs_er] = tACS_ER_makelist_ObjSize_xDivaStim(thePath)
% make lists for tacs encoding & retrieval (tacs_er) experiment
% This design:
% 1) One encoding block with flashing object stimuli.
% 2) Responses will be made to the stimulus itself, no cue conditions.
% 3) The decision will be a smaller or larger than a shoe box for each object. 
%
% Important variables
% subjNum       -> subject  ID
% thePath       -> path for all stimulus and data directories
%
% nEncStim      -> # of UNIQUE encoding stimuli
% nEncTrials    -> # of TOTAL trials (should be multiple of nEncStim)
% nEncCond      -> # of conditions at encoding (that are critical for retrieval)
% nEncBlocks    -> # of encoding blocks
%
% nRetTrials    -> # of TOTAL retrieval trials
% nFoilTrials   -> # of NEW trials (foils)
% nRetConds     -> # of retrieval conditions
% nRetBlocks    -> # of retrieval blocks
% nRetCondTrials-> VECTOR of # trials per Retrieal Condition
% maxNumConsecutiveOld -> MAX # of consecutive old trials
%
% RandStream    -> random stream seed
%
% EncStimType       -> VECTOR trial of indicating type (1  faces, 2 scene)
% EncStimTypeIDs    -> {face,scene}
% EncCondCodeIDs    -> {'Face0','Face90','Face180','Scn0','Scn90','Scn180'}
% EncCondTrialCode  -> VECTOR trial encoding indicating EncCondCodeID
% EncStimNames      -> VECTOR trial stimulus ID
% EncStimUniqueIDs  -> VECTOR trial with unique ID for each stim
% EncBlockID        -> VECTOR trial block ID
%
% RetCondIDs        -> {'OldFace','OldScene','NewFace','NewScene'};
% RetCondTrialCode  -> VECTOR trial 1:4, indicating RetCondID
% RetStimNames      -> VECTOR trial of retrieaval IDs
% RetBlockID        -> VECTOR of trial block ID
%
% Stimuli           -> key, matrix map of stimuli's names to image matrix.

%------------------------------------------------------------------------%
% Author:       Alex Gonzalez
% Created:      July 6, 2016
% LastUpdate:   July 6, 2016
%------------------------------------------------------------------------%

%% Set-up

% parameter list
nTotalStim          = 600;
nEncLargeStim       = 200;
nEncSmallStim       = 200;
nRetNewSmallStim    = 100;
nRetNewLargeStim    = 100;
assert((nEncLargeStim+nEncSmallStim+nRetNewLargeStim+nRetNewSmallStim==nTotalStim),'incorrect numbers of stims')
nSmallStim = nEncSmallStim + nRetNewSmallStim;
nLargeStim = nEncLargeStim + nRetNewLargeStim;

nEncStim    = nEncLargeStim + nEncSmallStim;
nEncPhases = 5;                 % constrained by xdiva design.
nEncConds  = 2*nEncPhases;      % small/large x phase (5 different phases)
EncPhases  = 0:(360/nEncPhases):359;

nEncTrials = nEncStim;
nRetTrials = nTotalStim;   % encoding trials + retrieval trials
nFoilTrials = nRetTrials-nEncStim;
nRetConds  = 4;             % old and new conditons per stim type (small/large)
maxNumConsecutiveOld = 8;   % maximum number of old trials in a row

% Set Random Seed based on SubjectNum
s = RandStream.create('mt19937ar','seed',thePath.subjNum);
RandStream.setGlobalStream(s)

%% load data names
% large
cd(fullfile(thePath.stim,'/Bigger'));
temp = dir('*.jpeg');
count = 1;
LargeObjMat     = cell(nLargeStim,1);
LargeObjNames   = cell(nLargeStim,1);
nFiles = size(temp,1);
randFileNums = datasample(1:nFiles,nFiles,'replace',false);
for n = randFileNums
    try
        x = imread(temp(n).name); %check if readable
        LargeObjMat{count} = imresize(x(:,:,1),2); % only resizing case for xdiva
        LargeObjNames(count) = temp(n);
        count = count + 1;
    catch
        fprintf(['\n' temp(n).name ' unreadable\n']);
    end
    if numel(LargeObjNames)>=nLargeStim
        break
    end
end
%get the names in a cell array and shuffle
LargeObjNames           = {LargeObjNames.name}';
[LargeObjNames,index]   = Shuffle(LargeObjNames);
LargeObjMat             = LargeObjMat(index);

% Small
cd(fullfile(thePath.stim,'/Smaller'));
temp = dir('*.jpeg');
count = 1;
SmallObjMat     = cell(nSmallStim,1);
SmallObjNames   = cell(nSmallStim,1);
nFiles = size(temp,1);
randFileNums = datasample(1:nFiles,nFiles,'replace',false);
for n = 1:randFileNums
    try
        x=imread(temp(n).name);  %check if readable
        SmallObjMat(count) = imresize(x(:,:,1),2); % only resizing case for xdiva
        SmallObjNames(count) = temp(n);
        count = count + 1;
    catch
        fprintf(['\n' temp(n).name ' unreadable\n']);
    end
    if numel(SmallObjNames)>=nSmallStim
        break
    end
end
%get the names in a cell array and shuffle
SmallObjNames = {SmallObjNames.name}';
[SmallObjNames,index] = Shuffle(SmallObjNames);
SmallObjMat = SmallObjMat(index);

StimObj = containers.Map( [SmallObjNames; LargeObjNames], [SmallObjMat; LargeObjMat]);

%% Encoding
% Equal numbers of small and large obj stims: #N encoding trials/2 per type
%
% EncStimType -> 1 for small, 2 for large
% EncStimNames -> name of the stimuli as it apperas on the database

% order of scenes and faces per block
EncStimTypeIDs  = {'Small','Large'};
EncStimType     = zeros(nEncTrials,1); 
EncStimNames    = cell(nEncTrials,1);

% Get Encoding Stimuli
temp        = Shuffle(SmallObjNames);
EncSmallObjNames = temp(1:nEncSmallStim);

temp        = Shuffle(LargeObjNames);
EncLargeObjNames = temp(1:nEncLargeStim);
EncStimUniqueIDs = zeros(nEncTrials,2); EncStimUniqueIDs(:,1)=1:nEncTrials;

% assign small/large trials
TrialID         = 1:nEncTrials;
SmallEncTrials  = datasample(TrialID,nEncSmallStim,'replace',false);
LargeEncTrials  = setdiff(TrialID,SmallEncTrials);
EncStimType(SmallEncTrials) = 1;
EncStimType(LargeEncTrials) = 2;

% assigning images to trials
EncStimNames(SmallEncTrials)    = EncSmallObjNames;
EncStimNames(LargeEncTrials)    = EncLargeObjNames;

%% assign encoding conditions
%  max 10 condtions for a 2 x 5 
% condition codes:
%  1:5 5 phases for stim type 1;
%  6:10 5 phases for stim type 2;
nEncCondTrials      = nEncTrials/nEncConds;
nEncCondsByStimType = nEncConds/2; % 5 for small; 5 for large
EncCondTrialCode    = zeros(nEncTrials,1); % vector of condition code for each trial

EncCondCodeIDs      = cell(nEncConds,1);
EncCondCodeIDs(smallObjConds) = strcat('Small',cellfun(@num2str,num2cell(EncPhases'),'UniformOutput',false));
EncCondCodeIDs(largeObjConds) = strcat('Large',cellfun(@num2str,num2cell(EncPhases'),'UniformOutput',false));

cnt = 1;
for jj = 1:2 % small/large
    availableTrials = find(EncStimType==jj);
    for ii = 1:nEncCondsByStimType
        condiTrials=datasample(availableTrials,nEncCondTrials,'replace',false);
        EncCondTrialCode(condiTrials) = cnt;
        availableTrials = setdiff(availableTrials,condiTrials);
        cnt = cnt+1;
    end
end

% check for equal number of conditions throughout
assert(sum(histc(EncCondTrialCode,1:nEncConds)==nEncCondTrials)==nEncConds,'unequal number of trials produced')

%% assign retrieval conditions
% Condtions
% 1 -> old small
% 2 -> old large
% 3 -> new small
% 4 -> new large

RetCondIDs          = {'OldSmall','OldLarge','NewSmal','NewLarge'};
nRetCondTrials      = [nEncSmallStim,nEncLargeStim,nRetNewSmallStim,nRetNewLargeStim];
RetCondTrialCode    = zeros(nRetTrials,1);    
% assign condition such that each is equally likely per block.
% also retry the sampling such that there are not more than maxNumConsecutiveOld
% i.e., no that many consecutive old trials.
counter = 1;
while true    
    availableTrials = 1:nRetTrials;
    for cc = 1:nRetConds
        condTrials  = datasample(availableTrials,nRetCondTrials(cc),'replace',false);
        RetCondTrialCode(condTrials)=cc;
        availableTrials  = setdiff(availableTrials,condTrials);
    end
    
    % quick way to the number of consecutive old trials
    temp=conv(ones(maxNumConsecutiveOld,1),double((RetCondTrialCode<=2)));
    if sum(temp>=maxNumConsecutiveOld)==0
        break
    end
    counter=counter+1;
    if counter >10000
        error('could not arrive at a stable sequence in 10000 attempts')
    end
end
% make sure that each condition is equally represented
assert(sum(histc(RetCondTrialCode,1:nRetConds)==nRetCondTrials)==nRetConds,'unequal number of trials produced')

%% assign stimuli based on retrieval condition
RetStimNames = cell(nRetTrials,1);
temp                = setdiff(SmallObjNames,EncSmallObjNames);
NewSmallObjNames    = Shuffle(temp(1:nRetCondTrials(3)));
temp         = setdiff(LargeObjNames,EncLargeObjNames);
NewLargeObjNames  = Shuffle(temp(1:nRetCondTrials(4)));

RetStimNames(RetCondTrialCode==1) = Shuffle(EncSmallObjNames); % old small
RetStimNames(RetCondTrialCode==2) = Shuffle(EncLargeObjNames); % old large
RetStimNames(RetCondTrialCode==3) = Shuffle(NewSmallObjNames); % new small
RetStimNames(RetCondTrialCode==4) = Shuffle(NewLargeObjNames); % new large

%% store necessary variables
tacs_er = [];
tacs_er.subjNum     = thePath.subjNum;
tacs_er.exptType    = thePath.exptType;
tacs_er.thePath     = thePath;

% parameterts for making stimuli list
tacs_er.nEncStim    = nEncStim;     % # of encoding stimuli
tacs_er.nEncTrials  = nEncTrials;   % # of encoding trials (same as nEncStim if # blocks=1)
tacs_er.nEncConds   = nEncConds;    % # of conditions
tacs_er.nEncPhases  = nEncPhases;   % # of tacs stimulation phases
tacs_er.EncPhases   = EncPhases;    % tacs stimulation phases

tacs_er.nRetTrials  = nRetTrials;   % # of retrieval trials
tacs_er.nFoilTrials = nFoilTrials;  % # of novel items at retrieval
tacs_er.nRetConds   = nRetConds;    % # of retrieval conditions (old/new * number of stimulus type)
tacs_er.maxNumConsecutiveOld = maxNumConsecutiveOld; % parameter that constrains the maximum # of old items in a row

tacs_er.RandStream = s; % random seed 

% encoding
tacs_er.EncStimType     = EncStimType;      % categorical vector of stimulus types (length = # of trials)
tacs_er.EncStimTypeIDs  = EncStimTypeIDs;   % id of the stim types
tacs_er.EncCondTrialCode= EncCondTrialCode; % categorical vector of condition codes  (length = # of trials)
tacs_er.EncCondCodeIDs  = EncCondCodeIDs;   % id of the condition codes
tacs_er.EncStimNames    = EncStimNames;     % cell vector with the names of the image for each trial
tacs_er.EncStimUniqueIDs=EncStimUniqueIDs;  % unique ID for each image
x = tacs_er.EncCondTrialCode;
condNums = x;
condNums(x>nEncPhases)=x(x>nEncPhases)-nEncPhases;
tacs_er.EncTrialPhaseConds  = condNums;     % categorical vector of tacs phase (length = # of trials) 
tacs_er.EncTrialPhase   = tacs_er.EncPhases(condNums); % actual phase for each trial

% retrieval
tacs_er.RetCondIDs      = RetCondIDs;       % Names of retrieval conditions
tacs_er.RetCondTrialCode= RetCondTrialCode; % codes of retrieval conditions by trial
tacs_er.RetStimNames    = RetStimNames;     % stimulus image names for each trial
tacs_er.nRetCondTrials  = nRetCondTrials;   % # of trials by condition

% sanity check: verify that old labels match stimuli
assert(sum(ismember(tacs_er.RetStimNames,tacs_er.EncStimNames) == ...
    (tacs_er.RetCondTrialCode<=2))==tacs_er.nRetTrials,'Labels do not match stimuli')

% encodign conditions at retrieval
[~,i]=ismember(tacs_er.RetStimNames,tacs_er.EncStimNames);
tacs_er.EncCondAtRet = zeros(tacs_er.nRetTrials,1);
tacs_er.EncCondAtRet(i>0) = tacs_er.EncCondTrialCode(i(i>0));

% Stimulus Object
tacs_er.Stimuli = StimObj;
cd(thePath.main)

%% Stimulus Presentation Parameters
PresParams                      = [];
PresParams.VideoFrameRate       = 60;
PresParams.nXDivaFrames         = 90;
PresParams.nXDivaPreludeFrames  = 10; % leave first 10 spots blank
PresParams.stimFrequency        = 6;
PresParams.stimDurationInCycles  = 0.5;
PresParams.totalTrialDuration   = 60; % in seconds;
PresParams.stimDurationInSecs   = 1/PresParams.stimFrequency*PresParams.stimDurationInCycles;
PresParams.nCycles              = PresParams.stimFrequency*PresParams.totalTrialDuration;
PresParams.stimPresCycles       = ones(PresParams.nCycles,1); % stimuli will be presented every cycle
PresParams.DegreesPerFrame      = 360/(PresParams.VideoFrameRate/PresParams.stimFrequency);
PresParams.VideoFrameDurSecs    = 1/PresParams.VideoFrameRate;
PresParams.stimDurationInFrames = PresParams.VideoFrameRate*PresParams.stimDurationInCycles/PresParams.stimFrequency;

%% Save individual images for xDiva presentation and prepare sequence
PresParams.FixCrossFrameIDs           = 1; % on the first frame
PresParams.FixCrossMinNFrames         = 20;% corresponding to 333ms or 2 cycles of 6Hz
PresParams.InstructionsLength         = 2*60*PresParams.VideoFrameRate;
PresParams.InstructionSet             = mod(thePath.subjectNum,2)==0 +1;
tacs_er.EncFaceRespID                 = mod(thePath.subjectNum,2)==0 +1;
tacs_er.EncSceneRespID                = mod(thePath.subjectNum,2)==1 +1;

PresParams.tACSPreTaskWarmUpDur       = 3*60*PresParams.VideoFrameRate;     
tacs_er.PresParams = PresParams;
% Save into subjects path
fileName = 'tacs_er_xdiva.task.mat';
if ~exist([thePath.subjectPath '/' fileName],'file')
    save([thePath.subjectPath '/' fileName],'tacs_er')
    mkdir([thePath.subjectPath '/xdiva/'])
    overwriteFlag=1;
else
    warning('tacs task for this subject already created')
    while 1
        s = input('Overwrite? (Y/N)','s');
        if strcmp(s,'Y')
            overwriteFlag=1;
            save([thePath.subjectPath '/' fileName],'tacs_er')
            break
        elseif strcmp(s,'N')
            break
        end
    end
end

if overwriteFlag
    ZeroPhaseImgFrameIDs       = (PresParams.FixCrossMinNFrames+1):...
        PresParams.stimDurationInFrames*2:PresParams.nXDivaFrames;
    ZeroPhaseImgFrameIDs(PresParams.stimFrequency+1:end)=[];
    
    ZeroPhaseBlankFrameIDs       = ZeroPhaseImgFrameIDs+PresParams.stimDurationInFrames;
    
    images = zeros([stimSize,1,3], 'uint8');
    images(:,:,1,2) = 128*ones(stimSize,'uint8');
    images(:,:,1,3) = fixCrossImage(stimSize(1),round(stimSize(1)/20));
    
    
    % Save individual stims into subjects path
    fileName = 'tacs_enc_xdiva_';
    for tt = 1:nEncTrials
        
        images(:,:,1,1) = StimObj(EncStimNames{tt});
        
        condNum = tacs_er.EncTrialPhaseConds(tt);
        % n phases go from 1 to nEncPhases, in jumps of 2 frames=72deg for 6Hz
        % stimulation
        imageSequence = zeros(PresParams.nXDivaFrames,1,'uint32');
        imageSequence(1)=3;
        imageSequence(ZeroPhaseImgFrameIDs+(condNum-1)*2)=1;
        imageSequence(ZeroPhaseBlankFrameIDs+(condNum-1)*2)=2;
        
        save([thePath.subjectPath '/xdiva/' fileName num2str(tt) '.mat'],'images','imageSequence')
    end
    
    % save instructions (changes by subj number)
    if PresParams.InstructionSet==1
        img  = rgb2gray(imread([thePath.stim '/EncInstructions_xDiva1.png']));
    else
        img  = rgb2gray(imread([thePath.stim '/EncInstructions_xDiva2.png']));
    end
    images = zeros([size(img),1,1], 'uint8');
    images(:,:,1,1) = img;        
    
    % 2 mins worth of instructions.
    imageSequence = zeros(PresParams.InstructionsLength,1,'uint32');
    imageSequence(1)=1;    
    save([thePath.subjectPath '/xdiva/Instructions_1.mat'],'images','imageSequence')
    
    % save pre-task warm up sequence    
    images = zeros([200,200,1,1], 'uint8');
    images(:,:,1,1) = 128*ones([200 200],'uint8');
    
    imageSequence = zeros(PresParams.tACSPreTaskWarmUpDur,1,'uint32');
    imageSequence(1)=1;    
    save([thePath.subjectPath '/xdiva/tACSPreTask_1.mat'],'images','imageSequence')
    
end

end

function im=fixCrossImage(imgSize,fixCrossSize)
% all in pixels
% only supports square imgages (for now).
center = round(imgSize(1)/2);
horizontalMarkIDs = (center-fixCrossSize):(center+fixCrossSize);
verticalMarkIDs   = horizontalMarkIDs;

im = 128*ones(imgSize(1),imgSize(1),'uint8');
im((-1:1)+center,horizontalMarkIDs) = 255;
im(verticalMarkIDs,(-1:1)+center)   = 255;

return

end

