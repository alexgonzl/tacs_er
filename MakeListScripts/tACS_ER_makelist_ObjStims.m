function [tacs_er] = tACS_ER_makelist_ObjStims(thePath)
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
% Created:      Sept 13, 2016
% LastUpdate:   Sept 13, 2016
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
temp = dir('*.jpg');
count = 1;
LargeObjMat     = cell(nLargeStim,1);
LargeObjNames   = cell(nLargeStim,1);
nFiles = size(temp,1);
randFileNums = datasample(1:nFiles,nFiles,'replace',false);
for n = randFileNums
    try
        x = imread(temp(n).name); %check if readable
        LargeObjMat{count} = x; % only resizing case for xdiva
        LargeObjNames{count} = temp(n).name;
        count = count + 1;
    catch
        fprintf(['\n' temp(n).name ' unreadable\n']);
    end
    if count>nLargeStim
        break
    end
end
%get the names in a cell array and shuffle
[LargeObjNames,index]   = Shuffle(LargeObjNames);
LargeObjMat             = LargeObjMat(index);

% Small
cd(fullfile(thePath.stim,'/Smaller'));
temp = dir('*.jpg');
count = 1;
SmallObjMat     = cell(nSmallStim,1);
SmallObjNames   = cell(nSmallStim,1);
nFiles = size(temp,1);
randFileNums = datasample(1:nFiles,nFiles,'replace',false);
for n = randFileNums
    try
        x=imread(temp(n).name);  %check if readable
        SmallObjMat{count} = x; % only resizing case for xdiva
        SmallObjNames{count} = temp(n).name;
        count = count + 1;
    catch
        fprintf(['\n' temp(n).name ' unreadable\n']);
    end
    if count>nSmallStim
        break
    end
end
%get the names in a cell array and shuffle
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
smallObjConds       = 1:nEncCondsByStimType;
largeObjConds       = (nEncCondsByStimType+1):nEncConds;
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
nRetCondTrials      = [nEncSmallStim,nEncLargeStim,nRetNewSmallStim,nRetNewLargeStim]';
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
    x = RetCondTrialCode<=2;
    temp=diff([0 find( ~ (x > 0))' numel(x)+1])-1;
    if max(temp)<=maxNumConsecutiveOld
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
PresParams.stimFrequency        = 6;
PresParams.stimDurationInCycles  = 0.5; % half a theta cycle
PresParams.totalTrialDuration   = 1.5; % in seconds;
PresParams.stimDurationInSecs   = 1/PresParams.stimFrequency*PresParams.stimDurationInCycles;
PresParams.DegreesPerFrame      = 360/(PresParams.VideoFrameRate/PresParams.stimFrequency);
PresParams.VideoFrameDurSecs    = 1/PresParams.VideoFrameRate;
PresParams.stimDurationInFrames = PresParams.VideoFrameRate*PresParams.stimDurationInCycles/PresParams.stimFrequency;

%% Save individual images for xDiva presentation and prepare sequence
PresParams.FixCrossFrameIDs           = 1; % on the first frame
PresParams.FixCrossMinNFrames         = 20;% corresponding to 333ms or 2 cycles of 6Hz
PresParams.InstructionsLength         = 2*60*PresParams.VideoFrameRate;
PresParams.InstructionSet             = mod(thePath.subjNum,2)==0 +1;
tacs_er.EncFaceRespID                 = mod(thePath.subjNum,2)==0 +1;
tacs_er.EncSceneRespID                = mod(thePath.subjNum,2)==1 +1;

PresParams.tACSPreTaskWarmUpDur       = 3*60*PresParams.VideoFrameRate;
tacs_er.PresParams = PresParams;

% Save into subjects path
fileName = 'tacs_er_objstim.task.mat';
if ~exist([thePath.subjectPath '/' fileName],'file')
    save([thePath.subjectPath '/' fileName],'tacs_er')
    mkdir([thePath.subjectPath '/xdiva/'])
else
    warning('tacs task for this subject already created')
    while 1
        s = input('Overwrite? (Y/N)','s');
        if strcmp(s,'Y')
            save([thePath.subjectPath '/' fileName],'tacs_er')
            break
        elseif strcmp(s,'N')
            break
        end
    end
end


end

function [Y,index] = Shuffle(X)

[~,index] = sort(rand(size(X)));
[n,m] = size(X);
Y = zeros(size(X));
if n == 1 || m == 1
    Y = X(index);
else
    for j = 1:m
        Y(:,j)  = X(index(:,j),j);
    end
end
end
