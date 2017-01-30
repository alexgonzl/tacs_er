function [tacs_er] = tACS_ER_makelist_Expt3b(thePath)
% make lists for tacs encoding & retrieval (tacs_er) experiment
% This design:
% 1) One encoding block with flashing scene stimuli.
% 2) Responses will be made to the stimulus itself, no cue conditions.
% 3) The decision will be man-made vs natural images
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
% EncStimType       -> VECTOR trial of indicating type (1  man-made, 2 natural)
% EncStimTypeIDs    -> {Mm,Na}; man made or natural
% EncCondCodeIDs    -> {'Mm000','Mm036','Mm072','Mm108'...,'Na000',...}; 
% EncCondTrialCode  -> VECTOR trial encoding indicating EncCondCodeID
% EncStimNames      -> VECTOR trial stimulus ID
% EncStimUniqueIDs  -> VECTOR trial with unique ID for each stim
% EncBlockID        -> VECTOR trial block ID
%
% RetCondIDs        -> {'OldMm','OldMm','NewNa','NewNa'};
% RetCondTrialCode  -> VECTOR trial 1:4, indicating RetCondID
% RetStimNames      -> VECTOR trial of retrieaval IDs
% RetBlockID        -> VECTOR of trial block ID
%
% Stimuli           -> key, matrix map of stimuli's names to image matrix.

%------------------------------------------------------------------------%
% Author:       Alex Gonzalez
% Created:      Jan 30, 2017
% LastUpdate:   Jan 30, 2017
%------------------------------------------------------------------------%

%% Set-up

% parameter list
nTotalStim          = 600;
nEncStim1           = 200;
nEncStim2           = 200;
nRetNewStim1        = 100;
nRetNewStim2        = 100;
stimSize            = [256 256];
assert((nEncStim1+nEncStim2+nRetNewStim1+nRetNewStim2==nTotalStim),'incorrect numbers of stims')
nStim1 = nEncStim1 + nRetNewStim1;
nStim2 = nEncStim2 + nRetNewStim2;

nEncStim    = nEncStim1 + nEncStim2;
nEncPhases  = 10;                 % constrained by xdiva design.
nEncConds   = 2*nEncPhases;       % stim1 vs stim2 x phase (10 different phases)
EncPhases   = 0:(360/nEncPhases):359;

nEncTrials  = nEncStim;
nRetTrials  = nTotalStim;   % encoding trials + retrieval trials
nFoilTrials = nRetTrials-nEncStim;
nRetConds   = 4;            % old and new conditons per stim type (small/large)
maxNumConsecutiveOld = 8;   % maximum number of old trials in a row

% Set Random Seed based on SubjectNum
s = RandStream.create('mt19937ar','seed',thePath.subjNum);
RandStream.setGlobalStream(s)
if strcmp(thePath.exptType,'behav_v14')
    save_xdiva_flag = 0;
else
    error('')
end
%% load stimulii

% stim 1 category: Indoor
load([thePath.stim,'/scene_categories/imageRanksv2.mat']);
load([thePath.stim,'/scene_categories/selInOutImgs_v14.mat']);

Stim1Mat     = cell(nStim1,1);
Stim1Names   = cell(nStim1,1);
for n = 1:nStim1
    Stim1Mat{n}     = squeeze(ImgMat(1,n,:,:));
    Stim1Names{n}   = selImgNames{1,n};
end
rand1FileNums = datasample(1:nStim1,nStim1,'replace',false);
Stim1Mat    = Stim1Mat(rand1FileNums);
Stim1Names  = Stim1Names(rand1FileNums);

Stim2Mat     = cell(nStim2,1);
Stim2Names   = cell(nStim2,1);
for n = 1:nStim2
    Stim2Mat{n}     = squeeze(ImgMat(2,n,:,:));
    Stim2Names{n}   = selImgNames{2,n};
end
rand2FileNums = datasample(1:nStim2,nStim2,'replace',false);
Stim2Mat    = Stim2Mat(rand2FileNums);
Stim2Names  = Stim2Names(rand2FileNums);

StimObj = containers.Map( [Stim1Names ; Stim2Names], [Stim1Mat; Stim2Mat]);

%% Encoding
% Equal numbers of stim1 and stim2: #N encoding trials/2 per type
%
% EncStimType -> 1 for man-made, 2 for natural 
% EncStimNames -> name of the stimuli as it apperas on the database

% order of man-made and natural per block
EncStimTypeIDs  = {'Mm','Na'};
EncStimType     = zeros(nEncTrials,1); 
EncStimNames    = cell(nEncTrials,1);

% Get Encoding man-made
temp        = Shuffle(Stim1Names);
EncStim1Names = temp(1:nEncStim1);

temp        = Shuffle(Stim2Names);
EncStim2Names = temp(1:nEncStim2);
EncStimUniqueIDs = zeros(nEncTrials,2); EncStimUniqueIDs(:,1)=1:nEncTrials;

% assign mm/na trials at encoding
TrialID         = 1:nEncTrials;
EncStim1Trials  = datasample(TrialID,nEncStim1,'replace',false);
EncStim2Trials  = setdiff(TrialID,EncStim1Trials);
EncStimType(EncStim1Trials) = 1;
EncStimType(EncStim2Trials) = 2;

% assigning images to trials
EncStimNames(EncStim1Trials)    = EncStim1Names;
EncStimNames(EncStim2Trials)    = EncStim2Names;

%% assign encoding conditions
%  max 2- condtions for a 2 x 10
% condition codes:
%  1:10 10 phases for stim type 1;
%  11:20 10 phases for stim type 2;
nEncCondTrials      = nEncTrials/nEncConds;
nEncCondsByStimType = nEncConds/2; % 10 for stim1; 10 for stim2
EncCondTrialCode    = zeros(nEncTrials,1); % vector of condition code for each trial
stim1ObjConds       = 1:nEncCondsByStimType;
stim2ObjConds       = (nEncCondsByStimType+1):nEncConds;
EncCondCodeIDs      = cell(nEncConds,1);
EncCondCodeIDs(stim1ObjConds) = strcat('mm',cellfun(@num2str,num2cell(EncPhases'),'UniformOutput',false));
EncCondCodeIDs(stim2ObjConds) = strcat('na',cellfun(@num2str,num2cell(EncPhases'),'UniformOutput',false));

cnt = 1;
for jj = 1:2 % mm/na
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
% 1 -> old stim1
% 2 -> old stim2
% 3 -> new stim1
% 4 -> new stim2

RetCondIDs          = {'OldStim1','OldStim2','NewStim1','NewStim2'};
nRetCondTrials      = [nEncStim1,nEncStim1,nRetNewStim1,nRetNewStim2]';
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
RetStimNames        = cell(nRetTrials,1);
temp                = setdiff(Stim1Names,EncStim1Names);
NewStim1Names       = Shuffle(temp(1:nRetCondTrials(3)));
temp                = setdiff(Stim2Names,EncStim2Names);
NewStim2Names       = Shuffle(temp(1:nRetCondTrials(4)));

RetStimNames(RetCondTrialCode==1) = Shuffle(EncStim1Names); % old small
RetStimNames(RetCondTrialCode==2) = Shuffle(EncStim2Names); % old large
RetStimNames(RetCondTrialCode==3) = Shuffle(NewStim1Names); % new small
RetStimNames(RetCondTrialCode==4) = Shuffle(NewStim2Names); % new large

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
PresParams.InstructionSet             = mod(thePath.subjNum,4)==0 +1;
tacs_er.EncStim1RespID                = mod(thePath.subjNum,2)==0 +1;
tacs_er.EncStim2RespID                = mod(thePath.subjNum,2)==1 +1;

PresParams.tACSPreTaskWarmUpDur       = 3*60*PresParams.VideoFrameRate;     
tacs_er.PresParams = PresParams;
% Save into subjects path
fileName = 'tacs_er3_xdiva.task.mat';
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

if save_xdiva_flag
if overwriteFlag
    ZeroPhaseImgFrameIDs        = PresParams.FixCrossMinNFrames+1;    
    ZeroPhaseBlankFrameIDs      = ZeroPhaseImgFrameIDs+PresParams.stimDurationInFrames;
    
   
    % Save individual stims into subjects path
    fileName = 'tacs_er3_xdiva_';
    
    % pseudo trials before the task //** this is an xdiva fix for trouble
    % going from condition to condition in run mode.
    images = zeros([10 10,1,3], 'uint8');
    images(:,:,1,1) = 128*ones([10 10],'uint8');
    imageSequence = zeros(PresParams.nXDivaFrames,1,'uint32');
    imageSequence(1)=1;    
    nPreTaskPseudoTrials = 50;
    for ii = 1:nPreTaskPseudoTrials
        save([thePath.subjectPath '/xdiva/' fileName num2str(ii) '.mat'],'images','imageSequence')
    end
    images = zeros([stimSize,1,3], 'uint8');
    images(:,:,1,2) = 128*ones(stimSize,'uint8');
    images(:,:,1,3) = fixCrossImage(stimSize(1),round(stimSize(1)/20));
    
    for tt = 1:nEncTrials        
        images(:,:,1,1) = StimObj(EncStimNames{tt});        
        condNum = tacs_er.EncTrialPhaseConds(tt);        
        % n phases go from 1 to nEncPhases, in jumps of 1 frame=36deg for 6Hz
        % stimulation
        imageSequence = zeros(PresParams.nXDivaFrames,1,'uint32');
        imageSequence(1)=3;
        imageSequence(ZeroPhaseImgFrameIDs+(condNum-1))=1;
        imageSequence(ZeroPhaseBlankFrameIDs+(condNum-1))=2;        
        save([thePath.subjectPath '/xdiva/' fileName num2str(tt+nPreTaskPseudoTrials) '.mat'],'images','imageSequence')
    end
    
    % save instructions (changes by subj number)
    if PresParams.InstructionSet==1
        img  = rgb2gray(imread([thePath.stim '/tacs_er3_ins1.png']));
    else
        img  = rgb2gray(imread([thePath.stim '/tacs_er3_ins2.png']));
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

function [Y,index] = Shuffle(X)

[null,index] = sort(rand(size(X)));
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
 
