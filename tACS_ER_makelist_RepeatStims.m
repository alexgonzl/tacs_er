
function [tacs_er] = tACS_ER_makelist_RepeatStims(subjNum, thePath)

% make lists for tacs encoding & retrieval (tacs_er) experiment
% This design:
% 1) all stimuli is presented twice, once per encoding block
% 2) responses are made to cue. counter balance cue responses per stimuli
% -> tom cruise: block 1. yes, block 2 no.

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
% EncStimCue        -> VECTOR trial code for cue ID
% EncBlockID        -> VECTOR trial block ID
%
% RetCondIDs        -> {'OldFace','OldScene','NewFace','NewScene'};
% RetCondTrialCode  -> VECTOR trial 1:4, indicating RetCondID
% RetStimNames      -> VECTOR trial of retrieaval IDs
% RetBlockID        -> VECTOR of trial block ID
%
% Stimuli           -> key, matrix map of stimuli's names to image matrix.

%% Set-up 

% parameter list
nUniqueFaceStimuli = 120;
nUniqueSceneStimuli = 120;
nEncStim    = nUniqueSceneStimuli+ nUniqueFaceStimuli;
nEncBlocks = 2;
nEncTrials = nEncBlocks*nEncStim;
nEncConds  = 6;     % faces, scenes x phase (3 different phases)
nCueTypes  = 2;
stimSize   = [225 225];

nRetTrials = 360;   % encoding trials + retrieval trials
nFoilTrials =nRetTrials-nEncStim;
nRetConds  = 4;     % old and new conditons per stim type (Face/scene)
nRetBlocks = nEncBlocks;    % blocks
maxNumConsecutiveOld = 8;   % maximum number of old trials in a row

% reset stream (to avoid duplicate lists)
s = RandStream.create('mt19937ar','seed',sum(100*clock));
if strfind(version,'R2014b')>0
    RandStream.setGlobalStream(s)
else
    RandStream.setDefaultStream(s);
end

%% load data names
cd(fullfile(thePath.stim,'landmarks'));
temp = dir('*.jpg');
count = 1;
ScenesMat = cell(numel(temp),1);
for n = 1:size(temp,1)
    ScenesMat{n} = zeros(stimSize(1),stimSize(2),'uint8');
    try        
        x=imread(temp(n).name); %check if readable        
        ScenesMat{n} = x(:,:,1);
        SceneNames(count) = temp(n);
        count = count + 1;
    catch
        fprintf(['\n' temp(n).name ' unreadable\n']);
    end
end

%get the names in a cell array and shuffle
SceneNames = {SceneNames.name}';
[SceneNames,index] = Shuffle(SceneNames);
ScenesMat = ScenesMat(index);

cd(fullfile(thePath.stim,'people'));
temp = dir('*.jpg');
count = 1;
FacesMat = cell(numel(temp),1);
for n = 1:size(temp,1)
    try
        x=imread(temp(n).name);  %check if readable
        FacesMat{n} = x(:,:,1);
        FaceNames(count) = temp(n);
        count = count + 1;
    catch
        fprintf(['\n' temp(n).name ' unreadable\n']);
    end
end
%get the names in a cell array and shuffle
FaceNames = {FaceNames.name}';
[FaceNames,index] = Shuffle(FaceNames);
FacesMat = FacesMat(index);

StimObj = containers.Map( [FaceNames; SceneNames], [FacesMat; ScenesMat]);
%% Encoding
% Equal numbers of faces and scenes: #N encoding trials/2 per type
% Each block will have one of each of stimuli, for a total of # of blocks x
% # of unique stimuli.
%
% Association Cue Form:
% Three options for associated cue (say red/blue)
% 1) Associated cue alternates by block
% 2) Associated cue 50% chance of one or the other
% 3) Associated cue stays the same per stimuli.
%
% EncStimType -> 1 for faces, 2 for scenes
% EncStimCue  -> 1 or 2
% EncStimNames -> name of the stimuli as it apperas on the database
% EncBlockID   -> block ID

% order of scenes and faces per block
EncStimTypeIDs = {'Face','Scene'};
EncStimType = zeros(nEncTrials,1);
EncBlockID  = zeros(nEncTrials,1);
EncStimCue  = zeros(nEncTrials,1);
EncStimNames = cell(nEncTrials,1);
nEncTrialsBlock = nEncTrials/nEncBlocks;
nFacesEncBlock  = nEncTrialsBlock/2;
nScenesEncBlock = nFacesEncBlock;

% Get Encoding Stimuli
temp = Shuffle(FaceNames);
EncFaces = temp(1:nFacesEncBlock);
temp = Shuffle(SceneNames);
EncScenes = temp(1:nScenesEncBlock);

for bb = 1:nEncBlocks
    TrialBlockID = ((bb-1)*nEncTrialsBlock+1):bb*nEncTrialsBlock;
    FaceTrials  = datasample(TrialBlockID,nFacesEncBlock,'replace',false);
    SceneTrials = setdiff(TrialBlockID,FaceTrials);
    EncBlockID(TrialBlockID) = bb;
    EncStimType(FaceTrials)  = 1;
    EncStimType(SceneTrials) = 2;
    
    % assign individial stimuli to trials
    EncStimNames(FaceTrials) = Shuffle(EncFaces);
    EncStimNames(SceneTrials) = Shuffle(EncScenes);
    
    % associate a cue
    if bb ==1
        % associated cue alternates by block; first one is random
        cues = Shuffle([ones(nFacesEncBlock/2,1); 2*ones(nFacesEncBlock/2,1)]);
        EncStimCue(FaceTrials) = Shuffle(cues);
        EncStimCue(SceneTrials) = Shuffle(cues);
        blockCues = EncStimCue(TrialBlockID);
    else
        % index of new order of stimuli
        [~,i]=ismember(EncStimNames(EncBlockID==bb),EncStimNames(EncBlockID==1));
        
        % alternate associated cue
        if mod(bb,2)==0
            EncStimCue(TrialBlockID)=mod(blockCues(i),2)+1;
        else
            EncStimCue(TrialBlockID)=blockCues(i);
        end
    end
end

% assign encoding conditions
% six conditions 2 x 3 design
% 1) faces at 0deg
% 2) faces at 90deg
% 3) faces at 180deg
% 4) scenes at 0deg
% 5) scenes at 90deg
% 6) scenes at 180deg
%
% ** Note, not counterbalancing by stimuli at the moment, that would have 
% to be done across subjects. **

nEncCondTrials      = nEncTrials/nEncConds;
nEncCondBlockTrials = nEncCondTrials/nEncBlocks;
EncCondTrialCode    = zeros(nEncTrials,1);
EncCondCodeIDs      = {'Face0','Face90','Face180','Scn0','Scn90','Scn180'};
faceConds = 1:3; sceneConds = 4:6;

if mod(nEncCondBlockTrials,1)~=0
    error('uneven encoding condition trials')
end

for bb = 1:nEncBlocks
    TrialBlockID = ((bb-1)*nEncTrialsBlock+1):bb*nEncTrialsBlock;
    
    % for face trials
    nFaceConds=numel(faceConds);
    for jj=1:nCueTypes
        availableFaceTrials = intersect(find(EncStimType==1 & EncStimCue==jj),TrialBlockID);
        for ii = 1:nFaceConds
            condiTrials=datasample(availableFaceTrials,nEncCondBlockTrials/nCueTypes,'replace',false);
            EncCondTrialCode(condiTrials) = faceConds(ii);
            availableFaceTrials = setdiff(availableFaceTrials,condiTrials);
        end
        
        % for scenes trials
        nSceneConds=numel(sceneConds);
        availableSceneTrials = intersect(find(EncStimType==2 & EncStimCue==jj),TrialBlockID);
        for ii = 1:nSceneConds
            condiTrials=datasample(availableSceneTrials,nEncCondBlockTrials/nCueTypes,'replace',false);
            EncCondTrialCode(condiTrials) = sceneConds(ii);
            availableSceneTrials = setdiff(availableSceneTrials,condiTrials);
        end
    end
end

% check for equal number of conditions throughout
assert(sum(histc(EncCondTrialCode,1:nEncConds)==nEncCondTrials)==nEncConds,'unequal number of trials produced')

%% assign retrieval conditions
% Condtions
% 1 -> old face
% 2 -> old scene
% 3 -> new face
% 4 -> new scene

RetCondIDs        = {'OldFace','OldScene','NewFace','NewScene'};
nRetTrialsBlock   = nRetTrials/nRetBlocks;
nRetTrialTypesBlock = [nUniqueFaceStimuli nUniqueSceneStimuli nFoilTrials/2 nFoilTrials/2]'/nRetBlocks;
nRetCondTrials      = nRetTrialTypesBlock*nRetBlocks;

% assign condition such that each is equally likely per block.
% also retry the sampling such that there are not more than maxNumConsecutiveOld
% i.e., no that many consecutive old trials.
counter = 1;
while true
    RetCondTrialCode = zeros(nRetTrials,1);
    RetBlockID = zeros(nRetTrials,1);
    for rr = 1:nRetBlocks        
        TrialBlockID = ((rr-1)*nRetTrialsBlock+1):rr*nRetTrialsBlock;
        availableTrials = TrialBlockID;
        for cc = 1:nRetConds
            condTrials  = datasample(availableTrials,nRetTrialTypesBlock(cc),'replace',false);
            RetCondTrialCode(condTrials)=cc;
            availableTrials  = setdiff(availableTrials,condTrials);
        end
        RetBlockID(TrialBlockID)=rr;
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
temp         = setdiff(FaceNames,EncFaces);
FoilFaces    = Shuffle(temp(1:nRetCondTrials(3)));
temp         = setdiff(SceneNames,EncScenes);
FoilScenes   = Shuffle(temp(1:nRetCondTrials(4)));

RetStimNames(RetCondTrialCode==1) = Shuffle(EncFaces); % old faces
RetStimNames(RetCondTrialCode==2) = Shuffle(EncScenes); % old scenes
RetStimNames(RetCondTrialCode==3) = Shuffle(FoilFaces); % new faces
RetStimNames(RetCondTrialCode==4) = Shuffle(FoilScenes); % new scenes


%% store necessary variables
tacs_er = [];
tacs_er.subjNum = subjNum;
tacs_er.thePath = thePath;

% parameterts for making stimuli list
tacs_er.nEncStim = nEncStim;
tacs_er.nEncTrials = nEncTrials;
tacs_er.nEncConds  = nEncConds;
tacs_er.nEncBlocks = nEncBlocks;

tacs_er.nRetTrials = nRetTrials;
tacs_er.nFoilTrials = nFoilTrials;
tacs_er.nRetConds  = nRetConds;
tacs_er.nRetBlocks = nRetBlocks;
tacs_er.maxNumConsecutiveOld = maxNumConsecutiveOld;

tacs_er.RandStream = s;

% encoding
tacs_er.EncStimType     = EncStimType;
tacs_er.EncStimTypeIDs  = EncStimTypeIDs;
tacs_er.EncCondTrialCode= EncCondTrialCode;
tacs_er.EncCondCodeIDs  = EncCondCodeIDs;
tacs_er.EncStimNames    = EncStimNames;
tacs_er.EncStimCue      = EncStimCue;
tacs_er.EncBlockID      = EncBlockID;

% retrieval
tacs_er.RetCondIDs      = RetCondIDs;
tacs_er.RetCondTrialCode= RetCondTrialCode;
tacs_er.RetStimNames    = RetStimNames;
tacs_er.RetBlockID      = RetBlockID;
tacs_er.nRetCondTrials  = nRetCondTrials;


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

end

