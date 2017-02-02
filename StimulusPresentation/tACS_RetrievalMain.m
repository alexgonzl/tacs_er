%function [ret_out,msg]=tACS_RetrievalMain(thePath)
% tACS experiment image recognition memory presentation script
% This script takes the 'tacs_er' structure that contains the stimuli and
% order that shall be used at retrieval. These images are then loaded and
% saved prior to presentation. The task is recognition memory.
% Subjects can respond to images as 'old', 'new', or 'unsure'. If the
% subjects respond old or new, they are then given a confidence scale to
% indicate their confidence in the decision.
%
%------------------------------------------------------------------------%
% Author:       Alex Gonzalez
% Created:      Aug 20th, 2015
% LastUpdate:   Jan 30, 2017
%------------------------------------------------------------------------%

% clear all the screens
close all;
sca;

% load the task
%fileName = strcat(thePath.subjectPath,'/', thePath.exptType,'.mat');
%fileName = strcat(thePath.subjectPath,'/','tacs_er_objstim.task.mat');
if any(strcmp(thePath.exptType,{'behav_v13','behav_v14','behav_v15'}))
    fileName = strcat(thePath.subjectPath,'/','tacs_er3_xdiva.task.mat');
end
if exist(fileName,'file')
    load(fileName);
else
    error('no task created, must run encoding task first!')
end

%PsychDebugWindowConfiguration;

% Presentation Parameters
PresParams  = [];
PresParams.stimDurationInSecs   = 3;
PresParams.ITI_Range            = [0.75 1.25]; % variable ITI in secs
PresParams.SaveEvNtrials        = 20; % save progress every X# of trials.
PresParams.lineWidthPix         = 5;       % Set the line width for our fixation cross
PresParams.dotColor             = [1 1 1];
PresParams.fixCrossColor        = [1 1 1];
PresParams.textColor            = [1 1 1];
PresParams.ConfidenceScale      = 0; % confidence scale flag
PresParams.UnsureButtonOpt      = 0;
PresParams.ConfBarColor         = [0.2 0.1385 1];
PresParams.SelfPaceFlag         = 0;

PresParams.preStartTime         = 1;
PresParams.MaxResponseTime      = 3;       % maximum to make recognition decision
PresParams.MaxConfDecInSecs     = 3;       % max time to make confidence decision
PresParams.TotalTrialDur        = 3.5;       %

if  any(strcmp(thePath.exptType,{'tacs_enc','tacs_enc_xdiva','tacs_enc_xdiva_obj','tacs_er_objstim'}))
    PresParams.StarStimEEG      = 1;
else
    PresParams.StarStimEEG      = 0;
end

% determine numbers for recognition decision
% depending on subject number and active Keyboard.
laptopResponseKeys = ['j','k','l'];
keypadResponseKeys = ['4','5','6'];
ConfidenceKeys     = ['1','2','3'];
RespConds         = {'old','unsure','new'};

if mod(thePath.subjNum,2)
    responseMap = [1,2,3];
else
    responseMap = [3,2,1];
end
laptopResponseKeys = laptopResponseKeys(responseMap);
keypadResponseKeys = keypadResponseKeys(responseMap);

PresParams.RespConds = RespConds;

% get the IDs of the trials
stimNames = tacs_er.RetStimNames;
nTrials   = numel(stimNames);


% EEG recording parameters
if PresParams.StarStimEEG
    PresParams.StarStimIP           = '10.0.0.42';
    PresParams.preStartTime         = 20; % in seconds.
    
    addpath(genpath('../MatNIC_v2.5/'))
    try
        [ret, status, socket] = MatNICConnect(PresParams.StarStimIP);
        if ret<0
            error(status)
        end
        [ret,LSLOutlet]=MatNICMarkerConnectLSL('alexLSL');
        if ret<0
            error('could not connect to streaming layer')
        end
        
    catch msg
        MatNICMarkerCloseLSL(LSLOutlet);
        MatNICStopEEG(socket);
        sca
        keyboard;
    end
end
%%

% Initialize trial timing structure
TimingInfo = [];
TimingInfo.preStimFixFlip   = cell(nTrials,1);
TimingInfo.stimPresFlip     = cell(nTrials,1);
TimingInfo.postStimFlip     = cell(nTrials,1);
TimingInfo.trialRT          = nan(nTrials,1);
TimingInfo.trialKeyPress    = cell(nTrials,1);
TimingInfo.CondResp         = cell(nTrials,1);
TimingInfo.Confidence       = nan(nTrials,1);
TimingInfo.ConfResp         = cell(nTrials,1);

try
    
    %---------------------------------------------------------------------%
    % Screen and additional presentation parameters
    %---------------------------------------------------------------------%
    % Get keyboard number
    [activeKeyboardID, laptopKeyboardID, pauseKey, resumeKey] = getKeyboardOr10key;
    % initialize Keyboard Queue
    KbQueueCreate(activeKeyboardID);
    % Start keyboard queue
    KbQueueStart(activeKeyboardID);
    
    % get correct mapping to keyboard
    if laptopKeyboardID==activeKeyboardID
        PresParams.RespButtons  = laptopResponseKeys;
        RespButtons             = laptopResponseKeys;
    else
        PresParams.RespButtons = keypadResponseKeys;
        RespButtons            = keypadResponseKeys;
    end
    % initialie window
    [window, windowRect] = initializeScreen;
    screenXpixels = windowRect(3);
    screenYpixels = windowRect(4);
    
    % Get the centre coordinate of the window
    [xCenter, yCenter] = RectCenter(windowRect);
    
    % Get coordinates for fixation cross
    fixCrossCoords = fixCross(xCenter, yCenter,screenXpixels,screenYpixels);
    
    % Query the frame duration
    ifi = Screen('GetFlipInterval', window);
    
    % Get the fixation time for every trial (random ITI)
    % in seconds.
    ITIsSecs        = randi(round(PresParams.ITI_Range/ifi),nTrials,1)*ifi;
    
    % post-stim max confidence response period duration
    MaxConfDescFrames  = round(PresParams.MaxConfDecInSecs /ifi);
    
    % pre-make image textures
    imgTextures = cell(nTrials,1);
    for ii = 1:nTrials
        imgTextures{ii}=Screen('MakeTexture', window, tacs_er.Stimuli(stimNames{ii}));
        loadStr = sprintf('Loading Stimuli  %g %%',floor(ii/nTrials*100));
        DrawFormattedText(window,loadStr,'center','center',255,50);
        Screen('Flip',window);
    end
    
    %if PresParams.ConfidenceScale
    % Get coordindates for confidence bar
    [ConfidenceBarCoords] = ConfBarParams(xCenter, yCenter,screenXpixels,screenYpixels);
    RightExtent = max(ConfidenceBarCoords(1,:));
    LeftExtent  = min(ConfidenceBarCoords(1,:));
    TopExtent   = max(ConfidenceBarCoords(2,:));
    BottomExtent= min(ConfidenceBarCoords(2,:));
    CenterYpos  = ConfidenceBarCoords(2,1);
    %end
    %---------------------------------------------------------------------%
    % Participant Instructions
    %---------------------------------------------------------------------%
    if PresParams.ConfidenceScale
        InstString = ['Instructions\n\n' ...
            'You will be presented with images that you might recognized from the previous experiment. '...
            'Your task is to indentify which images were presented before and which ones are new '...
            'by pressing a button. For ' RespConds{1} ' images you will be pressing the '...
            PresParams.RespButtons(1) ' key, for ' RespConds{3} ' images you will press the '...
            PresParams.RespButtons(3) ' key. If you are unsure, you can press the '...
            PresParams.RespButtons(2) ' key. You will have ' num2str(PresParams.MaxResponseTime) ...
            ' seconds to respond, and please do so as quickly and as accurately as possible. '...
            'If you identify the image as old or new, you will also be indicating your '...
            'confidence in that decision by clicking on a scale with the mouse. \n\n'...
            'If no questions, \n'...'
            'Press ''' resumeKey ''' to begin the experiment.'];
    elseif PresParams.UnsureButtonOpt
        InstString = ['Instructions\n\n' ...
            'You will be presented with images that you might recognized from the previous experiment. '...
            'Your task is to indentify which images were presented before and which ones are new '...
            'by pressing a button. For ' RespConds{1} ' images you will be pressing the '...
            PresParams.RespButtons(1) ' key, for ' RespConds{3} ' images you will press the '...
            PresParams.RespButtons(3) ' key. If you are unsure, you can press the '...
            PresParams.RespButtons(2) ' key. You will have ' num2str(PresParams.MaxResponseTime) ...
            ' seconds to respond, and please do so as quickly and as accurately as possible. '...
            'If you identify the image as old or new, you will also be indicating your '...
            'confidence by pressing ' PresParams.RespButtons(1) ', ' PresParams.RespButtons(2) ' and ' ...
            PresParams.RespButtons(3) ' for low, medium and high confidence, respectively. \n\n'...
            'If there are no questions, \n'...'
            'Press ''' resumeKey ''' to begin the experiment.'];
    else
        InstString = ['Instructions\n\n' ...
            'You will be presented with images that you might recognized from the previous task. '...
            'Your task in this part is to indentify which images were presented before and which ones are new '...
            'by pressing a button. \n'...
            'For ' RespConds{1} ' images you press the ''' PresParams.RespButtons(1) ''' key\n'...
            'For ' RespConds{3} ' images you press the ''' PresParams.RespButtons(3) ''' key\n\n' ...
            'After making the old/new judgment on the image, you will also be indicating your '...
            'confidence by pressing:\n'...
            '''' ConfidenceKeys(1) ''' press for low confidence \n'...
            '''' ConfidenceKeys(2) ''' press for mid confidence \n'...
            '''' ConfidenceKeys(3) ''' press for high confidence \n\n'...
            'Please make your responses as quickly and as accurate as possible.\n'...
            'If there are no questions, \n'...'
            'Press ''' resumeKey ''' to begin the experiment.'];
    end
    
    NoRespText = 'No response recorded, please answer quicker!';
    
    if PresParams.UnsureButtonOpt
        WrongKeyText = ['Please use only the following keys: \n \n' ...
            PresParams.RespButtons(1) ' for ' RespConds{1} '\n' ...
            PresParams.RespButtons(2) ' for ' RespConds{2} '\n' ...
            PresParams.RespButtons(3) ' for ' RespConds{3} '\n\n'...
            'Press ''' resumeKey ''' to continue'];
    else
        WrongKeyText = ['Please use only the following keys: \n \n' ...
            PresParams.RespButtons(1) ' for ' RespConds{1} '\n' ...
            PresParams.RespButtons(3) ' for ' RespConds{3} '\n\n'...
            'Press ''' resumeKey ''' to continue'];
    end
    
    WrongConfKeyText = ['Please use only the following keys: \n \n' ...
        ConfidenceKeys(1) ' for low confidence \n' ...
        ConfidenceKeys(2) ' for mid confidence \n' ...
        ConfidenceKeys(3) ' for high confidence \n\n'...
        'Press ''' resumeKey ''' to continue'];
    
    DrawFormattedText(window,InstString, 'wrapat', 'center', 255, ...
        75, [],[],[],[],[xCenter*0.2,0,screenXpixels*0.6,screenYpixels]);
    Screen('Flip',window);
    
    % resume if Resume Key is pressed
    WaitTillResumeKey(resumeKey,activeKeyboardID)
    %%
    %---------------------------------------------------------------------%
    % Trials
    %---------------------------------------------------------------------%
    
    % Maximum priority level
    topPriorityLevel = MaxPriority(window);
    Priority(topPriorityLevel);
    
    if PresParams.StarStimEEG
        [ret] = MatNICStartEEG (['s' num2str(thePath.subjNum)], true, false, socket);
        if ret<0
            error('could not start eeg')
        end
        % start of experiment marker
        sendMarker(1,LSLOutlet)
    end
    
    % Draw blank for a bit
    Screen('Flip', window);
    WaitSecs(PresParams.preStartTime);
    
    % iterate through trials
    for tt = 1:nTrials
        
        % empty flip var
        flip     = [];
        
        % Pre-stimulus fixation (variable ITI).
        Screen('DrawLines', window, fixCrossCoords,PresParams.lineWidthPix, PresParams.fixCrossColor, [0 0], 2);
        [flip.VBLTimestamp, flip.StimulusOnsetTime, flip.FlipTimestamp, flip.Missed, flip.Beampos,] ...
            = Screen('Flip', window);
        TimingInfo.preStimFixFlip{tt}=flip;
        vbl = flip.VBLTimestamp;
        
        % Re-draw for last frame of ITI, taking into account the previous
        % presentation
        itiDur = vbl - 1.5*ifi + ITIsSecs(tt);
        Screen('DrawLines', window, fixCrossCoords,PresParams.lineWidthPix, PresParams.fixCrossColor, [0 0], 2);
        vbl = Screen('Flip', window, itiDur);
        
        % Checks if the Pause Key has been pressed.
        CheckForPauseKey(pauseKey,resumeKey,activeKeyboardID)
        KbQueueFlush(activeKeyboardID);
        
        % Draw Stimulus
        Screen('DrawTexture', window, imgTextures{tt}, [], [], 0);
        [flip.VBLTimestamp, flip.StimulusOnsetTime, flip.FlipTimestamp, flip.Missed, flip.Beampos,] ...
            = Screen('Flip', window, vbl + 0.5*ifi);
        TimingInfo.stimPresFlip{tt}=flip;
        trialTime = GetSecs;
        vbl = flip.VBLTimestamp;
        
        if PresParams.StarStimEEG
            % true retrieval condition codes marker
            sendMarker((tacs_er.RetCondTrialCode(tt)<3)+6,LSLOutlet)
        end
        
        % Wait for Response
        [secs,key]=KbQueueWait2(activeKeyboardID,PresParams.stimDurationInSecs-2*ifi);
        if secs<inf
            if numel(key)>=1
                key = key(1);
            end
                TimingInfo.trialKeyPress{tt} = key;
                TimingInfo.trialRT(tt) = secs-trialTime;
                
                if PresParams.UnsureButtonOpt
                    switch key
                        case PresParams.RespButtons(1)
                            TimingInfo.CondResp{tt} = RespConds{1};
                        case PresParams.RespButtons(2)
                            TimingInfo.CondResp{tt} = RespConds{2};
                        case PresParams.RespButtons(3)
                            TimingInfo.CondResp{tt} = RespConds{3};
                        otherwise
                            TimingInfo.CondResp{tt} = 'wrongkey';
                            DrawFormattedText(window, WrongKeyText, 'center' , 'center');
                            vbl=Screen('Flip', window, vbl + 0.5*ifi);
                            WaitTillResumeKey(resumeKey,activeKeyboardID)
                    end
                else
                    switch key
                        case PresParams.RespButtons(1)
                            TimingInfo.CondResp{tt} = RespConds{1};
                            if PresParams.StarStimEEG
                                % old response marker
                                sendMarker(6,LSLOutlet)
                            end
                        case PresParams.RespButtons(3)
                            TimingInfo.CondResp{tt} = RespConds{3};
                            if PresParams.StarStimEEG
                                % new response marker
                                sendMarker(7,LSLOutlet)
                            end
                        otherwise
                            TimingInfo.CondResp{tt} = 'wrongkey';
                            DrawFormattedText(window, WrongKeyText, 'center' , 'center');
                            vbl=Screen('Flip', window, vbl + 0.5*ifi);
                            WaitTillResumeKey(resumeKey,activeKeyboardID)
                    end
                end            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else
            DrawFormattedText(window,NoRespText, 'center' , 'center');
            vbl=Screen('Flip', window, vbl + 0.5*ifi);
            WaitSecs(1);
        end
        
        if strcmp(TimingInfo.CondResp{tt},'old') || strcmp(TimingInfo.CondResp{tt},'new')
            if PresParams.ConfidenceScale
                % If response is old or new, probe confidence.
                confResp = 0; % confidence response flag.
                SetMouse(xCenter,CenterYpos,window);
                for ii = 1:(MaxConfDescFrames)
                    % Draw Confidence Bar
                    Screen('DrawLines', window, ConfidenceBarCoords,PresParams.lineWidthPix, PresParams.ConfBarColor, [0 0], 2);
                    DrawFormattedText(window, [' Confidence for ' TimingInfo.CondResp{tt}], 'center' , TopExtent-0.1*screenYpixels, PresParams.textColor);
                    DrawFormattedText(window, ' 0 ', LeftExtent-20, TopExtent-0.1*screenYpixels, PresParams.textColor);
                    DrawFormattedText(window, '100' , RightExtent-20, TopExtent-0.1*screenYpixels, PresParams.textColor);
                    
                    % Get the current position of the mouse
                    [mx, ~, buttons] = GetMouse(window);
                    
                    % draw dot indicating position of mouse
                    if mx>= RightExtent
                        vbl=Screen('DrawDots', window, [RightExtent CenterYpos], 15, PresParams.dotColor, [], 2);
                        mx = RightExtent;
                    elseif mx<= LeftExtent
                        vbl=Screen('DrawDots', window, [LeftExtent CenterYpos], 15, PresParams.dotColor, [], 2);
                        mx = LeftExtent;
                    else
                        vbl=Screen('DrawDots', window, [mx CenterYpos], 15, PresParams.dotColor, [], 2);
                    end
                    
                    % If there was a click, record. else continue to draw
                    if sum(buttons)
                        confResp = 1;
                        Pct = (mx-LeftExtent)/(RightExtent-LeftExtent);
                        TimingInfo.Confidence(tt) = Pct;
                        
                        buttonPress = sprintf('Response: %.2g',Pct);
                        DrawFormattedText(window,buttonPress,'center', BottomExtent+100, PresParams.textColor);
                        HideCursor();
                        vbl=Screen('Flip', window, vbl + 0.5* ifi);
                        SetMouse(xCenter,CenterYpos,window);
                        WaitSecs(0.5);
                        break
                    else
                        vbl  = Screen('Flip', window, vbl +  0.5*ifi);
                    end
                end
                if ~confResp
                    DrawFormattedText(window,NoRespText, 'center' , 'center');
                    vbl=Screen('Flip', window, vbl + 0.5*ifi);
                    WaitSecs(1);
                end
            else
                DrawFormattedText(window, [' Confidence for ' TimingInfo.CondResp{tt} '?'], 'center' , TopExtent-0.1*screenYpixels, PresParams.textColor);
                vbl=Screen('Flip', window, vbl + 0.5* ifi);
                if PresParams.StarStimEEG
                    % confidence probe
                    sendMarker(8,LSLOutlet)
                end
                % Wait for Response
                [secs,key]=KbQueueWait2(activeKeyboardID,PresParams.MaxConfDecInSecs-2*ifi);
                if secs<inf && ischar(key)
                    if numel(key)>1
                        key = key(1);
                    end
                    switch key
                        case ConfidenceKeys(1)
                            TimingInfo.ConfResp{tt} = 'low';
                            TimingInfo.Confidence(tt)   = 1;
                        case ConfidenceKeys(2)
                            TimingInfo.ConfResp{tt} = 'mid';
                            TimingInfo.Confidence(tt)   = 2;
                        case ConfidenceKeys(3)
                            TimingInfo.ConfResp{tt} = 'high';
                            TimingInfo.Confidence(tt)   = 3;
                        otherwise
                            TimingInfo.ConfResp{tt} = 'wrongkey';
                            DrawFormattedText(window, WrongConfKeyText, 'center' , 'center');
                            vbl=Screen('Flip', window, vbl + 0.5*ifi);
                            WaitTillResumeKey(resumeKey,activeKeyboardID)
                    end
                    if PresParams.StarStimEEG
                        % confidence resp
                        sendMarker(8,LSLOutlet)
                    end
                else
                    DrawFormattedText(window,NoRespText, 'center' , 'center');
                    Screen('Flip', window, vbl + 0.5*ifi);
                    WaitSecs(1);
                end
            end
        end
        
        if ~PresParams.SelfPaceFlag
            % Draw Post-Stim blank for remainder of trial (and at least one
            % frame)
            currentTime = GetSecs;
            currentTrialDur = currentTime-trialTime;
            postStimFrames = round((PresParams.TotalTrialDur-(currentTrialDur))/ifi);
            
            [flip.VBLTimestamp, flip.StimulusOnsetTime, flip.FlipTimestamp, flip.Missed, flip.Beampos,] ...
                = Screen('Flip', window, vbl + 0.5*ifi);
            TimingInfo.postStimFlip{tt}=flip;
            vbl = flip.VBLTimestamp;
            if postStimFrames>2
                for ii = 1:(postStimFrames-1)
                    vbl  = Screen('Flip', window,vbl + 0.5*ifi);
                end
            end
        end
        
        % save every PresParams.SaveEvNtrials
        if mod(tt,PresParams.SaveEvNtrials)==0
            tempName = sprintf('/tacs_er.s%i.test.%s.mat', thePath.subjNum, datestr(now,'dd.mm.yyyy.HH.MM'));
            save([thePath.subjectPath,tempName],'TimingInfo');
        end
        
        % Discard used image texture
        Screen('Close', imgTextures{tt})
        
    end
    %---------------------------------------------------------------------%
    % End of Experiment
    %---------------------------------------------------------------------%
    % store additional outputs
    ret_out = [];
    ret_out.PresParams  = PresParams;
    tacs_er.Stimuli = []; % don't re-store stimuli
    ret_out.expInfo     = tacs_er;
    ret_out.TimingInfo  = TimingInfo;
    if PresParams.StarStimEEG
        % end of experiment
        sendMarker(5,LSLOutlet)
        MatNICMarkerCloseLSL(LSLOutlet);
        MatNICStopEEG(socket);
    end
    
    % save
    fileName = [thePath.exptType '.test'];
    cnt = 0;
    while 1
        savePath = strcat(thePath.subjectPath,'/',fileName ,'.mat');
        if ~exist(savePath,'file')
            save(savePath,'ret_out')
            break
        else
            cnt = cnt+1;
            warning(strcat(fileName,' already existed.'))
            fileName = strcat(thePath.exptType,'.test','-',num2str(cnt),'.mat');
            warning(strcat('saving as ', fileName))
        end
    end
    
    InstString = ['End of Experiment.\n \n' ...
        'Press ''' resumeKey ''' to exit.'];
    
    DrawFormattedText(window,InstString, 'center', 'center', 255, 40);
    Screen('Flip',window);
    WaitTillResumeKey(resumeKey,activeKeyboardID)
    KbQueueStop(activeKeyboardID);
    
    msg='allGood';
catch msg
    if PresParams.StarStimEEG
        MatNICMarkerCloseLSL(LSLOutlet);
        MatNICStopEEG(socket);
    end
    sca
    ShowCursor
    keyboard
end

% Clear the screen
Priority(0);
sca;
Screen('CloseAll');
ShowCursor;

%end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% auxiliary functions and definitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %-------------------------------------------------------------------------%
% % fixCrossCoords
% % Set Fixation Cross Coordinates
% %-------------------------------------------------------------------------%
% function fixCrossCoords = fixCross(xCenter, yCenter,screenXpixels,screenYpixels)
%
% fixCrossXlength = max(0.02*screenXpixels,0.02*screenYpixels); % max of 2% screen dims
% fixCrossYlength = fixCrossXlength;
%
% LeftExtent  = xCenter-fixCrossXlength/2;
% RightExtent = xCenter+fixCrossXlength/2 ;
% BottomExtent = yCenter+fixCrossYlength/2 ;
% TopExtent   =  yCenter- fixCrossYlength/2 ;
%
% fixCrossXCoords   = [LeftExtent RightExtent; yCenter yCenter];
% fixCrossYCoords   = [xCenter xCenter; BottomExtent TopExtent];
%
% fixCrossCoords       = [fixCrossXCoords fixCrossYCoords];
%
% end
%
% %-------------------------------------------------------------------------%
% % WaitTillResumeKey
% % Wait until Resume Key is pressed on the keyboard
% %-------------------------------------------------------------------------%
% function WaitTillResumeKey(resumeKey,activeKeyboardID)
%
% KbQueueFlush(activeKeyboardID);
% while 1
%     [pressed,firstPress] = KbQueueCheck(activeKeyboardID);
%     if pressed
%         if strcmp(resumeKey,KbName(firstPress));
%             break
%         end
%     end
%     WaitSecs(0.1);
% end
% KbQueueFlush(activeKeyboardID);
% end

% %-------------------------------------------------------------------------%
% % CheckForPauseKey
% % Check if the resume key has been pressed, and pause exection until resume
% % key is pressed.
% %-------------------------------------------------------------------------%
% function CheckForPauseKey(pauseKey,resumeKey,activeKeyboardID)
%
% [pressed,firstPress] = KbQueueCheck(activeKeyboardID);
% if pressed
%     if strcmp(pauseKey,KbName(firstPress));
%         WaitTillResumeKey(resumeKey,activeKeyboardID)
%     end
% end
% end
%
% %-------------------------------------------------------------------------%
% % ConfBarParams
% % Confidence Bar Coordinates
% %-------------------------------------------------------------------------%
% function [ConfidenceBarCoords] = ConfBarParams(xCenter, yCenter,screenXpixels,screenYpixels)
%
% % Here we set the size of our confidence bar
% BarLength = 0.5*screenXpixels; % 50% of the width screen
% HeightOfBarWhisks = 0.05*screenYpixels; % 5% of the height of screen
%
% LeftExtent  = xCenter-BarLength/2;
% MidLeftExtent = xCenter-BarLength/4;
% RightExtent = xCenter+BarLength/2 ;
% MidRightExtent = xCenter+BarLength/4;
% BottomExtent = yCenter+HeightOfBarWhisks/2 ;
% TopExtent   =  yCenter- HeightOfBarWhisks/2 ;
%
% HorizontalBarCoords   = [LeftExtent RightExtent; yCenter yCenter];
%
% LeftWhisk             = [LeftExtent LeftExtent ; BottomExtent TopExtent];
% RightWhisk            = [RightExtent RightExtent ; BottomExtent TopExtent];
% MidLeftWhisk          = [MidLeftExtent MidLeftExtent ; BottomExtent TopExtent];
% MidRightWhisk         = [MidRightExtent MidRightExtent ; BottomExtent TopExtent];
% CenterWhisk           = [xCenter xCenter ; BottomExtent TopExtent];
%
% ConfidenceBarCoords   = [HorizontalBarCoords LeftWhisk RightWhisk ...
%     MidLeftWhisk MidRightWhisk CenterWhisk];
% end
