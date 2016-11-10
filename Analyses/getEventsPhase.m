function out = getEventsPhase(subj, expt)
out = [];

if strcmp(expt,'tacs_er')
    out.info = subjFileInfo(subj);
elseif strcmp(expt,'tacs_enc_xdiva')
    out.info = subjFileInfo_xdiva(subj);
elseif strcmp(expt,'tacs_er_objstim')
    out.info = subjFileInfo(subj,expt);
end

stimData=load([out.info.dataPath out.info.stimulationFileName]);
eegData =load([out.info.dataPath out.info.eegEncodingFileName]);

%% Timing and event markers

stimSamples = stimData(:,9)-stimData(1,9);
stimTime    = stimSamples/1000;  % in seconds
eegSamples  = eegData(:,13)-stimData(1,9);
eegTime     = eegSamples/1000;  % in seconds
eventCodesSamps = eegSamples(eegData(:,12)~=0);
eventCodes  = eegData(eegData(:,12)~=0,12);

%%  get the phase of the samples

% Fz-FC2
stimElec = (stimData(:,1)-stimData(:,2))/1e3; % in mA
stimElecAngle = angle(hilbert(stimElec));

% 6Hz 2nd order IIR Resonator Filter
 a = [0.0124         0   -0.0124];
 b =  [1.0000   -1.9696    0.9752];
 
%Pz Signal
x = eegData(:,7)/1e6; %in mV
y = filtfilt(a,b,x);
eegElecAngle = angle(hilbert(y));
 
%stimElecAngle3Quant = quant(stimElecAngle,pi/3);

%% get the phase of each event

out.EventCods           = eventCodes;
out.EventSamples        = eventCodesSamps;
out.FixationIdx         = find(eventCodes==2); % fixations
out.SmallIdx            = find(eventCodes==3); % Small stims
out.BigIdx              = find(eventCodes==4); % Big stims
out.AllStimIdx          = find(eventCodes==3|eventCodes==4);

out.FzFC2Phases         = stimElecAngle(eventCodesSamps(out.AllStimIdx));
out.PzPhases            = eegElecAngle( ismember(eegSamples,eventCodesSamps(out.AllStimIdx)));

save([out.info.dataPath 'EventsPhase'],'out')
