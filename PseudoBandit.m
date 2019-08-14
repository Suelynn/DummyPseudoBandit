function PseudoBandit()
% Learning to Nose Poke side ports, max number of trials ~1400

global BpodSystem
global TaskParameters
global nidaq

%% Task parameters
TaskParameters = BpodSystem.ProtocolSettings;
if isempty(fieldnames(TaskParameters))
    %general
    TaskParameters.GUI.Ports_LMR = '234';
    TaskParameters.GUI.FI = 1; % (s)
    TaskParameters.GUI.VI = false;
    TaskParameters.GUIMeta.VI.Style = 'checkbox';
    TaskParameters.GUI.ChoiceDeadline = 10;
    TaskParameters.GUI.LightGuided = 0;
    TaskParameters.GUIMeta.LightGuided.Style = 'checkbox';
    TaskParameters.GUIPanels.General = {'Ports_LMR','FI','VI','ChoiceDeadline','LightGuided'};
    
    %"stimulus"
    TaskParameters.GUI.PlayStimulus = 1;
    TaskParameters.GUIMeta.PlayStimulus.Style = 'popupmenu';
    TaskParameters.GUIMeta.PlayStimulus.String = {'No stim.','Click stim.','Freq. stim.'};
    TaskParameters.GUI.MinSampleTime = 0.2;
    TaskParameters.GUI.MaxSampleTime = 0.2;
    TaskParameters.GUI.AutoIncrSample = false;
    TaskParameters.GUIMeta.AutoIncrSample.Style = 'checkbox';
    TaskParameters.GUI.MinSampleIncr = 0.01;
    TaskParameters.GUI.MinSampleDecr = 0.005;
    TaskParameters.GUI.EarlyWithdrawalTimeOut = 1;
    TaskParameters.GUI.EarlyWithdrawalNoise = true;
    TaskParameters.GUIMeta.EarlyWithdrawalNoise.Style='checkbox';
    TaskParameters.GUI.GracePeriod = 0;
    TaskParameters.GUI.SampleTime = TaskParameters.GUI.MinSampleTime;
    TaskParameters.GUIMeta.SampleTime.Style = 'text';
    TaskParameters.GUIPanels.Sampling = {'PlayStimulus','MinSampleTime','MaxSampleTime','AutoIncrSample','MinSampleIncr','MinSampleDecr','EarlyWithdrawalTimeOut','EarlyWithdrawalNoise','GracePeriod','SampleTime'};
    
    %Reward
    TaskParameters.GUI.rewardAmount = 35;
    TaskParameters.GUI.banditSigma = 2.8; %reward noise
    TaskParameters.GUI.TrialsPerBout = 80; %average trial length of each bout
    TaskParameters.GUI.CenterPortRewAmount = 0.5;
    TaskParameters.GUI.CenterPortProb = 1;
    TaskParameters.GUI.RewardProb = 1;
    TaskParameters.GUI.Deplete = false;
    TaskParameters.GUIMeta.Deplete.Style = 'checkbox';
    TaskParameters.GUI.DepleteRate = 0.8;
    TaskParameters.GUI.Jackpot = 1;
    TaskParameters.GUIMeta.Jackpot.Style = 'popupmenu';
    TaskParameters.GUIMeta.Jackpot.String = {'No Jackpot','Fixed Jackpot','Decremental Jackpot','RewardCenterPort'};
    TaskParameters.GUI.JackpotMin = 1;
    TaskParameters.GUI.JackpotTime = 1;
    TaskParameters.GUIMeta.JackpotTime.Style = 'text';
    TaskParameters.GUIPanels.Reward = {'rewardAmount','banditSigma', 'TrialsPerBout','CenterPortRewAmount','CenterPortProb','RewardProb','Deplete','DepleteRate','Jackpot','JackpotMin','JackpotTime'};
        
    %Reward Dealy
    TaskParameters.GUI.DelayMean = 0;
    TaskParameters.GUI.DelaySigma=0;
    TaskParameters.GUI.DelayGracePeriod=0;
    TaskParameters.GUIPanels.RewardDelay = {'DelayMean','DelaySigma','DelayGracePeriod'};
    
    %% Photometry
    %photometry general
    TaskParameters.GUI.Photometry=1;
    TaskParameters.GUIMeta.Photometry.Style='checkbox';
    TaskParameters.GUI.DbleFibers=0;
    TaskParameters.GUIMeta.DbleFibers.Style='checkbox';
    TaskParameters.GUIMeta.DbleFibers.String='Auto';
    TaskParameters.GUI.Isobestic405=0;
    TaskParameters.GUIMeta.Isobestic405.Style='checkbox';
    TaskParameters.GUIMeta.Isobestic405.String='Auto';
    TaskParameters.GUI.RedChannel=1;
    TaskParameters.GUIMeta.RedChannel.Style='checkbox';
    TaskParameters.GUIMeta.RedChannel.String='Auto';    
    TaskParameters.GUIPanels.PhotometryRecording={'Photometry','DbleFibers','Isobestic405','RedChannel'};
    
    %plot photometry
    TaskParameters.GUI.TimeMin=-1;
    TaskParameters.GUI.TimeMax=15;
    TaskParameters.GUI.NidaqMin=-5;
    TaskParameters.GUI.NidaqMax=10;
    TaskParameters.GUI.SidePokeIn=1;
	TaskParameters.GUIMeta.SidePokeIn.Style='checkbox';
    TaskParameters.GUI.SidePokeLeave=1;
	TaskParameters.GUIMeta.SidePokeLeave.Style='checkbox';
    TaskParameters.GUI.RewardDelivery=1;
	TaskParameters.GUIMeta.RewardDelivery.Style='checkbox';    
     TaskParameters.GUI.BaselineBegin=0.1;
    TaskParameters.GUI.BaselineEnd=1.1;
    TaskParameters.GUIPanels.PhotometryPlot={'TimeMin','TimeMax','NidaqMin','NidaqMax','SidePokeIn','SidePokeLeave','RewardDelivery','BaselineBegin','BaselineEnd'};
    
    %% Nidaq and Photometry
    TaskParameters.GUI.PhotometryVersion=1;
    TaskParameters.GUI.Modulation=1;
    TaskParameters.GUIMeta.Modulation.Style='checkbox';
    TaskParameters.GUIMeta.Modulation.String='Auto';
	TaskParameters.GUI.NidaqDuration=4;
    TaskParameters.GUI.NidaqSamplingRate=6100;
    TaskParameters.GUI.DecimateFactor=610;
    TaskParameters.GUI.LED1_Name='Fiber1 470-A1';
    TaskParameters.GUIMeta.LED1_Name.Style='edittext';
    TaskParameters.GUI.LED1_Amp=1;
    TaskParameters.GUI.LED1_Freq=211;
    TaskParameters.GUI.LED2_Name='Fiber1 405 / 565';
    TaskParameters.GUIMeta.LED2_Name.Style='edittext';
    TaskParameters.GUI.LED2_Amp=5;
    TaskParameters.GUI.LED2_Freq=531;
    TaskParameters.GUI.LED1b_Name='Fiber2 470-mPFC';
    TaskParameters.GUIMeta.LED1b_Name.Style='edittext';
    TaskParameters.GUI.LED1b_Amp=2;
    TaskParameters.GUI.LED1b_Freq=531;

    TaskParameters.GUIPanels.PhotometryNidaq={'PhotometryVersion','Modulation','NidaqDuration',...
                            'NidaqSamplingRate','DecimateFactor',...
                            'LED1_Name','LED1_Amp','LED1_Freq',...
                            'LED2_Name','LED2_Amp','LED2_Freq',...
                            'LED1b_Name','LED1b_Amp','LED1b_Freq'};
                        
    %% rig-specific
        TaskParameters.GUI.nidaqDev='Dev2';
        TaskParameters.GUIMeta.nidaqDev.Style='edittext';
        
        TaskParameters.GUIPanels.PhotometryRig={'nidaqDev'};
    
    TaskParameters.GUITabs.General = {'General','Sampling','Reward','RewardDelay'};
    TaskParameters.GUITabs.Photometry = {'PhotometryRecording','PhotometryNidaq','PhotometryPlot','PhotometryRig'};
    
        
    TaskParameters.GUI = orderfields(TaskParameters.GUI);
    TaskParameters.Figures.OutcomePlot.Position = [200, 200, 1000, 400];
end
BpodParameterGUI('init', TaskParameters);

%% Initialize reward trains for left and right ports
left_rewardTrain=generate_rewardTrain();
right_rewardTrain=generate_rewardTrain();

n=min(numel(left_rewardTrain),numel(right_rewardTrain));
BpodSystem.Data.Custom.RewardMagnitude=horzcat(left_rewardTrain(1:n),right_rewardTrain(1:n));

%% Initializing everything else
BpodSystem.Data.Custom.ChoiceLeft = NaN;
BpodSystem.Data.Custom.SampleTime(1) = TaskParameters.GUI.MinSampleTime;
BpodSystem.Data.Custom.EarlyWithdrawal(1) = false;
BpodSystem.Data.Custom.Jackpot(1) = false;
%BpodSystem.Data.Custom.RewardMagnitude = [TaskParameters.GUI.rewardAmount,TaskParameters.GUI.rewardAmount];
BpodSystem.Data.Custom.CenterPortRewAmount =TaskParameters.GUI.CenterPortRewAmount;
BpodSystem.Data.Custom.Rewarded = false;
BpodSystem.Data.Custom.CenterPortRewarded = false;
BpodSystem.Data.Custom.GracePeriod = 0;
BpodSystem.Data.Custom.LightLeft = rand(1,1)<0.5;
BpodSystem.Data.Custom.RewardAvailable = rand(1,1)<TaskParameters.GUI.RewardProb;
BpodSystem.Data.Custom.RewardDelay = TaskParameters.GUI.DelayMean;
BpodSystem.Data.Custom = orderfields(BpodSystem.Data.Custom);
%server data
[~,BpodSystem.Data.Custom.Rig] = system('hostname');
[~,BpodSystem.Data.Custom.Subject] = fileparts(fileparts(fileparts(fileparts(BpodSystem.DataPath))));
BpodSystem.Data.Custom.PsychtoolboxStartup=false;
BpodSystem.Data.Custom.MaxSampleTime = 1; %only relevant for max stimulus length
[BpodSystem.Data.Custom.RightClickTrain,BpodSystem.Data.Custom.LeftClickTrain] = getClickStimulus(BpodSystem.Data.Custom.MaxSampleTime);
BpodSystem.Data.Custom.FreqStimulus = getFreqStimulus(BpodSystem.Data.Custom.MaxSampleTime);

BpodSystem.SoftCodeHandlerFunction = 'SoftCodeHandler';

%% Configuring PulsePal
load PulsePalParamStimulus.mat
load PulsePalParamFeedback.mat
BpodSystem.Data.Custom.PulsePalParamStimulus=PulsePalParamStimulus;
BpodSystem.Data.Custom.PulsePalParamFeedback=PulsePalParamFeedback;
clear PulsePalParamFeedback PulsePalParamStimulus
if ~BpodSystem.EmulatorMode
    ProgramPulsePal(BpodSystem.Data.Custom.PulsePalParamStimulus);
    SendCustomPulseTrain(1, BpodSystem.Data.Custom.RightClickTrain, ones(1,length(BpodSystem.Data.Custom.RightClickTrain))*5);
    SendCustomPulseTrain(2, BpodSystem.Data.Custom.LeftClickTrain, ones(1,length(BpodSystem.Data.Custom.LeftClickTrain))*5); 
    if TaskParameters.GUI.PlayStimulus == 3
        InitiatePsychtoolbox();
        PsychToolboxSoundServer('Load', 1, BpodSystem.Data.Custom.FreqStimulus);
    end
end

%% Initialize plots
BpodSystem.ProtocolFigures.SideOutcomePlotFig = figure('Position', TaskParameters.Figures.OutcomePlot.Position,'name','Outcome plot','numbertitle','off', 'MenuBar', 'none', 'Resize', 'off');
BpodSystem.GUIHandles.OutcomePlot.HandleOutcome = axes('Position',    [  .055            .15 .91 .3]);
BpodSystem.GUIHandles.OutcomePlot.HandleGracePeriod = axes('Position',  [1*.05           .6  .1  .3], 'Visible', 'off');
BpodSystem.GUIHandles.OutcomePlot.HandleTrialRate = axes('Position',    [3*.05 + 2*.08   .6  .1  .3], 'Visible', 'off');
BpodSystem.GUIHandles.OutcomePlot.HandleST = axes('Position',           [5*.05 + 4*.08   .6  .1  .3], 'Visible', 'off');
BpodSystem.GUIHandles.OutcomePlot.HandleMT = axes('Position',           [6*.05 + 6*.08   .6  .1  .3], 'Visible', 'off');
Bandit_PlotSideOutcome(BpodSystem.GUIHandles.OutcomePlot,'init');


%% NIDAQ Initialization and Plots
if TaskParameters.GUI.Photometry
if (TaskParameters.GUI.DbleFibers+TaskParameters.GUI.Isobestic405+TaskParameters.GUI.RedChannel)*TaskParameters.GUI.Photometry >1
    disp('Error - Incorrect photometry recording parameters')
    return
end

Nidaq_photometry('ini');

FigNidaq1=Online_NidaqPlot('ini','470');
if TaskParameters.GUI.DbleFibers || TaskParameters.GUI.Isobestic405 || TaskParameters.GUI.RedChannel
    FigNidaq2=Online_NidaqPlot('ini','channel2');
end
end

%% Main loop
RunSession = true;
iTrial = 1;

while RunSession
    TaskParameters = BpodParameterGUI('sync', TaskParameters);
    InitiatePsychtoolbox();
    
    sma = stateMatrix(iTrial);
    SendStateMatrix(sma);
    
    %% NIDAQ Get nidaq ready to start
    if TaskParameters.GUI.Photometry
    Nidaq_photometry('WaitToStart');
    end
    
    %% Run Trial
    RawEvents = RunStateMatrix;
    
    %% NIDAQ Stop acquisition and save data in bpod structure
    if TaskParameters.GUI.Photometry
    Nidaq_photometry('Stop');
    [PhotoData,Photo2Data]=Nidaq_photometry('Save');
    BpodSystem.Data.NidaqData{iTrial}=PhotoData;
    if TaskParameters.GUI.DbleFibers || TaskParameters.GUI.RedChannel
        BpodSystem.Data.Nidaq2Data{iTrial}=Photo2Data;
    end
    end
    
    %% Bpod save
    if ~isempty(fieldnames(RawEvents))
        BpodSystem.Data = AddTrialEvents(BpodSystem.Data,RawEvents);
        SaveBpodSessionData;
    end
    HandlePauseCondition; % Checks to see if the protocol is paused. If so, waits until user resumes.
    if BpodSystem.BeingUsed == 0
        return
    end
    
    %% update fields
    updateCustomDataFields(iTrial)
    
    %% update figures
    NosePoke_PlotSideOutcome(BpodSystem.GUIHandles.OutcomePlot,'update',iTrial);

    % plot photometry data
    if TaskParameters.GUI.Photometry
            
        Alignments = {[],[],[]};
        if TaskParameters.GUI.SidePokeIn && ~BpodSystem.Data.Custom.EarlyWithdrawal(iTrial)
            Alignments{1} = 'wait_Sin';
        end
        if TaskParameters.GUI.SidePokeLeave && ~BpodSystem.Data.Custom.EarlyWithdrawal(iTrial) && BpodSystem.Data.Custom.Rewarded(iTrial)==0
            Alignments{2} = 'ITI';
        end
        if TaskParameters.GUI.RewardDelivery && BpodSystem.Data.Custom.Rewarded(iTrial)==1
            Alignments{3} = 'water_';
        end
        
        for k =1:length(Alignments)
             align = Alignments{k};
             if ~isempty(align)
            [currentNidaq1, rawNidaq1]=Online_NidaqDemod(PhotoData(:,1),nidaq.LED1,TaskParameters.GUI.LED1_Freq,TaskParameters.GUI.LED1_Amp,align);
            FigNidaq1=Online_NidaqPlot('update',[],FigNidaq1,currentNidaq1,rawNidaq1,k);
            
            if TaskParameters.GUI.Isobestic405 || TaskParameters.GUI.DbleFibers || TaskParameters.GUI.RedChannel
                if TaskParameters.GUI.Isobestic405
                    [currentNidaq2, rawNidaq2]=Online_NidaqDemod(PhotoData(:,1),nidaq.LED2,TaskParameters.GUI.LED2_Freq,TaskParameters.GUI.LED2_Amp,align);
                elseif TaskParameters.GUI.RedChannel
                    [currentNidaq2, rawNidaq2]=Online_NidaqDemod(Photo2Data(:,1),nidaq.LED2,TaskParameters.GUI.LED2_Freq,TaskParameters.GUI.LED2_Amp,align);
                elseif TaskParameters.GUI.DbleFibers
                    [currentNidaq2, rawNidaq2]=Online_NidaqDemod(Photo2Data(:,1),nidaq.LED2,TaskParameters.GUI.LED1b_Freq,TaskParameters.GUI.LED1b_Amp,align);
                end
                FigNidaq2=Online_NidaqPlot('update',[],FigNidaq2,currentNidaq2,rawNidaq2,k);
            end
             end%if non-empty align
        end%alignment loop
    end%if photometry
    
    iTrial = iTrial + 1;    
end

%% photometry check
if TaskParameters.GUI.Photometry
    thismax=max(PhotoData(TaskParameters.GUI.NidaqSamplingRate:TaskParameters.GUI.NidaqSamplingRate*2,1))
    if thismax>4 || thismax<0.3
        disp('WARNING - Something is wrong with fiber #1 - run check-up! - unpause to ignore')
        BpodSystem.Pause=1;
        HandlePauseCondition;
    end
    if TaskParameters.GUI.DbleFibers
    thismax=max(Photo2Data(TaskParameters.GUI.NidaqSamplingRate:TaskParameters.GUI.NidaqSamplingRate*2,1))
    if thismax>4 || thismax<0.3
        disp('WARNING - Something is wrong with fiber #2 - run check-up! - unpause to ignore')
        BpodSystem.Pause=1;
        HandlePauseCondition;
    end
    end
end

end

function sma = stateMatrix(iTrial)
global BpodSystem
global TaskParameters
%% Define ports
LeftPort = floor(mod(TaskParameters.GUI.Ports_LMR/100,10));
CenterPort = floor(mod(TaskParameters.GUI.Ports_LMR/10,10));
RightPort = mod(TaskParameters.GUI.Ports_LMR,10);
LeftPortOut = strcat('Port',num2str(LeftPort),'Out');
CenterPortOut = strcat('Port',num2str(CenterPort),'Out');
RightPortOut = strcat('Port',num2str(RightPort),'Out');
LeftPortIn = strcat('Port',num2str(LeftPort),'In');
CenterPortIn = strcat('Port',num2str(CenterPort),'In');
RightPortIn = strcat('Port',num2str(RightPort),'In');

LeftValve = 2^(LeftPort-1);
CenterValve = 2^(CenterPort-1);
RightValve = 2^(RightPort-1);

LeftValveTime  = GetValveTimes(BpodSystem.Data.Custom.RewardMagnitude(iTrial,1), LeftPort);
if rand(1,1) <= TaskParameters.GUI.CenterPortProb && TaskParameters.GUI.Jackpot == 4
    CenterValveTime  = GetValveTimes(BpodSystem.Data.Custom.CenterPortRewAmount(iTrial), CenterPort);
else
    CenterValveTime=0;
end
RightValveTime  = GetValveTimes(BpodSystem.Data.Custom.RewardMagnitude(iTrial,2), RightPort);

if TaskParameters.GUI.Jackpot == 3 % Decremental Jackpot reward
    JackpotFactor = max(2,10 - sum(BpodSystem.Data.Custom.Jackpot)); 
else 
    JackpotFactor = 2; % Fixed Jackpot reward
end
LeftValveTimeJackpot  = JackpotFactor*GetValveTimes(BpodSystem.Data.Custom.RewardMagnitude(iTrial,1), LeftPort);
RightValveTimeJackpot  = JackpotFactor*GetValveTimes(BpodSystem.Data.Custom.RewardMagnitude(iTrial,2), RightPort);

if TaskParameters.GUI.PlayStimulus == 1 %no
    StimStartOutput = {};
    StimStart2Output = {};
    StimStopOutput = {};
elseif TaskParameters.GUI.PlayStimulus == 2 %click
    StimStartOutput = {'BNCState',1};
    StimStart2Output = {'BNCState',1};
    StimStopOutput = {'BNCState',0};
elseif TaskParameters.GUI.PlayStimulus == 3 %freq
    StimStartOutput = {'SoftCode',21};
    StimStopOutput = {'SoftCode',22};
    StimStart2Output = {};
end

if TaskParameters.GUI.EarlyWithdrawalNoise
    PunishSoundAction=11;
else
    PunishSoundAction=0;
end

%light guided task
if TaskParameters.GUI.LightGuided 
    if BpodSystem.Data.Custom.LightLeft(iTrial)
        LeftLight=255;
        RightLight = 0;
    elseif ~BpodSystem.Data.Custom.LightLeft(iTrial)
        LeftLight=0;
        RightLight=255;
    else
        error('Light guided state matrix error');
    end
else
    LeftLight=255;
    RightLight=255;
end

% reward available?
RightWaitAction = 'ITI';
LeftWaitAction = 'ITI';
if BpodSystem.Data.Custom.RewardAvailable(iTrial)
    DelayTime = BpodSystem.Data.Custom.RewardDelay(iTrial);
    if TaskParameters.GUI.LightGuided && BpodSystem.Data.Custom.LightLeft(iTrial)
        LeftWaitAction = 'water_L';
    elseif TaskParameters.GUI.LightGuided && ~BpodSystem.Data.Custom.LightLeft(iTrial)
        RightWaitAction = 'water_R';
    else
        LeftWaitAction = 'water_L';
        RightWaitAction = 'water_R';
    end
else
    DelayTime = 30;
end

    
    
sma = NewStateMatrix();
sma = SetGlobalTimer(sma,1,TaskParameters.GUI.SampleTime);
sma = SetGlobalTimer(sma,2,DelayTime);
sma = AddState(sma, 'Name', 'state_0',...
    'Timer', 0,...
    'StateChangeConditions', {'Tup', 'wait_Cin'},...
    'OutputActions', {});
sma = AddState(sma, 'Name', 'wait_Cin',...
    'Timer', 0,...
    'StateChangeConditions', {CenterPortIn, 'StartSampling'},...
    'OutputActions', {strcat('PWM',num2str(CenterPort)),255});
sma = AddState(sma, 'Name', 'StartSampling',...
    'Timer', 0.01,...
    'StateChangeConditions', {'Tup', 'Sampling'},...S
    'OutputActions', {'GlobalTimerTrig',1});
sma = AddState(sma, 'Name', 'Sampling',...
    'Timer', TaskParameters.GUI.SampleTime,...
    'StateChangeConditions', {CenterPortOut, 'GracePeriod','Tup','stillSampling','GlobalTimer1_End','stillSampling'},...
    'OutputActions', StimStartOutput);
sma = AddState(sma, 'Name', 'GracePeriod',...
    'Timer', TaskParameters.GUI.GracePeriod,...
    'StateChangeConditions', {CenterPortIn, 'Sampling','Tup','EarlyWithdrawal','GlobalTimer1_End','EarlyWithdrawal',LeftPortIn,'EarlyWithdrawal',RightPortIn,'EarlyWithdrawal'},...
    'OutputActions',{});
if TaskParameters.GUI.Jackpot == 1 % Jackpot == none
    sma = AddState(sma, 'Name', 'stillSampling',...
        'Timer', TaskParameters.GUI.ChoiceDeadline,...
        'StateChangeConditions', {CenterPortOut, 'stop_stim','Tup','stop_stim'},...
        'OutputActions', StimStart2Output);
    sma = AddState(sma, 'Name', 'stop_stim',...
        'Timer',0.001,...
        'StateChangeConditions', {'Tup','wait_Sin'},...
        'OutputActions',[StimStopOutput {strcat('PWM',num2str(LeftPort)),255,strcat('PWM',num2str(RightPort)),255}]);
elseif TaskParameters.GUI.Jackpot == 2 || TaskParameters.GUI.Jackpot == 3 % Jackpot activated (either Fixed or Decremental)
    sma = AddState(sma, 'Name', 'stillSampling',...
        'Timer', TaskParameters.GUI.JackpotTime-TaskParameters.GUI.SampleTime,...
        'StateChangeConditions', {CenterPortOut, 'stop_stim','Tup','stillSamplingJackpot'},...
        'OutputActions', StimStart2Output);
    sma = AddState(sma, 'Name', 'stillSamplingJackpot',...
        'Timer', TaskParameters.GUI.ChoiceDeadline-TaskParameters.GUI.JackpotTime-TaskParameters.GUI.SampleTime,...
        'StateChangeConditions', {CenterPortOut, 'stop_stim_jackpot','Tup','ITI'},...
        'OutputActions', StimStart2Output);
    sma = AddState(sma, 'Name', 'stop_stim_jackpot',...
        'Timer',0.001,...
        'StateChangeConditions', {'Tup','wait_SinJackpot'},...
        'OutputActions',[StimStopOutput {strcat('PWM',num2str(LeftPort)),255,strcat('PWM',num2str(RightPort)),255}]);
    sma = AddState(sma, 'Name', 'stop_stim',...
        'Timer',0.001,...
        'StateChangeConditions', {'Tup','wait_Sin'},...
        'OutputActions',[StimStopOutput {strcat('PWM',num2str(LeftPort)),255,strcat('PWM',num2str(RightPort)),255}]);
elseif TaskParameters.GUI.Jackpot ==4 % 
    sma = AddState(sma, 'Name', 'stillSampling',...
        'Timer', CenterValveTime,...
        'StateChangeConditions', {'Tup','lat_Go_signal'},...
        'OutputActions', [StimStopOutput {'ValveState', CenterValve}]);
    sma = AddState(sma, 'Name', 'lat_Go_signal',...
        'Timer',0,...
        'StateChangeConditions', {CenterPortOut,'wait_Sin'},...
        'OutputActions',{strcat('PWM',num2str(LeftPort)),255,strcat('PWM',num2str(RightPort)),255});
end
sma = AddState(sma, 'Name', 'wait_Sin',...
    'Timer',TaskParameters.GUI.ChoiceDeadline,...
    'StateChangeConditions', {LeftPortIn,'wait_L_start',RightPortIn,'wait_R_start','Tup','ITI'},...
    'OutputActions',{strcat('PWM',num2str(LeftPort)),LeftLight,strcat('PWM',num2str(RightPort)),RightLight});
sma = AddState(sma, 'Name', 'wait_L_start',...
    'Timer',0,...
    'StateChangeConditions', {'Tup','wait_L'},...
    'OutputActions',{'GlobalTimerTrig',2});
sma = AddState(sma, 'Name', 'wait_L',...
    'Timer',DelayTime,...
    'StateChangeConditions', {'Tup',LeftWaitAction,'GlobalTimer2_End',LeftWaitAction,LeftPortOut,'wait_L_grace'},...
    'OutputActions',{strcat('PWM',num2str(LeftPort)),0});
sma = AddState(sma, 'Name', 'wait_L_grace',...
    'Timer',TaskParameters.GUI.DelayGracePeriod,...
    'StateChangeConditions', {'Tup','ITI','GlobalTimer2_End',LeftWaitAction,LeftPortIn,'wait_L'},...
    'OutputActions',{});
sma = AddState(sma, 'Name', 'wait_R_start',...
    'Timer',0,...
    'StateChangeConditions', {'Tup','wait_R'},...
    'OutputActions',{'GlobalTimerTrig',2});
sma = AddState(sma, 'Name', 'wait_R',...
    'Timer',DelayTime,...
    'StateChangeConditions', {'Tup',RightWaitAction,'GlobalTimer2_End',RightWaitAction,RightPortOut,'wait_R_grace'},...
    'OutputActions',{strcat('PWM',num2str(RightPort)),0});
sma = AddState(sma, 'Name', 'wait_R_grace',...
    'Timer',TaskParameters.GUI.DelayGracePeriod,...
    'StateChangeConditions', {'Tup','ITI','GlobalTimer2_End',RightWaitAction,RightPortIn,'wait_R'},...
    'OutputActions',{});
sma = AddState(sma, 'Name', 'water_L',...
    'Timer', LeftValveTime,...
    'StateChangeConditions', {'Tup','ITI'},...
    'OutputActions', {'ValveState', LeftValve});
sma = AddState(sma, 'Name', 'water_R',...
    'Timer', RightValveTime,...
    'StateChangeConditions', {'Tup','ITI'},...
    'OutputActions', {'ValveState', RightValve});
sma = AddState(sma, 'Name', 'wait_SinJackpot',...
    'Timer',TaskParameters.GUI.ChoiceDeadline,...
    'StateChangeConditions', {LeftPortIn,'water_LJackpot',RightPortIn,'water_RJackpot','Tup','ITI'},...
    'OutputActions',{strcat('PWM',num2str(LeftPort)),255,strcat('PWM',num2str(RightPort)),255});
sma = AddState(sma, 'Name', 'water_LJackpot',...
    'Timer', LeftValveTimeJackpot,...
    'StateChangeConditions', {'Tup','ITI'},...
    'OutputActions', {'ValveState', LeftValve});
sma = AddState(sma, 'Name', 'water_RJackpot',...
    'Timer', RightValveTimeJackpot,...
    'StateChangeConditions', {'Tup','ITI'},...
    'OutputActions', {'ValveState', RightValve});
sma = AddState(sma, 'Name', 'EarlyWithdrawal',...
    'Timer', TaskParameters.GUI.EarlyWithdrawalTimeOut,...
    'StateChangeConditions', {'Tup','ITI'},...
    'OutputActions', {'SoftCode',PunishSoundAction});
if TaskParameters.GUI.VI
    sma = AddState(sma, 'Name', 'ITI',...
        'Timer',exprnd(TaskParameters.GUI.FI),...
        'StateChangeConditions',{'Tup','exit'},...
        'OutputActions',{});
else
    sma = AddState(sma, 'Name', 'ITI',...
        'Timer',TaskParameters.GUI.FI,...
        'StateChangeConditions',{'Tup','exit'},...
        'OutputActions',{});
end

end

function updateCustomDataFields(iTrial)
global BpodSystem
global TaskParameters
BpodSystem.Data.TrialTypes(iTrial)=1;
%% OutcomeRecord
statesThisTrial = BpodSystem.Data.RawData.OriginalStateNamesByNumber{iTrial}(BpodSystem.Data.RawData.OriginalStateData{iTrial});
BpodSystem.Data.Custom.ST(iTrial) = NaN;
BpodSystem.Data.Custom.MT(iTrial) = NaN;
BpodSystem.Data.Custom.DT(iTrial) = NaN;
BpodSystem.Data.Custom.GracePeriod(1:50,iTrial) = NaN(50,1);
if any(strcmp('Sampling',statesThisTrial))
    if any(strcmp('stillSampling',statesThisTrial)) && any(strcmp('lat_Go_signal',statesThisTrial))==0
        if any(strcmp('stillSamplingJackpot',statesThisTrial))
            BpodSystem.Data.Custom.ST(iTrial) = BpodSystem.Data.RawEvents.Trial{iTrial}.States.stillSamplingJackpot(1,2) - BpodSystem.Data.RawEvents.Trial{iTrial}.States.StartSampling(1,1);
        else
            BpodSystem.Data.Custom.ST(iTrial) = BpodSystem.Data.RawEvents.Trial{iTrial}.States.stillSampling(1,2) - BpodSystem.Data.RawEvents.Trial{iTrial}.States.StartSampling(1,1);
        end
    else
            BpodSystem.Data.Custom.ST(iTrial) = BpodSystem.Data.RawEvents.Trial{iTrial}.States.Sampling(1,end) - BpodSystem.Data.RawEvents.Trial{iTrial}.States.StartSampling(1,1); 
    end
end

% Compute grace period:
if any(strcmp('GracePeriod',statesThisTrial))
    for nb_graceperiod =  1: size(BpodSystem.Data.RawEvents.Trial{iTrial}.States.GracePeriod,1)
        BpodSystem.Data.Custom.GracePeriod(nb_graceperiod,iTrial) = (BpodSystem.Data.RawEvents.Trial{iTrial}.States.GracePeriod(nb_graceperiod,2)...
            -BpodSystem.Data.RawEvents.Trial{iTrial}.States.GracePeriod(nb_graceperiod,1));
    end
end  

if any(strncmp('wait_L',statesThisTrial,6))
    BpodSystem.Data.Custom.ChoiceLeft(iTrial) = 1;
    BpodSystem.Data.Custom.MT(iTrial) = BpodSystem.Data.RawEvents.Trial{iTrial}.States.wait_L_start(1,2) - BpodSystem.Data.RawEvents.Trial{iTrial}.States.wait_Sin(1,1);
    FeedbackPortTimes = BpodSystem.Data.RawEvents.Trial{end}.States.wait_L_start;
    BpodSystem.Data.Custom.DT(iTrial) = FeedbackPortTimes(end,end)-FeedbackPortTimes(1,1);
elseif any(strncmp('wait_R',statesThisTrial,6))
    BpodSystem.Data.Custom.ChoiceLeft(iTrial) = 0;
    BpodSystem.Data.Custom.MT(iTrial) = BpodSystem.Data.RawEvents.Trial{iTrial}.States.wait_R_start(1,2) - BpodSystem.Data.RawEvents.Trial{iTrial}.States.wait_Sin(1,1);
    FeedbackPortTimes = BpodSystem.Data.RawEvents.Trial{end}.States.wait_R_start;
    BpodSystem.Data.Custom.DT(iTrial) = FeedbackPortTimes(end,end)-FeedbackPortTimes(1,1);
elseif any(strcmp('EarlyWithdrawal',statesThisTrial))
    BpodSystem.Data.Custom.EarlyWithdrawal(iTrial) = true;
end



if any(strncmp('water_L',statesThisTrial,7))
    BpodSystem.Data.Custom.Rewarded(iTrial) = true;
elseif any(strncmp('water_R',statesThisTrial,7)) 
    BpodSystem.Data.Custom.Rewarded(iTrial) = true;
end

if any(strcmp('water_LJackpot',statesThisTrial)) || any(strcmp('water_RJackpot',statesThisTrial))
    BpodSystem.Data.Custom.Jackpot(iTrial) = true;
    BpodSystem.Data.Custom.Rewarded(iTrial) = true;
    if any(strcmp('water_LJackpot',statesThisTrial))
        BpodSystem.Data.Custom.MT(iTrial) = BpodSystem.Data.RawEvents.Trial{iTrial}.States.water_LJackpot(1,2) - BpodSystem.Data.RawEvents.Trial{iTrial}.States.wait_SinJackpot(1,1);
    elseif any(strcmp('water_LJackpot',statesThisTrial))
        BpodSystem.Data.Custom.MT(iTrial) = BpodSystem.Data.RawEvents.Trial{iTrial}.States.water_RJackpot(1,2) - BpodSystem.Data.RawEvents.Trial{iTrial}.States.wait_SinJackpot(1,1);
    end
end

if any(strcmp('lat_Go_signal',statesThisTrial))
    BpodSystem.Data.Custom.CenterPortRewarded(iTrial) = true;
end

% correct/error?
BpodSystem.Data.Custom.Correct(iTrial) = NaN;
if TaskParameters.GUI.LightGuided
if BpodSystem.Data.Custom.LightLeft(iTrial) && BpodSystem.Data.Custom.ChoiceLeft(iTrial)==1
    BpodSystem.Data.Custom.Correct(iTrial) = true;
elseif ~BpodSystem.Data.Custom.LightLeft(iTrial) && BpodSystem.Data.Custom.ChoiceLeft(iTrial)==0
    BpodSystem.Data.Custom.Correct(iTrial) = true;
elseif BpodSystem.Data.Custom.LightLeft(iTrial) && BpodSystem.Data.Custom.ChoiceLeft(iTrial)==0
    BpodSystem.Data.Custom.Correct(iTrial) = false;
elseif ~BpodSystem.Data.Custom.LightLeft(iTrial) && BpodSystem.Data.Custom.ChoiceLeft(iTrial)==1
    BpodSystem.Data.Custom.Correct(iTrial) = false;
end
else
    if ~isnan(BpodSystem.Data.Custom.ChoiceLeft)
        BpodSystem.Data.Custom.Correct(iTrial) = true;
    end
end

%% initialize next trial values
BpodSystem.Data.Custom.ChoiceLeft(iTrial+1) = NaN;
BpodSystem.Data.Custom.EarlyWithdrawal(iTrial+1) = false;
BpodSystem.Data.Custom.Jackpot(iTrial+1) = false;
BpodSystem.Data.Custom.ST(iTrial+1) = NaN;
BpodSystem.Data.Custom.MT(iTrial+1) = NaN;
BpodSystem.Data.Custom.DT(iTrial+1) = NaN;
BpodSystem.Data.Custom.Rewarded(iTrial+1) = false;
BpodSystem.Data.Custom.CenterPortRewarded(iTrial+1) = false;
BpodSystem.Data.Custom.GracePeriod(1:50,iTrial+1) = NaN(50,1);
BpodSystem.Data.Custom.LightLeft(iTrial+1) = rand(1,1)<0.5;
BpodSystem.Data.Custom.RewardAvailable(iTrial+1) = rand(1,1)<TaskParameters.GUI.RewardProb;
BpodSystem.Data.Custom.RewardDelay(iTrial+1) = abs( randn(1,1)*TaskParameters.GUI.DelaySigma+TaskParameters.GUI.DelayMean);

%stimuli
if ~BpodSystem.EmulatorMode
    if TaskParameters.GUI.PlayStimulus == 2
        [BpodSystem.Data.Custom.RightClickTrain,BpodSystem.Data.Custom.LeftClickTrain] = getClickStimulus(BpodSystem.Data.Custom.MaxSampleTime);
        SendCustomPulseTrain(1, BpodSystem.Data.Custom.RightClickTrain, ones(1,length(BpodSystem.Data.Custom.RightClickTrain))*5);
        SendCustomPulseTrain(2, BpodSystem.Data.Custom.LeftClickTrain, ones(1,length(BpodSystem.Data.Custom.LeftClickTrain))*5);
    elseif TaskParameters.GUI.PlayStimulus == 3
        InitiatePsychtoolbox();
        BpodSystem.Data.Custom.FreqStimulus = getFreqStimulus(BpodSystem.Data.Custom.MaxSampleTime);
        PsychToolboxSoundServer('Load', 1, BpodSystem.Data.Custom.FreqStimulus);
    end
end

%jackpot time
if  TaskParameters.GUI.Jackpot ==2 || TaskParameters.GUI.Jackpot ==3
    if sum(~isnan(BpodSystem.Data.Custom.ChoiceLeft(1:iTrial)))>10
        TaskParameters.GUI.JackpotTime = max(TaskParameters.GUI.JackpotMin,quantile(BpodSystem.Data.Custom.ST,0.95));
    else
        TaskParameters.GUI.JackpotTime = TaskParameters.GUI.JackpotMin;
    end
end

%reward depletion
% if BpodSystem.Data.Custom.ChoiceLeft(iTrial) == 1 && TaskParameters.GUI.Deplete
%     BpodSystem.Data.Custom.RewardMagnitude(iTrial+1,1) = BpodSystem.Data.Custom.RewardMagnitude(iTrial,1)*TaskParameters.GUI.DepleteRate;
%     BpodSystem.Data.Custom.RewardMagnitude(iTrial+1,2) = TaskParameters.GUI.rewardAmount;
% elseif BpodSystem.Data.Custom.ChoiceLeft(iTrial) == 0 && TaskParameters.GUI.Deplete
%     BpodSystem.Data.Custom.RewardMagnitude(iTrial+1,2) = BpodSystem.Data.Custom.RewardMagnitude(iTrial,2)*TaskParameters.GUI.DepleteRate;
%     BpodSystem.Data.Custom.RewardMagnitude(iTrial+1,1) = TaskParameters.GUI.rewardAmount;
% elseif isnan(BpodSystem.Data.Custom.ChoiceLeft(iTrial)) && TaskParameters.GUI.Deplete
%     BpodSystem.Data.Custom.RewardMagnitude(iTrial+1,:) = BpodSystem.Data.Custom.RewardMagnitude(iTrial,:);
% else
%     BpodSystem.Data.Custom.RewardMagnitude(iTrial+1,:) = [TaskParameters.GUI.rewardAmount,TaskParameters.GUI.rewardAmount];
% end

%center port reward amount
BpodSystem.Data.Custom.CenterPortRewAmount(iTrial+1) =TaskParameters.GUI.CenterPortRewAmount;

%increase sample time
if TaskParameters.GUI.AutoIncrSample
    History = 50; % Rat: History = 50
    Crit = 0.8; % Rat: Crit = 0.8
    if iTrial<5
        ConsiderTrials = iTrial;
    else
        ConsiderTrials = max(1,iTrial-History):1:iTrial;
    end
    ConsiderTrials = ConsiderTrials(~isnan(BpodSystem.Data.Custom.ChoiceLeft(ConsiderTrials))|BpodSystem.Data.Custom.EarlyWithdrawal(ConsiderTrials));
    if sum(~BpodSystem.Data.Custom.EarlyWithdrawal(ConsiderTrials))/length(ConsiderTrials) > Crit % If SuccessRate > crit (80%)
        if ~BpodSystem.Data.Custom.EarlyWithdrawal(iTrial) % If last trial is not EWD
            BpodSystem.Data.Custom.SampleTime(iTrial+1) = min(TaskParameters.GUI.MaxSampleTime,max(TaskParameters.GUI.MinSampleTime,BpodSystem.Data.Custom.SampleTime(iTrial) + TaskParameters.GUI.MinSampleIncr)); % SampleTime increased
        else % If last trial = EWD
            BpodSystem.Data.Custom.SampleTime(iTrial+1) = min(TaskParameters.GUI.MaxSampleTime,max(TaskParameters.GUI.MinSampleTime,BpodSystem.Data.Custom.SampleTime(iTrial))); % SampleTime = max(MinSampleTime or SampleTime)
        end
    elseif sum(~BpodSystem.Data.Custom.EarlyWithdrawal(ConsiderTrials))/length(ConsiderTrials) < Crit/2  % If SuccessRate < crit/2 (40%)
        if BpodSystem.Data.Custom.EarlyWithdrawal(iTrial) % If last trial = EWD
            BpodSystem.Data.Custom.SampleTime(iTrial+1) = max(TaskParameters.GUI.MinSampleTime,min(TaskParameters.GUI.MaxSampleTime,BpodSystem.Data.Custom.SampleTime(iTrial) - TaskParameters.GUI.MinSampleDecr)); % SampleTime decreased
        else
            BpodSystem.Data.Custom.SampleTime(iTrial+1) = min(TaskParameters.GUI.MaxSampleTime,max(TaskParameters.GUI.MinSampleTime,BpodSystem.Data.Custom.SampleTime(iTrial))); % SampleTime = max(MinSampleTime or SampleTime)
        end
    else % If crit/2 < SuccessRate < crit
        BpodSystem.Data.Custom.SampleTime(iTrial+1) =  BpodSystem.Data.Custom.SampleTime(iTrial); % SampleTime unchanged
    end
else
    BpodSystem.Data.Custom.SampleTime(iTrial+1) = TaskParameters.GUI.MinSampleTime;
end
if BpodSystem.Data.Custom.Jackpot(iTrial) % If last trial is Jackpottrial
    BpodSystem.Data.Custom.SampleTime(iTrial+1) = BpodSystem.Data.Custom.SampleTime(iTrial+1)+0.05*TaskParameters.GUI.JackpotTime; % SampleTime = SampleTime + 5% JackpotTime
end
TaskParameters.GUI.SampleTime = BpodSystem.Data.Custom.SampleTime(iTrial+1); % update SampleTime

%send bpod status to server
% try
% script = 'receivebpodstatus.php';
% %create a common "outcome" vector
% outcome = BpodSystem.Data.Custom.ChoiceLeft(1:iTrial); %1=left, 0=right
% outcome(BpodSystem.Data.Custom.EarlyWithdrawal(1:iTrial))=3; %early withdrawal=3
% outcome(BpodSystem.Data.Custom.Jackpot(1:iTrial))=4;%jackpot=4
% SendTrialStatusToServer(script,BpodSystem.Data.Custom.Rig,outcome,BpodSystem.Data.Custom.Subject,BpodSystem.CurrentProtocolName);
% catch
% end

end

function [RightClickTrain,LeftClickTrain]=getClickStimulus(time)
rr = rand(1,1)*0.6+0.2;
l = ceil(rr*100);
r=100-l;
RightClickTrain=GeneratePoissonClickTrain(r,time);
LeftClickTrain=GeneratePoissonClickTrain(l,time);
end

function Sound = getFreqStimulus(time)
StimulusSettings=struct();
            StimulusSettings.SamplingRate = 192000; % Sound card sampling rate;
            StimulusSettings.ramp = 0.003;
            StimulusSettings.nFreq = 18; % Number of different frequencies to sample from
            StimulusSettings.ToneOverlap = 0.6667;
            StimulusSettings.ToneDuration = 0.03;
            StimulusSettings.Noevidence=0;
            StimulusSettings.minFreq = 5000 ;
            StimulusSettings.maxFreq = 40000 ;
            StimulusSettings.UseMiddleOctave=0;
            StimulusSettings.Volume=50;
            StimulusSettings.nTones = floor((time-StimulusSettings.ToneDuration*StimulusSettings.ToneOverlap)/(StimulusSettings.ToneDuration*(1-StimulusSettings.ToneOverlap))); %number of tones
            newFracHigh = rand(1,1);
            [Sound, ~, ~] = GenerateToneCloudDual(newFracHigh, StimulusSettings);
end
function trialReward=generate_rewardTrain()
    global TaskParameters
 
    waterReward=TaskParameters.GUI.rewardAmount;

    boutLen=num2cell(round(normrnd(TaskParameters.GUI.TrialsPerBout,15,18,1))); %generate n bouts per session of trial len 80 ~gaussian distrib
    trialSplit=cellfun(@(x)  ones(x,1), boutLen, 'UniformOutput', false); %generate separate trials

    idx = rand(length(boutLen),1); %randomly assign splits to low, med, high reward

       for i = 1:length(idx)
            if idx(i)<=0.3
                rewardSize=waterReward*0.50;
                trialSplit(i)=cellfun(@(x) x*rewardSize,trialSplit(i),'un',0);
            elseif idx(i)<=0.6
                rewardSize=waterReward*0.75;
                trialSplit(i)=cellfun(@(x) x*rewardSize,trialSplit(i),'un',0);
            else
                rewardSize=waterReward;
                trialSplit(i)=cellfun(@(x) x*rewardSize,trialSplit(i),'un',0);
            end
       end

    trialSplit=vertcat(trialSplit{:});

    decay=0.9836; %add noise centered around the average reward
    noise=normrnd(0,TaskParameters.GUI.banditSigma,[length(trialSplit),1]);
    trialReward=trialSplit*decay+(1-decay)*mean(trialSplit)+noise;
end