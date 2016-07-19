% Display All Channels
function displayAllChannelsSRC(subjectName,expDate,protocolName,folderSourceString,gridType,capType)

if ~exist('folderSourceString','var');   folderSourceString='E:';       end
if ~exist('gridType','var');             gridType='EEG';                end

% this is for running EEG datasets;For microelectrode-this argument need not be passed
if ~exist('capType','var');              capType='actiCap64';           end

folderName = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName);

% Get folders
folderExtract = fullfile(folderName,'extractedData');
folderSegment = fullfile(folderName,'segmentedData');
folderLFP = fullfile(folderSegment,'LFP');
folderSpikes = fullfile(folderSegment,'Spikes');
folderMontage='C:\'; % Have to fix a location to store the montages folder that contains the channel location files for topoplots.

% load LFP Information
[analogChannelsStored,timeVals] = loadlfpInfo(folderLFP);
[neuralChannelsStored,SourceUnitID] = loadspikeInfo(folderSpikes);

% Get Combinations
[parameterCombinations,cIndexValsUnique,tIndexValsUnique,eIndexValsUnique,...
    aIndexValsUnique,sIndexValsUnique] = loadParameterCombinations(folderExtract);

% Get properties of the Stimulus
stimResults = loadStimResults(folderExtract);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display main options
fontSizeSmall = 10; fontSizeMedium = 12; fontSizeLarge = 16;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display layout is modified if the grid type is EEG :
if strcmpi(gridType,'EEG') % the panel gets divided into two halves
    panel_factor=0.5;
    eegPanelWidth=0.25*panel_factor; eegStartPos=0.150;
else
    panel_factor=1;
end
% Make Panels
panelHeight = 0.25; panelStartHeight = 0.725;
staticPanelWidth = 0.25*panel_factor; staticStartPos = 0.025;
dynamicPanelWidth = 0.25; dynamicStartPos = 0.275;
timingPanelWidth = 0.25; timingStartPos = 0.525;
plotOptionsPanelWidth = 0.2; plotOptionsStartPos = 0.775;
backgroundColor = 'w';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Static Panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% staticTitle = [subjectName '_' expDate '_' protocolName];
% hStaticPanel = uipanel('Title','Information','fontSize', fontSizeLarge, ...
%    'Unit','Normalized','Position',[staticStartPos panelStartHeight staticPanelWidth panelHeight]);
%
% if isfield(stimResults,'orientation')
%     oriStr = num2str(stimResults.orientation);
% else
%     oriStr = [num2str(stimResults.baseOrientation0) ',' num2str(stimResults.baseOrientation1)];
% end
%
% staticText = [{ '   '};
%     {['Monkey Name: ' subjectName]}; ...
%     {['Date: ' expDate]}; ...
%     {['Protocol Name: ' protocolName]}; ...
%     {'   '}
%     {['Orientation  (Deg): ' oriStr]}; ...
%     {['Spatial Freq (CPD): ' num2str(stimResults.spatialFrequency)]}; ...
%     {['Azimuth      (Deg): ' num2str(stimResults.azimuth)]}; ...
%     {['Elevation    (Deg): ' num2str(stimResults.elevation)]}; ...
%     {['Sigma        (Deg): ' num2str(stimResults.sigma)]}; ...
%     {['Radius       (Deg): ' num2str(stimResults.radius)]}; ...
%     ];
%
% disp(staticText);
%
% uicontrol('Parent',hStaticPanel,'Unit','Normalized', ...
%     'Position',[0 0 1 1], 'Style','text','String',staticText,'FontSize',fontSizeSmall);

electrodeGridPos = [staticStartPos panelStartHeight staticPanelWidth panelHeight];
hElectrodes = showElectrodeLocations(electrodeGridPos,[],[],[],1,0,gridType);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Dynamic panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dynamicHeight = 0.08; dynamicGap=0.02; dynamicTextWidth = 0.6;
hDynamicPanel = uipanel('Title','Parameters','fontSize', fontSizeLarge, ...
    'Unit','Normalized','Position',[dynamicStartPos panelStartHeight dynamicPanelWidth panelHeight]);

% Contrast
contrastString = getContrastString(cIndexValsUnique,stimResults);
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight],...
    'Style','text','String','Contrast (%)','FontSize',fontSizeSmall);
hContrast = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth 1-(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
    'Style','popup','String',contrastString,'FontSize',fontSizeSmall);

% Temporal Frequency
temporalFreqString = getTemporalFreqString(tIndexValsUnique,stimResults);
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-2*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight], ...
    'Style','text','String','Temporal Freq (Hz)','FontSize',fontSizeSmall);
hTemporalFreq = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position',...
    [dynamicTextWidth 1-2*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
    'Style','popup','String',temporalFreqString,'FontSize',fontSizeSmall);

% EOT Codes
EOTCodeString = getEOTCodeString(eIndexValsUnique);
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-3*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight], ...
    'Style','text','String','EOT Code','FontSize',fontSizeSmall);
hEOTCode = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth 1-3*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
    'Style','popup','String',EOTCodeString,'FontSize',fontSizeSmall);

% Attend Loc
attendLocString = getAttendLocString(aIndexValsUnique);
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-4*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight], ...
    'Style','text','String','Attended Location','FontSize',fontSizeSmall);
hAttendLoc = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth 1-4*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
    'Style','popup','String',attendLocString,'FontSize',fontSizeSmall);

% Stimulus Type
stimTypeString = getStimTypeString(sIndexValsUnique);
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-5*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight], ...
    'Style','text','String','Stimulus Type','FontSize',fontSizeSmall);
hStimType = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth 1-5*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
    'Style','popup','String',stimTypeString,'FontSize',fontSizeSmall);

% Analysis Type
analysisTypeString = 'ERP|Firing Rate|FFT|delta FFT';
uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'Position',[0 1-6*(dynamicHeight+dynamicGap) dynamicTextWidth dynamicHeight], ...
    'Style','text','String','Analysis Type','FontSize',fontSizeSmall);
hAnalysisType = uicontrol('Parent',hDynamicPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, 'Position', ...
    [dynamicTextWidth 1-6*(dynamicHeight+dynamicGap) 1-dynamicTextWidth dynamicHeight], ...
    'Style','popup','String',analysisTypeString,'FontSize',fontSizeSmall);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Timing panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timingHeight = 0.1; timingTextWidth = 0.5; timingBoxWidth = 0.25;
hTimingPanel = uipanel('Title','Timing','fontSize', fontSizeLarge, ...
    'Unit','Normalized','Position',[timingStartPos panelStartHeight timingPanelWidth panelHeight]);

signalRange = [timeVals(1) timeVals(end)];
fftRange = [0 250];
baseline = [-0.2 0];
stimPeriod = [0 0.2];

% Signal Range
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','Parameter','FontSize',fontSizeMedium);

uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[timingTextWidth 1-timingHeight timingBoxWidth timingHeight], ...
    'Style','text','String','Min','FontSize',fontSizeMedium);

uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[timingTextWidth+timingBoxWidth 1-timingHeight timingBoxWidth timingHeight], ...
    'Style','text','String','Max','FontSize',fontSizeMedium);

% Stim Range
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-3*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','Stim Range (s)','FontSize',fontSizeSmall);
hStimMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth 1-3*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(signalRange(1)),'FontSize',fontSizeSmall);
hStimMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth 1-3*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(signalRange(2)),'FontSize',fontSizeSmall);

% FFT Range
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-5*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','FFT Range (Hz)','FontSize',fontSizeSmall);
hFFTMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth 1-5*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(fftRange(1)),'FontSize',fontSizeSmall);
hFFTMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth 1-5*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(fftRange(2)),'FontSize',fontSizeSmall);

% Baseline
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-6*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','Basline (s)','FontSize',fontSizeSmall);
hBaselineMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth 1-6*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(baseline(1)),'FontSize',fontSizeSmall);
hBaselineMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth 1-6*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(baseline(2)),'FontSize',fontSizeSmall);

% Stim Period
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-7*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','Stim period (s)','FontSize',fontSizeSmall);
hStimPeriodMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth 1-7*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(stimPeriod(1)),'FontSize',fontSizeSmall);
hStimPeriodMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth 1-7*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String',num2str(stimPeriod(2)),'FontSize',fontSizeSmall);

% Y Range
uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'Position',[0 1-8*timingHeight timingTextWidth timingHeight], ...
    'Style','text','String','Y Range','FontSize',fontSizeSmall);
hYMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth 1-8*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String','0','FontSize',fontSizeSmall);
hYMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[timingTextWidth+timingBoxWidth 1-8*timingHeight timingBoxWidth timingHeight], ...
    'Style','edit','String','1','FontSize',fontSizeSmall);

% include a ERP computation range- to make ERP calculation independant of TF plots
% also to calculate ERP for very small time ranges
if strcmp(gridType,'EEG')
    ERPPeriod=[0 0.250];
    uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
        'Position',[0 1-9*timingHeight timingTextWidth timingHeight], ...
        'Style','text','String','ERP period (s)','FontSize',fontSizeSmall);
    hERPPeriodMin = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
        'BackgroundColor', backgroundColor, ...
        'Position',[timingTextWidth 1-9*timingHeight timingBoxWidth timingHeight], ...
        'Style','edit','String',num2str(ERPPeriod(1)),'FontSize',fontSizeSmall);
    hERPPeriodMax = uicontrol('Parent',hTimingPanel,'Unit','Normalized', ...
        'BackgroundColor', backgroundColor, ...
        'Position',[timingTextWidth+timingBoxWidth 1-9*timingHeight timingBoxWidth timingHeight], ...
        'Style','edit','String',num2str(ERPPeriod(2)),'FontSize',fontSizeSmall);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotOptionsHeight = 0.1;
hPlotOptionsPanel = uipanel('Title','Plotting Options','fontSize', fontSizeLarge, ...
    'Unit','Normalized','Position',[plotOptionsStartPos panelStartHeight plotOptionsPanelWidth panelHeight]);

% Button for Plotting
[colorString, colorNames] = getColorString;
uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 1-plotOptionsHeight 0.6 plotOptionsHeight], ...
    'Style','text','String','Color','FontSize',fontSizeSmall);
hChooseColor = uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'BackgroundColor', backgroundColor, ...
    'Position',[0.6 1-plotOptionsHeight 0.4 plotOptionsHeight], ...
    'Style','popup','String',colorString,'FontSize',fontSizeSmall);

uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 4*plotOptionsHeight 1 plotOptionsHeight], ...
    'Style','pushbutton','String','cla','FontSize',fontSizeMedium, ...
    'Callback',{@cla_Callback});

hHoldOn = uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 3*plotOptionsHeight 1 plotOptionsHeight], ...
    'Style','togglebutton','String','hold on','FontSize',fontSizeMedium, ...
    'Callback',{@holdOn_Callback});

uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 2*plotOptionsHeight 1 plotOptionsHeight], ...
    'Style','pushbutton','String','rescale Y','FontSize',fontSizeMedium, ...
    'Callback',{@rescaleY_Callback});

uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 plotOptionsHeight 1 plotOptionsHeight], ...
    'Style','pushbutton','String','rescale X','FontSize',fontSizeMedium, ...
    'Callback',{@rescaleData_Callback});

uicontrol('Parent',hPlotOptionsPanel,'Unit','Normalized', ...
    'Position',[0 0 1 plotOptionsHeight], ...
    'Style','pushbutton','String','plot','FontSize',fontSizeMedium, ...
    'Callback',{@plotData_Callback});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% EEG Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the layout is changed only if the gridType is EEG;
if strcmpi(gridType,'EEG')
    eegHeight = 0.08; eegGap=0.02; eegTextWidth = 0.6;
    heegPanel = uipanel('Title','EEG Plots','fontSize', fontSizeLarge, ...
        'Unit','Normalized','Position',[eegStartPos panelStartHeight eegPanelWidth panelHeight]);
    
    %  To choose Channel for ERP and pooled ERP (if pooled elecs=0) and  TF Plots
    analogChannelStringList = getStringFromValues(analogChannelsStored,1); % using the function used in displySingleChannelSRC
    uicontrol('Parent',heegPanel,'Unit','Normalized', ...
        'Position',[0 1-(eegHeight+eegGap) eegTextWidth eegHeight],...
        'Style','text','String','Analog Channel','FontSize',fontSizeSmall);
    hAnalogChannel = uicontrol('Parent',heegPanel,'Unit','Normalized', ...
        'BackgroundColor', backgroundColor, 'Position', ...
        [eegTextWidth 1-(eegHeight+eegGap) 1-eegTextWidth eegHeight], ...
        'Style','popup','String',analogChannelStringList,'FontSize',fontSizeSmall);
    
    %Referencing- whether single,bipolar,average,with a specific channel
    uicontrol('Parent',heegPanel,'Unit','Normalized', ...
        'Position',[0 1-2*(eegHeight+eegGap) eegTextWidth eegHeight],...
        'Style','text','String','Ref Channel','FontSize',fontSizeSmall);
    hRefChannel = uicontrol('Parent',heegPanel,'Unit','Normalized', ...
        'BackgroundColor', backgroundColor, 'Position', ...
        [eegTextWidth 1-2*(eegHeight+eegGap) 1-eegTextWidth eegHeight], ...
        'Style','popup','String',(['Single Wire|Hemisphere|Average|Bipolar|' analogChannelStringList]),'FontSize',fontSizeSmall,'Callback',{@resetRef_Callback});
    
    % Plot and clear all buttons
    uicontrol('Parent',heegPanel,'Unit','Normalized', ...
        'Position',[0 1-4*(eegHeight+eegGap) eegTextWidth+0.4 eegHeight],...
        'Style','pushbutton','String','Plot ERP-topoplot','FontSize',fontSizeSmall,'Callback',{@plot_Callback});
    
    uicontrol('Parent',heegPanel,'Unit','Normalized', ...
        'Position',[0 1-5*(eegHeight+eegGap) eegTextWidth+0.4 eegHeight],...
        'Style','pushbutton','String','Clear all','FontSize',fontSizeSmall,'Callback',{@cla_eegplots_Callback});
    
    % Elecs for pooling: ERP
    % Plot button for electrode pooling
    uicontrol('Parent',heegPanel,'Unit','Normalized', ...
        'Position',[0 1-7*eegHeight eegTextWidth+0.4 eegHeight],...
        'Style','pushbutton','String','Pool elecs-ERP/TF ','FontSize',fontSizeSmall,'Callback',{@plotPoolElecs_Callback});
    hERPElecPool = uicontrol('Parent',heegPanel,'Unit','Normalized', ...
        'BackgroundColor', backgroundColor,...
        'Position',[ 0 1-8*eegHeight eegTextWidth+0.4 eegHeight],...
        'Style','edit','FontSize',fontSizeSmall);
    
    % Topoplots for power spectrum
    % topolots for the entire scalp -requires two inputs- the frequency band & type of data (single/bipolar/average/hemisphere)
    % Plot button for power spectrum topolot
    alphaRange=[8 13];
    freqRange=[30 80];
    uicontrol('Parent',heegPanel,'Unit','Normalized', ...
        'Position',[0 1-9*eegHeight (eegTextWidth+0.4)/2 eegHeight],...
        'Style','pushbutton','String','Plot TF-topoplot ','FontSize',fontSizeSmall,'Callback',{@plotTFtopo_Callback});
    uicontrol('Parent',heegPanel,'Unit','Normalized', ...
        'Position',[eegTextWidth 1-9*eegHeight (eegTextWidth+0.2)/2 eegHeight],...
        'Style','pushbutton','String','Plot TF ','FontSize',fontSizeSmall,'Callback',{@plotTF_Callback});
    uicontrol('Parent',heegPanel,'Unit','Normalized', ...
        'Position',[0 1-11*eegHeight eegTextWidth eegHeight], ...
        'Style','text','String','Alpha range (Hz)','FontSize',fontSizeSmall);
    halphaRangeMin = uicontrol('Parent',heegPanel,'Unit','Normalized', ...
        'BackgroundColor', backgroundColor, ...
        'Position',[eegTextWidth 1-11*eegHeight eegGap+0.15 eegHeight], ...
        'Style','edit','String',num2str(alphaRange(1)),'FontSize',fontSizeSmall);
    halphaRangeMax = uicontrol('Parent',heegPanel,'Unit','Normalized', ...
        'BackgroundColor', backgroundColor, ...
        'Position',[eegTextWidth+eegGap+0.2 1-11*eegHeight eegGap+0.15 eegHeight], ...
        'Style','edit','String',num2str(alphaRange(2)),'FontSize',fontSizeSmall);
    
    uicontrol('Parent',heegPanel,'Unit','Normalized', ...
        'Position',[0 1-12*eegHeight eegTextWidth eegHeight], ...
        'Style','text','String','Freq range (Hz)','FontSize',fontSizeSmall);
    hfreqRangeMin = uicontrol('Parent',heegPanel,'Unit','Normalized', ...
        'BackgroundColor', backgroundColor, ...
        'Position',[eegTextWidth 1-12*eegHeight eegGap+0.15 eegHeight], ...
        'Style','edit','String',num2str(freqRange(1)),'FontSize',fontSizeSmall);
    hfreqRangeMax = uicontrol('Parent',heegPanel,'Unit','Normalized', ...
        'BackgroundColor', backgroundColor, ...
        'Position',[eegTextWidth+eegGap+0.2 1-12*eegHeight eegGap+0.15 eegHeight], ...
        'Style','edit','String',num2str(freqRange(2)),'FontSize',fontSizeSmall);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get main plot and message handles
if strcmpi(gridType,'ECoG')
    plotHandles = getPlotHandles(8,10,[0.05 0.05 0.9 0.625]);
elseif strcmpi(gridType,'EEG')
    % Panel gets divided into two to house the topoplots
    % additional variables- defining the topoplot locations on the right side of the panel
    plotHandles = getPlotHandles(9,11,[0.05 0.05 0.9*panel_factor 0.625]);
    capStartPos=0.515;
    capStartHeight=0.05;
    capBoxWidth=0.45;
    capBoxHeight=0.55;
    capGap = 0.04;
    % topoplot locations -ERP plots
    electrodeCapPosERP=[capStartPos capStartHeight+capBoxHeight/2+capGap capBoxWidth/2 capBoxHeight/2];
    capERPHandle = subplot('Position',electrodeCapPosERP); axis off;
    electrodeCapPosERPRef=[capStartPos+(capBoxWidth/2)+capGap capStartHeight+capBoxHeight/2+capGap capBoxWidth/2 capBoxHeight/2];
    capERPRefHandle = subplot('Position',electrodeCapPosERPRef); axis off;
    % topoplot locations -  TF plots
    electrodeCapPosTFalpha=[capStartPos capStartHeight capBoxWidth/2 capBoxHeight/2];
    capTFalphaHandle=subplot('Position',electrodeCapPosTFalpha); axis off;
    electrodeCapPosTF=[capStartPos+(capBoxWidth/2)+capGap capStartHeight capBoxWidth/2 capBoxHeight/2];
    capTFHandle=subplot('Position',electrodeCapPosTF); axis off;
else
    plotHandles = getPlotHandles;
end

hMessage = uicontrol('Unit','Normalized','Position',[0 0.975 1 0.025],...
    'Style','text','String','DisplayAllChannelsSRC','FontSize',fontSizeLarge);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% functions
    function plotData_Callback(~,~)
        c=get(hContrast,'val');
        t=get(hTemporalFreq,'val');
        e=get(hEOTCode,'val');
        a=get(hAttendLoc,'val');
        s=get(hStimType,'val');
        analysisType = get(hAnalysisType,'val');
        plotColor = colorNames(get(hChooseColor,'val'));
        
        blRange = [str2double(get(hBaselineMin,'String')) str2double(get(hBaselineMax,'String'))];
        stRange = [str2double(get(hStimPeriodMin,'String')) str2double(get(hStimPeriodMax,'String'))];
        
        goodPos = parameterCombinations{c,t,e,a,s};
        set(hMessage,'String',[num2str(length(goodPos)) ' stimuli found' ]);
        
        if ~isempty(goodPos)
            if analysisType == 2
                holdOnState = get(hHoldOn,'val');
                channelsStored = neuralChannelsStored;
                [baselineFiringRate,stimulusFiringRate] = plotSpikeData(plotHandles,channelsStored,goodPos,folderSpikes,...
                    timeVals,plotColor,SourceUnitID,holdOnState,blRange,stRange,gridType);
                
                responsiveElectrodes = (channelsStored(stimulusFiringRate>=5));
                rareElectrodes = setdiff((channelsStored(stimulusFiringRate>=2)),responsiveElectrodes);
                inhibitedElectrodes = (channelsStored(intersect(find(baselineFiringRate>stimulusFiringRate),find(baselineFiringRate>2))));
                disp(['rareUnits: ' num2str(rareElectrodes)]);
                disp(['responsive: ' num2str(responsiveElectrodes)]);
                disp(['inhibited : ' num2str(inhibitedElectrodes)]);
                
                showElectrodeLocations(electrodeGridPos,responsiveElectrodes,'b',hElectrodes,0,0,gridType,subjectName);
                showElectrodeLocations(electrodeGridPos,inhibitedElectrodes,'g',hElectrodes,1,0,gridType,subjectName);
            else
                channelsStored = analogChannelsStored;
                plotLFPData(plotHandles,channelsStored,goodPos,folderLFP,...
                    analysisType,timeVals,plotColor,blRange,stRange,gridType,subjectName);
            end
            
            if analysisType<=2 % ERP or spikes
                xRange = [str2double(get(hStimMin,'String')) str2double(get(hStimMax,'String'))];
            else
                xRange = [str2double(get(hFFTMin,'String')) str2double(get(hFFTMax,'String'))];
            end
            
            yRange = getYLims(plotHandles,channelsStored,gridType,subjectName);
            rescaleData(plotHandles,channelsStored,[xRange yRange],gridType,subjectName);
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function rescaleY_Callback(~,~)
        
        analysisType = get(hAnalysisType,'val');
        
        if analysisType == 2
            channelsStored = neuralChannelsStored;
        else
            channelsStored = analogChannelsStored;
        end
        
        if analysisType<=2 % ERP or spikes
            xRange = [str2double(get(hStimMin,'String')) str2double(get(hStimMax,'String'))];
        else
            xRange = [str2double(get(hFFTMin,'String')) str2double(get(hFFTMax,'String'))];
        end
        
        yRange = [str2double(get(hYMin,'String')) str2double(get(hYMax,'String'))];
        rescaleData(plotHandles,channelsStored,[xRange yRange],gridType,subjectName);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function rescaleData_Callback(~,~)
        
        analysisType = get(hAnalysisType,'val');
        
        if analysisType == 2
            channelsStored = neuralChannelsStored;
        else
            channelsStored = analogChannelsStored;
        end
        
        if analysisType<=2 % ERP or spikes
            xRange = [str2double(get(hStimMin,'String')) str2double(get(hStimMax,'String'))];
        else
            xRange = [str2double(get(hFFTMin,'String')) str2double(get(hFFTMax,'String'))];
        end
        
        yRange = getYLims(plotHandles,channelsStored,gridType,subjectName);
        rescaleData(plotHandles,channelsStored,[xRange yRange],gridType,subjectName);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function holdOn_Callback(source,~)
        holdOnState = get(source,'Value');
        
        [numRow,numCol] = size(plotHandles);
        
        if holdOnState
            for i=1:numRow
                for j=1:numCol
                    set(plotHandles(i,j),'Nextplot','add');
                end
            end
        else
            for i=1:numRow
                for j=1:numCol
                    set(plotHandles(i,j),'Nextplot','replace');
                end
            end
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function cla_Callback(~,~)
        [numRow,numCol] = size(plotHandles);
        for i=1:numRow
            for j=1:numCol
                cla(plotHandles(i,j));
            end
        end
        
    end
%__________________________EEG Specific functions__________________________
%this function plots the ERP topoplots for single and chosen reference scheme
    function plot_Callback(~,~)
        % Intitialise
        c=get(hContrast,'val');
        t=get(hTemporalFreq,'val');
        e=get(hEOTCode,'val');
        a=get(hAttendLoc,'val');
        s=get(hStimType,'val');
        refChanIndex = get(hRefChannel,'val'); % reference channel-single/bipolar/avg
        
        % ERP  Variables
        ERPMin = str2double(get(hERPPeriodMin,'string'));
        ERPMax = str2double(get(hERPPeriodMax,'string'));
        tERP = (timeVals>=ERPMin) & (timeVals<=ERPMax);
        blPeriod = (timeVals<=0);
        if ~blPeriod
            blPeriod = tERP;
        end
        % Get Data
        [plotData] = getDataSRC(c,t,e,a,s,folderName,folderLFP,analogChannelsStored);
        
        % plot data
        refChanIndex=1; % by default one topoplot is single wire referenced
        [~,chanlocs] = bipolarRef(refChanIndex,plotData,folderMontage);
        erpTopo = getERPData(plotData,tERP,blPeriod);
        subplot(capERPHandle); topoplot(erpTopo,chanlocs,'electrodes','numbers','style','both','drawaxis','off');
        colorbar; title(['RMS ERP: ' num2str(ERPMin) ' to ' num2str(ERPMax) ' s' ' (SingleWire)']);
        clear Data;
        
        %Plotting topoplot based on the reference scheme-right side plot
        refChanIndex = get(hRefChannel,'val');
        if refChanIndex~=1 % plot the rereferenced topoplot only for bipolar/average/hemisphere/other chan reference
            [Data,chanlocs] = bipolarRef(refChanIndex,plotData,folderMontage);
            
            if refChanIndex ~= 4 % anything apart from bipolar
                erpTopoRef = getERPData(Data,tERP,blPeriod);
                subplot(capERPRefHandle); topoplot(erpTopoRef,chanlocs,'electrodes','numbers','style','both','drawaxis','off');
                colorbar; title(['Rereferenced']);
            else % bipolar reference
                erpTopoRef = (rms(squeeze(mean((Data(:,:,tERP)-repmat(mean(Data(:,:,blPeriod),3),1,1,size(Data(:,:,tERP),3))),2))'))'...
                    -(rms(squeeze(mean((Data(:,:,blPeriod)-repmat(mean(Data(:,:,blPeriod),3),1,1,size(Data(:,:,blPeriod),3))),2))'))';
                subplot(capERPRefHandle); topoplot(erpTopoRef,chanlocs,'electrodes','numbers','style','both','drawaxis','off','nosedir','-Y');
                colorbar; title('Re-referenced');
            end
        else
            disp('Rereferencing was not chosen');
        end
    end

% function for Plotting TF topoplot- alpha and specified freq plots____
    function plotTFtopo_Callback(~,~)
        % Intitialise
        a=get(hAttendLoc,'val');
        e=get(hEOTCode,'val');
        s=get(hStimType,'val');
        c=get(hContrast,'val');
        t=get(hTemporalFreq,'val');
        refChanIndex = get(hRefChannel,'val');
        
        %get the baseline and stimulus epoch ranges
        BLMin = str2double(get(hBaselineMin,'string'));
        BLMax = str2double(get(hBaselineMax,'string'));
        STMin = str2double(get(hStimPeriodMin,'string'));
        STMax = str2double(get(hStimPeriodMax,'string'));
        
        %get the ranges for the frequency bands
        alphaBandMin=str2double(get(halphaRangeMin,'string'));
        alphaBandMax=str2double(get(halphaRangeMax,'string'));
        fBandMin = str2double(get(hfreqRangeMin,'string'));
        fBandMax = str2double(get(hfreqRangeMax,'string'));
        
        % Get Data
        [plotData] = getDataSRC(c,t,e,a,s,folderName,folderLFP,analogChannelsStored);
        
        % checking the state of Data re-referencing option chosen
        [Data,chanlocs] = bipolarRef(refChanIndex,plotData,folderMontage);
        
        % Set MT parameters
        fMax=100;
        params.tapers   = [1 1];
        params.pad      = -1;
        params.Fs       = 1/abs(timeVals(1)-timeVals(2));
        params.fpass    = [0 fMax];
        params.trialave = 1;
        movingWin = [0.25 0.025];
        
        %get the spectral estimate for each electrode
        for i=1:size(Data,1)
            dataTF=squeeze(Data(i,:,:));
            [~,dS1,t2,f2] = getSTFT(dataTF,movingWin,params,timeVals,BLMin,BLMax);
            dSPower(i,:,:) = dS1;
        end
        
        % for any specified frequency band
        meanTF=[];
        for i=1:size(Data,1)
            tST =  (t2>=STMin) & (t2<=STMax);
            fPL=   (f2>=fBandMin) & (f2<=fBandMax);
            dS1 = squeeze(dSPower(i,:,:));
            S1ST = dS1(tST,fPL);
            meanTF(i)=mean(mean(S1ST,2));
        end
        
        % for alpha band
        meanAlphaTF=[];
        for j=1:size(Data,1)
            tST1 =  (t2>=STMin) & (t2<=STMax);
            fPL1=   (f2>=alphaBandMin) & (f2<=alphaBandMax);
            dS11 = squeeze(dSPower(j,:,:));
            S2ST = dS11(tST1,fPL1);
            meanAlphaTF(j)=mean(mean(S2ST,2));
        end
        
        %Plotting the topoplots
        if refChanIndex ~=4
            subplot(capTFalphaHandle);cla(gca,'reset');axis off;
            subplot(capTFHandle);cla(gca,'reset');axis off;
            subplot(capTFalphaHandle); topoplot(meanAlphaTF,chanlocs,'electrodes','numbers','style','both','drawaxis','off'); colorbar; title('Change in alpha power');
            subplot(capTFHandle); topoplot(meanTF,chanlocs,'electrodes','numbers','style','both','drawaxis','off'); colorbar; title('Change in chosen band power');
        else
            subplot(capTFalphaHandle);cla(gca,'reset');axis off;
            subplot(capTFHandle);cla(gca,'reset');axis off;
            subplot(capTFalphaHandle); topoplot(meanAlphaTF,chanlocs,'electrodes','numbers','style','both','drawaxis','off','nosedir','-Y'); colorbar; title('Change in alpha power');
            subplot(capTFHandle); topoplot(meanTF,chanlocs,'electrodes','numbers','style','both','drawaxis','off','nosedir','-Y'); colorbar; title('Change in  chosen band power');
        end
    end

% function for plotting ERP waveforms from pooled electrodes
    function plotPoolElecs_Callback(~,~)
        a=get(hAttendLoc,'val');
        e=get(hEOTCode,'val');
        s=get(hStimType,'val');
        c=get(hContrast,'val');
        t=get(hTemporalFreq,'val');
        refChanIndex = get(hRefChannel,'val');
        epochMin = timeVals(1);
        epochMax = timeVals(end);
        ERPMin = str2double(get(hStimMin,'string'));
        ERPMax = str2double(get(hStimMax,'string'));
        tERP = (timeVals>=epochMin) & (timeVals<=epochMax);
        blPeriod = (timeVals<=0);
        if ~blPeriod
            blPeriod = tERP;
        end
        % get  electrodes list
        allElecData = [];
        ERPPoolElecs = str2num(get(hERPElecPool,'string'));
        if isempty(ERPPoolElecs)
            disp('Electrodes to pool not specified...');
            ERPPoolElecs = get(hAnalogChannel,'val');
        else
            disp(['Computing pooled ERP for electrodes : ',num2str(ERPPoolElecs)]);
        end
        % Get Data
        try
            if refChanIndex==1
                [plotData] = getDataSRC(c,t,e,a,s,folderName,folderLFP,ERPPoolElecs);
                Data=plotData;
                for iP = 1:length(ERPPoolElecs)
                    elecData = squeeze(Data(iP,:,:));
                    allElecData = [allElecData;elecData];
                end
            else % data re referencing
                %Data for all the electrodes : this is to ensure the re referencing if any
                [plotData] = getDataSRC(c,t,e,a,s,folderName,folderLFP,analogChannelsStored);
                
                % checking the state of Data re-referencing option chosen
                [Data]=bipolarRef(refChanIndex,plotData,folderMontage);
                
                % pick out the data for the set of electrodes to be pooled
                for iP = 1:length(ERPPoolElecs)
                    elecData = squeeze(Data(ERPPoolElecs(iP),:,:));
                    allElecData = [allElecData;elecData];
                end
            end
            clear erp erpData erpDataBLCor;
            erpData = allElecData(:,tERP);
            erpDataBLCor = (erpData - repmat(mean(erpData(:,blPeriod),2),1,size(erpData,2))); % Correct for DC Shift (baseline correction)
            erp = mean(erpDataBLCor,1);
            
            % pooled ERP for the specified electrodes-plots in a new figure
            figure('name','Pooled ERP plot');
            plot(timeVals(tERP),erp);
            xlim([epochMin epochMax]);
            [refType]=getRefschemeName(refChanIndex);
            text(0.1,0.9,['ERP; n = ' num2str(size(erpData,1))],'unit','normalized','fontsize',9);
            text(0.3,0.9,['Electrodes pooled: ' num2str(ERPPoolElecs)],'unit','normalized','fontsize',9);
            text(0.6,0.9,['Reference Scheme: ',refType],'unit','normalized','fontsize',9);
        catch
            disp('Please enter electrodes according to the reference option chosed; This electrode no does not exist for this ref scheme');
        end
    end

% this function is used to plot the TF plots for the chosen electrode/group of electrodes
    function plotTF_Callback(~,~)
        % Intitialise
        a=get(hAttendLoc,'val');
        e=get(hEOTCode,'val');
        s=get(hStimType,'val');
        c=get(hContrast,'val');
        t=get(hTemporalFreq,'val');
        refChanIndex = get(hRefChannel,'val');
        
        %get the baseline and stimulus epoch ranges
        BLMin = str2double(get(hBaselineMin,'string'));
        BLMax = str2double(get(hBaselineMax,'string'));
        STMin = str2double(get(hStimPeriodMin,'string'));
        STMax = str2double(get(hStimPeriodMax,'string'));
        
        %channel selection : single/pooled
        ERPPoolElecs = str2num(get(hERPElecPool,'string'));
        if isempty(ERPPoolElecs)
            disp('Electrodes to pool not specified...');
            ERPPoolElecs = get(hAnalogChannel,'val');
        else
            disp(['Computing pooled Time Frequency for electrodes : ',num2str(ERPPoolElecs)]);
        end
        if refChanIndex==1
            [plotData] = getDataSRC(c,t,e,a,s,folderName,folderLFP,ERPPoolElecs);
            Data=plotData;
        else
            % Get Data for all the electrodes : this is to ensure the re
            [plotData] = getDataSRC(c,t,e,a,s,folderName,folderLFP,analogChannelsStored);
            % checking the state of Data re-referencing option chosen
            [DataAll] = bipolarRef(refChanIndex,plotData,folderMontage);
            % pick out the data for the set of electrodes to be pooled
            for e=1:length(ERPPoolElecs)
                Data(e,:,:)=DataAll(ERPPoolElecs(e),:,:);
            end
        end
        % Set MT parameters
        fMax=100;
        params.tapers   = [1 1];
        params.pad      = -1;
        params.Fs       = 1/abs(timeVals(1)-timeVals(2));
        params.fpass    = [0 fMax];
        params.trialave = 1;
        movingWin = [0.25 0.025];
        cLims=[-5 5];
        %get the spectral estimate for each electrode
        for i=1:size(Data,1)
            dataTF=squeeze(Data(i,:,:));
            [~,dS1,t2,f2] = getSTFT(dataTF,movingWin,params,timeVals,BLMin,BLMax);
            dSPower(i,:,:) = dS1;
        end
        disp(['Computing pooled Time Frequency for electrodes : ',num2str(ERPPoolElecs)]);
        tfData=squeeze(mean(dSPower,1)); % averaging the log spectrum across different electrodes
        
        %plotting the TF freq plot for that electrode/pooled electrode
        figure('name',['TF plot for electrode :',num2str(ERPPoolElecs)]);
        pcolor(t2,f2,tfData');
        colormap('jet')
        shading('interp');
        caxis(cLims);
        axis([-0.2 0.7 0 fMax]);
        [refType]=getRefschemeName(refChanIndex);
        title(refType);   
    end

% data rereference function- gets the scheme for referencing and changing
% the list of channels accordingly.
    function resetRef_Callback(~,~)
        refChanIndex = get(hRefChannel,'val');
        if refChanIndex == 4 % if the scheme is bipolar , form the bipolar sets
            chanloc = loadChanLocs(capType,refChanIndex,folderMontage);
            bipAnalogChannelList = getBipAnalogChannelList(chanloc);
            set(hAnalogChannel,'String',bipAnalogChannelList);
        else
            analogChannelStringList = getStringFromValues(analogChannelsStored,1);
            set(hAnalogChannel,'String',analogChannelStringList);
        end
    end

% function for clearing the EEG related plots
    function cla_eegplots_Callback(~,~)
        subplot(capERPHandle); cla(gca,'reset'); axis off;
        subplot(capERPRefHandle); cla(gca,'reset'); axis off;
        subplot(capTFalphaHandle);cla(gca,'reset');axis off;
        subplot(capTFHandle);cla(gca,'reset');axis off;
        hElectrodes = showElectrodeLocations(electrodeGridPos,[],[],[],1,0,gridType);
    end
% this function returns the electrode wise data according to the reference
% option chosen; right now common bad trials are used.
    function [Data,chanlocs] = bipolarRef(refChanIndex,plotData,folderMontage)
        
        if refChanIndex == 1 % Single wire referencing
            disp(' working on single wire reference data..');
            Data = plotData;
            chanlocs = loadChanLocs(capType,refChanIndex,folderMontage);%load the channel locations based on the captype
        elseif refChanIndex == 2 % hemisphere referencing
            [chanlocs,hemBipolarLocs] = loadChanLocs(capType,refChanIndex,folderMontage);
            disp(' working on hemisphere referenced data..');
            for iH = 1:size(plotData,1)
                Data(iH,:,:) = plotData(hemBipolarLocs(iH,1),:,:) - plotData(hemBipolarLocs(iH,2),:,:);
            end
            
        elseif refChanIndex == 3 % average referencing
            chanlocs = loadChanLocs(capType,refChanIndex,folderMontage);
            disp(' working on average referenced data..');
            aveData = mean(plotData,1);
            for iH = 1:size(plotData,1)
                Data(iH,:,:) = plotData(iH,:,:) - aveData;
            end
        elseif refChanIndex == 4 % bipolar referencing
            [chanlocs,~,bipolarLocs] = loadChanLocs(capType,refChanIndex,folderMontage);
            disp(' working on bipolar referenced data..');
            maxChanKnown = 96; % default set by MD while creating bipolar montage; this might be different for different montages!!
            for iH = 1:size(bipolarLocs,1)
                chan1 = bipolarLocs(iH,1);
                chan2 = bipolarLocs(iH,2);
                if chan1<(maxChanKnown+1)
                    unipolarChan1 = plotData(chan1,:,:);
                else
                    unipolarChan1 = Data(chan1,:,:);
                end
                if chan2<(maxChanKnown+1)
                    unipolarChan2 = plotData(chan2,:,:);
                else
                    unipolarChan2 = Data(chan2,:,:);
                end
                Data(iH,:,:) = unipolarChan1 - unipolarChan2;
            end
        else
            refChanIndex = refChanIndex-4;
            chanlocs = loadChanLocs(capType,refChanIndex,folderMontage);
            disp(['Rereferencing data wrto electrode :' num2str(refChanIndex)]);
            for iR = 1:size(plotData,1)
                Data(iR,:,:) = plotData(iR,:,:) - plotData(refChanIndex,:,:);
            end
        end
    end

% this function gets the channel list for bipolar electrode reference scheme
    function outString = getBipAnalogChannelList(chanloc)
        chanSize = size(chanloc,2);
        outString='';
        for iC = 1:chanSize
            outArray{iC} = ['elec' num2str(iC)];
            outString = cat(2,outString,[outArray{iC} '|']);
        end
    end

% This function is used to get the  time averaged ERP data for all the channels
    function erp = getERPData(Data,tERP,blPeriod)
        erp = zeros(size(Data,1),1);
        for iE = 1:size(Data,1)
            for goodPos=1:size(Data,2) % for every trial compute the baseline factor
                dataBL = mean(mean(squeeze(Data(iE,goodPos,blPeriod)),1),2);
                dataERP(goodPos,:)=squeeze(Data(iE,goodPos,tERP))-dataBL;
            end
            erp(iE,1) = mean(mean(dataERP,2),1); % mean ERP value for each trial (averaged across time first) and mean across trials
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main function that plots the data
function plotLFPData(plotHandles, channelsStored, goodPos, ...
    folderData, analysisType, timeVals, plotColor,blRange,stRange,gridType,subjectName)

if isempty(goodPos)
    disp('No entries for this combination..')
else
    if analysisType>2 % FFT
        Fs = round(1/(timeVals(2)-timeVals(1)));
        blPos = find(timeVals>=blRange(1),1)+ (1:diff(blRange)*Fs);
        stPos = find(timeVals>=stRange(1),1)+ (1:diff(stRange)*Fs);
        
        xsBL = 0:1/(diff(blRange)):Fs-1/(diff(blRange));
        xsST = 0:1/(diff(stRange)):Fs-1/(diff(stRange));
    end
    
    for i=1:length(channelsStored)
        disp(i)
        channelNum = channelsStored(i);
        
        % get position
        [row,column] = electrodePositionOnGrid(channelNum,gridType,subjectName);
        
        if analysisType == 1        % compute ERP
            clear signal analogData
            load(fullfile(folderData ,['elec' num2str(channelNum)]));
            erp = mean(analogData(goodPos,:),1); %#ok<*NODEF>
            
            %Plot
            plot(plotHandles(row,column),timeVals,erp,'color',plotColor);
            
        elseif analysisType == 2    % compute Firing rates
            disp('Use plotSpikeData instead of plotLFPData...');
        else
            clear signal analogData
            load(fullfile(folderData,['elec' num2str(channelNum)]));
            
            fftBL = abs(fft(analogData(goodPos,blPos),[],2));
            fftST = abs(fft(analogData(goodPos,stPos),[],2));
            
            if analysisType == 3
                plot(plotHandles(row,column),xsBL,log10(mean(fftBL)),'g');
                hold(plotHandles(row,column),'on');
                plot(plotHandles(row,column),xsST,log10(mean(fftST)),plotColor);
            end
            
            if analysisType == 4
                if xsBL == xsST %#ok<BDSCI>
                    plot(plotHandles(row,column),xsBL,log10(mean(fftST))-log10(mean(fftBL)),'color',plotColor);
                else
                    disp('Choose same baseline and stimulus periods..');
                end
            end
        end
    end
end
end
function [baselineFiringRate,stimulusFiringRate] = plotSpikeData(plotHandles,channelsStored,goodPos, ...
    folderData, timeVals, plotColor, SourceUnitID,holdOnState,blRange,stRange,gridType)

unitColors = ['r','m','y','c','g'];
binWidthMS = 10;

if isempty(goodPos)
    disp('No entries for this combination..')
else
    numChannels = length(channelsStored);
    baselineFiringRate = zeros(1,numChannels);
    stimulusFiringRate = zeros(1,numChannels);
    
    for i=1:numChannels
        %disp(i)
        channelNum = channelsStored(i);
        
        % get position
        [row,column] = electrodePositionOnGrid(channelNum,gridType,subjectName);
        
        clear neuralInfo spikeData
        load([folderData 'elec' num2str(channelNum) '_SID' num2str(SourceUnitID(i))]);
        [psthVals,xs] = psth_SR(spikeData(goodPos),binWidthMS,timeVals(1),timeVals(end));
        
        % Compute the mean firing rates
        blPos = find(xs>=blRange(1),1)+ (1:(diff(blRange))/(binWidthMS/1000));
        stPos = find(xs>=stRange(1),1)+ (1:(diff(stRange))/(binWidthMS/1000));
        
        baselineFiringRate(i) = mean(psthVals(blPos));
        stimulusFiringRate(i) = mean(psthVals(stPos));
        
        if SourceUnitID(i)==0
            plot(plotHandles(row,column),xs,psthVals,'color',plotColor);
        elseif SourceUnitID(i)> 5
            disp('Only plotting upto 6 single units per electrode...')
        else
            plot(plotHandles(row,column),xs,smooth(psthVals),'color',unitColors(SourceUnitID(i)));
        end
        
        if (i<length(channelsStored))
            if channelsStored(i) == channelsStored(i+1)
                disp('hold on...')
                set(plotHandles(row,column),'Nextplot','add');
            else
                if ~holdOnState
                    set(plotHandles(row,column),'Nextplot','replace');
                end
            end
        end
    end
end
end
function yRange = getYLims(plotHandles,channelsStored,gridType,subjectName)

% Initialize
yMin = inf;
yMax = -inf;

for i=1:length(channelsStored)
    channelNum = channelsStored(i);
    % get position
    [row,column] = electrodePositionOnGrid(channelNum,gridType,subjectName);
    
    axis(plotHandles(row,column),'tight');
    tmpAxisVals = axis(plotHandles(row,column));
    if tmpAxisVals(3) < yMin
        yMin = tmpAxisVals(3);
    end
    if tmpAxisVals(4) > yMax
        yMax = tmpAxisVals(4);
    end
end
yRange = [yMin yMax];
end
function rescaleData(plotHandles,channelsStored,axisLims,gridType,subjectName)

[numRows,numCols] = size(plotHandles);
labelSize=12;
for i=1:length(channelsStored)
    channelNum = channelsStored(i);
    % get position
    [row,column] = electrodePositionOnGrid(channelNum,gridType,subjectName);
    
    axis(plotHandles(row,column),axisLims);
    if (row==numRows && rem(column,2)==1)
        if column~=1
            set(plotHandles(row,column),'YTickLabel',[],'fontSize',labelSize);
        end
    elseif (rem(row,2)==0 && column==1)
        set(plotHandles(row,column),'XTickLabel',[],'fontSize',labelSize);
    else
        set(plotHandles(row,column),'XTickLabel',[],'YTickLabel',[],'fontSize',labelSize);
    end
end

% Remove Labels on the four corners
set(plotHandles(1,1),'XTickLabel',[],'YTickLabel',[]);
set(plotHandles(1,numCols),'XTickLabel',[],'YTickLabel',[]);
set(plotHandles(numRows,1),'XTickLabel',[],'YTickLabel',[]);
set(plotHandles(numRows,numCols),'XTickLabel',[],'YTickLabel',[]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function contrastString = getContrastString(cIndexValsUnique,stimResults)
if isfield(stimResults,'contrast0PC')
    [cVals0Unique,cVals1Unique] = getValsFromIndex(cIndexValsUnique,stimResults,'contrast');
    if length(cVals0Unique)==1
        contrastString = [num2str(cVals0Unique) ',' num2str(cVals1Unique)];
    else
        contrastString = '';
        for i=1:length(cVals0Unique)
            contrastString = cat(2,contrastString,[num2str(cVals0Unique(i)) ',' num2str(cVals1Unique(i)) '|']);
        end
        contrastString = [contrastString 'all'];
    end
    
else % Work with indices
    if length(cIndexValsUnique)==1
        if cIndexValsUnique ==0
            contrastString = '0';
        else
            contrastString = num2str(100/2^(7-cIndexValsUnique));
        end
        
    else
        contrastString = '';
        for i=1:length(cIndexValsUnique)
            if cIndexValsUnique(i) == 0
                contrastString = cat(2,contrastString,[ '0|']); %#ok<*NBRAK>
            else
                contrastString = cat(2,contrastString,[num2str(100/2^(7-cIndexValsUnique(i))) '|']);
            end
        end
        contrastString = [contrastString 'all'];
    end
end
end
function temporalFreqString = getTemporalFreqString(tIndexValsUnique,stimResults)

if isfield(stimResults,'temporalFreq0Hz')
    [tVals0Unique,tVals1Unique] = getValsFromIndex(tIndexValsUnique,stimResults,'temporalFreq');
    if length(tIndexValsUnique)==1
        temporalFreqString = [num2str(tVals0Unique) ',' num2str(tVals1Unique)];
    else
        temporalFreqString = '';
        for i=1:length(tIndexValsUnique)
            temporalFreqString = cat(2,temporalFreqString,[num2str(tVals0Unique(i)) ',' num2str(tVals1Unique(i)) '|']);
        end
        temporalFreqString = [temporalFreqString 'all'];
    end
else
    if length(tIndexValsUnique)==1
        if tIndexValsUnique ==0
            temporalFreqString = '0';
        else
            temporalFreqString = num2str(min(50,80/2^(7-tIndexValsUnique)));
        end
        
    else
        temporalFreqString = '';
        for i=1:length(tIndexValsUnique)
            if tIndexValsUnique(i) == 0
                temporalFreqString = cat(2,temporalFreqString,[ '0|']);
            else
                temporalFreqString = cat(2,temporalFreqString,[num2str(min(50,80/2^(7-tIndexValsUnique(i)))) '|']);
            end
        end
        temporalFreqString = [temporalFreqString 'all'];
    end
end
end
function EOTCodeString = getEOTCodeString(eIndexValsUnique)

if length(eIndexValsUnique)==1
    if eIndexValsUnique == 0
        EOTCodeString = 'Correct';
    elseif eIndexValsUnique == 1
        EOTCodeString = 'Wrong';
    elseif eIndexValsUnique == 2
        EOTCodeString = 'Failed';
    elseif eIndexValsUnique == 3
        EOTCodeString = 'Broke';
    elseif eIndexValsUnique == 4
        EOTCodeString = 'Ignored';
    elseif eIndexValsUnique == 5
        EOTCodeString = 'False Alarm';
    elseif eIndexValsUnique == 6
        EOTCodeString = 'Distracted';
    elseif eIndexValsUnique == 7
        EOTCodeString = 'Force Quit';
    else
        disp('Unknown EOT Code');
    end
else
    EOTCodeString = '';
    for i=1:length(eIndexValsUnique)
        if eIndexValsUnique(i) == 0
            EOTCodeString = cat(2,EOTCodeString,[ 'Correct|']);
        elseif eIndexValsUnique(i) == 1
            EOTCodeString = cat(2,EOTCodeString,[ 'Wrong|']);
        elseif eIndexValsUnique(i) == 2
            EOTCodeString = cat(2,EOTCodeString,[ 'Failed|']);
        elseif eIndexValsUnique(i) == 3
            EOTCodeString = cat(2,EOTCodeString,[ 'Broke|']);
        elseif eIndexValsUnique(i) == 4
            EOTCodeString = cat(2,EOTCodeString,[ 'Ignored|']);
        elseif eIndexValsUnique(i) == 5
            EOTCodeString = cat(2,EOTCodeString,[ 'False Alarm|']);
        elseif eIndexValsUnique(i) == 6
            EOTCodeString = cat(2,EOTCodeString,[ 'Distracted|']);
        elseif eIndexValsUnique(i) == 7
            EOTCodeString = cat(2,EOTCodeString,[ 'Force Quit|']);
        else
            disp('Unknown EOT Code');
        end
    end
    EOTCodeString = [EOTCodeString 'all'];
end
end
function attendLocString = getAttendLocString(aIndexValsUnique)

if length(aIndexValsUnique)==1
    if aIndexValsUnique == 0
        attendLocString = '0 (right)';
    elseif aIndexValsUnique == 1
        attendLocString = '1 (left)';
    else
        disp('Unknown attended location');
    end
else
    attendLocString = '';
    for i=1:length(aIndexValsUnique)
        if aIndexValsUnique(i) == 0
            attendLocString = cat(2,attendLocString,[ '0 (right)|']);
        elseif aIndexValsUnique(i) == 1
            attendLocString = cat(2,attendLocString,[ '1 (left)|']);
        else
            disp('Unknown attended location');
        end
    end
    attendLocString = [attendLocString 'Both'];
end
end
function stimTypeString = getStimTypeString(sIndexValsUnique)

if length(sIndexValsUnique)==1
    if sIndexValsUnique == 0
        stimTypeString = 'Null';
    elseif sIndexValsUnique == 1
        stimTypeString = 'Correct';
    elseif sIndexValsUnique == 2
        stimTypeString = 'Target';
    elseif sIndexValsUnique == 3
        stimTypeString = 'FrontPad';
    elseif sIndexValsUnique == 4
        stimTypeString = 'BackPad';
    else
        disp('Unknown Stimulus Type');
    end
else
    stimTypeString = '';
    for i=1:length(sIndexValsUnique)
        if sIndexValsUnique(i) == 0
            stimTypeString = cat(2,stimTypeString,['Null|']);
        elseif sIndexValsUnique(i) == 1
            stimTypeString = cat(2,stimTypeString,['Correct|']);
        elseif sIndexValsUnique(i) == 2
            stimTypeString = cat(2,stimTypeString,['Target|']);
        elseif sIndexValsUnique(i) == 3
            stimTypeString = cat(2,stimTypeString,['FrontPad|']);
        elseif sIndexValsUnique(i) == 4
            stimTypeString = cat(2,stimTypeString,['BackPad|']);
        else
            disp('Unknown Stimulus Type');
        end
    end
    stimTypeString = [stimTypeString 'all'];
end

end
function [colorString, colorNames] = getColorString

colorNames = 'brkgcmy';
colorString = 'blue|red|black|green|cyan|magenta|yellow';

end
function [valList0Unique,valList1Unique] = getValsFromIndex(indexListUnique,stimResults,fieldName)
if isfield(stimResults,[fieldName 'Index'])
    
    indexList = getfield(stimResults,[fieldName 'Index']); %#ok<*GFLD>
    if strcmpi(fieldName,'contrast')
        valList0 = getfield(stimResults,[fieldName '0PC']);
        valList1 = getfield(stimResults,[fieldName '1PC']);
    else
        valList0 = getfield(stimResults,[fieldName '0Hz']);
        valList1 = getfield(stimResults,[fieldName '1Hz']);
    end
    
    numList = length(indexListUnique);
    valList0Unique = zeros(1,numList);
    valList1Unique = zeros(1,numList);
    for i=1:numList
        valList0Unique(i) = unique(valList0(indexListUnique(i)==indexList));
        valList1Unique(i) = unique(valList1(indexListUnique(i)==indexList));
    end
else
    valList0Unique = indexListUnique;
    valList1Unique = indexListUnique;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load Data
function [analogChannelsStored,timeVals,goodStimPos] = loadlfpInfo(folderLFP) %#ok<*STOUT>
load(fullfile(folderLFP,'lfpInfo.mat'));
end
function [neuralChannelsStored,SourceUnitID] = loadspikeInfo(folderSpikes)
fileName = fullfile(folderSpikes,'spikeInfo.mat');
if exist(fileName,'file')
    load(fileName);
else
    neuralChannelsStored=[];
    SourceUnitID=[];
end
end
function [parameterCombinations,cValsUnique,tValsUnique,eValsUnique,...
    aValsUnique,sValsUnique] = loadParameterCombinations(folderExtract)

load(fullfile(folderExtract,'parameterCombinations.mat'));
end
function stimResults = loadStimResults(folderExtract)
load (fullfile(folderExtract,'stimResults'));
end

% added associate functions- for EEG analysis
% same function as used in displaySingleChannelSRC- To get list of  electrodes
function outString = getStringFromValues(valsUnique,decimationFactor)

if length(valsUnique)==1
    outString = convertNumToStr(valsUnique(1),decimationFactor);
else
    outString='';
    for i=1:length(valsUnique)
        outString = cat(2,outString,[convertNumToStr(valsUnique(i),decimationFactor) '|']);
    end
    outString = [outString 'all'];
end

    function str = convertNumToStr(num,f)
        if num > 16384
            num=num-32768;
        end
        str = num2str(num/f);
    end
end

% to get data from all electrodes for the specified condition- for the topolots
function [Data]=getDataSRC(c,t,e,a,s,folderName,folderLFP,analogChannels)

% Load Trial Numbers for  given Parameter Combinations
folderExtract = fullfile(folderName,'extractedData');
folderSegment = fullfile(folderName,'segmentedData');
[parameterCombinations] = loadParameterCombinations(folderExtract);
trialNums = cell2mat(parameterCombinations(c,t,e,a,s));
try
    load(fullfile(folderSegment,'badTrials.mat'));
catch
    disp('No Bad Trials file found...');
    badTrials=[];
end
goodPos=setdiff(trialNums,badTrials);
% Extraction- for the specified channels
for iC = 1:size(analogChannels,2)
    analogData = load(fullfile(folderLFP,['elec' num2str(analogChannels(iC)) '.mat']));
    Data(iC,:,:)=analogData.analogData(goodPos,:);
    clear analogData
end
end

% this function loads the channel location files to be passed to
% topoplots;these locs files are stored in a folder called montages;
% the montages folder is currently stored in the folderSourceString
function [chanlocs,hemBipolarLocs,bipolarLocs] = loadChanLocs(capType,refType,folderMontage)

if nargin<2;    refType = 1; end

if strcmp(capType,'actiCap64')
    if (refType == 4) %bipolar
        load(fullfile(folderMontage,'Montages','bipolarChanlocsActiCap64.mat'));
        chanlocs = eloc;
    else
        load(fullfile(folderMontage,'Montages','actiCap64.mat'));
    end
    load(fullfile(folderMontage,'Montages','hemBipChInfoActiCap_64.mat'));
    load(fullfile(folderMontage,'Montages','bipChInfoActiCap64.mat'));
    
else % passive cap
    if (refType == 4)
        % have to generate the bipolar montage file for passive cap
        load(fullfile(folderMontage,'Montages','bipolarChanlocsActiCap64.mat'));
        chanlocs = eloc;
    else
        load(fullfile(folderMontage,'Montages','brainCap64.mat'));
    end
end
end

% this function is used to compute the STFT for the time frequency
% topoplots and single channel TF plots.
function [SRaw,SChange,t,f] = getSTFT(analogData,movingWin,params,timeVals,BLMin,BLMax)

[SRaw,t,f]=mtspecgramc(analogData',movingWin,params);
t = t + timeVals(1); % shift the t values to the actual time
tBL = intersect(find(t>=BLMin),find(t<=BLMax)); % baseline time indices

SRawBL = SRaw(tBL,:,:);

mlogSRawBL = conv2Log(mean(mean(SRawBL,1),3));
SChange = 10*(conv2Log(mean(SRaw,3)) - repmat(mlogSRawBL,size(SRaw,1),1));
SRaw = 10*(conv2Log(mean(SRaw,3)));

end

% to get titles for the figures- reference type used
function [refType]=getRefschemeName(refChanIndex)
switch(refChanIndex)
    case 1
        refType='Single wire';
    case 2
         refType='Hemisphere';
    case 3
        refType='Average';
    case 4
       refType='Bipolar';
    otherwise
        refType=strcat('Referenced to Electrode- ',num2str(refChanIndex-4));
end
end