% This is the main program used to find bad trials in EEG data. 

% Note: This program was built on top of findBadTrialsEEG_GAV_v2 to _v5.
% This program was used for finding bad trials for 350 subjects who were
% part of ADGammaProject. This program will be modified in future commits
% to be compatible with the data format used here.

function [badTrials,allBadTrials,badTrialsUnique,badElecs,totalTrials,slopeValsVsFreq] = findBadTrialsWithEEG(eegDataAllElecs,timeVals,eegChannelsStored,elecNames,elecImpedance,capType,eyeDataDeg,eyeRangeMS,FsEye,saveFolder)   

    % Initializations
    highPassCutOff = 1.6; % Hz   
    checkPeriod = [-0.500 0.750]; % s
    ImpedanceCutOff = 25; % KOhm
    time_threshold  = 6;
    psd_threshold = 6;
    badTrialThreshold = 30; % Percentage
    
    tapersPSD = 1; % No. of tapers used for computation of slopes    
    slopeRange = {[56 86]}; % Hz, slope range used to compute slopes
    freqsToAvoid = {[0 0] [8 12] [46 54] [96 104]}; % Hz

    Fs = 1/(timeVals(2) - timeVals(1)); %Hz
    checkBaselinePeriod = [-0.5 0]; % For computing slopes for artifact rejection
    numElectrodes = length(eegChannelsStored);

    % Setting MT parameters
    params.tapers   = [3 5];
    params.pad      = -1;
    params.Fs       = Fs;
    params.fpass    = [0 200];
    params.trialave = 0;

    % 1. Get bad trials from eye data
    if exist('FsEye','var') && ~isempty(FsEye)
        badEyeTrials = findBadTrialsFromEyeData_v2(eyeDataDeg,eyeRangeMS,FsEye,checkPeriod)'; % added by MD 10-09-2017; Modified by MD 03-09-2019
    else
        badEyeTrials = [];
    end
    originalTrialInds = 1:size(eegDataAllElecs,2);
    originalTrialInds(badEyeTrials) = [];
    clear eyeDataDeg
    
    totalTrials = size(eegDataAllElecs,2);     
    badTrialsUnique.badEyeTrials = badEyeTrials;
    
    % 2. Get electrode impedances for rejecting noisy electrodes (impedance > 25k)
    EEGelectrodeLabels = load('EEGelectrodeLabels.mat','EEGelectrodeLabels');
    EEGelectrodeLabels = EEGelectrodeLabels.EEGelectrodeLabels;    
    electrodeLabelsList = EEGelectrodeLabels(2:end,(ismember(EEGelectrodeLabels(1,:),capType))); electrodeLabelsList(cellfun(@isempty,electrodeLabelsList))=[];
    clear elecInds; elecInds = NaN(1,length(electrodeLabelsList));
    for iML = 1:length(electrodeLabelsList)
        elecInds(iML) = find(strcmp(electrodeLabelsList(iML),elecNames));
    end
    elecImpedance = elecImpedance(elecInds); % Remap the electrodes according to the standard montage workspace
    GoodElec_Z = elecImpedance<ImpedanceCutOff;
    nBadElecs{1} = ~GoodElec_Z;
    
    % 3. Analysis for each trial and each electrode        
    if exist('highPassCutOff','var') || ~isempty(highPassCutOff) % Defining filter
        d1 = designfilt('highpassiir','FilterOrder',8, ...
            'PassbandFrequency',highPassCutOff,'PassbandRipple',0.2, ...
            'SampleRate',Fs);
    end
    
    allBadTrials = cell(1,numElectrodes);
    hW1 = waitbar(0,'Processing electrodes...');
    for iElec=1:numElectrodes

        waitbar((iElec-1)/numElectrodes,hW1,['Processing electrode: ' num2str(iElec)]);
        if ~GoodElec_Z(iElec); allBadTrials{iElec} = NaN; continue; end % Analyzing only those electrodes with impedance < 25k        
        analogData = squeeze(eegDataAllElecs(iElec,:,:));
        analogData(badEyeTrials,:) = [];
        
        % determine indices corresponding to the check period
        checkPeriodIndices = timeVals>=checkPeriod(1) & timeVals<checkPeriod(2);
        analogData = analogData(:,checkPeriodIndices);
        
        clear analogDataSegment; analogDataSegment = analogData;

        if exist('highPassCutOff','var') || ~isempty(highPassCutOff)    % high pass filter            
            clear analogData; analogData = filtfilt(d1,analogDataSegment')';
        end       

        % subtract dc
        analogData = analogData - repmat(mean(analogData,2),1,size(analogData,2));

        % Check time-domain waveforms
        numTrials = size(analogData,1);                            % excluding bad eye trials
        meanTrialData = nanmean(analogData,1);                     % mean trial trace
        stdTrialData = nanstd(analogData,[],1);                    % std across trials

        tDplus = (meanTrialData + (time_threshold)*stdTrialData);    % upper boundary/criterion
        tDminus = (meanTrialData - (time_threshold)*stdTrialData);   % lower boundary/criterion

        tBoolTrials = sum((analogData > ones(numTrials,1)*tDplus) | (analogData < ones(numTrials,1)*tDminus),2);

        clear badTrialsTimeThres
        badTrialsTimeThres = find(tBoolTrials>0);

        % Check PSD
        clear powerVsFreq;
        [powerVsFreq,~] = mtspectrumc(analogDataSegment',params);
        powerVsFreq = powerVsFreq';

        clear meanTrialData stdTrialData tDplus
        meanTrialData = nanmean(powerVsFreq(setdiff(1:size(powerVsFreq,1),badTrialsTimeThres),:),1);                     % mean trial trace
        stdTrialData = nanstd(powerVsFreq(setdiff(1:size(powerVsFreq,1),badTrialsTimeThres),:),[],1);                    % std across trials

        tDplus = (meanTrialData + (psd_threshold)*stdTrialData);    % upper boundary/criterion
        clear tBoolTrials; tBoolTrials = sum((powerVsFreq > ones(numTrials,1)*tDplus),2);
        clear badTrialsFreqThres; badTrialsFreqThres = find(tBoolTrials>0);

        tmpBadTrialsAll = unique([badTrialsTimeThres;badTrialsFreqThres]);

        % Remap bad trial indices to original indices        
        allBadTrials{iElec} = originalTrialInds(tmpBadTrialsAll);
        
        % Calculate number of unique bad trials for each thresholding criterion
        badTrialsUnique.timeThres{iElec} = originalTrialInds(badTrialsTimeThres);
        badTrialsUnique.freqThres{iElec} = originalTrialInds(setdiff(badTrialsFreqThres,badTrialsTimeThres));
        
    end
    close(hW1);

    % 4. Remove electrodes containing more than x% bad trials
    badTrialUL = (badTrialThreshold/100)*numTrials;
    badTrialLength=cellfun(@length,allBadTrials);
    badTrialLength(nBadElecs{1})=NaN; % Removing the bad impedance electrodes
    nBadElecs{2} = logical(badTrialLength>badTrialUL)';
    allBadTrials(nBadElecs{2}) = {NaN};
    
    % 5. Find common bad trials across all electrodes subject to conditions
    commonBadTrialsAllElecs = trimBadTrials(allBadTrials);
    
    % 6. Find common bad trials across visual electrodes
    electrodeListVisGamma = load('electrodeListVisGamma.mat','electrodeListVisGamma');
    electrodeListVisGamma = electrodeListVisGamma.electrodeListVisGamma;
    electrodeList = electrodeListVisGamma.(capType);
    commonBadTrialsVisElecs=[];
    for iElec=1:length(electrodeList)
        if ~isnan(allBadTrials{1,electrodeList(iElec)}); commonBadTrialsVisElecs=union(commonBadTrialsVisElecs,allBadTrials{electrodeList(iElec)}); end
    end

    badTrialsUnique.commonBadTrialsAllElecs = commonBadTrialsAllElecs;
    badTrialsUnique.commonBadTrialsVisElecs = commonBadTrialsVisElecs;
    
    badTrials = union(commonBadTrialsVisElecs,commonBadTrialsAllElecs);
    
    % 6. PSD Slope calculation across baseline period
    checkPeriodIndicesPSD = timeVals>=checkBaselinePeriod(1) & timeVals<checkBaselinePeriod(2);    
    params.tapers   = [(tapersPSD+1)/2 tapersPSD];
    slopeValsVsFreq = cell(1,numElectrodes);
    
    eegDataAllElecs = eegDataAllElecs(:,setdiff(originalTrialInds,badTrials),checkPeriodIndicesPSD);
    for iElec=1:numElectrodes
        if isnan(allBadTrials{1,iElec}); slopeValsVsFreq{iElec} = {NaN,NaN}; goodSlopeFlag(iElec) = false; continue; end %#ok<AGROW>
        
        % Computing slopes
        analogDataPSD = squeeze(eegDataAllElecs(iElec,:,:));
%         analogDataPSD = analogDataPSD - repmat(mean(analogDataPSD,2),1,size(analogDataPSD,2));
        
        clear powerVsFreq freqVals
        [powerVsFreq,freqVals] = mtspectrumc(analogDataPSD',params);
        slopeValsVsFreq{iElec} = getSlopesPSDBaseline_v2((log10(mean(powerVsFreq,2)))',freqVals,slopeRange,[],freqsToAvoid);
        goodSlopeFlag(iElec) = slopeValsVsFreq{iElec}{2}>0; %#ok<AGROW>
    end
    
    nanElecs = find(cell2mat(cellfun(@(x)any(isnan(x)),allBadTrials,'UniformOutput',false))); % MD: 09-09-2019
    
    badElecs.elecImpedance = elecImpedance;
    badElecs.badImpedanceElecs = find(nBadElecs{1});
    badElecs.noisyElecs = find(nBadElecs{2});
    badElecs.flatPSDElecs = setdiff(find(~goodSlopeFlag),nanElecs)'; % storing bad electrode labels  
    
    if exist('saveFolder','var') && ~isempty(saveFolder)
        disp(['Saving ' num2str(length(union(badTrialsUnique.badEyeTrials,badTrials))) ' bad trials']);
        badTrialsFileName = fullfile(saveFolder,'badTrials_v5.mat');
        if exist(badTrialsFileName,'file'); delete(fullfile(saveFolder,'badTrials_v5.mat')); end;
        save(badTrialsFileName,'badTrials','allBadTrials','badTrialsUnique','badElecs','totalTrials','slopeValsVsFreq');
    else
        disp('Bad trials will not be saved..');
    end
end

function [newBadTrials] =  trimBadTrials(allBadTrials)
    badElecThreshold = 10; % Percentage
    
    % a. Taking union across bad electrodes for conditions 1 and 2
    newBadTrials=[];
    numElectrodes = length(allBadTrials);
    for iElec=1:numElectrodes
        if ~isnan(allBadTrials{1,iElec}); newBadTrials=union(newBadTrials,allBadTrials{iElec}); end
    end

    % b. Co-occurence condition - Counting the trials which occurs in more than x% of the electrodes
    badTrialElecs = zeros(1,length(newBadTrials));
    for iTrial = 1:length(newBadTrials)
        for iElec = 1:numElectrodes
            if isnan(allBadTrials{1,iElec}); continue; end % Discarding the electrodes where the bad trials are NaN because of this NaN entries in badTrials have zero in 'badTrialElecs'
            if find(newBadTrials(iTrial)==allBadTrials{1,iElec})
                badTrialElecs(iTrial) = badTrialElecs(iTrial)+1;
            end
        end
    end
    newBadTrials(badTrialElecs<(badElecThreshold/100.*numElectrodes))=[];
end
