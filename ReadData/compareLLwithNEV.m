% NEV stores information of only one of the Gabors (activeSide)
% 0 - Mapping 0
% 1 - Mapping 1
% 2 - Task Gabor

function matchingParameters=compareLLwithNEV(folderExtract,activeSide,showResults)

load(fullfile(folderExtract,'LL.mat'));
load(fullfile(folderExtract,'stimResults.mat'));

% Use trialResults instead of digitalEvents, as it has the same data nicely
% stored. Moreover, it has pruned data in case of NEV digital code loss.
load(fullfile(folderExtract,'trialResults.mat'));

% Need this flag to determine how to compare NEV and LL data
if ~exist('digitalCodeLoss','var'); digitalCodeLoss = 0; end

% Compare basic properties
% NEV stores information of only one of the Gabors (activeSide)

if activeSide==3 % Plaid
    validMap1 = find(LL.stimType1==5);
    validMap2 = find(LL.stimType2==5);
    if validMap1 ~= validMap2
        error('Mapping0/1 components for Plaid do not match');
    end
    aziLL = alternateCombineArrays(LL.azimuthDeg1(validMap1),LL.azimuthDeg2(validMap2));
    eleLL = alternateCombineArrays(LL.elevationDeg1(validMap1),LL.elevationDeg2(validMap2));
    sigmaLL = alternateCombineArrays(LL.sigmaDeg1(validMap1),LL.sigmaDeg2(validMap2));
    
    if isfield(LL,'radiusDeg1')
        radiusExists = 1;
        radiusLL = alternateCombineArrays(LL.radiusDeg1(validMap1),LL.radiusDeg2(validMap2));
    else
        radiusExists = 0;
    end
    sfLL = alternateCombineArrays(LL.spatialFreqCPD1(validMap1),LL.spatialFreqCPD2(validMap2));
    oriLL = alternateCombineArrays(LL.orientationDeg1(validMap1),LL.orientationDeg2(validMap2));
    conLL = alternateCombineArrays(LL.contrastPC1(validMap1),LL.contrastPC2(validMap2)); 
    tfLL = alternateCombineArrays(LL.temporalFreqHz1(validMap1),LL.temporalFreqHz2(validMap2)); 
    timeLL = alternateCombineArrays(LL.time1(validMap1),LL.time2(validMap2));
    
elseif activeSide==2 % SRC protocol
    aziLL = LL.azimuthDeg;
    eleLL = LL.elevationDeg;
    sigmaLL = LL.sigmaDeg;
    radiusExists=0;
    sfLL = LL.spatialFreqCPD;
    oriLL = LL.orientationDeg;
    
elseif activeSide==0 % Map0
    validMap = find(LL.stimType1==1);
    aziLL = LL.azimuthDeg1(validMap);
    eleLL = LL.elevationDeg1(validMap);
    sigmaLL = LL.sigmaDeg1(validMap);
    
    if isfield(LL,'radiusDeg1')
        radiusExists = 1;
        radiusLL = LL.radiusDeg1(validMap);
    else
        radiusExists = 0;
    end
    sfLL = LL.spatialFreqCPD1(validMap);
    oriLL = LL.orientationDeg1(validMap);
    conLL = LL.contrastPC1(validMap); 
    tfLL = LL.temporalFreqHz1(validMap); 
    timeLL = LL.time1(validMap);
    
elseif activeSide==1 % Map2
    
    validMap = find(LL.stimType2==1);
    aziLL = LL.azimuthDeg2(validMap);
    eleLL = LL.elevationDeg2(validMap);
    sigmaLL = LL.sigmaDeg2(validMap);
    if isfield(LL,'radiusDeg2')
        radiusExists = 1;
        radiusLL = LL.radiusDeg2(validMap);
    else
        radiusExists = 0;
    end
    sfLL = LL.spatialFreqCPD2(validMap);
    oriLL = LL.orientationDeg2(validMap);
    conLL = LL.contrastPC2(validMap); 
    tfLL = LL.temporalFreqHz2(validMap);
    timeLL = LL.time2(validMap);
end

% Compare

if ~isfield(stimResults,'azimuth') %#ok<NODEF>
    disp('Stimulus parameter values not present in the digital data stream, taking from the LL file...');
    
    stimResults.azimuth = aziLL;
    stimResults.elevation = eleLL;
    stimResults.contrast = conLL;
    stimResults.temporalFrequency = tfLL;
    if radiusExists
        stimResults.radius = radiusLL;
    end
    stimResults.sigma = sigmaLL;
    stimResults.orientation = oriLL;
    stimResults.spatialFrequency = sfLL;
    matchingParameters = [];
    
    % Saving updated stimResults 
    disp('Saving stimResults after taking values from LL file...');
    save(fullfile(folderExtract,'stimResults.mat'),'stimResults');
else

    % In case of NEV digital code loss, we need to only compare the
    % corresponding leading information from LL.
    if digitalCodeLoss
        warning(['NEV digital code loss case detected! ' ...
            'LL/NEV comparison may be misleading...']); %#ok<*UNRCH>
        disp(['numStimCodesLL: ' num2str(length(aziLL)) ...
            ', numStimCodesNEV: ' num2str(length(stimResults.azimuth))]);
        aziLL = aziLL(1:length(stimResults.azimuth));
        eleLL = eleLL(1:length(stimResults.elevation));
        sigmaLL = sigmaLL(1:length(stimResults.sigma));
        if radiusExists
            radiusLL = radiusLL((1:length(stimResults.radius)));
        end
        sfLL = sfLL(1:length(stimResults.spatialFrequency));
        oriLL = oriLL(1:length(stimResults.orientation));
        conLL = conLL(1:length(stimResults.contrast));
        tfLL = tfLL(1:length(stimResults.temporalFrequency));
        if activeSide == 3
            % For Plaid, Map0+Map1 onset times
            timeLL = timeLL(1:2*length(stimResults.time));
        else
            timeLL = timeLL(1:length(stimResults.time));
        end
    end
    
    if compareValues(aziLL,stimResults.azimuth)
        matchingParameters.azimuth=1;
        disp('Azimuths match.');
    else
        matchingParameters.azimuth=0; %#ok<*STRNU>
        disp('*****************Azimuths do not match!!');
    end
    
    if compareValues(eleLL,stimResults.elevation)
        matchingParameters.elevation=1;
        disp('Elevations match.');
    else
        matchingParameters.elevation=0;
        disp('*****************Elevations do not match!!');
    end
    
    if compareValues(sigmaLL,stimResults.sigma)
        matchingParameters.sigma=1;
        disp('Sigmas match.');
    else
        matchingParameters.sigma=0;
        disp('*****************Sigmas do not match!!');
    end
    
    if radiusExists
        if compareValues(radiusLL,stimResults.radius)
            matchingParameters.radius=1;
            disp('Radii match.');
        else
            matchingParameters.radius=0;
            disp('******************Radii do not match!!');
        end
    else
        disp('Radius parameter does not exist');
    end
    
    if compareValues(sfLL,stimResults.spatialFrequency)
        matchingParameters.spatialFrequency=1;
        disp('Spatial Freq match.');
    else
        matchingParameters.spatialFrequency=0;
        disp('**************Spatial Freq do not match!!');
    end
    
    if compareValues(oriLL,stimResults.orientation)
        matchingParameters.orientation=1;
        disp('Orientations match.');
    else
        matchingParameters.orientation=0;
        disp('***************Orientations do not match!!');
    end
    
    if activeSide~=2
        if compareValues(conLL,stimResults.contrast)
            matchingParameters.contrast=1;
            disp('Contrasts match.');
        else
            matchingParameters.contrast=0;
            disp('*******************Contrasts do not match!!');
        end
        
        if compareValues(tfLL,stimResults.temporalFrequency)
            matchingParameters.temporalFrequency=1;
            disp('Temporal frequencies match.');
        else
            matchingParameters.temporalFrequency=0;
            disp('*******************Temporal frequencies do not match!!');
        end
    end
end

if showResults
    % always plot in a new figure, helpful when extracting data in a batch
    figure;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Match start Times

    clear xD xL lxD lxL
    
    % Get NEV start times
    endPos = find(strcmp(trialEvents,'TS'));
    if isempty(endPos)
        error('No trialEvent named TS');
    end
    xD = trialResults(endPos).times;
    xD = diff(xD); lxD = length(xD); xD=xD(:);
    
    xL = LL.startTime;    
    xL = diff(xL); lxL = length(xL); xL=xL(:);
    
    if lxD == lxL
        disp(['Number of startTrials: ' num2str(lxD+1)]);
        subplot(231)
        plot(xD,'b.'); hold on; plot(xL,'ro'); hold off;
        ylabel('Difference in Start Times (s)');
        axis tight
        legend('Dig','LL','Location','SouthEast');
        
        subplot(234)
        plot(1000*(xD-xL),'b');
        ylabel('Digital-Lablib times (ms)');
        xlabel('Trial Number');
        axis tight
        
        matchingParameters.maxChangeTrialsPercent = 100*max(abs(xD-xL) ./ xD);
    else
        if ~digitalCodeLoss
            warning('Number of trials different in Digital vs Lablib data even though no digital codes were lost!');
        end
        disp(['Num of startTrials: digital: ' num2str(lxD+1) ' , LL: ' num2str(lxL+1)]);
        mlx = min(lxD,lxL);
        subplot(231)
        plot(xD(1:mlx),'b.'); hold on; plot(xL(1:mlx),'ro'); hold off;
        ylabel('Start Times (s)');
        axis tight
        legend('Dig','LL','Location','SouthEast');
        
        subplot(234)
        plot(1000*(xD(1:mlx)-xL(1:mlx)),'b');
        ylabel('Difference in start times (ms)');
        xlabel('Trial Number');
        axis tight
        
        matchingParameters.maxChangeTrialsPercent = 100*max(abs(xD(1:mlx)-xL(1:mlx)) ./ xD(1:mlx));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if activeSide~=2
        % Match stimulus Frame positions
        clear xD xL lxD lxL
        xD = stimResults.time;
        if activeSide == 3
            % For Plaid, take Map1 onset times
            xL = timeLL(2:2:end)/1000;
        else
            xL = timeLL/1000;
        end
        
        xD = diff(xD); lxD = length(xD); xD=xD(:);
        xL = diff(xL); lxL = length(xL); xL=xL(:);
        
        if lxD == lxL
            disp(['Number of stimOnset: ' num2str(lxD+1)]);
            subplot(232)
            plot(xD,'b.'); hold on; plot(xL,'ro'); hold off;
            ylabel('Difference in StimOn time (s)');
            axis tight; 
            %legend('Dig','LL','Location','SouthEast');
            
            subplot(235)
            plot(1000*(xD-xL),'b');
            ylabel('Digital-Lablib times (ms)');
            xlabel('Stim Number');
            axis tight
            
            matchingParameters.maxChangeStimOnPercent = 100*max(abs(xD-xL) ./ xD);
        else
            if ~digitalCodeLoss
                warning('Number of stimuli different in Digital vs Lablib data even though no digital codes were lost!');
            end
            disp(['Number of stimOnset: digital: ' num2str(lxD+1) ' , LL: ' num2str(lxL+1)]);
            mlx = min(lxD,lxL);
            subplot(232)
            plot(xD(1:mlx),'b.'); hold on; plot(xL(1:mlx),'ro'); hold off;
            ylabel('StartOn Frame');
            axis tight
            %legend('Dig','LL','Location','SouthEast');
            
            subplot(235)
            plot(1000*(xD(1:mlx)-xL(1:mlx)),'b');
            ylabel('Digital-Lablib times (ms)');
            xlabel('Stim Number');
            axis tight
            
            matchingParameters.maxChangeStimOnPercent = 100*max(abs(xD(1:mlx)-xL(1:mlx)) ./ xD(1:mlx));
        end
    end
    
    %%%% Match EOTCodes
    
    clear xD xL lxD lxL
    % Get NEV start times
    clear endPos

    endPos = find(strcmp(trialEvents,'TE'));
    if isempty(endPos)
        error('No trialEvent named TE');
    end
    xD = trialResults(endPos).value;
    lxD = length(xD); xD = xD(:);
    
    xL = LL.eotCode;
    lxL = length(xL); xL = xL(:);
    
    if lxD == lxL
        disp(['Number of eotCodes: ' num2str(lxD)]);
        subplot(233)
        plot(xD,'b.'); hold on; plot(xL,'ro'); hold off;
        ylabel('eotCode number');
        axis tight
        %legend('Dig','LL','Location','SouthEast');
        
        subplot(236)
        plot(xD-xL,'b');
        ylabel('\delta eotCode number');
        xlabel('Trial Number');
        axis tight
        
    else
        if ~digitalCodeLoss
            warning('Number of eotCodes different in Digital vs Lablib data even though no digital codes were lost!');
        end
        disp(['Number of eotCodes: digital: ' num2str(lxD) ' , LL: ' num2str(lxL)]);
        mlx = min(lxD,lxL);
        subplot(233)
        plot(xD(1:mlx),'b.'); hold on; plot(xL(1:mlx),'ro'); hold off;
        ylabel('eotCode number');
        axis tight
        %legend('Dig','LL','Location','SouthEast');
        
        subplot(236)
        plot(xD(1:mlx)-xL(1:mlx),'b');
        ylabel('\Delta eotCode number');
        xlabel('Trial Number');
        axis ([0 mlx -7 7]);
    end
end

end
function result = compareValues(x,y)
thres = 10^(-2);
if max(x-y) < thres;
    result=1;
else
    result=0;
end
end
function z = alternateCombineArrays(x,y)
z = [x' y']';
z = z(:)';
end
