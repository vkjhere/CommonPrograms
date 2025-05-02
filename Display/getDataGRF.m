% This function is used to get all relevant spike and LFP data for GRF
% protocols

function dataOut = getDataGRF(dataIn,a,e,s,f,o,c,t,blRange,stRange,removeERPFlag,tapers,movingWin)

if ~exist('removeERPFlag','var');      removeERPFlag = 0;               end
if ~exist('tapers','var');             tapers = [1 1];                  end
if ~exist('movingWin','var');          movingWin = [0.25 0.025];        end

goodPos = dataIn.parameterCombinations{a,e,s,f,o,c,t};
goodPos = setdiff(goodPos,dataIn.badTrials);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LFP measures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
signal = dataIn.analogData(goodPos,:);
timeVals = dataIn.timeVals;

% ERP
erp = mean(signal,1);
dataOut.erp = erp;
dataOut.timeVals = timeVals;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PSDs and TFs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fs = round(1/(timeVals(2)-timeVals(1)));
if removeERPFlag
    signal = signal - repmat(erp,size(signal,1),1);
end

% Setup multitaper
params.tapers   = tapers;
params.pad      = -1;
params.Fs       = Fs;
params.fpass    = [0 250];
params.trialave = 1;  % Averaging across trials

% ranges
range = blRange;
rangePos = round(diff(range)*Fs);
blPos = find(timeVals>=blRange(1),1)+ (1:rangePos);
stPos = find(timeVals>=stRange(1),1)+ (1:rangePos);

% PSDs
[dataOut.SBL,dataOut.freqBL] = mtspectrumc(signal(:,blPos)',params);
[dataOut.SST,dataOut.freqST] = mtspectrumc(signal(:,stPos)',params);

% Time-frequency data
[dataOut.STF,timeTF,dataOut.freqTF] = mtspecgramc(signal',movingWin,params);
dataOut.timeTF = timeTF+timeVals(1)-1/Fs;

% Spike data
dataOut.raster = dataIn.spikeData(goodPos);
[dataOut.frVals,dataOut.frTimeVals] = getPSTH(dataIn.spikeData(goodPos),10,[timeVals(1) timeVals(end)]);

% Information about the ranges
dataOut.blRange = blRange;
dataOut.stRange = stRange;
end