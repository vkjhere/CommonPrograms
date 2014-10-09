% This program computes the CCG between two spike trains computed using the
% equation from Kohn and Smith, 2005.

% Inputs
% spike1, spike2: cell arrays of spike times in seconds. They must be of
% the same length 
% tRangeS = [tMinS tMaxS] is the time range (in seconds)
% dMS: The time resolution at which spikes are binned (default: 1 ms)

function [ccg,xs,rSC,ccgShift] = getCCG(spike1,spike2,tRangeS,dMS)

if ~exist('dMS','var');                    dMS=1;                        end

% Convert spike times to binned analog data of zeros and ones
analogSpikeTrain1 = convertSpikeTimes2Bins(spike1,tRangeS,dMS);
analogSpikeTrain2 = convertSpikeTimes2Bins(spike2,tRangeS,dMS);

for i=1:size(analogSpikeTrain1,2)
    R12(i,:) = xcorr(analogSpikeTrain1(:,i),analogSpikeTrain2(:,i),'unbiased')/(dMS/1000); %#ok<*AGROW>
end

% Shift predictor
for i=1:size(analogSpikeTrain1,2)-1
    R12Shift(i,:) = xcorr(analogSpikeTrain1(:,i),analogSpikeTrain2(:,i+1),'unbiased')/(dMS/1000);
end

h1 = getSpikeCounts(spike1,tRangeS);
h2 = getSpikeCounts(spike2,tRangeS);

% mean fr
dT = diff(tRangeS);
meanFR1 = mean(h1)/dT;
meanFR2 = mean(h2)/dT;

% CCG
ccg = mean(R12)/sqrt(meanFR1*meanFR2);
ccgShift = mean(R12Shift)/sqrt(meanFR1*meanFR2);

dS=dMS/1000;
xs = -dT+dS:dS:dT-dS;

% rSC
rSC = (mean(h1.*h2) - mean(h1)*mean(h2))/(std(h1)*std(h2));
end