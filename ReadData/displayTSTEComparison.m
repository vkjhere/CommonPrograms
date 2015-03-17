% Compares TrialStart (TS) and TrialEnd (TE) obtained from different data
% streams and uses LL data to constrcut stimResults

function maxDiffTrialStartTimePercent=displayTSTEComparison(folderExtract)

load(fullfile(folderExtract,'LL.mat'));
load(fullfile(folderExtract,'digitalEvents.mat'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% trialStart
clear xD xL lxD lxL
[~,xD] = getDigitalData(digitalCodeInfo,'TS');
xL = LL.startTime;

xD = diff(xD); lxD = length(xD); xD=xD(:);
xL = diff(xL); lxL = length(xL); xL=xL(:);

if lxD == lxL
    disp(['Number of startTrials: ' num2str(lxD)]);
    subplot(221)
    plot(xD,'b.'); hold on; plot(xL,'ro'); hold off;
    ylabel('Difference in Start Times (s)');
    legend('Dig','LL','Location','SouthEast');
    
    subplot(223)
    plot(1000*(xD-xL),'b');
    ylabel('Digital-Lablib times (ms)');
    xlabel('Trial Number');
    
    maxDiffTrialStartTimePercent = 100*max(abs(xD-xL) ./ xD);
else
    error(['Num of startTrials: digital: ' num2str(lxD+1) ' , LL: ' num2str(lxL+1)]);
%     mlx = min(lxD,lxL);
%     subplot(231)
%     plot(xD,'b.'); hold on; plot(xL,'ro'); hold off;
%     ylabel('Start Times (s)');
%     legend('Dig','LL','Location','SouthEast');
%     
%     subplot(234)
%     plot(1000*(xD(1:mlx)-xL(1:mlx)),'b');
%     ylabel('Difference in start times (ms)');
%     xlabel('Trial Number');
%     
%     matchingParameters.maxChangeTrialsPercent = 100*max(abs(xD(1:mlx)-xL(1:mlx)) ./ xD(1:mlx));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% TrialEnd: Match EOTCodes
xD = getDigitalData(digitalCodeInfo,'TE');
xL=LL.eotCode;

lxD = length(xD); xD=xD(:);
lxL = length(xL); xL=xL(:);

if lxD == lxL
    disp(['Number of eotCodes: ' num2str(lxD)]);
    subplot(222)
    plot(xD,'b.'); hold on; plot(xL,'ro'); hold off;
    ylabel('eotCode number');
    %legend('Dig','LL','Location','SouthEast');
    
    subplot(224)
    plot(xD-xL,'b');
    ylabel('\delta eotCode number');
    xlabel('Trial Number');
    
    if max(abs(xD-xL))>0
        error('EOT codes do not match');
    end
else
    error(['Number of stimOnset: digital: ' num2str(lxD) ' , LL: ' num2str(lxL)]);
%     mlx = min(lxD,lxL);
%     subplot(233)
%     plot(xD,'b.'); hold on; plot(xL,'ro'); hold off;
%     ylabel('eotCode number');
%     axis tight
%     %legend('Dig','LL','Location','SouthEast');
%     
%     subplot(236)
%     plot(xD(1:mlx)-xL(1:mlx),'b');
%     ylabel('\Delta eotCode number');
%     xlabel('Trial Number');
%     axis ([0 mlx -7 7]);
end
end
function [value,time] = getDigitalData(digitalCodeInfo,codeName)

for i=1:length(digitalCodeInfo)
    if strcmp(digitalCodeInfo(i).codeName,codeName)
        endPos=i;
        break;
    end
end

if ~exist('endPos','var')
    error(['No trialEvent named ' codeName]);
else
    value=digitalCodeInfo(endPos).value;
    time=digitalCodeInfo(endPos).time;
end

end