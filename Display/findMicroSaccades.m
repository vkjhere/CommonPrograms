% This program finds the microsaccades in the 2D array X, which has the eye
% speeds. MS is a cell array of size(X,1), each cell entry gives the
% position when the data exceeds the cutOff. If the data exceeds the cutoff
% for more than 1 consecutive points, only the first one is considered.

function [MS,numMSInRange] = findMicroSaccades(X,cutoff,xs,timeLims)

if ~exist('timeLims','var');     timeLims=[xs(1) xs(end)];               end
N = size(X,1);

MS = cell(1,N);
numMSInRange = zeros(1,N);
for i=1:N
    y = X(i,:);
    ms =  sort([find(y>cutoff) find(y<-cutoff)]);
    
    if ~isempty(ms)
        msI = 1:length(ms);
        dms = [0 diff(ms)];
        msI(dms==1)=[];
        MS{i} = xs(ms(msI));
    else
        MS{i}=[];
    end
    
    numMSInRange(i) = length(intersect(find(MS{i}>=timeLims(1)),find(MS{i}<timeLims(2))));
end
end