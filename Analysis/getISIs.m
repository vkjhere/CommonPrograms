% H =getISIs(X,tRange)
% This function takes all spikes in the time range tRange and returns an
% array of inter-spike interval (ISI) times

% Supratim Ray 07/01/14

function Y = getISIs(X,tRange)

if isempty(X)
    Y=[];
else
    Y = [];
    for i=1:length(X)
        spk = X{i};
        spkShort = spk(intersect(find(spk>=tRange(1)),find(spk<tRange(2))));
        if length(spkShort)>1
            d = diff(spkShort);
            Y = [Y;d(:)]; %#ok<AGROW>
        end
    end
end