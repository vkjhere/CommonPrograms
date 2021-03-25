function OK = isStimTypeOk(stimTypeString)
OK = 0;
% length check
if isempty(stimTypeString); return; end
% format check
if ~prod(double(ismember(stimTypeString,'0':'9'))); return; end
% value check
stimType = str2num(stimTypeString);
if isempty(stimType) || stimType < 1 || stimType > 9 || ...
        stimType ~= uint32(stimType); return; end
OK = 1;
end
