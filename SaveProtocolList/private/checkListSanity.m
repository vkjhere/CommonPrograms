function status = checkListSanity(expDates,protocolNames,stimTypes)

status = 1;

% make sure all are cell arrays...
if ~iscell(expDates) || ~iscell(protocolNames) || ~iscell(stimTypes)
    return;
end

% ...and of the same length
if length(expDates) ~= length(protocolNames) || length(protocolNames) ~= length(stimTypes)
    return;
end

% check whether all entries are as expected
for i=1:length(expDates)
    if ~isDateOk(expDates{i}) || ~isProtocolOk(protocolNames{i}) || ...
            ~isStimTypeOk(num2str(stimTypes{i}))
        return;
    end
end

status = 0;

end

