function trimProtocolList(subjectName,gridType)

% check whether file is accessible
[protocolListFileName,protocolListPath] = getProtocolListPath(subjectName,gridType);
if isempty(protocolListPath)
    explainString = ['The protocol list database file is not in your MATLAB path. ' ...
        'Make sure the file exists and is added to the path, and try again.'];
    uiwait(errordlg(explainString,'Error finding protocol list','modal'));
    return;
end

protocolListFile = fullfile(protocolListPath,protocolListFileName);

% try to guess whether file is local and warn user
[fileRoot,nextDir,~]=fileparts(protocolListFile);
while ~isempty(nextDir)
    mountDir = nextDir; %#ok<*NASGU>
    [fileRoot,nextDir,~]=fileparts(fileRoot);
end
explainString = {'The protocol list file', protocolListFile, ...
    'appears to be a file on the local drive and not network-mapped.', '', ...
    'Are you sure you want to continue?'};
if ispc % Windows
    % guess based on drive letter
    if ~exist('localFileCreateResponse','var') && fileRoot(1) < 'L'
        localFileConfirmResponse = questdlg(explainString,'Confirm file modification','OK','Cancel','Cancel');
        if strcmp(localFileConfirmResponse,'Cancel')
            return;
        end
    end
else % Linux/Mac
    % TODO: guess based on mount point
end

% pre-update sanity check
try
    [expDates,protocolNames,stimType] = getAllProtocols(subjectName,gridType);
catch protocolsReadException
    uiwait(errordlg(protocolsReadException.message, ...
        'Error reading pre-update protocol list','modal'));
    return;
end
if checkListSanity(expDates,protocolNames,stimType)
    uiwait(errordlg(['The protocol list is not in correct format. ' ...
        'Please verify protocol list database sanity offline by hand.'], ...
        'Pre-update sanity check failed','modal'));
    return;
end
    
% present dialog to gather values for current entry
validInputs = 0;
expDateString = '';
protocolString = '';
while ~validInputs
    prompt = {'Experiment date string (ddmmyy):', ...
        'Protocol name string (like GRF_001):'};
    dlg_title = ['Updating ' protocolListFile];
    num_lines = repmat([1 length(protocolListFile)+10],2,1);
    defaultans = {expDateString,protocolString};
    answer = inputdlg(prompt,dlg_title,num_lines,defaultans,'on');
    if isempty(answer)
        return;
    end

    % update default strings
    expDateString = answer{1};
    protocolString = answer{2};
    
    % do some sanity check on the inputs
    if ~isDateOk(expDateString)
        uiwait(msgbox('Invalid date string. Please correct and try again.',mfilename,'modal'));
        continue;
    end
    if ~isProtocolOk(protocolString)
        uiwait(msgbox('Invalid protocol string. Please correct and try again.',mfilename,'modal'));
        continue;
    end
    validInputs = 1;
end

% check whether entry exists
try
    [stimType,comment] = searchProtocolList(protocolListFile,expDateString,protocolString);
catch protocolListSearchException
    uiwait(errordlg(protocolListSearchException.message, ...
        'Error searching protocol list','modal'));
    return;
end

% return if protocol entry doesn't exist
if isempty(stimType)
    explainString = ['This entry does not exist in the protocol list ' ...
        'database file. Please check the input values and try again.'];
    uiwait(errordlg(explainString,'Error finding protocol entry','modal'));
    return;
end

% otherwise confirm deletion
explainString = {'You are about to delete the following entry from ' ...
    'the protocol list database:','', ...
    ['Experiment date: ''' expDateString ''''], ...
    ['Protocol name: ''' strrep(protocolString,'_','\_') ''''], ...
    ['Stimulus type: ''' num2str(stimType) ''''], ...
    ['Comment: ''' comment ''''],'', ...
    'Are you sure you want to continue?'};
options.Interpreter = 'tex';
options.Default = 'Cancel';
entryDeleteConfirmResponse = questdlg(explainString,'Confirm file modification','OK','Cancel',options);
if strcmp(entryDeleteConfirmResponse,'Cancel')
    return;
end

% delete entry from the list!
h = waitbar(1,'Updating protocol list...','Name',mfilename);
try
    deleteFromProtocolList(protocolListFile,expDateString,protocolString);
catch protocolDeleteException
    uiwait(errordlg(protocolDeleteException.message, ...
        'Error updating protocol list','modal'));
    close(h);
    return;
end

% check whether everything was successful
waitbar(1,h,'Performing post-update sanity check...');
try
    [expDatesUpdated,protocolNamesUpdated,stimTypeUpdated,commentsUpdated] = getAllProtocols(subjectName,gridType);
catch protocolsReadException
    uiwait(errordlg(protocolsReadException.message, ...
        'Error reading post-update protocol list','modal'));
    close(h);
    return;
end

% post-update sanity check
if checkListSanity(expDatesUpdated,protocolNamesUpdated,stimTypeUpdated)
    uiwait(errordlg(['The protocol list is not in correct format. ' ...
        'Please verify protocol list database sanity offline by hand.'], ...
        'Post-update sanity check failed','modal'));
    close(h);
    return;
end

% update successful: generate the human-readable file
waitbar(1,h,'Re-generating human-readable protocol list...');
generateProtocolList(subjectName,gridType,...
    expDatesUpdated,protocolNamesUpdated,stimTypeUpdated,commentsUpdated,'m');
close(h);

end
