function updateProtocolList(subjectName,gridType)

% check whether file is accessible and if not, confirm new file creation
[protocolListFileName,protocolListPath] = getProtocolListPath(subjectName,gridType);
if isempty(protocolListPath)
    explainString = ['The protocol list database file is not in your MATLAB path. ' ...
        'A new file will be created in the current directory. Are you sure you want to continue?'];
    localFileCreateResponse = questdlg(explainString,'Confirm file creation','OK','Cancel','Cancel');
    if strcmp(localFileCreateResponse,'Cancel')
        return;
    end
    protocolListPath = pwd;
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

% if file exists, guess protocol prefix/suffix and stimulus type using
% previous history
expDateString = dateStringToday();
protocolString = 'GRF_001';
stimTypeString = '1';
if ~exist('localFileCreateResponse','var')
    % file exists and is mapped
    try
        [expDates,protocolNames,stimType] = getAllProtocols(subjectName,gridType);
    catch protocolsReadException
        uiwait(errordlg(protocolsReadException.message, ...
            'Error reading pre-update protocol list','modal'));
        return;
    end
    % pre-update sanity check
    if checkListSanity(expDates,protocolNames,stimType)
        uiwait(errordlg(['The protocol list is not in correct format. ' ...
            'Please verify protocol list database sanity offline by hand.'], ...
            'Pre-update sanity check failed','modal'));
        return;
    end
    
    if ~isempty(expDates) && strcmp(expDates{end},expDateString) % another protocol already added today
        % increment protocol number but retain stim type
        lastProtocol = protocolNames{end};
        upos = strfind(lastProtocol,'_');
        if ~isempty(upos)
            pnum = str2num(lastProtocol(upos+1:end)); %#ok<*ST2NM>
            protocolString = [lastProtocol(1:upos) sprintf(['%0' num2str(length(lastProtocol)-upos) 'd'],pnum+1)];
        end
        stimTypeString = num2str(stimType{end});
    end
else
    % empty protocol list
    expDates = {}; protocolNames = {}; stimType = {};
end

% present dialog to gather values for current entry
validInputs = 0;
while ~validInputs
    prompt = {'Experiment date string (ddmmyy):', ...
        'Protocol name string (like GRF_001):', ...
        'Stimulus type number (1 - RFMap, 2 - Gamma, 3 - TF, 4 - GammaDiscont):', ...
        'Comment (optional):'};
    dlg_title = ['Updating ' protocolListFile];
    num_lines = repmat([1 length(protocolListFile)+10],4,1);
    defaultans = {expDateString,protocolString,stimTypeString,''};
    answer = inputdlg(prompt,dlg_title,num_lines,defaultans,'on');
    if isempty(answer)
        return;
    end

    % update default strings
    expDateString = answer{1};
    protocolString = answer{2};
    stimTypeString = answer{3};
    commentString = answer{4};
    
    % do some sanity check on the inputs
    if ~isDateOk(expDateString)
        uiwait(msgbox('Invalid date string. Please correct and try again.',mfilename,'modal'));
        continue;
    end
    if ~isProtocolOk(protocolString)
        uiwait(msgbox('Invalid protocol string. Please correct and try again.',mfilename,'modal'));
        continue;
    end
    if ~isStimTypeOk(stimTypeString)
        uiwait(msgbox('Invalid stimulus type. Please correct and try again.',mfilename,'modal'));
        continue;
    end
    validInputs = 1;
end

% add entry to the list!
h = waitbar(1,'Updating protocol list...','Name',mfilename);
if exist('localFileCreateResponse','var')
    try
        createProtocolList(protocolListFile);
    catch protocolListCreateException
        uiwait(errordlg(protocolListCreateException.message, ...
            'Error initializing protocol list','modal'));
        close(h);
        return;
    end
end
try
    insertIntoProtocolList(protocolListFile,expDateString,protocolString,str2num(stimTypeString),commentString);
catch protocolInsertException
    uiwait(errordlg(protocolInsertException.message, ...
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
waitbar(1,h,'Generating human-readable protocol list...');
generateProtocolList(subjectName,gridType,...
    expDatesUpdated,protocolNamesUpdated,stimTypeUpdated,commentsUpdated,'m');
close(h);

if exist('localFileCreateResponse','var')
    helpdlg(['The protocol list files have been created in the current ' ...
        'working directory. You should copy them to a shared network ' ...
        'location and map that to your MATLAB path.'],mfilename);
end

end

function d = dateStringToday()
t = datevec(now);
d = sprintf('%02d%02d%02d',t(3),t(2),rem(t(1),100));
end
