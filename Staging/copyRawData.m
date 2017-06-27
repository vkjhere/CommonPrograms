function copyRawData(subjectName,gridType,useTheseIndices,folderSourceString,folderDestinationString)
% Copies raw data from a specified source to a specified destination. Uses
% a few checks to skip copying unnecessary files and data.
%
% Siddhesh Salelkar     14-Sep-2016

%%%%%%%%%%%%%% Get subject details %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(subjectName,'abu') || strcmpi(subjectName,'rafiki') || strcmpi(subjectName,'alpa')
    [allExpDates,allProtocolNames] = eval(['allProtocols' ...
        upper(subjectName(1)) subjectName(2:end) upper(gridType(1)) gridType(2:end)]);
else % thereafter, we use the automated scripts
    [allExpDates,allProtocolNames] = getAllProtocols(subjectName,gridType);
end

%%%%%%%%%%%% Copy all data from source to destination %%%%%%%%%%%%%%
for i=1:length(useTheseIndices)
    expDate = allExpDates{useTheseIndices(i)};
    protocolName = allProtocolNames{useTheseIndices(i)};
    
    %%%%%%%%%%%% Check for input %%%%%%%%%%%%%
    protocolPrefix = [subjectName expDate];
    protocolFilePrefix = [protocolPrefix protocolName];
    if exist(fullfile(folderSourceString,'Data','rawData',protocolPrefix),'dir')
        protocolSourceData = fullfile(folderSourceString,'Data','rawData',protocolPrefix);
    elseif exist(fullfile(folderSourceString,'data','rawData',protocolPrefix),'dir')
        protocolSourceData = fullfile(folderSourceString,'data','rawData',protocolPrefix);
    elseif exist(fullfile(folderSourceString,'Data','rawData',subjectName,protocolPrefix),'dir')
        protocolSourceData = fullfile(folderSourceString,'Data','rawData',subjectName,protocolPrefix);
    elseif exist(fullfile(folderSourceString,'data','rawData',subjectName,protocolPrefix),'dir')
        protocolSourceData = fullfile(folderSourceString,'data','rawData',subjectName,protocolPrefix);
    else
        error('Unable to find data folder in source path');
    end
    
    %%%%%%%%%%%% Prepare for output %%%%%%%%%%%%%
    folderOut = fullfile(folderDestinationString,'data','rawData',protocolPrefix);
    makeDirectory(folderOut);

    fileListToCopy = {};
    folderInContents = dir(protocolSourceData);
    fileListInFolder = {folderInContents.name};
    
    % compile list of files to copy
    dataFileTypes = {['^' protocolFilePrefix]}; % protocol data files
    impedanceFileTypes = {'impedance'}; % impedance files
    crosstalkFileTypes = {'crosstalk'}; % crosstalk files
    otherFileTypes = {'^screen shot' 'screenshot' 'ssvep' 'alpha' 'noise'}; % screenshots
    fileTypesToCopy = [dataFileTypes impedanceFileTypes crosstalkFileTypes otherFileTypes];
    for j=1:length(fileTypesToCopy)
        idx = ~cellfun('isempty', ...
            cellfun(@(x)regexpi(x,fileTypesToCopy{j}),fileListInFolder, ...
            'UniformOutput',false));
        fileListToCopy = union(fileListToCopy,fileListInFolder(idx));
    end
    
    %%%%%%%%%%%% Copy everything %%%%%%%%%%%%%
    disp([num2str(useTheseIndices(i)) ': ' protocolFilePrefix]); 
    for j=1:length(fileListToCopy)
        fileIn = fullfile(protocolSourceData,fileListToCopy{j});
        fileCopy(fileIn,folderOut);
    end
end
end

function fileCopy(sourceFile,destFolder)

% Change this for your PC
teraCopyExists = 1;
teraCopyPath = '"C:\Program Files\TeraCopy\TeraCopy.exe"';

[~,sourceFileName,sourceFileExt] = fileparts(sourceFile);
fprintf(['=> ' sourceFileName sourceFileExt ' ... ']);

% try to use TeraCopy on Windows machines that have it
if ispc && teraCopyExists && exist(teraCopyPath,'file')
    system([teraCopyPath ' Copy "' sourceFile '" "' destFolder '" /OverwriteOlder /Close']);
    fprintf('TeraCopy''ed.\n');
    return;
end

% try to check if destination file exists and is newer than source file
useMatlabCopy = 0;
destFile = fullfile(destFolder,[sourceFileName sourceFileExt]);
if ~exist(destFile,'file')
    useMatlabCopy = 1;
else
    sourceFileInfo = dir(sourceFile);
    destFileInfo = dir(destFile);
    if isempty(sourceFileInfo.bytes) || isempty(sourceFileInfo.datenum) || ...
            isempty(destFileInfo.bytes) || isempty(destFileInfo.datenum)
        fprintf('Not copied.\n');
        error('Error querying source/destination file.');
    end
    if sourceFileInfo.bytes ~= destFileInfo.bytes || ...
            sourceFileInfo.datenum > destFileInfo.datenum
        useMatlabCopy = 1;
    end
end

if useMatlabCopy
    [status,errMsg] = copyfile(sourceFile,destFolder,'f');
    if status == 0
        fprintf('Not copied.\n');
        error(['Error copying file:' errMsg]);
    else
        fprintf('Copied.\n');
    end
else
    fprintf('Skipped.\n');
end
end