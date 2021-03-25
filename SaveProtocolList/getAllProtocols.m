function [expDates,protocolNames,stimTypes,comments] = getAllProtocols(subjectName,gridType)

% try to see whether protocol list is mapped to the path
[protocolListFileName,protocolListPath] = getProtocolListPath(subjectName,gridType);
if isempty(protocolListPath)
    throw(MException([mfilename ':fileNotFound'], ...
        'Protocol list file not found. Please check your MATLAB path.'));
end

% read data
protocolListFile = fullfile(protocolListPath,protocolListFileName);
try
    [expDates,protocolNames,stimTypes,comments] = readProtocolList(protocolListFile);
catch readException
    throw(MException([mfilename ':protocolListError'], readException.message));
end

end
