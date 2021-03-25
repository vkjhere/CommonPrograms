function [protocolListFileName,protocolListPath] = getProtocolListPath(subjectName,gridType)

% construct file name of database file
protocolListFileName = ['allProtocols' ...
    upper(subjectName(1)) subjectName(2:end) ...
    upper(gridType(1)) gridType(2:end) '.sqlite'];

% try to locate database file: if not mapped, protocolListPath is empty
[protocolListPath,~,~] = fileparts(which(protocolListFileName));

end
