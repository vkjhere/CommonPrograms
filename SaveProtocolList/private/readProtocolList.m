function [expDates,protocolNames,stimTypes,comments] = readProtocolList(dbFileName)

try
    % open protocol list database file
    dbID = openProtocolList(dbFileName,'ro');
    % read entries
    list = mksqlite(dbID, ...
        ['SELECT EXPDATE,PROTOCOLNAME,STIMTYPE,COMMENT ' ...
         'FROM PROTOCOLLIST ' ...
         'ORDER BY EXPDATE,PROTOCOLNAME ASC;']);
    % close database
    closeProtocolList(dbID,dbFileName);
catch sqliteException
    if exist('dbID','var'); closeProtocolList(dbID,dbFileName); end
    rethrow(sqliteException);
end

expDates = {}; protocolNames = {}; stimTypes = {}; comments = {};
for i=1:length(list)
    % restore ddmmyy format
    expDates{i} = [list(i).EXPDATE(5:6) list(i).EXPDATE(3:4) ...
        list(i).EXPDATE(1:2)]; %#ok<*AGROW>
    protocolNames{i} = list(i).PROTOCOLNAME;
    stimTypes{i} = list(i).STIMTYPE;
    comments{i} = list(i).COMMENT;
end

end
