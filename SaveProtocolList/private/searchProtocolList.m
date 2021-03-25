function [stimType,comment] = searchProtocolList(dbFileName,expDate,protocolName)

% convert from ddmmyy to yymmdd
expDate = [expDate(5:6) expDate(3:4) expDate(1:2)];

try
    % open protocol list database file
    dbID = openProtocolList(dbFileName,'ro');
    % search protocol list table in database
    entry = mksqlite(dbID, ...
        ['SELECT STIMTYPE,COMMENT FROM PROTOCOLLIST ' ...
        'WHERE EXPDATE=''' expDate ''' AND PROTOCOLNAME=''' protocolName ''';']);
    if isempty(entry)
        stimType = '';
        comment = '';
    else
        stimType = entry.STIMTYPE;
        comment = entry.COMMENT;
    end
    % close database
    closeProtocolList(dbID,dbFileName);
catch sqliteException
    if exist('dbID','var'); closeProtocolList(dbID,dbFileName); end
    rethrow(sqliteException);
end

end
