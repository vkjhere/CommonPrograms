function deleteFromProtocolList(dbFileName,expDate,protocolName)

% convert from ddmmyy to yymmdd
expDate = [expDate(5:6) expDate(3:4) expDate(1:2)];

try
    % open protocol list database file
    dbID = openProtocolList(dbFileName,'rw');
    % delete from protocol list table in database
    mksqlite(dbID, ...
        ['DELETE FROM PROTOCOLLIST ' ...
        'WHERE EXPDATE=''' expDate ''' AND PROTOCOLNAME=''' protocolName ''';']);
    % close database
    closeProtocolList(dbID,dbFileName);
catch sqliteException
    if exist('dbID','var'); closeProtocolList(dbID,dbFileName); end
    rethrow(sqliteException);
end

end
