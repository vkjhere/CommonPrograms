% This function is used to save the digital data recorded from MonkeyLogic

function saveDigitalDataML(digitalEvents,digitalTimeStamps,folderExtract)

digitalEvents = getValue(digitalEvents);

% Write the digitalCodes
makeDirectory(folderExtract);
save(fullfile(folderExtract,'digitalEvents.mat'),'digitalTimeStamps','digitalEvents');
end
function outNumList = getValue(numList)
% Only 8-bits are used. The remaining ones are usually 1 but could be anything
binStr = dec2bin(numList);
outNumList = bin2dec(binStr(:,9:16));
end