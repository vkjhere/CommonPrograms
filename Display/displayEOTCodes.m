function [eotFraction,colorNames,numTrials] = displayEOTCodes(gridPosition,eotCodes,protocolName,showEOTbar,plotHandle,removeBreaks,removeIgnores)

if ~exist('showEOTbar','var');          showEOTbar=1;                   end
if ~exist('removeBreaks','var');        removeBreaks=0;                 end
if ~exist('removeIgnores','var');       removeIgnores=0;                end

%%%%% Some protocol specific definitions
if strncmpi(protocolName,'SRC',3)
    % Definition in SRCContrast
    %{kEOTCorrect = 0, kEOTWrong, kEOTFailed, kEOTBroke, kEOTIgnored, kEOTQuit, kEOTTypes}
    %#define		kEOTFAlarm			(kEOTIgnored + 1) 
    %#define		kEOTDistracted		(kEOTIgnored + 2)
    %#define		kEOTTotal			(kEOTIgnored + 3)
    %#define		kEOTForceQuit		kEOTQuit + 2
    
    eotCodeMeaning{1} = 'correct';              colorNames{1} = 'g';
    eotCodeMeaning{2} = 'wrong';                colorNames{2} = 'r';
    eotCodeMeaning{3} = 'failed';               colorNames{3} = [0.6 0.2 0]; % brown
    eotCodeMeaning{4} = 'broke';                colorNames{4} = 'b';
    eotCodeMeaning{5} = 'ignored';              colorNames{5} = [0.5 0.5 0.5]; % gray
    eotCodeMeaning{6} = 'falseAlarm';           colorNames{6} = [1 0.6 0]; %orange
    eotCodeMeaning{7} = 'distracted';           colorNames{7} = 'm';
    eotCodeMeaning{8} = 'quit';                 colorNames{8} = 'k';

else
   
    eotCodeMeaning{1} = 'correct';              colorNames{1} = 'g';
    eotCodeMeaning{2} = 'wrong';                colorNames{2} = [1 0.6 0]; %orange, because this is like a false alarm
    eotCodeMeaning{3} = 'failed';               colorNames{3} = [0.6 0.2 0]; % brown
    eotCodeMeaning{4} = 'broke';                colorNames{4} = 'b';
    eotCodeMeaning{5} = 'ignored';              colorNames{5} = [0.5 0.5 0.5]; % gray
    eotCodeMeaning{6} = 'quit';                 colorNames{6} = 'k';
end

if removeBreaks
    eotCodes(eotCodes==3)=[];
end
if removeIgnores
    eotCodes(eotCodes==4)=[];
end

numTrials = length(eotCodes);

eotFraction =zeros(1,length(eotCodeMeaning));
for i=1:length(eotCodeMeaning)
    eotFraction(i) = length(find(eotCodes==i-1))/numTrials;
end
    
if showEOTbar
    if ~exist('plotHandle','var') || isempty(plotHandle)
        plotHandle = subplot('Position',gridPosition,'XTickLabel',[],'YTickLabel',[],'box','on');
        cla(plotHandle);
    end
    axes(plotHandle); %#ok<*MAXES>
    plotEOTbar(eotFraction,colorNames);
    set(gca,'XTickLabel',[],'YTickLabel',[]);
    axis([0 1 0 1]);
end
end

function plotEOTbar(eotFraction,colorNames)

startYPos = 0;
for i=1:length(eotFraction)
    fill([0 0 1 1],[startYPos startYPos+eotFraction(i) ...
        startYPos+eotFraction(i) startYPos],colorNames{i});
    
    hold on;
    
    if eotFraction(i)>0.02
        %text(0.5,startYPos+eotFraction(i)/2,[num2str(round(100*eotFraction(i))) '%'],'HorizontalAlignment','center');
        text(0.5,startYPos+eotFraction(i)/2,num2str(round(100*eotFraction(i))),'HorizontalAlignment','center');
    end
    startYPos=startYPos+eotFraction(i);
end
end