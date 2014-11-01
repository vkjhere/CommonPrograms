% displays the contents of entries, which is a 2D cell array.
function makeTable(tablePos,entries,clearTable,textColor)

if ~exist('clearTable','var');              clearTable=0;               end

% Make the plot
if length(tablePos)==4     % instead of tablePos, a handle can also be passed
    h=subplot('Position',tablePos);
else
    h=tablePos;
end

if clearTable
    cla(h);
end
axes(h);

% size of the table
[numRows,numCols] = size(entries);

if ~exist('textColor','var')
    for i=1:numRows
        textColor{i} = 'k';            
    end
end

dX=1/numCols;
dY=1/numRows;

% Make lines
for i=1:numRows
    lineXRow(:,i) = [0 1]; lineYRow(:,i) = [i*dY i*dY]; %#ok<*AGROW>
end
for i=1:numCols
    lineXCol(:,i) = [i*dX i*dX]; lineYCol(:,i) = [0 1];
end

line(lineXRow,lineYRow,'color','k'); hold on;
line(lineXCol,lineYCol,'color','k');


% Show the entries
for i=1:numRows
    for j=1:numCols
        textPosX = dX/2+(j-1)*dX;
        textPosY = dY/2+(numRows-i)*dY;
        text(textPosX,textPosY,entries{i,j},'Units','Normalized','HorizontalAlignment','center','color',textColor{i},'Parent',h);
    end
end

set(h,'box','on','XTick',[],'YTick',[],'XTickLabel',[],'YTickLabel',[]);