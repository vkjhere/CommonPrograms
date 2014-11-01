function barPlotForCellArrays(Y,xVals,colorList,plotHandle,plotSEM)

L=length(Y);
if ~exist('xVals','var');            xVals=1:L;                         end
if ~exist('colorList','var');        colorList=jet(L);                  end
if ~exist('plotHandle','var');       plotHandle=gca;                    end
if ~exist('plotSEM','var');          plotSEM=0;                         end

axes(plotHandle); %#ok<*MAXES> % makes this axis current

for i=1:L
    bar(xVals(i),mean(Y{i}),'FaceColor',colorList(i,:)); hold on;
    
    if plotSEM
        errorbar(xVals(i),mean(Y{i}),std(Y{i})/sqrt(length(Y{i})),'k');
    else
        errorbar(xVals(i),mean(Y{i}),std(Y{i}),'k');
    end
end
axis tight;
end