% This is the same program as compareAndDisplayData.m in the parent folder.
% The program compareAndDisplayData.m in this folder is an extension that
% allows a comparison of more than two datasets.

% X - A cell array of length L, with each cell being a 2D array of size
% NixK (i=1...L). Ni is the number of stimuli. xs must be of length K.
% testMethod - Anova or KruskalWallis

function p=compareAndDisplayData_AnovaOrKW(X,xs,hPlot,hSig,testMethod,colorX,colorSig,displayPlot,showStdErr)

if ~exist('hPlot','var');        hPlot = subplot(211);                  end
if ~exist('hSig','var');         hSig  = subplot(212);                  end
if ~exist('testMethod','var');   testMethod = 'anova';                  end
if ~exist('colorSig','var');     colorSig = 'k';                        end
if ~exist('displayPlot','var');  displayPlot=1;                         end
if ~exist('showStdErr','var');   showStdErr = 0;                        end

color1 = 'g'; % 0.01
color2  = 'r'; % Bonferroni

L = length(X);
if ~exist('colorX','var');     colorX = jet(L);                          end

N=zeros(1,L);
Ktmp=zeros(1,L);
for i=1:L
    [N(i),Ktmp(i)] = size(X{i});
end

if strcmpi(testMethod,'anova')
    disp('Significance computed using Anova');
elseif strcmpi(testMethod,'kruskalWallis') || strcmpi(testMethod,'KW')
    disp('Significance computed using KruskalWallis');
else
    error('Unknown test name');
end
    

% Get significance
if max(abs(diff(Ktmp)))>0
    error(['The y dimension of the ' num2str(L) ' matrices must be the same']);
else
    K=Ktmp(1);
    p=zeros(1,K);
    for i=1:K  % for each time bin
        
        comparisonData = []; groupValue = [];
        for j=1:L
            comparisonData = cat(1,comparisonData,X{j}(:,i));
            groupValue     = cat(1,groupValue,j+zeros(N(j),1));
        end
        
        if strcmpi(testMethod,'anova')
            p(i) = anova1(comparisonData,groupValue,'off');
        elseif strcmpi(testMethod,'kruskalWallis') || strcmpi(testMethod,'KW')
            p(i) = kruskalwallis(comparisonData,groupValue,'off');
        end
        
    end
end

% Plot data
if ~exist('xs','var');          xs = 1:K;                               end

% Plot data
if displayPlot
    set(hPlot,'Nextplot','replace');
    for i=1:L
        plot(hPlot,xs,mean(X{i}),'Color',colorX(i,:));
        set(hPlot,'Nextplot','add');
        
        if showStdErr
            plot(hPlot,xs,mean(X{i})+std(X{i})/sqrt(N(i)),'Color',colorX(i,:),'LineStyle','--');
            plot(hPlot,xs,mean(X{i})-std(X{i})/sqrt(N(i)),'Color',colorX(i,:),'LineStyle','--');
        end
        
    end
    set(hPlot,'Nextplot','replace');
    axis(hPlot,'tight');
      
    %Plot Significance
    set(hSig,'Nextplot','replace');
    plot(hSig,xs,log10(p),'Color',colorSig);
    set(hSig,'Nextplot','add');
    plot(hSig,xs,log10(0.01+zeros(1,K)),'Color',color1);
    plot(hSig,xs,log10(0.05/K+zeros(1,K)),'Color',color2);
    
    % Significant points
    axis(hSig,'tight');
    axisLims = axis(hSig);
    axis(hSig,[axisLims(1:2) axisLims(3) 0]);
    
    sig1 = find(p<0.01);
    sig2 = find(p<0.05/K);
    
    if ~isempty(sig1)
        plot(hSig,xs(sig1),0,'Color',color1,'marker','*');
    end
    
    if ~isempty(sig2)
        plot(hSig,xs(sig2),0,'Color',color2,'marker','*');
    end
    set(hSig,'Nextplot','replace');
end
end