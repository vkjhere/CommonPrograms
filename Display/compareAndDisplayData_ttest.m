% This is the same program as compareAndDisplayData.m in the parent folder.
% The program compareAndDisplayData.m in this folder is an extension that
% allows a comparison of more than two datasets.

% X, Y - 2 arrays of size MxK and NxK. M,N are the number of trials.
% p(i) = p value of a t-test on X(:,i) and Y(i,:)

function p=compareAndDisplayData_ttest(X,Y,xs,hPlot,hSig,colorX,colorY,colorSig,displayPlot,showStdErr)

if ~exist('hPlot','var');       hPlot = subplot(211);                   end
if ~exist('hSig','var');        hSig  = subplot(212);                   end
if ~exist('colorX','var');      colorX = 'b';                           end
if ~exist('colorY','var');      colorY = 'r';                           end
if ~exist('colorSig','var');    colorSig = 'k';                         end
if ~exist('displayPlot','var'); displayPlot=1;                          end

color1 = 'g'; % 0.01
color2  = 'r'; % Bonferroni

% Get significance
[M,K1] = size(X);
[N,K2] = size(Y);

if K1~=K2
    error('The y dimension of the two matrices must be the same');
else
    K=K1;
    [~,p] = ttest2(X,Y,0.05,'both','unequal');
end

% Plot data
if ~exist('xs','var');          xs = 1:K;                               end

% Plot data
if displayPlot
    set(hPlot,'Nextplot','replace');
    plot(hPlot,xs,mean(X),'Color',colorX);
    set(hPlot,'Nextplot','add');
    plot(hPlot,xs,mean(Y),'Color',colorY);
    
    if showStdErr
        plot(hPlot,xs,mean(X)+std(X)/sqrt(M),'Color',colorX,'LineStyle','--');
        plot(hPlot,xs,mean(X)-std(X)/sqrt(M),'Color',colorX,'LineStyle','--');
        plot(hPlot,xs,mean(Y)+std(Y)/sqrt(N),'Color',colorY,'LineStyle','--');
        plot(hPlot,xs,mean(Y)-std(Y)/sqrt(N),'Color',colorY,'LineStyle','--');
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