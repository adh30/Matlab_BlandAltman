% ba1
%% Bland Altman analysis 
% version 1alpha
% 01/04/24
% ADH

%%
clear
close all
%% select datafile
[file,path] = uigetfile('*.xlsx');
selectedfile = fullfile(path,file);
% import data
data = importfile(selectedfile, "Results", [2, Inf]);

%% calculations
diffval = data{:,1}- data{:,2};
meanval = (data{:,1} + data{:,2})/2;
mn=(mean(data{:,1})+mean(data{:,1}))/2;
meandiff = mean(diffval);
sddiff = std(diffval);
sx=std(data{:,1});
sy=std(data{:,2});
sxy   = sqrt((sx^2+sy^2)/2);
n=length(diffval); % number of obs
% ICC
icc = sum((data{:,1}-mn).*(data{:,2}-mn)) / sxy^2 / (n-1); % https://warwick.ac.uk/fac/sci/statistics/staff/academic-research/nichols/scripts/matlab/icc.m

% Calculate LOA
uloa = meandiff+1.96*sddiff;
lloa = meandiff-1.96*sddiff;
% calculate 95% CI for LOA
seloa=1.71*sddiff/sqrt(n); % Barlett 2008
uluci = uloa + 1.96*seloa;
ullci = uloa - 1.96*seloa;
lluci = lloa + 1.96*seloa;
lllci = lloa - 1.96*seloa;
% Bradley-Blackwood Test tests the joint hypothesis (Hj) that : mean1=mean2  
% and variance1=variance2 in a Bland Altman test. 
% Bradley, E.L., Blackwood, L.G. (1989). Comparing paired data: 
% A simultaneous test for means and variances.Am. Stat.43:234–235
% This function also tests if there is evidence of a correlation
% between the difference and the average (confusingly also termed a 
% Bradley Blackwood test by some but Pitman's test by others).
% x1 and x2 are two sets of paired measurements using different
% methods (or repeated measures using the same method). 
% See Macro for Pairwise Reliability Indicators, Carl M. Russell, 
% Medical College of Georgia for original data used to test/validate code
% ADH 01/06/21
% Calculate Bradley Blackwood
regmdl = fitlm(meanval,diffval);                            % Regress
Pitman_p = regmdl.Coefficients.pValue(2);                   % P value slope
[rho,~] = corr(meanval, diffval);                           % r for Pitman
Fvalue=((sum(diffval.^2)-regmdl.SSE)/2)/(regmdl.SSE/(n-2)); % F value (BB)
Brad_Black_p=1-fcdf(Fvalue, 2, n-2);                        % P value (BB) 

% is there statistically significant bias
[~,Bias_p,Bias_ci,stats] = ttest(diffval);

% Calculate Coefficient of variation


%% plot
plot1 = subplot(2,1,1);
plot(plot1,meanval, diffval, 'o','MarkerFaceColor', 'k', 'MarkerEdgeColor','k', 'DisplayName','')
hold on; 
ax = gca;
ax.XAxisLocation = "origin";
h = legend('Location','northeastoutside'); 
% add lines
yline(meandiff,'b--')
yline([uloa lloa],'r--')
% 
xaxl = xlim; % query axis limits

% add bias CI
z = patch([xaxl(1) xaxl(2) xaxl(2) xaxl(1)], [Bias_ci(1) Bias_ci(1) ...
    Bias_ci(2) Bias_ci(2)], 'r');
z.FaceAlpha = 0.1;
z.EdgeColor = "none";

% make shaded areas for 95%CI of LOA

yaxlu = [uluci ullci]; 
yaxll = [lluci lllci];
a=patch([xaxl(1) xaxl(2) xaxl(2) xaxl(1)], [yaxlu(1) yaxlu(1) yaxlu(2) ...
    yaxlu(2)], 'g');
b=patch([xaxl(1) xaxl(2) xaxl(2) xaxl(1)], [yaxll(1) yaxll(1) yaxll(2)...
    yaxll(2)], 'g');
a.FaceAlpha = 0.05;
a.EdgeColor = "none";
b.FaceAlpha = 0.05;
b.EdgeColor = "none";
% yline([uluci ullci lluci lllci],'g--')

% if pValueslope <0.05 draw regression line on plot
if Pitman_p<0.05 
    plot (meanval, regmdl.Fitted,'b-')
end

% add labels
xlabel('Mean');
ylabel('Difference');
legend off
% add legend
% lgd = legend("","Mean difference","Upper LOA","Lower LOA");

% add results to subplot 2
subplot(2,1,2)
axis off
% Results strings line by line
str1=strcat("Mean difference = " + num2str(meandiff, 3) + " LOA (" + ...
    num2str(uloa,3) + ", " + num2str(lloa, 3) +")");
str2 = strcat("Bradley Blackwood p = " + num2str(Brad_Black_p,2) + "; r = "+ ...
    num2str(rho, 2) + "; Pitman p = " + num2str(Pitman_p, 2) + ...
    "; p (bias) = " + num2str(Bias_p,2)); 
str3 = strcat("ICC = " + num2str(icc,2));
str= {str1,str2,str3}; 
text(min(xlim), max(ylim), str, 'Horiz','left', 'Vert','top')

