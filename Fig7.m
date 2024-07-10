% Initial code used to reproduce Figure 7 of the paper `Unsupervised image
% velocimetry for automated computation of river flow velocities'. Data
% required to produced this plot can be accessed at: 
% https://data.ncl.ac.uk/articles/dataset/User input files/19762027.
% The user also needs to download readtext.m and parsing_html.m from the 
% GitHub pages and have these files accessible on the MATLAB path.

clear all; close all; clc;

% The user must specify the location where the files downloaded from the
% above link are saved to on the hard drive of the PC by modifying the
% `pathIn' variable below, note the trailing backslash at the end:
pathIn = 'D:\OneDrive - Newcastle University\Documents - TENDERLY Archive [Geog]\General\Archive_Dart\testing\';

% Specify the output location where the Figure will be saved to:
pathOut = 'D:\OneDrive - Newcastle University\Documents - TENDERLY Archive [Geog]\General\Archive_Dart\testing\';

validationDataIn     = [pathIn '2018.10.25 Austins_Gaugings.csv']; % location of the flow gauging data
xsIn                 = [pathIn 'transect_2.474max.csv']; % location of the cross-section survey data
addpath(genpath(pathIn));

% if you want to extract the data from each of the html discharge report
% files, set parseData as 1, otherwise leave as 0 to simply use the summary
% csv file (quicker)
parseData   = 0;
if parseData == 1 
    [obsLevel, obsQ, obsVelocity, obdDT] = parsingHTML(pathIn, q_path);
else
    KLTin        = [pathIn 'Combined_out.csv'];
    dataIn_      = readtext(KLTin,',', '','','textual');
    obsLevel     = str2double(dataIn_(2:end,2));
    obsQ         = str2double(dataIn_(2:end,4:6));
    obsVelocity  = str2double(dataIn_(2:end,3));
    obsDT        = dataIn_(2:end,1);
end

% Load the reference flow gauging data
validationIn            = readtext(validationDataIn,',', '','','textual'); % read in the csv file
compiledSampleDates     = strcat(validationIn(17:end,1), {' ' }, validationIn(17:end,2));
validationDischarge     = str2double(validationIn(17:end,3));
validationLevel_t       = str2double(validationIn(17:end,5));
validationLevel         = 4997.916 + validationLevel_t;


% Load the cross-section data
xs_raw      = readtext(xsIn,',', '','','textual'); % read in the csv file
xi          = str2double(xs_raw(2:end,1));
dist        = xi-abs(min(xi));
zi          = str2double(xs_raw(2:end,3));

for a = 1:length(validationLevel)
    zi_t            = zi-validationLevel(a,1);
    idx             = zi_t<0;
    totalArea(a,1)  = trapz(dist(idx),abs(zi_t(idx)));
end

% Load the KLT-IV generated outputs for the calibration period
cal_period   = 4242; % limit to only the calibration period
dataIn_      = readtext(KLTin,',', '','','textual');
obsLevel     = str2double(dataIn_(2:cal_period,2));
obsQ         = str2double(dataIn_(2:cal_period,4:6));
obsVelocity  = str2double(dataIn_(2:cal_period,3));
obsDT        = dataIn_(2:cal_period,1);

for a = 1:length(validationLevel)
    diff1                   = abs(obsLevel - validationLevel(a,1));
    idx                     = find(diff1<0.01);
    compV(a,1:length(idx))  = obsVelocity(idx);
end

comp_dist_v     = transpose(nanmedian(replace_num(compV',0,NaN))); % distributed velocity
yAxisData       = replace_num(compV,0,NaN)'; 
xAxisData       = validationDischarge./totalArea;
keeper1         = find(~isnan(nanmedian(yAxisData)));
yAxisData       = yAxisData(:,keeper1);
xAxisData       = xAxisData(keeper1);
crossSection    = totalArea(keeper1);
stageIn         = validationLevel(keeper1);


f0=figure(1); hold on;
f0.Units='pixels';
set(f0,'Position',[100, 0, 2480./2, 3508./2]); % A4 aspect ratio
f0.Units='normalized';

sub = subplot(3,1,1);
hold on;
pbaspect([1 1 1])
ax1 = gca();
set(ax1,'fontsize', 12)
axis tight


% x = reference; y = image vel
[s1,s2]         = sort(nanmedian(yAxisData));

b_plot = scatter(nanmedian(yAxisData(:,s2)),xAxisData(s2),...
    'MarkerEdgeColor', [0.5,0.5,0.5],...
    'MarkerFaceColor', [0.5,0.5,0.5],...
    'SizeData', 10);

xLims       = [0 2.8];
yLims       = [0 2.1];
set(ax1, 'XLim', xLims)
set(ax1, 'YLim',yLims )
set(ax1,'XTick',(0.5:0.5:xLims(2)),...
    'XTickLabel',(0.5:0.5:xLims(2)))
set(ax1,'YTick',(0:0.5:yLims(2)),...
    'YTickLabel',(0:0.5:yLims(2)))

intervals   = 0:0.01:nanmax(nanmedian(yAxisData(:,s2)));
plot(intervals,intervals,'k--'); 

% Do the statistical analysis - linear regression
[xData,yData]       = prepareCurveData(nanmedian(yAxisData(:,s2)),xAxisData(s2));
t_out               = length(xData);
xData(1:t_out,2)    = 1; % Prep the data for regression
[b,~,~,~,stats]     = regress(yData,xData); % Run the regression
Original_r2         = stats(1,1); % Pull out the R2 value
Original_p          = stats(1,3); % Pull out the p value

dlm                     = fitlm(xData(:,1),yData(:,1),'Intercept',true);
interceptCoefficient    = table2array(dlm.Coefficients(1,1)); % intercept
slopeCoefficient        = table2array(dlm.Coefficients(2,1)); % slope
Pvalue                  = table2array(dlm.Coefficients(2,4)); % p-value
R2                      = round(dlm.Rsquared.Ordinary,2);
plot(intervals,interceptCoefficient + (slopeCoefficient.*intervals),'k-') % Plot the best-fit

annotationText1     = ['$U_{a} = ' num2str(interceptCoefficient) '+'  num2str(slopeCoefficient) ...
    ' \textnormal{(} \tilde{U}_{s} \textnormal{)}$'];
annotationText2     = ['$\mathrm{R^{2} =' num2str(R2) '}$'];


s1 = [0.37 0.885646401771864 0.12399641122175 0.0315107493246396];
annotation(f0,'textbox',...
    s1,...
    'String',{annotationText1, annotationText2} ,...
    'FitBoxToText','on',...
    'Interpreter','LaTex',...
    'LineStyle', 'none',...
    'fontsize', 12);

annotation(gcf,'textbox',...
    s1 - [0.05 -0.01 0 0],...
    'String',{'[A]'},...
    'FitBoxToText','on',...
    'Interpreter','LaTex',...
    'FontSize', 12,...
    'EdgeColor','none')

set(ax1,'TickLabelInterpreter','latex')
ylabel('$U_{a} \ \textnormal{[m s\textsuperscript{-1}]}$','Interpreter','LaTex');
set(ax1,'fontsize', 12)
%xlabel('$\tilde{U}_{s} \ \textnormal{[m s\textsuperscript{-1}]}$','Interpreter','LaTex');
set(ax1,'Box','on')
grid on;

%%

% Load the KLT-IV generated outputs for the validation period
cal_period   = 4243:length(dataIn_); % limit to only the validation period
dataIn_      = readtext(KLTin,',', '','','textual');
obsLevel     = str2double(dataIn_(cal_period,2));
obsQ         = str2double(dataIn_(cal_period,4:6));
obsVelocity  = str2double(dataIn_(cal_period,3));
obsDT        = dataIn_(cal_period,1);

for a = 1:length(validationLevel)
    diff1                   = abs(obsLevel - validationLevel(a,1));
    idx                     = find(diff1<0.01);
    compV(a,1:length(idx))  = obsVelocity(idx);
end

comp_dist_v     = transpose(nanmedian(replace_num(compV',0,NaN))); % distributed velocity
yAxisData       = replace_num(compV,0,NaN)'; 
yAxisData       = interceptCoefficient + slopeCoefficient .* (yAxisData);
xAxisData       = validationDischarge./totalArea;
keeper1         = find(~isnan(nanmedian(yAxisData)));
yAxisData       = yAxisData(:,keeper1);
xAxisData       = xAxisData(keeper1);
crossSection    = totalArea(keeper1);
stageIn         = validationLevel(keeper1);

ax2                 = subplot(3,1,2); 
set(ax2,'Position',[0.13 0.470899227292163 0.775 0.215735294117647]);
set(ax2,'fontsize', 12)
hold on
pbaspect([1 1 1])
axis equal
set(ax2,'TickLabelInterpreter','latex')

% x = reference; y = image vel
[s1,s2]         = sort(nanmedian(yAxisData));

b_plot = scatter(nanmedian(yAxisData(:,s2)),xAxisData(s2),...
    'MarkerEdgeColor', [0.5,0.5,0.5],...
    'MarkerFaceColor', [0.5,0.5,0.5],...
    'SizeData', 10);

xLims       = [0 2.4];
yLims       = [0 2.4];
set(ax2, 'XLim', xLims)
set(ax2, 'YLim',yLims )
set(ax2,'XTick',(0.5:0.5:xLims(2)),...
    'XTickLabel',(0.5:0.5:xLims(2)))
set(ax2,'YTick',(0:0.5:yLims(2)),...
    'YTickLabel',(0:0.5:yLims(2)))

intervals   = 0:0.01:nanmax(nanmedian(yAxisData(:,s2)));
plot(intervals,intervals,'k--'); 

% Do the statistical analysis - linear regression
[xData,yData]       = prepareCurveData(nanmedian(yAxisData(:,s2)),xAxisData(s2));
t_out               = length(xData);
xData(1:t_out,2)    = 1; % Prep the data for regression
[b,~,~,~,stats]     = regress(yData,xData); % Run the regression
Original_r2         = stats(1,1); % Pull out the R2 value
Original_p          = stats(1,3); % Pull out the p value

dlm                     = fitlm(xData(:,1),yData(:,1),'Intercept',true);
interceptCoefficient    = table2array(dlm.Coefficients(1,1)); % intercept
slopeCoefficient        = table2array(dlm.Coefficients(2,1)); % slope
Pvalue                  = table2array(dlm.Coefficients(2,4)); % p-value
R2                      = round(dlm.Rsquared.Ordinary,2);
plot(intervals,interceptCoefficient + (slopeCoefficient.*intervals),'k-') % Plot the best-fit

annotationText1     = ['$U_{a} = ' num2str(interceptCoefficient) '+'  num2str(slopeCoefficient) ...
    ' \textnormal{(} \tilde{U}_{s} \textnormal{)}$'];
annotationText2     = ['$\mathrm{R^{2} =' num2str(R2) '}$'];

s1 = [0.37 0.646286652431097 0.126390047546268 0.0339302952919482];

annotation(gcf,'textbox',...
    s1,... 
    'String',{annotationText1, annotationText2} ,...
    'FitBoxToText','on',...
    'Interpreter','LaTex',...
    'FontSize', 12,...
    'EdgeColor','none');

annotation(gcf,'textbox',...
    s1 - [0.05 -0.01 0 0],...
    'String',{'[B]'},...
    'FitBoxToText','on',...
    'Interpreter','LaTex',...
    'FontSize', 12,...
    'EdgeColor','none')

set(ax2,'TickLabelInterpreter','latex')
ylabel('$U_{a} \ \textnormal{[m s\textsuperscript{-1}]}$','Interpreter','LaTex');
set(ax2,'fontsize', 12)
xlabel('$\tilde{U}_{s} \ \textnormal{[m s\textsuperscript{-1}]}$','Interpreter','LaTex');
set(ax2,'Box','on')
grid on;

%% Export the outputs
outDir    = [pathOut 'fig7.svg'];
saveas(f0,outDir,'svg');




