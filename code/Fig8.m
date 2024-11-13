% Initial code used to reproduce Figure 8 of the paper `Unsupervised image
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

% You must also specify the location where the unzipped Q Reports html
% files are saved to:
q_path = [pathIn 'Q Reports\'];

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
    totalArea_init(a,1)  = trapz(dist(idx),abs(zi_t(idx)));
end

for a = 1:length(obsLevel)
    zi_t            = zi-obsLevel(a,1);
    idx             = zi_t<0;
    totalArea(a,1)  = trapz(dist(idx),abs(zi_t(idx)));
end

for a = 1:length(validationLevel)
    diff1                       = abs(obsLevel - validationLevel(a,1));
    idx                         = find(diff1<0.01);
    compV(a,1:length(idx))      = obsVelocity(idx);
    totalArea2(a,1:length(idx)) = totalArea(idx);

    for b = 1:3 % three Q estimates
        FieldName                           = ['q' num2str(b)];
        compQ.(FieldName)(a,1:length(idx))  = obsQ(idx,b);
    end

    for c = 1
        FieldName                           = ['time' num2str(c)];
        timeQ.(FieldName)(a,1:length(idx))  = obsDT(idx,c);
    end

end

totalArea2      = replace_num(totalArea2,0,NaN);
unpackStruct    = @(s) cellfun(@(name) assignin('base',name,getfield(s,name)),fieldnames(s));
unpackStruct(compQ);
unpackStruct    = @(s) cellfun(@(name) assignin('base',name,getfield(s,name)),fieldnames(s));
unpackStruct(timeQ);


comp_dist_v   = transpose(nanmedian(replace_num(compV',0,NaN))); % distributed velocity
comp_q_v(:,1) = transpose(nanmedian(replace_num(q1./totalArea2,0,NaN)'));
comp_q_v(:,2) = transpose(nanmedian(replace_num(q2./totalArea2,0,NaN)'));
comp_q_v(:,3) = transpose(nanmedian(replace_num(q3./totalArea2,0,NaN)'));

comp_q_v_all{1} = transpose(replace_num(q1./totalArea2,0,NaN)');
comp_q_v_all{2} = transpose(replace_num(q2./totalArea2,0,NaN)');
comp_q_v_all{3} = transpose(replace_num(q3./totalArea2,0,NaN)');

comp_q(:,1) = transpose(nanmedian(replace_num(q1,0,NaN)'));
comp_q(:,2) = transpose(nanmedian(replace_num(q2,0,NaN)'));
comp_q(:,3) = transpose(nanmedian(replace_num(q3,0,NaN)'));


xData = validationDischarge./totalArea_init;
yData = replace_num(compV,0,NaN);

%%

f0=figure(1); hold on;
f0.Units='pixels';
set(f0,'Position',[100, 0, 2480./2, 3508./2]); % A4 aspect ratio
f0.Units='normalized';
%ax1 = gca();
ax1 = subplot(1,2,1);
set(ax1,'fontsize', 12)
pbaspect([1 1 1])
axis tight
set(ax1,'TickLabelInterpreter','latex')


% Using level versus reference q
comparison(:,1) = (validationLevel_t);
comparison(:,2) = validationDischarge;
[xData,yData]   = prepareCurveData(comparison(:,1),comparison(:,2));
scatter(yData,xData,...
    'MarkerEdgeColor', [0.5,0.5,0.5],...
    'MarkerFaceColor', [0.5,0.5,0.5],...
    'SizeData', 20);

% Using level vs Uxs q
comparison(:,1) = (validationLevel_t);
comparison(:,2) = (comp_q_v(:,3).*totalArea_init); % using the Fr output
perDiff(:,1)    = ((comparison(:,2) - validationDischarge)./validationDischarge).*100;
[xData,yData]   = prepareCurveData(comparison(:,1),comparison(:,2));
scatter(yData,xData,...
    'MarkerEdgeColor', [252,141,89]./255,...
    'MarkerFaceColor', [252,141,89]./255,...
    'SizeData', 20); hold on;

%Calculate the summary stats for when using the Fr output
simData = [1:length(validationDischarge);...
    (comp_q_v(:,3).*totalArea_init)']';
fr_use = simData;
obsData = [1:length(validationDischarge);...
    validationDischarge']';
idx_use = ~isnan(simData(:,2));

% Calculate NSE
[fr_n_out] = nashsutcliffe(obsData(idx_use,1:2), simData(idx_use,1:2)); % not appropriate

% Calculate r2 instead
mdl = fitlm(obsData(idx_use,2), simData(idx_use,2));
fr_r2  = mdl.Rsquared.Ordinary;

% Predicted values assuming intercept = 0 and slope = 1
y_pred = obsData(idx_use,2);

% Sum of squared residuals
SS_res = sum((simData(idx_use,2) - y_pred).^2);

% Total sum of squares
y_mean = mean(simData(idx_use,2));
SS_tot = sum((simData(idx_use,2) - y_mean).^2);

% R-squared calculation
fr_r2_constrained = 1 - (SS_res / SS_tot);

% Calculate RMSE
simulatedData  = simData(idx_use,2) ;
experimentalData = obsData(idx_use,2);
RMSE_fr = sqrt(mean((simulatedData - experimentalData).^2));

% Calculate PBIAS
sum_observed = sum(experimentalData);
sum_difference = sum(experimentalData - simulatedData);
pbias_fr = 100 * (sum_difference / sum_observed);

% Using level vs distributed q
comparison(:,1) = (validationLevel_t);
temp1 = replace_num(compV,0,NaN);
temp2 = nanmedian(temp1')';
comparison(:,2) = (0.017652+0.84506.*temp2).*totalArea_init;
perDiff(:,2)    = ((comparison(:,2) - validationDischarge)./validationDischarge).*100;
[xData,yData]   = prepareCurveData(comparison(:,1),comparison(:,2));

scatter(yData,xData,...
    'MarkerEdgeColor', [145,191,219]./255,...
    'MarkerFaceColor', [145,191,219]./255,...
    'SizeData', 20); hold on;

set(gca,'xlim',[0 162]);
set(gca,'ylim',[0 2.5]);
daspect([max(xlim)/max(ylim) 1 1]);

set(ax1,'TickLabelInterpreter','latex')
ylabel('$\textnormal{Stage [m]}$','Interpreter','LaTex');
set(ax1,'fontsize', 12)
xlabel('$\textnormal{Discharge [m\textsuperscript{3} s\textsuperscript{-1}]}$','Interpreter','LaTex');
set(ax1,'Box','on')

l1 = legend('Reference measurements', 'Constant Froude', 'Index approach',...
    'location', 'northwest',...
    'Interpreter','latex',...
    'fontsize', 11);
legend boxoff
set(l1,'Position', [get(l1,'Position')] - [0.01 0 0 0]  )

annotationText1 =  {'[A]'};
annotation(f0,'textbox',...
    [0.085677419354839 0.623147092360321 0.0264193548387096 0.0193842645381985],...
    'String',{char(annotationText1)} ,...
    'FitBoxToText','on',...
    'Interpreter','LaTex',...
    'LineStyle', 'none',...
    'fontsize', 12);



%Calculate the summary stats for when using the Distributed index output
simData = [1:length(validationDischarge);...
    comparison(:,2)']';
dist_use = simData;
obsData = [1:length(validationDischarge);...
    validationDischarge']';
idx_use = ~isnan(comparison(:,2));

% Calculate NSE
[dist_n_out] = nashsutcliffe(obsData(idx_use,1:2), simData(idx_use,1:2)); % not appropriate

% Calculate r2 instead
mdl      = fitlm(obsData(idx_use,2), simData(idx_use,2));
dist_r2  = mdl.Rsquared.Ordinary;

% Predicted values assuming intercept = 0 and slope = 1
y_pred = obsData(idx_use,2);

% Sum of squared residuals
SS_res = sum((simData(idx_use,2) - y_pred).^2);

% Total sum of squares
y_mean = mean(simData(idx_use,2));
SS_tot = sum((simData(idx_use,2) - y_mean).^2);

% R-squared calculation
dist_r2_constrained = 1 - (SS_res / SS_tot);


% Calculate RMSE
simulatedData  = simData(idx_use,2) ;
experimentalData = obsData(idx_use,2);
RMSE_dist = sqrt(mean((simulatedData - experimentalData).^2));

% Calculate PBIAS
sum_observed = sum(experimentalData);
sum_difference = sum(experimentalData - simulatedData);
pbias_dist = 100 * (sum_difference / sum_observed);


ax2 = subplot(1,2,2);
set(ax2,'fontsize', 12);
pbaspect([1 1 1]);
axis tight
set(ax2,'TickLabelInterpreter','latex');
temp1 = get(ax1,'xlim');
set(ax2,'xlim', temp1);


scatter(obsData(:,2),fr_use(:,2),...
    'MarkerEdgeColor', [252,141,89]./255,...
    'MarkerFaceColor',  [252,141,89]./255,...
    'SizeData', 20); hold on;


scatter(obsData(:,2),dist_use(:,2),...
    'MarkerEdgeColor', [145,191,219]./255,...
    'MarkerFaceColor', [145,191,219]./255,...
    'SizeData', 20); hold on;

plot(0:max(obsData(:,2)), 0:max(obsData(:,2)),...
    '--',...
    'Color', [0.5,0.5,0.5]),...
    hold on;

axis equal
set(gca,'xlim',[0 162]);
set(gca,'ylim',[0 162]);

set(ax2,'TickLabelInterpreter','latex')
xlabel('$\textnormal{Reference Discharge [m\textsuperscript{3} s\textsuperscript{-1}]}$','Interpreter','LaTex');
ylabel('$\textnormal{Image-based Discharge [m\textsuperscript{3} s\textsuperscript{-1}]}$','Interpreter','LaTex');
set(ax2,'fontsize', 12)
set(ax2,'Box','on')
yticks(xticks);

pos = get(ax2, 'Position'); % Get current position
set(ax2, 'Position', [pos(1)-0.03, pos(2), pos(3), pos(4)]); % Adjust position

annotationText2 =  {'[B]'};
annotation(f0,'textbox',...
    [0.489709677419355 0.623147092360321 0.0264193548387096 0.0193842645381985],...
    'String',{char(annotationText2)} ,...
    'FitBoxToText','on',...
    'Interpreter','LaTex',...
    'LineStyle', 'none',...
    'fontsize', 12);




%% Export the outputs
outDir    = [pathOut 'figure8.svg'];
saveas(f0,outDir,'svg');

