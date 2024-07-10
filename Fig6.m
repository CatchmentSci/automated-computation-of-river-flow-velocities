% Initial code used to reproduce Figure 6 of the paper `Unsupervised image
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

%% Select one PIV/ADCP measurement pair 100,000 times and calculate the conversion coefficient 
iter            = 100000;
user1           = find(~isnan(yData(:,1))==1);
xData           = xData(user1,:);
yData           = yData(user1,:);
msize           = numel(xData);

for a = 1:50%length(user1)
    for b = 1:iter
        warning('off')
        idx                      = randperm(msize, msize); % randomize order of xData collection
        xData_temp               = xData(idx,:);
        yData_temp               = yData(idx,:);
        
        unique_xData             = xData_temp;
        [~,l1]                   = size(time1);
        xData_temp               = repmat(xData_temp,1,l1);
        xData_temp               = xData_temp(:);
        yData_temp               = yData_temp(:);
        
        [xData_temp, yData_temp] = prepareCurveData (xData_temp, yData_temp);
        
        matched_idx              = find(xData_temp == unique_xData(a)); % all idx for a given val measurement
        y_out                    = yData_temp(matched_idx(randi(length(matched_idx),1,1)'));
        x_out                    = xData_temp(matched_idx(randi(length(matched_idx),1,1)'));
        temp1(b,a)               = y_out - x_out;
        clear xData_temp yData_temp
    end
    differ(:,a)         = (temp1(:,a) ./ x_out) .* 100;
    collatedDiffs(:,a)  = sum(differ(:,1:a),2)./a;
end
    
f0=figure(1); hold on;
f0.Units='pixels';
%set(f0,'Position',[3162, 1162, 1718, 1314]); % A4 aspect ratio
%f0.Units='normalized';
hold on
pbaspect([1 1 1])
ax = gca();
set(ax,'fontsize', 14)
axis tight


boxplot(collatedDiffs,...
    'plotstyle','compact',.... % compact style
    'colors', [0.5,0.5,0.5],... % grey color
    'symbol',''); % remove outliers
%clear collatedDiffs

set(ax,'TickLabelInterpreter','latex')
xticks(0:5:length(xData));
set(ax,'Xticklabel',[]);
xticklabels(0:5:length(xData));
set(ax, 'XLim', [0, get(ax, 'XLim') * [0; 1]]);
set(ax, 'YLim', [-30, 70]);

set(ax,'Box','on')
set(ax,'fontsize', 12)
xlabel('$\mathrm{Number \ of \ observations}$','Interpreter','LaTex');
ylabel('$\frac{1}{n} \sum_{i=1}^{n} \ (\bar{U}_{s_i} - U_{a_i}) / U_{a_i}  \cdot 100 \ [\%]$','Interpreter','LaTex') 



%% Export the outputs
outDir    = [pathOut 'figure6.svg'];
saveas(f0,outDir,'svg');
