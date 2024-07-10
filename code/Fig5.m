% Initial code used to reproduce Figure 4 of the paper `Unsupervised image
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

%% Produce the plots
for a = 1
    comparison(:,1) = (validationDischarge./totalArea_init);
    comparison(:,2) = (comp_q_v(:,a)); % 1D xs velocity)

    forensicCells = [...
        492, 493, 501, 502, 504, 500, 505, 506]; % rmeoved 3rd, 4th
    forensicCells = forensicCells - 205; % account for removal of data prior to 1980
    

    temp  = replace_num(comp_q_v_all{a},0,NaN);
    temp2 = replace_num(validationDischarge,0,NaN);
    refQ  = temp2(forensicCells);
    [ordering, order_idx] = sort(refQ,'descend');

    for k = 1:length(forensicCells)
        [~,idx]        = min(abs(comparison(forensicCells(k),2) - temp(forensicCells(k),:)));
        dp(k,1)        = sum(~isnan(temp(forensicCells(k),:)));
        minVal(k,1)    = temp(forensicCells(k),idx);
        video_analyse(k,1) = time1(forensicCells(k),idx);
    end

    video_analyse = video_analyse(order_idx);

    for k = 1:length(forensicCells)

        fileimportname{k,1} = ['devon_dart' datestr(datenum(video_analyse(k,1),'dd/mm/yyyy HH:MM'),'yyyymmdd_HHMM') '0_discharge_summary_report.html' ];


        % bring in the data from the html file
        htmlIn = [q_path...
            fileimportname{k,1}];

        htmlDataIn  = readtext(htmlIn,'>', '','','textual');
        textstr     = '1 </div';
        popCells    = find(strcmp(htmlDataIn, textstr)==1);
        popCells    = popCells(1);
        tVal        = htmlDataIn([popCells+1, popCells + (1:19)*28+1]);
        distanceVal = cellfun(@(x) str2double(x(1:end - 6)), tVal);
        distancecor = distanceVal - min(distanceVal) + ((distanceVal(2)-distanceVal(1))/2); %=G2-MIN(G:G)+(($G$3-$G$2)/2)
        tVal        = htmlDataIn([popCells+3, popCells + (1:19)*28+3]);
        depthVal    = cellfun(@(x) str2double(x(1:end - 6)), tVal);
        tVal        = htmlDataIn([popCells+5, popCells + (1:19)*28+5]);
        areaVal     = cellfun(@(x) str2double(x(1:end - 6)), tVal);
        tVal        = htmlDataIn([popCells+7, popCells + (1:19)*28+7]);
        cubic_Val   = cellfun(@(x) str2double(x(1:end - 6)), tVal);
        tVal        = htmlDataIn([popCells+9, popCells + (1:19)*28+9]);
        quad_Val    = cellfun(@(x) str2double(x(1:end - 6)), tVal);
        tVal        = htmlDataIn([popCells+11, popCells + (1:19)*28+11]);
        fr_Val      = cellfun(@(x) str2double(x(1:end - 6)), tVal);
        tVal        = htmlDataIn([popCells+13, popCells + (1:19)*28+13]);
        velMode     = cellfun(@(x) (x(1:end - 6)), tVal,'UniformOutput',0);
        velMode     = ~strcmp(velMode,'Estimated'); % zero is estimated


        % bring in the KLT discharge
        opts            = extractFileText([htmlIn]);
        str             = extractHTMLText(opts);
        newStr          = split(str,' ');
        obsQ_KLT(k,1:3)     = [str2num(newStr{19}(2:end-1)),...
        str2num(newStr{20}(1:end-1)),...
        str2num(newStr{21}(1:end-1))];


        obsLevel     = str2double(dataIn_(2:end,2));
        obsQ         = str2double(dataIn_(2:end,4:6));
        obsVelocity  = str2double(dataIn_(2:end,3));
        obsDT        = dataIn_(2:end,1);

        if k ==1
            fileIn = [pathIn 'compiled_data_v5_20240524.xlsx'];
            [~,sheet_name]      = xlsfinfo(fileIn);
            f0                  = figure(); hold on;
            figPos              = set(f0,'Position',[2,42,958,1074]);
        end

        data{k}       = readtable(fileIn,'Sheet', sheet_name{order_idx(k)}); % this is not reading the file correctly
        %ordering      = order_idx(k);
        ax{k}         = subplot(5,2,k);

        set(ax{k},'fontsize', 12)
        set(ax{k},'TickLabelInterpreter','latex')
        set(ax{k},'Box','on')
        pbaspect([2 1 1])
        axis tight
        hold on

        tidx1 = find(string(data{k}.Properties.VariableNames) == "DMGCorrected");
        tidy1 = find(string(data{k}.Properties.VariableNames) == "MeanSpeed_m_s_");

        if ismember(k,[1,2,4]) % if stationary

            ref_data{k} = scatter( table2array(data{k}(:,tidx1)),...
                table2array(data{k}(:,tidy1)),...
                'MarkerEdgeColor', [0.5,0.5,0.5],...
                'MarkerFaceColor', [0.5,0.5,0.5],...
                'SizeData', 10);
            hold on;

            tidx2 = distancecor;
            tidy2 = fr_Val;

            % plot measured
            iv_data{k} = scatter(tidx2(velMode),...
                tidy2(velMode),...
                'MarkerEdgeColor', [0.8,0.0,0.0],...
                'MarkerFaceColor', [0.8,0.0,0.0],...
                'SizeData', 10);

            % plot estimated
            iv_data{k} = scatter(tidx2(~velMode),...
                tidy2(~velMode),...
                'MarkerEdgeColor', [0.8,0.0,0.0],...
                'MarkerFaceColor', [1,1,1],...
                'SizeData', 10);

        else % if moving


            tidx2 = distancecor;
            tidy2 = fr_Val;

            % plot measured
            iv_data{k} = scatter(tidx2(velMode),...
                tidy2(velMode),...
                'MarkerEdgeColor', [0.8,0.0,0.0],...
                'MarkerFaceColor', [0.8,0.0,0.0],...
                'SizeData', 10);

            % plot estimated
            iv_data{k} = scatter(tidx2(~velMode),...
                tidy2(~velMode),...
                'MarkerEdgeColor', [0.8,0.0,0.0],...
                'MarkerFaceColor', [1,1,1],...
                'SizeData', 10);

            %tidx3 = find(string(data{k}.Properties.VariableNames) == "corrected_dist");
            %tidy3 = find(string(data{k}.Properties.VariableNames) == "FroudeVelocity");
            %distIn_adcp = table2array(data{k}(:,tidx3));
            distIn_adcp = distancecor;
            steps = 0:distIn_adcp(1).*2:nanmax(distIn_adcp)+distIn_adcp(1)*1.2;


            tidx1 = find(string(data{k}.Properties.VariableNames) == "DMGCorrected");
            tidy1 = find(string(data{k}.Properties.VariableNames) == "MeanSpeed_m_s_");
            distIn = table2array(data{k}(:,tidx1));
            [val,idx] = sort(distIn);
            distIn = val;
            velIn  = table2array(data{k}(:,tidy1));
            velIn = velIn(idx);

            for a = 1:length(steps)-1
                idx2 = distIn > steps(a) & distIn < steps(a)+1;
                adcp_vel(a) = nanmean(velIn(idx2));
                adcp_std(a) = std(velIn(idx2));
            end

            errorbar(distIn_adcp(1:length(adcp_std)),adcp_vel,adcp_std,...
                'Marker', 'o',...
                'LineStyle', 'none',...
                'Color',[0.5,0.5,0.5],...
                'MarkerSize',3,...
                'MarkerEdgeColor', [0.5,0.5,0.5],...
                'MarkerFaceColor', [0.5,0.5,0.5])
        end

        % ensure ylims are not negative
        ax_lim_tmp = get(ax{k},'ylim');
        set(ax{k},'ylim',[0,ax_lim_tmp(2)]);

        if ~mod(k,2)
            initpos = get(ax{k}, 'Position'); % move axes closer together
            initpos (1) =  initpos (1) - 0.1;
            set(ax{k},'Position', initpos)

            t1 = get(ax{k},'yticklabels');
            set(ax{k}, 'yticklabels', {});
            yyaxis right
            set(ax{k}.YAxis(2), 'Limits', get(ax{k}.YAxis(1),'Limits'));
            set(ax{k}.YAxis(2), 'TickValues', get(ax{k}.YAxis(1),'tickValues'));
            set(ax{k}.YAxis(2), 'TickLabels', get(ax{k}.YAxis(1),'tickLabels'));

            set(ax{k}, 'yticklabels', []);
            set(ax{k}, 'yticklabels', t1);
            ax{k}.YAxis(1).Color = 'k';
            ax{k}.YAxis(2).Color = 'k';

            s1 = s1 + [0.36 0 0 0];

        elseif k == 1
            s0 = get(ax{k},'position');
            s1 = s0 - [0.03   -0.1018    0.2885    0.0971];

        else
            s1 = s1 - [0.36 0 0 0];
            s1 = s1 - [0 0.173 0 0];
        end

        subLabels = {'[A]','[B]','[C]','[D]','[E]','[F]','[G]','[H]','[I]','[J]'};
        annotationText1 = char([subLabels(k)]);
        annotation(f0,'textbox',...
            [s1],...
            'String',{annotationText1} ,...
            'FitBoxToText','on',...
            'Interpreter','LaTex',...
            'LineStyle', 'none',...
            'fontsize', 12);

        % at the point where we are trying to write teh Q value on the
        % plot. This will be achieved below:
        annotationTemp = num2str(round(refQ(order_idx(k))*100)/100) ;
        annotationText1 = [annotationTemp 'm\textsuperscript{3} s\textsuperscript{-1}'];
        if ~mod(k,2)
            movFac = [0.04 -0.005];
        else
            movFac = [0.06 -0.005];
        end
        annotation(f0,'textbox',...
            [s1 + [movFac 0 0]],...
            'String',{annotationText1} ,...
            'FitBoxToText','on',...
            'Interpreter','LaTex',...
            'LineStyle', 'none',...
            'fontsize', 11);

        % display the difference in Q measurements
        diff(k) = ((obsQ_KLT(k,3)-refQ(order_idx(k))) ./ refQ(order_idx(k))) .* 100;
        annotationTemp2 = num2str(round(diff(k)*100)/100) ;
        annotationText2 = [annotationTemp2 '\%'];
        if ~mod(k,2)
            movFac = [0.04 -0.018];
        else
            movFac = [0.06 -0.018];
        end        
        annotation(f0,'textbox',...
            [s1 + [movFac 0 0]],...
            'String',{annotationText2} ,...
            'FitBoxToText','on',...
            'Interpreter','LaTex',...
            'LineStyle', 'none',...
            'fontsize', 11);

    end

end

han                 = axes(f0,'visible','off');
han.Title.Visible   = 'on';
han.XLabel.Visible  = 'on';
han.XLabel.FontSize = 17;
han.YLabel.Visible  = 'on';
han.YLabel.FontSize = 17;

xlabel(han,'$\textnormal{Distance from left bank \ [m]}$','Interpreter','LaTex');
set(han.XLabel,'Position',get(han.XLabel,'Position') - [0.07 -0.21 0]);
set(han,'TickLabelInterpreter','latex')

ylabel(han,'$\textnormal{Depth averaged velocity [m s\textsuperscript{-1}]}$','Interpreter','LaTex');
set(han.YLabel,'Position',get(han.YLabel,'Position') - [0 -0.10 0]);

outDir    = [pathOut 'figure5.svg'];
saveas(f0,outDir,'svg');



