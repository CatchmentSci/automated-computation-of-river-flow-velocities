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
    [obsLevel, obsQ, obsVelocity, obsDT] = parsingHTML(pathIn, q_path);
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
        FieldName1                           = ['time_num'];
        timeQ.(FieldName1)(a,1:length(idx))  = datenum(obsDT(idx,c), 'dd/mm/yyyy HH:MM');        
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

for a = 1:3
    clear comparison

    comparison(:,1) = (validationDischarge./totalArea_init);
    comparison(:,2) = (comp_q_v(:,a)); % 1D xs velocity)

    [xData,yData]   = prepareCurveData(comparison(:,2),comparison(:,1)); % These have been switched

    % Do the statistical analysis - linear regression
    xData(1:length(xData),2)    = 1; % Prep the data for regression
    [b,~,~,~,stats]             = regress(yData,xData); % Run the regression
    Original_r2                 = stats(1,1); % Pull out the R2 value
    Original_p                  = stats(1,3); % Pull out the p value
    aCoefficient                = b(2,1); % Pull out the a coefficient
    bCoefficient                = b(1,1); % Pull out the b coefficient
    intervals                   = min(xData(:,1)):0.01:max(xData(:,1));

    % Alternative method - Set up fittype and options.
    ft                  = fittype( 'a+b*(x)', 'independent', 'x', 'dependent', 'y' );
    opts                = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display        = 'Off';
    opts.StartPoint     = [0 1.00];
    [fitresult, gof]    = fit( xData(:,1), yData, ft, opts ); % Fit model to data.

    if a == 1
        f0=figure(1); hold on;
        f0.Units='pixels';
        set(f0,'Position',[100, 0, 2480./2, 3508./2]); % A4 aspect ratio
        f0.Units='normalized';
    end

    sub = subplot(3,1,a);
    hold on;
    pbaspect([1 1 1])
    ax = gca();
    set(ax,'fontsize', 12)
    axis tight

    scatter(xData(:,1),yData, ...
        'MarkerEdgeColor', [0.5,0.5,0.5],...
        'MarkerFaceColor', [0.5,0.5,0.5],...
        'SizeData', 10);

    forensicCells = [...
        492, 493, 497, 499, 501, 502, 504, 500, 505, 506];

    forensicCells = forensicCells - 205; % account for removal of data prior to 1980

    plot_red = [comparison(forensicCells,2),comparison(forensicCells,1)];

    temp  = replace_num(comp_q_v_all{a},0,NaN);
    temp2 = replace_num(validationDischarge,0,NaN);
    for aa = 1:length(forensicCells)
        [~,idx]         = min(abs(comparison(forensicCells(aa),2) - temp(forensicCells(aa),:)));
        dp(aa,1)        = sum(~isnan(temp(forensicCells(aa),:)));
        minVal(aa,1)    = temp(forensicCells(aa),idx);
        refQ(aa,1)      = temp2(forensicCells(aa));
        video_analyse(aa,1) = time1(forensicCells(aa),idx);

        if isa(video_analyse(aa,1),'double')
            fileimportname{aa,1} = ['devon_dart' datestr(video_analyse(aa,1),'yyyymmdd_HHMM') '0_discharge_summary_report.html' ];
        else
            fileimportname{aa,1} = ['devon_dart' datestr(datenum(video_analyse(aa,1),'dd/mm/yyyy HH:MM'),'yyyymmdd_HHMM') '0_discharge_summary_report.html' ];
        end

        % bring in the data from the html file
        htmlIn = [q_path ...
            fileimportname{aa,1}];

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


        %obsLevel     = str2double(dataIn_(2:end,2));
        %obsQ         = str2double(dataIn_(2:end,4:6));
        %obsVelocity  = str2double(dataIn_(2:end,3));
        %obsDT        = dataIn_(2:end,1);

    end


    scatter(minVal(:,1),plot_red(:,2), ...
        'MarkerEdgeColor', [0.8,0,0],...
        'MarkerFaceColor', [0.8,0,0],...
        'SizeData', 10);

    plot(intervals,intervals,'k--');

    plot(intervals,aCoefficient+bCoefficient.*intervals,'k-')
    set(ax,'TickLabelInterpreter','latex')
    if a == 3
        xlabel('$\mathrm{\tilde{U}_{xs} \ [m \ s^{-1}]}$','Interpreter','LaTex');
    end

    if a == 2
        ylabel('$\mathrm{U_{a} \ [m \ s^{-1}]}$','Interpreter','LaTex');
    end

    if a == 1
        s1 = [0.368967741935484 0.886915340595894 0.12399641122175 0.0315107493246396];
    elseif a == 2
        s1 = [0.368967741935484 0.646286652431097 0.126390047546268 0.0339302952919482];
    else
        s1 = [0.368967741935484 0.404985089736709 0.126390047546268 0.0339302952919482];
    end

    annotationText1 = ['$\mathrm{U_{a} = ' num2str(round(aCoefficient*1000)/1000) ' + ' num2str(round(bCoefficient*1000)/1000) ' \ \tilde{U}_{xs}' '}$'];
    annotationText2 = ['$\mathrm{R^{2} =' num2str(round(Original_r2*100)/100) '}$'];
    annotation(f0,'textbox',...
        s1,...
        'String',{annotationText1, annotationText2} ,...
        'FitBoxToText','on',...
        'Interpreter','LaTex',...
        'LineStyle', 'none',...
        'fontsize', 11);

    annotationText3 =  {'[A]', '[B]', '[C]'};
    s2 = s1 - [0.05 -0.01 0 0];
    annotation(f0,'textbox',...
        [s2],...
        'String',{char(annotationText3(a))} ,...
        'FitBoxToText','on',...
        'Interpreter','LaTex',...
        'LineStyle', 'none',...
        'fontsize', 12);

    xLims       = [0 2.55];
    yLims       = [0 2.55];
    set(ax, 'XLim', xLims)
    set(ax, 'YLim',yLims )

    if a < 3
        set(ax,'XTick',0.5:0.5:xLims(2), 'XTickLabel',[])
    else
        set(ax,'XTick',(0.5:0.5:xLims(2)),...
            'XTickLabel',(0.5:0.5:xLims(2)))
    end

    if a == 2
        set(ax,'Position',[0.13 0.470899227292163 0.775 0.215735294117647]);
    end

    if a == 3
        set(ax,'Position',[0.13 0.231495327102803 0.775 0.215735294117647]);
    end

    set(ax,'YTick',(0:0.5:yLims(2)),...
        'YTickLabel',(0:0.5:yLims(2)))
    %set(AxesHandle,'Position');

    grid on;

end

outDir    = [pathOut 'Fig4.svg'];
saveas(f0,outDir,'svg');

% this section deals with the temporal variability of reconstructions for a
% given flow stage. Four randomly selected stages are chosen and the
% Q estimates across the experimental period are plotted.
reviewer_test = 0;
if reviewer_test == 1

    f0 = figure();
    set(f0, 'Position', [ 1000         425         961         813]  )
    hold on;
    testing_idx = [102, 38, 76, 95];
    timeQ.time_num = replace_num(timeQ.time_num,0,NaN);

    for a = 1:4

        nexttile
        hold on;
        pbaspect([1 1 1])
        ax = gca();

        set(ax,'fontsize', 12)
        set(ax,'TickLabelInterpreter','latex')
        %axis tight

        scatter(timeQ.time_num(testing_idx(a),:), compQ.q3(testing_idx(a),:), '+k' );

        xlim([min(timeQ.time_num(:)) max(timeQ.time_num(:))])
        lims = xlim;
        datetick('x', 'mmm yy')
        xticks(lims(1) :(lims(2) - lims(1))./5: lims(2)    )
        xlabel('Date','Interpreter','LaTex');
        ylabel('$\mathrm{Q \ [m^{3} \ s^{-1}]}$','Interpreter','LaTex');


        annotation(f0,'textbox',...
            [0.0900915712799173 0.90218992791082 0.0432420794240754 0.0359163588116823],...
            'String','A)' ,...
            'FitBoxToText','on',...
            'Interpreter','LaTex',...
            'LineStyle', 'none',...
            'fontsize', 12);

        annotation(f0,'textbox',...
            [0.544826222684706 0.907109977111312 0.0432420794240754 0.0359163588116823],...
            'String','B)' ,...
            'FitBoxToText','on',...
            'Interpreter','LaTex',...
            'LineStyle', 'none',...
            'fontsize', 12);

        annotation(f0,'textbox',...
            [0.0900915712799173 0.42248513086285 0.0432420794240754 0.0359163588116823],...
            'String','C)' ,...
            'FitBoxToText','on',...
            'Interpreter','LaTex',...
            'LineStyle', 'none',...
            'fontsize', 12);

        annotation(f0,'textbox',...
            [0.544826222684706 0.42248513086285 0.0432420794240754 0.0359163588116823],...
            'String','D)' ,...
            'FitBoxToText','on',...
            'Interpreter','LaTex',...
            'LineStyle', 'none',...
            'fontsize', 12);
    end


    outDir    = ['D:\OneDrive - Newcastle University\Documents - TENDERLY Archive [Geog]\General\Archive_Dart\Paper Drafts\Reviewers_comments_2025_v1\discharge_test.png'];
    exportgraphics(f0,outDir,'Resolution',600, 'BackgroundColor','none')
end
