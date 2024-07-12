% Load the KLT-IV generated outputs in html format
function [obsLevel, obsQ, obsVelocity, obsDT] = parsingHTML(pathIn,q_path)

listing = dir(q_path);
for a = 1:length(listing)
    listingName(a,1) = cellstr(listing(a).name);
end
listingName = listingName(contains(listingName,'.html'));

% extract data from html format
for a = 1:length(listingName)
    opts            = extractFileText([q_path listingName{a,1}]);
    str             = extractHTMLText(opts);
    newStr          = split(str,' ');
    obsLevel(a,1)   = str2num(newStr(5));
    obsArea(a,1)    = str2num(newStr(8));
    obsQ(a,1:3)     = [str2num(newStr{19}(2:end-1)),...
        str2num(newStr{20}(1:end-1)),...
        str2num(newStr{21}(1:end-1))];
    obsDT(a,1)      = datenum(listingName{a}(11:23),'yyyymmdd_HHMM');
end

% match the distributed velocities from KLT csv to the correct video time
KLTin        = [pathIn 'cal_val_stage_discharge_combined.csv'];
dataIn_      = readtext(KLTin,',', '','','textual');
obsQ_        = str2double(dataIn_(2:end,4:6));
obsVelocity_ = str2double(dataIn_(2:end,3));
obsDT_       = dataIn_(2:end,1);

for a = 1:length(obsQ)
    temp            = find(sum(abs(obsQ(a,1:3) - obsQ_(:,1:3)),2) == 0);
    if length(temp) == 1
        diff(a,1)   = temp;
    else
        temp2       = obsDT(a);
        temp3       = datenum(obsDT_(temp),'dd/mm/yyyy HH:MM');
        [~, u]      = min(abs(temp3 - temp2));
        diff(a,1)   = temp(u);
    end
end

obsVelocity         = obsVelocity_(diff);