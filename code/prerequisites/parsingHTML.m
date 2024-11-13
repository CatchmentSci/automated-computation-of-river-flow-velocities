% Load the KLT-IV generated outputs in html format
function [obsLevel, obsQ, obsDT] = parsingHTML(q_path)

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
