% script to compare perimeter/area of flakes with or without "cleaning
% boarders"

clear all; close all;

path_cleaned = '/home/praz/Documents/MASC/sample_snow_20150620/sample_cool/PROCESSED_with_cleaning_v2/DATA/GOOD/';
cleaned_names = dir(fullfile(path_cleaned,'*.mat'));
cleaned_names = {cleaned_names.name};

path_not_cleaned = '/home/praz/Documents/MASC/sample_snow_20150620/sample_cool/PROCESSED_wo_cleaning/DATA/GOOD/'; 
uncleaned_names = dir(fullfile(path_not_cleaned,'*.mat'));
uncleaned_names = {uncleaned_names.name};

for i=1:length(uncleaned_names)
    if ~isempty(find(strcmp(cleaned_names,uncleaned_names{i})))
        roi1 = load(fullfile(path_not_cleaned,uncleaned_names{i}));
        roi2 = load(fullfile(path_cleaned,uncleaned_names{i}));
        
        diff_perim(i) = length(roi1.roi.x_perim) - length(roi2.roi.x_perim);  
        error_perim(i) = abs(diff_perim(i)/length(roi1.roi.x_perim))*100;
    
        diff_area(i) = roi1.roi.area - roi2.roi.area;
        error_area(i) = abs(diff_area(i)/roi1.roi.area)*100;
        
        if error_perim(i) > 5
            fprintf('%s \n',roi1.roi.name);
        end
    end
    
end

figure;
subplot(211);
histogram(diff_area);
title('diff area [pix]');
subplot(212);
histogram(error_area);
title('error area [%]');
figure;
subplot(211);
histogram(diff_perim);
title('diff perim [pix]');
subplot(212);
histogram(error_perim);
title('error perim [%]');

