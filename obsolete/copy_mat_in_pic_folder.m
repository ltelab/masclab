% small script to copy flake data .mat in the same folder as some selected
% pictures

%% loading data
clear all; close all;

dir_data = '/media/praz/Masc-Data/DATA_meteoswiss/PROCESSED_20151109/DATA/GOOD';
dir_img  = '/home/praz/Documents/MASC/sample_clustering';

subdir_img = {'capped_columns','dry_aggregates','graupels','pristines','rimed_aggregates','rimed_pristines','sleet'};

for i=1:length(subdir_img)
    
    path_img = fullfile(dir_img,subdir_img{i});
    img_list = dir(fullfile(path_img,'*.png'));
    img_list = {img_list.name};
    
    % find and copy mat files for each subdir
    for j=1:length(img_list)
        
        data_file = img_list{j};
        data_file = strcat(data_file(1:end-3),'mat');
        data_fullfile = fullfile(dir_data,data_file);
        if exist(data_fullfile,'file')
            copyfile(data_fullfile,fullfile(dir_img,subdir_img{i},data_file));
        else
            sprintf('error : impossible to find %s \n',data_fullfile);
        end
        
    end
    
end


% file_list = dir(fullfile(dir_data,'*.mat'));
% file_only_list = {file_list.name};
% file_list = fullfile(dir_data,file_only_list);
% 
% img_list = dir(fullfile(dir_img,'*.png'));
% img_only_list = {img_list.name};
% img_list = fullfile(dir_img,img_only_list);