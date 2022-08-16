% === Reorganize Masc folders hierarchy ==================================
%
% Initially the Masc saves the pictures in a slightly messy way creating
% a new folder every hour and a pointless subfolder called "masc_%d_Hr_%h"
% where the pictures and the text files are located.
%
% This script reorganizes the folders with respect to the following
% hierarchy : CAMPAIGN_NAME/DATE(yyyy.mm.dd)/hour(HH)/ which is the format
% Meteoswiss adopted.
% 
% Empty folders (only txt files but no pictures) are ignored, ie. not
% copied
%
% Author : Christophe Praz (christophe.praz@epfl.ch)
% Last update : September 2015
% ========================================================================

%script to reorganize MASC data in a less messy way
function reorganize_MASC_folders(path2dir,path2copy,type)

    if strcmp(type,'EPFL')

        %path2dir = '/home/praz/Documents/MASC/sample_rain_20150615_T_2.5';
        alldir_names = dir(fullfile(path2dir,'*BASE*'));
        alldir_names = {alldir_names.name};
        annoying_folder_name = 'masc_%d_Hr_%h';
        dir_list = {};

        for i=1:length(alldir_names)

            if ~isempty(dir(fullfile(path2dir,alldir_names{i},annoying_folder_name,'*.png*')))

                folder_date = alldir_names{i}(7:16);
                folder_hour = alldir_names{i}(21:22);

                if isempty(dir(fullfile(path2copy,folder_date,folder_hour)))

                    mkdir(fullfile(path2copy,folder_date,folder_hour));
                    %copyfile(fullfile(path2dir,alldir_names{i},annoying_folder_name,'*'),fullfile(path2dir,folder_date,folder_hour));

                    cd(fullfile(path2dir,alldir_names{i},annoying_folder_name));
                    system(sprintf('find -name "*" -exec cp -t %s {} \\;',fullfile(path2copy,folder_date,folder_hour)));
                    %system('find -name "*.m" -exec cp {} target \;')
                    %find -maxdepth 1 -name '*.prj' -exec mv -t ../prjshp {} +
                    fprintf('Copying content From : %s \n',strcat(alldir_names{i},'/',annoying_folder_name));
                    fprintf('To : %s \n',fullfile(path2copy,folder_date,folder_hour));
                    fprintf('********************************************\n');

                end

            end

        end
    
    elseif strcmp(type,'Vanderbilt')
        
        %path2dir = '/media/praz/MASC_Vanderbilt/RAW';
        %path2copy = '/media/praz/MASC_Vanderbilt/part2';
        alldir_names = dir(fullfile(path2dir,'MASC*'));
        alldir_names = {alldir_names.name};

        for i=1:length(alldir_names)

            if ~isempty(dir(fullfile(path2dir,alldir_names{i},'*.png*')))

                folder_date = alldir_names{i}(6:15);
                folder_hour = alldir_names{i}(17:18);

                if isempty(dir(fullfile(path2copy,folder_date,folder_hour)))

                    mkdir(fullfile(path2copy,folder_date,folder_hour));
                    cd(fullfile(path2dir,alldir_names{i}));
                    system(sprintf('find -name "*" -exec cp -t %s {} \\;',fullfile(path2copy,folder_date,folder_hour)));
                    fprintf('Copying content From : %s \n',fullfile(path2dir,alldir_names{i}));
                    fprintf('To : %s \n',fullfile(path2copy,folder_date,folder_hour));
                    fprintf('********************************************\n');

                end

            end

        end
        
    end

end