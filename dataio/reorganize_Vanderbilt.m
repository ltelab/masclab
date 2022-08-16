clear all; close all;

path2dir = '/media/praz/MASC_Vanderbilt/RAW';
path2copy = '/media/praz/MASC_Vanderbilt/part2';
alldir_names = dir(fullfile(path2dir,'MASC*'));
alldir_names = {alldir_names.name};

for i=1:100%length(alldir_names)

    if ~isempty(dir(fullfile(path2dir,alldir_names{i},'*.png*')))
    
        folder_date = alldir_names{i}(6:15);
        folder_hour = alldir_names{i}(17:18);

        if isempty(dir(fullfile(path2copy,folder_date,folder_hour)))
            
            mkdir(fullfile(path2copy,folder_date,folder_hour));
            cd(fullfile(path2dir,alldir_names{i}));
            system(sprintf('find -name "*" -exec cp -t %s {} \\;',fullfile(path2copy,folder_date,folder_hour)));
%                 %system('find -name "*.m" -exec cp {} target \;')
%                 %find -maxdepth 1 -name '*.prj' -exec mv -t ../prjshp {} +
            fprintf('Copying content From : %s \n',fullfile(path2dir,alldir_names{i}));
            fprintf('To : %s \n',fullfile(path2copy,folder_date,folder_hour));
            fprintf('********************************************\n');
            
        end
        
    end
    
end


