% === upload all Masc directories satisfying label_params ================
%
% This functions uploads directory names of all the MASC folders found in
% campaigndir which are within the interval defined by
% starthr_vec and endhr_vec
%
% This function is called in MASC_process.m
%
% Author : Christophe Praz (christophe.praz@epfl.ch)
%
% Last Update : April 2016
%
% ========================================================================
function dir_list = uploaddirs(campaigndir,starthr_vec,endhr_vec)



    alldir_names = dir(campaigndir);
    alldir_names = {alldir_names.name};
    dir_list = {};
    
    %loop around the folders
     for i=1:length(alldir_names)
      
         %check if folder format is yyyy.mm.dd        
         try 
             
             folder_datenum = datenum(alldir_names{i},'yyyy.mm.dd');     
             t_min = datenum(starthr_vec);
             t_max = datenum(endhr_vec);
             
             if folder_datenum >= floor(t_min) && folder_datenum <= t_max
                 
                 allsubdir_names = dir(fullfile(campaigndir,alldir_names{i}));
                 allsubdir_names = {allsubdir_names.name};
                 
                 for j=1:length(allsubdir_names)
                     
                     try 
                          
                         timevec = datevec(folder_datenum);
                         timevec(1,4) = str2num(allsubdir_names{j});
                         sub_folder_datenum = datenum(timevec);
                         
                         % t_min approximated rounded down at the hour
                         t_min_vec = datevec(t_min);
                         t_min_vec(:,5) = 0;
                         t_min_vec(:,6) = 0;
                         t_min = datenum(t_min_vec);
                         
                         if sub_folder_datenum >= t_min && sub_folder_datenum <= t_max
                             
                             dir_content = dir(fullfile(campaigndir,alldir_names{i},allsubdir_names{j},'*.png'));
                             if ~isempty(dir_content) % -> at least one image in the folder
                                dir_list{end+1} = fullfile(campaigndir,alldir_names{i},allsubdir_names{j});
                                disp(strcat(datestr(sub_folder_datenum),' added to the dir. list'));
                             end
                             
                         end
                         
                     catch err             
                         
                     end
                     
                 end
                 
             end
             
         catch err
               
         end
         
     end
     
end



% function dirlist = uploaddirs(dirall,dirlistall,starthr,endhr)
% %UPLOAD Uploads directory names of the MASC fallspeed and image filenames to be processed
% %   dirlist = uploaddirs(dirlistall,starthr,endhr) where inputs are specified in label_params 
% %   dirall is the main directory, dirlistall lists all directories in the main directory, and starthr and endhr are in
% %   the format [yyyy mm dd hh].
% %
% %   Copyright Tim Garrett, University of Utah. This code is freely available for
% %   non-commercial distribution and modification
% 
% disp('Uploading list of directories that contain flakes')
% disp('This may take some time...')
% %%% The following parses the names of the directories for timestamp content
% dirtime = zeros(length(dirlistall),1); %Directory hours
% flakes = zeros(length(dirlistall),1); %Directories that contain flakes
% 
% parfor i = 1:length(dirlistall)
%     temp = cell2mat(dirlistall(i));
%     flakedir = cell2mat(strcat(dirall,'/',dirlistall(i)));
%   
% 
% % Add to the following list as need be if getting an error in loaddata.m
% %%omit '.' and '..', '.mat' .ps .pdf .zip etc.
%     if strcmp(temp(1),'.') == 1 | ...
%        isempty(strfind(temp,'.mat')) ~= 1 | ...
%        isempty(strfind(temp,'.ps')) ~= 1 | ...
%        isempty(strfind(temp,'.pdf')) ~= 1 |...
%        isempty(strfind(temp,'.zip')) ~= 1;
%        continue
%     end
%     r = regexp(temp, '_', 'split' );
%     m = regexp( r{1}, '\.','split' );
%     d = regexp( r{2}, '\.','split' );
%     t = regexp( r{4}, '\.','split' );
%     
%     dirMASC = str2num(m{1});
%     diryr = str2num(d{1});
%     dirmonth = str2num(d{2});
%     dirday = str2num(d{3});
%     dirhr = str2num(t{1});
%     
%     dirtime(i) = datenum([diryr dirmonth dirday dirhr 0 0]);
%     
%         
%     %%%%% Delete extraneous files in each directory that begin with '.' %%%%%
%     dircon = dir(flakedir);
%     dircon = {dircon([dircon.isdir] == 0).name};
%     dotfiles = strcat(flakedir,'/',dircon(~cellfun(@isempty, regexpi(dircon, '.[^.}*'))));
%     if length(dotfiles) > 0
%         for j = 1:length(dotfiles)
%             delete(dotfiles(j));
%         end
%     end
%     
%     %Skip empty directories or those with insufficient images for a single
%     %flake
%     
%     dircon = dir(flakedir);
%     dircon = {dircon([dircon.isdir] == 0).name};
%     dotfiles = strcat(flakedir,'/',dircon(~cellfun(@isempty, regexpi(dircon, '.[^.}*'))));
%     
%     if length(dircon)<7 
%         flakes(i) = 0; %%% Empty directory
%     else
%         flakes(i) = 1; %% Directory with flakes
%     end
%        
%     
% end
% 
% disp('Proportion of hours that contain flakes')
% disp(length(find(flakes == 1))/length(dirlistall))
% %%%Find directories that are not empty
% 
%     %use only directories between starthr and endhr that contain flakes
%     good = find(dirtime >= datenum([starthr 0 0]) & dirtime <= datenum([endhr 0 0]) & flakes == 1);
%     dirlist = strcat(dirall,'/',dirlistall(good));
