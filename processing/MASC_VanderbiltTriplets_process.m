% upload_Vanderbilt_triplets
clear all; close all;

label = label_params;
process = process_params;
dirname = label.campaigndir;
pic_list.files = dir(fullfile(dirname,'FLAKE*'));
pic_list.files = {pic_list.files.name}';



for i=1:length(pic_list.files) 
   
    tmp_str = pic_list.files{i};
    date_str = tmp_str(7:14);
    hour_str = tmp_str(16:22);
    pic_list.time_str{i,1} = [date_str hour_str];
    pic_list.time_num(i,1) = datenum(strcat(pic_list.time_str(i),'00'),'yyyymmddHHMMSSFFF');
    pic_list.cam(i,1) = str2num(tmp_str(end-4));
    
    if i==1
        pic_list.id(i,1) = 1;
    elseif strcmp(pic_list.time_str(i),pic_list.time_str(i-1)) %might fail if 2 flakes have the exact same timestamp!
        pic_list.id(i,1) = pic_list.id(i-1);
    else
        pic_list.id(i,1) = pic_list.id(i-1)+1;
    end
    
end
            
            
if process.parallel
    n_workers = 16; % if you have less than 16 threads matlab will take as many as possible
else
    n_workers = 0;
end

N_bad = 0;
N_blurry = 0;
N_good = 0;

% Main loop on all the pictures located in dirname
parfor (j=1:length(pic_list.id),n_workers) % length(pic_list.id)
    
        if isfield(pic_list,'files')

            fprintf('Processing %s... \n',pic_list.files{j});

        end

        [tt,flag] = MASC_picture_process(j,dirname,pic_list,label,process);
        % flags for MASC_triplet_process:
        % -1 : something went wrong
        % 0 : no roi found on one image
        % 1 : no matching triplet found
        % 2 : GOOD (at least one triplet found and saved)

        if flag == 0 

            N_bad = N_bad + 1;

        elseif flag == 1 

            N_blurry = N_blurry + 1;

        elseif flag == 2

            N_good = N_good + 1;

        end
        
end
