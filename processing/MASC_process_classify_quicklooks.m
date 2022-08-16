% ==================== CORE MATLAB ROUTINE ====================
% Script used to process MASC images (identify and crop around the "best" snowflake triplet and calculate various geometrical and
% textural descriptors for each view), classify each image according to the classification scheme proposed in Praz et al., AMT 2017.
% Finally, the user can also choose to generate daily quicklooks of the classification output 
%
% The user must enter the user parameters (lines ~10-50) and provide a "process_params_file" text file containing MASC processing parameters.
%
% Author : Christophe Praz (christophe.praz@epfl.ch)
% Last update : March 2018
% =============================================================
try

    clearvars;

% USER-PARAMETER %%%%% ----do not modify this comment line

    % part 1 : image processing and cropping
    masc_process = true;
    masc_regen_from_start = true;
    process_params_file = 'userspecs/params_files/process_params_ICEGENESIS_2021.txt'; % -> path to file containing processing parameters. Examples are given in masclab/userspecs/params_files/*.txt
    % just use process_params_DEFAULT.txt if you do not know how to fill the parameters
    label_params.campaigndir = '/ltedata/ICEGENESIS_2021/MASC/Raw_data'; % -> path to directory where raw MASC data is located.
    label_params.output_dir = '/ltedata/ICEGENESIS_2021/MASC/Proc_data'; % -> path to directory where processed MASC data will be stored. 
    
    % time interval you want to process
    label_params.starthr_vec = [2021 01 01 00 00 00];
    label_params.endhr_vec   = [2021 01 04 00 00 00]; %datevec(datenum(date)+datenum(days(1))); 

    % part 2 : classification
    classif_process = true;
    classif_regen_from_start = true;
    classif_parallel_proc = true;

    % part 3 : quicklooks
    gen_quicklooks = true;
    regen_quicklooks = true;
    options.savefigs = true;
    options.savepath = '/ltedata/ICEGENESIS_2021/MASC/Quicklooks';
    options.pixres = 33.5/1000; % MASC pixel resolution in mm/pixel
    options.xi_thresh = 9; % quality parameter threshold
    options.Nmin_interval = 30; % time window (min)
    options.Nmin_shift = 10;
    options.Nclasses_masc = 6;
    options.MASC_classes = {'SP','CC','PC','AG','GR','CPC'};
    options.MASC_classes_desired = [1 2 3 4 5 6];
    options.N_MascSamples_min = 0;
    options.use_triplet = false; % whether we treat each image independently or by triplet
    options.use_quintuplet = false;
    options.verbose = 0;
    options.OR_180 = true;

% END-USER-PARAMETER -do not modify this comment line
    
    starthr_vec_ini = label_params.starthr_vec;

    % part 1 : image processing and cropping
    if masc_process

        label_params.outdir = label_params.output_dir;

        if ~masc_regen_from_start
            proc_pic_list = dir(fullfile(label_params.outdir,'**','*.png'));
            if ~isempty(proc_pic_list)
                proc_pic_list = {proc_pic_list.name};
                proc_last_pic = proc_pic_list{end};
                tnum_last_pic = datenum(proc_last_pic(1:19),'yyyy.mm.dd_HH.MM.SS');
                label_params.starthr_vec = datevec(tnum_last_pic); 
            end
        end

        fprintf('Part I : processing MASC images from %s to %s \n.',datestr(label_params.starthr_vec,'yyyy.mm.dd HH:MM'),datestr(label_params.endhr_vec,'yyyy.mm.dd HH:MM'));
        process_params = process_params(process_params_file);
        cam_params = cam_params();
        MASC_process(label_params, cam_params, process_params);

    end

    % part 2 : hydrometeor classification 
    if classif_process

        classif_dir = label_params.output_dir;
        if classif_regen_from_start
            classif_starthr_vec = [2010 01 01 00 00 00];
        else
            classif_starthr_vec = label_params.starthr_vec;
        end
        classif_endhr_vec= label_params.endhr_vec;
        classifiers.class = 'classification/logit_trained_models/logit_v.1.1_FINAL3_6classes_weighted_scheme.mat';
        classifiers.riming = 'classification/logit_trained_models/logit_v.1.1_FINAL3_riming_5classes_weighted_scheme.mat';
        classifiers.melting = 'classification/logit_trained_models/logit_v.1.1_FINAL3_melting_2classes_noweight.mat'; 
        make_predictions_for_campaign(classif_dir,classif_starthr_vec,classif_endhr_vec,classifiers,classif_parallel_proc);

    end

    % part 3 : quicklooks
    if gen_quicklooks

        % display
        disp.mode = 'Off';
        disp.now = true;
        disp.ht_props = false;
        disp.ht_piechart = false;
        disp.riming = false;
        disp.melting = false;
        disp.n_MASC = false;
        disp.overview_classif = true;
        disp.overview_microstruct = true;

        day_stop = datetime(label_params.endhr_vec);
        day_stop.TimeZone = 'UTC';
        %day_stop = datetime('now','TimeZone','UTC');
        day_stop_midnight_vec = datevec(day_stop);
        day_stop_midnight_vec = day_stop_midnight_vec(1:3);
        if regen_quicklooks
            day_ini = datetime(starthr_vec_ini(1:3));
        else
            day_ini =  datetime(day_stop_midnight_vec) - days(1);%datetime(label_params.starthr_vec(1:3));
        end
        day_ini.TimeZone = 'UTC';
        day_current = [day_ini:days(1):day_stop];

        for i=1:numel(day_current)

            data = merge_predictions_for_campaign(label_params.output_dir,datestr(day_current(i),'yyyymmddHHMMSS'),datestr(day_current(i)+days(1),'yyyymmddHHMMSS'));
            if options.savefigs
                options.savedir = fullfile(options.savepath,datestr(day_current(i),'yyyy_mm'));
                if ~exist(options.savedir,'dir')
                    mkdir(options.savedir);
                end      
            end
            make_masc_time_series(data,datevec(day_current(i)),datevec(day_current(i)+days(1)),options,disp);

        end

    end

catch err
    
    fprintf('Error in %s at line %u : %s\n',err.stack.file,err.stack.line,err.message);
    
    
    
end




