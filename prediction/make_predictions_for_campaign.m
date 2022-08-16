% ========== main script used to classify processed MASC images in batch (.mat files) ==========
%
% inputs :
% ========
% campaingdir  : path to the folder where the processed MASC .mat files are located
% t_vec_start  : start time vector [yyyy mm dd HH MM SS] (MASC images recorded before this will be ignored)
% t_vec_stop   : stop time vector  [yyyy mm dd HH MM SS] (MASC images recorded before this will be ignored)
% classifiers  : a structure containing the path the each of the 3 classifiers used for hydrometeor class, degree of riming and melting snow detection, respectively
%   by default, use the following classifiers :
%   -      classifiers.class = 'classification/logit_trained_models/logit_v.1.1_FINAL3_6classes_weighted_scheme.mat';
%   -      classifiers.riming = 'classification/logit_trained_models/logit_v.1.1_FINAL3_riming_5classes_weighted_scheme.mat';
%   -      classifiers.melting = 'classification/logit_trained_models/logit_v.1.1_FINAL3_melting_2classes_noweight.mat'; 
% parallel_proc : boolean value (0/1) indicating whether the user want to use Matlab parallel processing toolbox to gain computing time (not suitable for all computers)
%
% outputs :
% =========
% none          : the script will directly write the label, riming degree and indicator of melting snow in the corresponding .mat file. Probabilities are also stored.
%
% Author : Christophe Praz (christophe.praz@epfl.ch)
% Last update : March 2018
% =============================================================================================
function make_predictions_for_campaign(campaigndir,t_vec_start,t_vec_stop,classifiers,parallel_proc)

t_str_start = datestr(t_vec_start,'yyyymmddHHMMSS');
t_str_stop = datestr(t_vec_stop,'yyyymmddHHMMSS');
snowflakes_dirs = uploaddirs(campaigndir,t_vec_start,t_vec_stop);
params_file_dir = campaigndir;
save_results = true;
if parallel_proc
    n_workers = 16; % if you have less than 16 threads matlab will take as many as possible
else
    n_workers = 0;
end


tic;
parfor (i=1:length(snowflakes_dirs),n_workers)
    
    fprintf('%s ... \n',snowflakes_dirs{i});
    %process_new_descriptors(snowflakes_dirs{i},t_str_start,t_str_stop)
    predict_snowflakes_class(snowflakes_dirs{i},t_str_start,t_str_stop,classifiers.class,save_results,0);
    predict_snowflakes_riming(snowflakes_dirs{i},t_str_start,t_str_stop,classifiers.riming,save_results,0);
    predict_snowflakes_melting(snowflakes_dirs{i},t_str_start,t_str_stop,classifiers.melting,save_results,0);
      
end

if save_results
    
    fileID = fopen(fullfile(params_file_dir,'proc_params.txt'),'a');
    fprintf(fileID,'\n%s: Classification performed from %s to %s : \n',datestr(now,'dd.mm.yyyy HH:MM.SS'),datestr(t_vec_start,'dd.mm.yyyy HH:MM.SS'),datestr(t_vec_stop,'dd.mm.yyyy HH:MM.SS'));
    fprintf(fileID,'Class : %s \nRiming : %s \nMelting : %s \n',classifiers.class,classifiers.riming,classifiers.melting);
    fclose(fileID);

end




end

%%%
%classifiers.class = 'classification/logit_trained_models/logit_v.1.1_FINAL3_6classes_weighted_scheme.mat';
%%%'classification/logit_trained_models/logit_v.1.1_extended_dataset_6classes_weighted_scheme.mat'; % FINAL3 seems to underestimate crystals towards more aggregates
%classifiers.riming = 'classification/logit_trained_models/logit_v.1.1_FINAL3_riming_5classes_weighted_scheme.mat'; % FINAL3 def. better!
%classifiers.melting = 'classification/logit_trained_models/logit_v.1.1_FINAL3_melting_2classes_noweight.mat'; % tbd