% === Labels for Masc data processing ===============
%
% include the path to the images to process, the time interval to process
% and the output directory
%
% Adapted from Tim Garrett, University of Utah. 
% Author : Christophe Praz (christophe.praz@epfl.ch)
% Last update : August 2015
% ===================================================
function label_params = label_params()

    %input dir
    label_params.campaigndir = '/data2/MASC/Davos_test_LTE'; 
    label_params.output_dir = '/data2/MASC/processed_Davos_test';  
    % time interval to process
    label_params.starthr_vec = [2015 06 19 00 00 00];
    label_params.endhr_vec   = [2015 06 21 00 00 00];   
    % output dir
    folder_name = strcat('processed_',datestr(now,'yyyy-mm-dd_HH:MM'));
    label_params.outdir = fullfile(label_params.output_dir,folder_name);

end