%OUTDIRS Processed data and image output directories
%   [CROPCAM,UNCROPCAM,REJECTS,TRIPLETS, FID1, FID2] =
%   outdirs(DIRNAME,CROPCAMDIR,UNCROPCAMDIR,REJECTDIR,TRIPLETDIR) places
%   in the directory DIRNAME three directories for image output for cropped
%   images CROPCAMDIR uncropped images UNCROPCAMDIR and rejects
%   REJECTDIR. FID1 and FID2 are the ids of the files for image statistics 
%   and diagnostic statistics.
  
%   Copyright Tim Garrett, University of Utah. This code is freely available for
%   non-commercial distribution and modification
function outdirs(dirname,label_params)

    cropcam = fullfile(dirname,label_params.cropcamdir);
    uncropcam = fullfile(dirname,label_params.uncropcamdir);
    rejects = fullfile(dirname,label_params.rejectdir);
    triplets = fullfile(dirname,label_params.tripletdir);

    % Create subdirectories
    if ~exist(cropcam,'dir')
        mkdir(cropcam); mkdir(uncropcam); mkdir(rejects); mkdir(triplets);   
    end

    % Create analyzed flake statistics filename
    statsname = fullfile(dirname,'stats.txt');
    
    % Create analyzed flake diagnostics filename
    diagsname = fullfile(dirname,'diagnostics.txt');
    
    % Create analyzed flake diagnostics filename
    labelsname = fullfile(dirname,'labels.txt');
              
    % Clear directories and create statistics and diagnostics file    
    delete(fullfile(cropcam,'*'));
    delete(fullfile(uncropcam,'*'));
    delete(fullfile(rejects,'*'));
    delete(fullfile(triplets,'*'));

    % Delete preexisting output files
    if exist(statsname,'file') || exist(diagsname,'file') || exist(labelsname,'file'); 
        delete(statsname);
        delete(diagsname);
        delete(labelsname);
        disp('older output files detected : they are deleted...');
    end
    
end
 



% this is bullshit to open fid at this point    
%         fid0 = fopen(labelsname,'a');
%         fid1 = fopen(statsname,'a');
%         fid2 = fopen(diagsname,'a');
% 
%         header1 = ['id\tcam\tidcam\tyr\tmonth\tday\thr\tmins\tsec\tnflakes\tvel\tmaxdim\txsec\tperim\tpartarea\trangeintens\tstrucdens\tflakeang\tasprat\n'];
%         header2 = ['id\tcam\tidcam\tyr\tmonth\tday\thr\tmins\tsec\tacceptid\ttotalflakes\tnflakes\tintens\tminintens\tmaxintens\trangeintens\tstrucdens\tfocus\tmaxareafocus\theight\twidth\tbotloc\thorzloc\n'];
% 
%         fprintf(fid1,header1);
%         fprintf(fid2,header2);
% [cropcam,uncropcam,rejects,triplets,fid0, fid1,fid2]
    



