% === Processing and analysis parameters ===========
%
% load MASC processing parameters from a text file given as input to the function. 
% The file includes some image processing parameters, criterions whether to accept/reject snowflakes detection 
% and data saving properties.
%
% Adapted from Tim Garrett, University of Utah. 
% Author : Christophe Praz (christophe.praz@epfl.ch)
% Last update : March 2018
% ===================================================
function process_params = process_params(params_file)

    if nargin == 0
        
        fprintf('Error : please enter a valid process_params_***.txt file as input. The default parameters feature has been disabled \n');
        process_params = NaN;
        
    else
        
        fid = fopen(params_file);
        C = textscan(fid,'%u %s','Delimiter','\t');
        C = C{1};
        process_params.use_triplet_algo = C(1);
        process_params.backthresh = C(2);
        if C(3) == 1
            process_params.camera_order = [0 1 2];
        elseif C(3) == 2
            process_params.camera_order = [4 3 2 0 1];
        end
        process_params.discardmat = [C(4) C(5) C(6) C(7)];
        process_params.LED_dist_from_top = C(8);
        process_params.parallel = C(9);
        process_params.retrieve_T = C(10);
        process_params.sizemin = C(11);
        process_params.min_area = C(12);
        process_params.minbright = C(13);
        process_params.max_intensthresh = C(14); 
        process_params.matching_tol_pix = C(15); 
        process_params.matching_tol_percent = C(16);
        process_params.generate_figs = C(17);
        if C(18)
            process_params.display_figs = 'on';
        else
            process_params.display_figs = 'off';
        end
        process_params.flakebrighten  = C(19);  
        process_params.min_hole_area = C(20);
        process_params.saveresults = C(21);
        process_params.save_figs = C(22);
        fclose(fid);
      
    end
  
            
    end
   
%     
%     % "default parameters" if no params file is provided
%     if nargin == 0
% 
%         % ===== general parameters =====
% 
%         % whether you want to process masc images independently or by triplets
%         % (using a particles matching algorithm)
%         process_params.use_triplet_algo = 0;
% 
%         % background threshold
%         process_params.backthresh = 12; % 12 usually
% 
%         % camera order (IR on left, middle, IR on right, add camera
%         process_params.camera_order = [0 1 2];
%         % [0 1 2] for regular MASC; [4 3 2 0 1] for CSU MASCRAD
% 
%         % masking to limit field of view [top bottom left right]
%         process_params.discardmat = [0 0 0 0];
%         % [50 50 210 250] : EA, !replace cam 0 by 4 in masking.m!!!
%         % [400 550 200 250] : Vanderbilt
%         % [400 548 320 320] : APRES3 2015-16;
%         % [300 300 300 300] : APRES3 2016-17; top and bot are a bit arbitrary
%         %[400 400 320 410];%[400 400 320 320];
%         % [400 400 410 410] for Davos MCH/LTE test (?)
%         % 410 right is for Davos meteoswiss data
%         % Antarctica : maybe try more margin top and bottom 548 bottom
% 
%         % distance of the LED array from the top boarder
%         process_params.LED_dist_from_top = 0;
% 
%         % whether you want to run the code using the Parallel Comp. ToolBox
%         process_params.parallel = false;
% 
%         % whether try to retrieve local temperature from Thygan
%         process_params.retrieve_T = false;
% 
% 
% 
%         % ===== criterions to accept/reject snowflake detection =====
% 
%         % minimum acceptable average width in pixels (TJ : 200 microns)
%         process_params.sizemin = 0; % 10 usually
% 
%         % minimum acceptable size of the ROI area
%         process_params.min_area = 20; %20 usually
% 
%         % minimum acceptable maximum pixel brightness [0 1].
%         process_params.minbright = 0; %0.10 usually
% 
%         % minimum threshold for the "max_intens" parameter
%         process_params.max_intensthresh = 0; % 0.15 usually
% 
%         % maximum acceptable length for a flake to touch the image frame edge
%         % process_params.edgetouchlength = 500;
% 
%         % minimum threshold for the "focus" parameter, lower value corresponds
%         % to fewer "rejects"
%         % process_params.focusthresh = 0.015;
% 
%         % minimum threshold for the "areafocus" parameter    
%         % process_params.areafocusthresh = 20;
% 
%         % minimum threshold for the "range_intens" parameter (internal var.)
%         % process_params.range_intensthresh = 0.045;
% 
%         % minimum threshold for the "modified Laplacian operator" indicating
%         % blurry snowflake detection
%         % process_params.laplacian_threshold = 5;
% 
%         % process_params.area_laplacian_threshold = 15000;
% 
% 
%         % ===== criterions to match particles of the same triplet =====
%         % (if use_triplet_algo = true)
% 
%         % tolerance for matching algo (in pixels)
%         process_params.matching_tol_pix = 50; % 30 usually
% 
%         % 100% of the height of the particle
%         process_params.matching_tol_percent = 1; % 1 usually
% 
% 
%         % ===== post-processing, saving and illustration =====
% 
%         % generate some post-processed figure (decrease speed)
%         process_params.generate_figs = true;
% 
%         % display these figures? (decrease speed)
%         process_params.display_figs = 'on';
% 
%         % modify processed images for improved display by brightening?
%         process_params.flakebrighten  = true;     
% 
%        % minimal area in pixels to consider a dark (<thresh) part in a
%        % snowflake as a hole
%        process_params.min_hole_area = 10;
% 
%        % whether you want to save croped images and data structures in the
%        % campaign folder
%        process_params.saveresults = true;
% 
%        % save generated figures in a FIGURES folder beside the croped images
%        % IMAGES folder and the .mat structure DATA folder
%        process_params.save_figs = false;

