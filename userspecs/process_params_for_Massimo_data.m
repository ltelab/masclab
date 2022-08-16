% === Processing and analysis parameters ===========
%
% include some image processing parameters, criterions whether to
% accept/reject snowflakes detection and saving properties. Script adapted
% to work with Massimo data.
%
% Adapted from Tim Garrett, University of Utah. 
% Author : Christophe Praz (christophe.praz@epfl.ch)
% Last update : August 2015
% ===================================================
function process_params = process_params_for_Massimo_data(cam_params)

    % ===== general parameters =====

    % background threshold
    process_params.backthresh = 14;
    
    % masking to limit field of view [top bottom left right]
    process_params.discardmat = [-5 -5 -5 -5];%[400 400 320 320];
    % 410 right is for new meteoswiss data
    % Antarctica : maybe try more margin top and bottom 548 bottom
   
    % distance of the LED array from the top boarder
    process_params.LED_dist_from_top = 0;
    
    % whether you want to run the code using the Parallel Comp. ToolBox
    process_params.parallel = false;
    
    % whether try to retrieve local temperature from Thygan
    process_params.retrieve_T = false;
    
    
    
    % ===== criterions to accept/reject snowflake detection =====
    
    % minimum acceptable average width in pixels (TJ : 200 microns)
    process_params.sizemin = 10;
    
    % minimum acceptable maximum pixel brightness [0 1].
    process_params.minbright = 0;
    
    % minimum threshold for the "max_intens" parameter
    process_params.max_intensthresh = 0;
   
    % maximum acceptable length for a flake to touch the image frame edge
    % process_params.edgetouchlength = 500;
    
    % minimum threshold for the "focus" parameter, lower value corresponds
    % to fewer "rejects"
    % process_params.focusthresh = 0.015;
    
    % minimum threshold for the "areafocus" parameter    
    % process_params.areafocusthresh = 20;
    
    % minimum threshold for the "range_intens" parameter (internal var.)
    % process_params.range_intensthresh = 0.045;
    
    % minimum threshold for the "modified Laplacian operator" indicating
    % blurry snowflake detection
    % process_params.laplacian_threshold = 5;
    
    % process_params.area_laplacian_threshold = 15000;
  
    
    

    % ===== post-processing, saving and illustration =====
    
    % generate some post-processed figure (decrease speed)
    process_params.generate_figs = true;

    % display these figures? (decrease speed)
    process_params.display_figs = 'off';
            
    % modify processed images for improved display by brightening?
    process_params.flakebrighten  = true;     

   % minimal area in pixels to consider a dark (<thresh) part in a
   % snowflake as a hole
   process_params.min_hole_area = 4;
    
   % whether you want to save croped images and data structures in the
   % campaign folder
   process_params.saveresults = true;
   
   % save generated figures in a FIGURES folder beside the croped images
   % IMAGES folder and the .mat structure DATA folder
   process_params.save_figs = true;
  
            
end

