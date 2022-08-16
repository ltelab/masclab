% === Core function to process MASC data =====================================
%
% Process and crop around snowflakes in the pictures located in
% campaigndir. The code assumes that the images are organized with the
% hierarchy described in reorganize_MASC_folders.m
%
% Outputs are the cropped image and a structure containing process
% parameters and snowflake features.
%
% The data processing is done taking into account three lists of parameters
% 
% - label -> label_params()
% - cam -> cam_params()
% - process -> process_params()
%
% Check out these 3 functions for more information.
%
% In particular, in process_params() you can choose to run the code using
% Matlab parallel processing toolbox or to save additional images
% illustrating how the images were processed.
%
% Adapted from Tim Garrett, University of Utah. 
% Author : Christophe Praz (christophe.praz@epfl.ch)
% Last update : October 2017
% ========================================================================

function MASC_process(label, cam, process)

% Load IO paths
% label = label_params();
% Camera and lens details
% cam = cam_params();
% Image processing thresholds    
% process = process_params(); 

% Find all relevant snowflake directories
dir_list = uploaddirs(label.campaigndir,label.starthr_vec,label.endhr_vec);

% some global stat variables
N_tot = 0;
if process.use_triplet_algo
    N_triplet_full = 0;
    N_triplet_miss = 0;
    N_matched = 0;
    N_no_match = 0;
    N_processed_indep = 0;
else    
    N_good = 0;
    N_bad = 0;
    N_blurry = 0;
end

% initialization of some times
tt_uploading = 0;
tt_loading = 0;
tt_clutter = 0;
tt_edging = 0;
tt_roiying = 0;
tt_feature = 0;
tt_plotting = 0;
tt_saving = 0;
if process.use_triplet_algo
    tt_matching = 0;
end


% Upload all snowflake pictures, dir-by-dir
t_startprogram = tic;
for i=1:length(dir_list)

    t_uploading = tic;
    % Retrieve snowflakes info
    pic_list = upload(dir_list{i},label);
    pic_list.id_unique = unique(pic_list.id);
    
    % Create output sub-folders and fils
    % outdirs(dir_list{i},label);
    
    % Number of new flakes added
    n_id = length(pic_list.id_unique);
    N_tot = N_tot + n_id;
    
    % Copy some data in tmp variables to avoid matlab parfor complain
    current_dir_list = dir_list{i};
    
    % Discard triplets missing 1 or 2 images (it can happen)
%     if process.use_triplet_algo
%         pic_id_unique = unique(pic_list.id);
%         for j=1:length(pic_id_unique)
%             idx = find(pic_list.id==pic_id_unique(j));
%             if idx < 3
%                 pic_list.time_vec(idx) = [];
%                 pic_list.id(idx) = [];
%                 pic_list.cam(idx) = [];
%                 pic_list.files
%         
%         
%     end
    
    
    % to analyse an image particularly
    % idx_oi = find(pic_list.id == 47851 & pic_list.cam == 1);
    % idx_oi = find(pic_list.id < 48100)
    
    tt_uploading = tt_uploading + toc(t_uploading);
    
    % whether to run the code in parallel or in serial
    if process.parallel
        n_workers = 16; % if you have less than 16 threads matlab will take as many as possible
    else
        n_workers = 0;
    end
    
    
    if process.use_triplet_algo == false
        % Main loop on all the pictures located in dir_list{i}
        parfor (j=1:length(pic_list.id),n_workers) % length(pic_list.id)

            if isfield(pic_list,'files')

                fprintf('Processing %s... \n',pic_list.files{j});

            end

            [tt,flag] = MASC_picture_process(j,current_dir_list,pic_list,label,process);
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

            if ~isempty(tt)

                tt_loading = tt_loading + tt.loading;
                tt_clutter = tt_clutter + tt.clutter;
                tt_edging = tt_edging + tt.edging;
                tt_roiying = tt_roiying + tt.roiying;
                tt_feature = tt_feature + tt.feature;
                tt_plotting = tt_plotting + tt.plotting;
                tt_saving = tt_saving + tt.saving;

            end

        end
    
    % triplet matching algo    
    else
    
        parfor(j=1:length(pic_list.id_unique),n_workers)
            
            idx_pics = find(pic_list.id == pic_list.id_unique(j));
            
            if length(idx_pics) >= 3
            
                %fprintf('Processing ID %u taken on %s ... \n',pic_list.id_unique(j),pic_list.time_str(idx_pics(1),:));
                try
                    [tt,flag] = MASC_triplet_process(pic_list.id_unique(j),idx_pics,current_dir_list,pic_list,label,process);
                catch
                    flag = -1;
                    str1 = pic_list.files{idx_pics(1)}; str1 = str1(1:end-5);
                    str2 = '1-2-3.png';
                    fprintf('Unknwon error while processing triplet %s, trying now to process images independently... \n',strcat(str1,str2));
                end
                N_triplet_full = N_triplet_full + 1;
                if flag < 2
                    N_no_match = N_no_match + 1;
                    for k=1:length(idx_pics)
                        %fprintf('No matched triplet found for ID %u taken on %s, image(s) found are processed independantly ! \n',pic_list.id_unique(j),pic_list.time_str(idx_pics(1),:));
                        try
                            [tt,flag] = MASC_picture_process(idx_pics(k),current_dir_list,pic_list,label,process);
                        catch
                            flag = -1;
                            fprintf('Unknwon error while processing image %s, consider checking the image manually. \n',pic_list.files{idx_pics(k)});
                        end
                        if flag==2
                            N_processed_indep = N_processed_indep + 1;
                        end
                        if ~isempty(tt)
                            tt_loading = tt_loading + tt.loading;
                            tt_clutter = tt_clutter + tt.clutter;
                            tt_edging = tt_edging + tt.edging;
                            tt_roiying = tt_roiying + tt.roiying;
                            tt_feature = tt_feature + tt.feature;
                            tt_plotting = tt_plotting + tt.plotting;
                            tt_saving = tt_saving + tt.saving;
                        end  
                    end
                else
                    N_matched = N_matched + 1;
                    if ~isempty(tt)
                        tt_loading = tt_loading + tt.loading;
                        tt_clutter = tt_clutter + tt.clutter;
                        tt_edging = tt_edging + tt.edging;
                        tt_roiying = tt_roiying + tt.roiying;
                        tt_matching = tt_matching + tt.matching;
                        tt_feature = tt_feature + tt.feature;
                        tt_plotting = tt_plotting + tt.plotting;
                        tt_saving = tt_saving + tt.saving;
                    end  
                end
                
            else %if the triplet is incomplete, we process the images independantly
                
                %fprintf('Missing image(s) for ID %u taken on %s, image(s) found are processed independantly ! \n',pic_list.id_unique(j),pic_list.time_str(idx_pics(1),:));
                N_triplet_miss = N_triplet_miss + 1;
                for k=1:length(idx_pics)
                    try
                        [tt,flag] = MASC_picture_process(idx_pics(k),current_dir_list,pic_list,label,process);
                    catch
                        flag = -1;
                        fprintf('Unknwon error while processing image %s, consider checking the image manually. \n',pic_list.files{idx_pics(k)});
                    end
                    if flag==2    
                        N_processed_indep = N_processed_indep + 1;
                    end
                    
                    if ~isempty(tt)
                        tt_loading = tt_loading + tt.loading;
                        tt_clutter = tt_clutter + tt.clutter;
                        tt_edging = tt_edging + tt.edging;
                        tt_roiying = tt_roiying + tt.roiying;
                        tt_feature = tt_feature + tt.feature;
                        tt_plotting = tt_plotting + tt.plotting;
                        tt_saving = tt_saving + tt.saving;
                    end
                end
                
            end
        
        end
    
    end
             
end

tt_program = toc(t_startprogram);
tt_all = tt_uploading + tt_loading + tt_clutter + tt_edging + tt_feature + tt_roiying + tt_plotting + tt_saving;
if process.use_triplet_algo
    tt_all = tt_all + tt_matching;
end

% display stats in Matlab
fprintf('*********************************************************\n');
fprintf('*********************************************************\n');
fprintf('*********************************************************\n');
fprintf('***** Task finished ! Total Time spent    : %2.2f seconds \n',tt_program);
if process.use_triplet_algo
    fprintf('***** Number of triplets processed        : %2.2d \n',N_tot);
    fprintf('***** Number of full triplets             : %d %2.1f (%%) \n',N_triplet_full,N_triplet_full/N_tot * 100);
    fprintf('***** Number of incomplete triplets       : %d %2.1f (%%) \n',N_triplet_miss,N_triplet_miss/N_tot * 100);
    fprintf('***** Number of matched triplets          : %d %2.1f (%%) \n',N_matched,N_matched/N_tot * 100);
    fprintf('***** Number of unmatched triplets        : %d %2.1f (%%) \n',N_no_match,N_no_match/N_tot * 100); 
    fprintf('***** Number of imgs processed indep.     : %d %2.1f (%%) \n',N_processed_indep,N_processed_indep/(N_processed_indep+3*N_triplet_full) * 100);
else
    N_tot_im = N_good + N_bad + N_blurry;
    fprintf('***** Number of flakes found              : %2.2d \n',N_tot);
    fprintf('***** Number of pictures processed        : %2.2d \n',N_tot_im);
    fprintf('***** Number of good snowflakes           : %d %2.1f (%%) \n',N_good,N_good/N_tot_im * 100);
    fprintf('***** Number of blurry snowflakes         : %d %2.1f (%%) \n',N_blurry,N_blurry/N_tot_im * 100);
    fprintf('***** Number of no/bad detections         : %d %2.1f (%%) \n',N_bad,N_bad/N_tot_im * 100);
end
fprintf('***** Total time uploading directories    : %2.1f %%  \n',tt_uploading/tt_all * 100);
fprintf('***** Total time loading pictures         : %2.1f %%  \n',tt_loading/tt_all * 100);
fprintf('***** Total time removing clutter         : %2.1f %%  \n',tt_clutter/tt_all * 100);
fprintf('***** Total time edging pictures          : %2.1f %% \n',tt_edging/tt_all * 100);
fprintf('***** Total time computing feature        : %2.1f %% \n',tt_feature/tt_all * 100);
fprintf('***** Total time selecting flake RO       : %2.1f %%  \n',tt_roiying/tt_all * 100);
if process.use_triplet_algo
    fprintf('***** Total time matching particules      : %2.1f %%  \n',tt_matching/tt_all * 100);
end
fprintf('***** Total time plotting                 : %2.1f %%  \n',tt_plotting/tt_all * 100);
fprintf('***** Total time saving processed data    : %2.1f %% \n',tt_saving/tt_all * 100);
if process.use_triplet_algo
fprintf('***** Net average time per triplet         : %2.2f \n',tt_program/N_tot);
else
fprintf('***** Net average time per image           : %2.2f \n',tt_program/N_tot_im);
end
fprintf('*********************************************************\n');
fprintf('*********************************************************\n');
fprintf('*********************************************************\n');


% save stats in a text file
if process.saveresults
    
    filename = 'proc_stats.txt';
    pathname = label.outdir;
    fileID = fopen(fullfile(pathname,filename),'w');

    
    fprintf(fileID,'*********************************************************\n');
    fprintf(fileID,'*********************************************************\n');
    fprintf(fileID,'*********************************************************\n');
    fprintf(fileID,'***** Task finished ! Total Time spent    : %2.2f seconds \n',tt_program);
    if process.use_triplet_algo
        fprintf(fileID,'***** Number of triplets processed        : %2.2d \n',N_tot);
        fprintf(fileID,'***** Number of full triplets             : %d %2.1f (%%) \n',N_triplet_full,N_triplet_full/N_tot * 100);
        fprintf(fileID,'***** Number of incomplete triplets       : %d %2.1f (%%) \n',N_triplet_miss,N_triplet_miss/N_tot * 100);
        fprintf(fileID,'***** Number of matched triplets          : %d %2.1f (%%) \n',N_matched,N_matched/N_tot * 100);
        fprintf(fileID,'***** Number of unmatched triplets        : %d %2.1f (%%) \n',N_no_match,N_no_match/N_tot * 100); 
        fprintf(fileID,'***** Number of imgs processed indep.     : %d %2.1f (%%) \n',N_processed_indep,N_processed_indep/(N_processed_indep+3*N_triplet_full) * 100);
    else
        fprintf(fileID,'***** Number of flakes found              : %2.2d \n',N_tot);
        fprintf(fileID,'***** Number of pictures processed        : %2.2d \n',N_tot_im);
        fprintf(fileID,'***** Number of good snowflakes           : %d %2.1f (%%) \n',N_good,N_good/N_tot_im * 100);
        fprintf(fileID,'***** Number of blurry snowflakes         : %d %2.1f (%%) \n',N_blurry,N_blurry/N_tot_im * 100);
        fprintf(fileID,'***** Number of no/bad detections         : %d %2.1f (%%) \n',N_bad,N_bad/N_tot_im * 100);
    end
    fprintf(fileID,'***** Total time uploading directories    : %2.1f %%  \n',tt_uploading/tt_all * 100);
    fprintf(fileID,'***** Total time loading pictures         : %2.1f %%  \n',tt_loading/tt_all * 100);
    fprintf(fileID,'***** Total time removing clutter         : %2.1f %%  \n',tt_clutter/tt_all * 100);
    fprintf(fileID,'***** Total time edging pictures          : %2.1f %% \n',tt_edging/tt_all * 100);
    fprintf(fileID,'***** Total time computing feature        : %2.1f %% \n',tt_feature/tt_all * 100);
    fprintf(fileID,'***** Total time selecting flake RO       : %2.1f %%  \n',tt_roiying/tt_all * 100);
    if process.use_triplet_algo
        fprintf(fileID,'***** Total time matching particules      : %2.1f %%  \n',tt_matching/tt_all * 100);
    end
    fprintf(fileID,'***** Total time plotting                 : %2.1f %%  \n',tt_plotting/tt_all * 100);
    fprintf(fileID,'***** Total time saving processed data    : %2.1f %% \n',tt_saving/tt_all * 100);
    if process.use_triplet_algo
    fprintf(fileID,'***** Net average time per triplet         : %2.2f \n',tt_program/N_tot);
    else
    fprintf(fileID,'***** Net average time per image           : %2.2f \n',tt_program/N_tot_im);
    end
    fprintf(fileID,'*********************************************************\n');
    fprintf(fileID,'*********************************************************\n');
    fprintf(fileID,'*********************************************************\n');

    fclose(fileID);

end
    

%delete(gcp);








% fprintf('\n');fprintf('\n');fprintf('\n');
% fprintf('*********************************************************\n');
% fprintf('*********************************************************\n');
% fprintf('*********************************************************\n');
% fprintf('***** Task finished ! Total Time spent    : %2.2f seconds \n',tt_program);
% if process.use_triplet_algo
%     fprintf('***** Number of triplets processed        : %2.2d \n',N_tot);
%     fprintf('***** Number of full triplets             : %d %2.1f (%%) \n',N_triplet_full,N_triplet_full/N_tot * 100);
%     fprintf('***** Number of incomplete triplets       : %d %2.1f (%%) \n',N_triplet_miss,N_triplet_miss/N_tot * 100);
%     fprintf('***** Number of matched triplets          : %d %2.1f (%%) \n',N_matched,N_matched/N_tot * 100);
%     fprintf('***** Number of unmatched triplets        : %d %2.1f (%%) \n',N_no_match,N_no_match/N_tot * 100);    
% else
%     N_tot_im = N_good + N_bad + N_blurry;
%     fprintf('***** Number of flakes found              : %2.2d \n',N_tot);
%     fprintf('***** Number of pictures processed        : %2.2d \n',N_tot_im);
%     fprintf('***** Number of good snowflakes           : %d %2.1f (%%) \n',N_good,N_good/N_tot_im * 100);
%     fprintf('***** Number of blurry snowflakes         : %d %2.1f (%%) \n',N_blurry,N_blurry/N_tot_im * 100);
%     fprintf('***** Number of no/bad detections         : %d %2.1f (%%) \n',N_bad,N_bad/N_tot_im * 100);
% end
% fprintf('***** Total time uploading directories    : %2.1f %%  \n',tt_uploading/tt_all * 100);
% fprintf('***** Total time loading pictures         : %2.1f %%  \n',tt_loading/tt_all * 100);
% fprintf('***** Total time removing clutter         : %2.1f %%  \n',tt_clutter/tt_all * 100);
% fprintf('***** Total time edging pictures          : %2.1f %% \n',tt_edging/tt_all * 100);
% fprintf('***** Total time computing feature        : %2.1f %% \n',tt_feature/tt_all * 100);
% fprintf('***** Total time selecting flake RO       : %2.1f %%  \n',tt_roiying/tt_all * 100);
% if process.use_triplet_algo
%     fprintf('***** Total time matching particules      : %2.1f %%  \n',tt_matching/tt_all * 100);
% end
% fprintf('***** Total time plotting                 : %2.1f %%  \n',tt_plotting/tt_all * 100);
% fprintf('***** Total time saving processed data    : %2.1f %% \n',tt_saving/tt_all * 100);
% fprintf('***** Net average time per triplet        : %2.1f \n',tt_program/N_tot);
% fprintf('*********************************************************\n');
% fprintf('*********************************************************\n');
% fprintf('*********************************************************\n');
% fprintf('\n');fprintf('\n');fprintf('\n');











        %% ******************************************************************** 
%         % Load img if possible
%         t_loading = tic;
%         try
%             
%             flake_filepath = fullfile(current_dir_list,pic_list.files{j}); %#ok<PFBNS>
%             flake_filename = pic_list.files{j};
%             flake_imfinfo = imfinfo(flake_filepath);
%             flake_data = imread(flake_filepath);      
%             flake_id = pic_list.id(j);
%             flake_fallspeed = pic_list.fallspeed(pic_list.fallid == flake_id);
%             flake_cam = pic_list.cam(j);
%             
%             % resolution [um]
%             flake_XRes = cam.fovmat(flake_cam+1)/flake_imfinfo.Width*1000; %#ok<PFBNS>
%             flake_YRes = cam.fovmat(flake_cam+1)/flake_imfinfo.Height*1000;
%             
%             
%         catch err
%             
%             %fprintf('Warning : image %s is skipped (couldn''t load the image)\n',flake_filename);
%             continue
%             
%         end   
%         tt_loading = tt_loading + toc(t_loading);
%     
%         % Remove clutter
%         t_clutter = tic;
%         data = masking(flake_data,flake_cam,process);
%         tt_clutter = tt_clutter + toc(t_clutter);
%         
%         % Detecting edges
%         t_edging = tic;
%         data_eroded = edge_detection(data);
%         tt_edging = tt_edging + toc(t_edging);
%         
%         % Detecting the likeliest ROI
%         t_roiying = tic;
%         roi = roi_detection(data,data_eroded,process,flake_cam);
%         tt_roiying = tt_roiying + toc(t_roiying);
%         
%         % If no ROI found at all
%         if strcmp(roi.status,'no ROI detected')
%             
%             N_discarded = N_discarded + 1;
%             fprintf('WARNING : no ROI candidate found in %s', flake_filename);
%             continue;
%             
%         end
%                    
%         t_feature = tic;
%         % Some additional focus measurements
% %         roi.lap1 = fmeasure(roi.data,'LAPE',[]);
% %         
% %         roi.lap3 = fmeasure(roi.data,'LAPD',[]);
% %         roi.lap4 = fmeasure(roi.data,'LAPV',[]);
%                 
%         % Generate mask of the snowflake + holes
%         roi.bw_mask = logical(roi.data>0); % mask + holes
%         roi.bw_mask_filled = imfill(roi.bw_mask,'holes'); % mask filled        
%         bw_holes_mask = logical(roi.bw_mask_filled - roi.bw_mask);
%         holes_roi = regionprops(bw_holes_mask,'PixelIdxList','Perimeter','Area');
%         holes_area = [holes_roi.Area];
%         idx = find(holes_area > process.min_hole_area);
%         roi.holes_mask = false(size(roi.data,1),size(roi.data,2));
% 
%         if ~isempty(idx)
% 
%             holes_roi = holes_roi(idx);    
%             roi.nb_holes = length(holes_roi);
% 
% 
%             for m=1:length(holes_roi)
% 
%                 roi.holes_mask(holes_roi(m).PixelIdxList) = true;
% 
%             end
% 
%         else 
% 
%             roi.nb_holes = 0;
% 
%         end
% 
%         % after holes identification, cleaning of the leftover pixels around the flake... 
%         % /!\ I don't do this step anymore !!!
%         % [in_x,in_y] = find(roi.data > 0 & roi.bw_mask_filled == 1);
%         % roi.bw_mask_filled = bwselect(roi.bw_mask_filled,in_y,in_x);
%         roi.bw_mask = logical(roi.bw_mask_filled - roi.holes_mask);
% 
%         % the perimeter (used for the ellipse fitting
%         roi.bw_perim = bwperim(roi.bw_mask_filled);
%         [roi.y_perim,roi.x_perim] = find(roi.bw_perim);                  
% 
%         % the area
%         roi.area = sum(roi.bw_mask_filled(:));
%         %roi.area2 = guess_roi(idx_selected).Area;
%         roi.area_porous = sum(roi.bw_mask(:));
% 
%         % dimensions
%         roi.width = size(roi.bw_mask,2);
%         roi.height = size(roi.bw_mask,1);
% 
%         % adjust brightness if desired
%         if process.flakebrighten
%             roi.min_intens = double(min(roi.data(roi.bw_mask)))/255;
%             roi.data(roi.bw_mask_filled) = brightening(roi.data(roi.bw_mask_filled),roi.min_intens,roi.max_intens,process.limitintens);
%             roi.lap = fmeasure(roi.data,'LAPM',[]);
%             roi.grayvar = fmeasure(roi.data,'GLVN',[]);
%             roi.wavelet = fmeasure(roi.data,'WAVS',[]);    
%         end     
%                 
%         % Retrieve correspondant fallspeed
%         roi.fallspeed = flake_fallspeed;
%                 
%         % Retrieve otehr quantities from other ground instruments
%         roi.tnum = pic_list.time_num(j);
%         roi.thygan = load_wfj_thygan(datestr(roi.tnum,'yyyymmddHHMMSS'));
%         
%         % increment N_processed/discarded
%         if roi.good      
%             N_processed = N_processed + 1;   
%             roi.status = 'good detection.';
%         else
%             N_discarded = N_discarded + 1;
%         end
%         
%         tt_feature = tt_feature + toc(t_feature);
% 
%         
%         % Plot generation (illustration if desired)
%         t_plotting = tic;
%         
%         data_outlined = roi.data;
%         data_outlined(roi.bw_perim) = 255;
%                     
%         fig1=figure('Visible',process.displayaccept);
%         subplot(4,5,[1 2 3 4 6 7 8 9 11 12 13 14]);
%         imshow(flake_data);
%         title(sprintf('v_{f} = %1.2f [m/s]     T = %1.1f [C]     RH = %1.2f [%%]',roi.fallspeed,roi.thygan.T,roi.thygan.RH));
%         subplot(4,5,5);
%         imshow(roi.data);
%         subplot(4,5,10);
%         imshow(roi.bw_mask_filled);
%         subplot(4,5,15);
%         imshow(roi.bw_perim | roi.holes_mask);
%         subplot(4,5,20);
%         imshow(data_outlined);
%         h_text = subplot(4,5,[16 17 18 19]);
%         xl = xlim(h_text); 
%         xPos = xl(1) + diff(xl) / 2; 
%         yl = ylim(h_text); 
%         yPos = yl(1) + diff(yl) / 2; 
%         plot_text = text(xPos, yPos, sprintf('mean brightness : %1.3f \n range intens : %1.4f \n focus : %1.3f \n area * focus : %1.3f \n area focus ratio : %1.2f \n lap : %1.3f \n gray variance : %1.3f \n wavelet : %1.3f \n status : %s', ...
%             roi.mean_intens,roi.range_intens,roi.focus,roi.area_focus,roi.area_focus_ratio,roi.lap,roi.grayvar,roi.wavelet,roi.status), 'Parent', h_text);
%         set(plot_text, 'HorizontalAlignment', 'center','FontSize',12,'FontWeight','demi');
%         set(h_text,'visible','off');
% 
% 
%         fig2=figure('Visible',process.displayaccept);
%         hold on;
%         image(roi.data);
%         set(gca,'YDir','reverse');
%         plot(roi.x_perim,roi.y_perim,'r.');
%         roi.ellipse = fit_ellipse(roi.x_perim,roi.y_perim,gca);
%         
%         tt_plotting = tt_plotting + toc(t_plotting);
% 
%         % save flakes processed in /path2campaign/PROCESSED/
%         t_saving = tic;
%         
%         if process.saveresults
%         
%             path2save = label.outdir;
%             path2save_im = fullfile(path2save,'IMAGES');
%             path2save_data = fullfile(path2save,'DATA');
%             path2save_fig = fullfile(path2save,'FIGURES');
%             
%             if ~exist(path2save,'dir');
%                 
%                 mkdir(path2save);
%                 mkdir(path2save_im);
%                 mkdir(path2save_data);
%                 mkdir(path2save_fig);
%                 
%                 mkdir(fullfile(path2save_im,'GOOD'));
%                 mkdir(fullfile(path2save_im,'BAD'));
%                 mkdir(fullfile(path2save_data,'GOOD'));
%                 mkdir(fullfile(path2save_data,'BAD'));
%                 mkdir(fullfile(path2save_fig,'GOOD'));
%                 mkdir(fullfile(path2save_fig,'BAD'));
%                 
%             end
%             
%             % saving the image
%             fprintf('Saving %s...\n',strcat(flake_filename));
%             if roi.good
%                 imwrite(roi.data,fullfile(path2save_im,'GOOD',flake_filename),'png','BitDepth', 8);
%                 save_mat(fullfile(path2save_data,'GOOD'),strcat(flake_filename(1:end-4),'.mat'),roi);
%                 saveas(fig1,fullfile(path2save_fig,'GOOD',strcat('fig1_',flake_filename)));
%                 saveas(fig2,fullfile(path2save_fig,'GOOD',strcat('fig2_',flake_filename)));
%             elseif ~strcmp(roi.status,'no ROI detected')
%                 imwrite(roi.data,fullfile(path2save_im,'BAD',flake_filename),'png','BitDepth', 8);
%                 save_mat(fullfile(path2save_data,'BAD'),strcat(flake_filename(1:end-4),'.mat'),roi);
%                 saveas(fig1,fullfile(path2save_fig,'BAD',strcat('fig1_',flake_filename)));
%                 saveas(fig2,fullfile(path2save_fig,'BAD',strcat('fig2_',flake_filename)));
%             end    
%               
%             
% 
%             
%             
%         end
%         
%         close all;
%         %clear fig1 fig2 h_text plot_text roi data data_eroded data_outlined flake_data flake_imfinfo
%         % clearvars -except tt_uploading cam current_dir_list dir_list i idx_oi label N_discarded N_flake n_id N_processed pic_list process t_startprogram t_uploading tt_clutter tt_edging tt_feature tt_loading tt_plotting tt_roiying tt_saving
% 
%         
%         tt_saving = tt_saving + toc(t_saving);
        %%******************************************************************






























       
%         %%%%%%% Matlab magic: a bit black box. Mostly yanked from image toolbox
%         %%%%%%% doc sheets. Detects snowflake internal edges and creates a snowflake
%         %%%%%%% cross-section
%         se0 = strel('line', floor(1.5*linefill)/res, 0); %horz
%         se90 = strel('line', floor(1.5*linefill)/res, 90); %vert
%         
%         
%         BW = edge(flakebw,'Sobel',0.008); %edge detection algorithm
%         BWsdil = imdilate(BW, [se0 se90]); %dilates edges
%         BWdfill = imfill(BWsdil,'holes'); %fills in dilated edgy image
%         BWfinal = imerode(BWdfill,[se0 se90]); %filled cross-section
%         
%         %Isn't matlab amazing: everything one could want in one command for
%         %analyzing the image cross-section
%         statsall = regionprops(BWfinal,'Image','BoundingBox','SubarrayIdx',...
%             'PixelList','PixelIdxList',...
%             'Perimeter','Area','MajorAxisLength',...
%             'Orientation','MinorAxisLength');
%         
%         %local variability within the snowflake. Fun to look at: imshow(rangearry,[])
%         rangearry = rangefilt(flakebw);
%         totalflakes = length(statsall);
%         good = zeros(1,totalflakes);
%         
%         %%%%%%%%%%%%%%%% ASSESS THE QUALITIES AND NUMBER OF EACH DISTINCT
%         %%%%%%%%%%%%%%%% OBJECT IN THE FRAME
%             
%             %%% The next section determines whether images get sent to a
%             %%% rejects folder
%         
%         %Send to rejects folder if a black frame 
%         if length(statsall) == 0
%             %imwrite(flakebw,[rejects '/' filesname],'png',...
%             %    'Title', flakeImInfo.Title, 'Author', flakeImInfo.Author, 'Description', flakeImInfo.Description, ...
%             %    'Copyright', flakeImInfo.Copyright, 'CreationTime', flakeImInfo.CreationTime, 'Source', flakeImInfo.Source, ...
%             %    'BitDepth', 8);
%             imwrite(flakebw,[rejects '/' filesname],'png');
%             acceptid = 0;
%             labelsdata = flakefile;
%             fprintf(fid0,'%s\n',labelsdata'); % save label data 
%             diagnosticdata = [id cam idcam date time acceptid totalflakes NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN];
%             fprintf(fid2,'%6d %4d %8.1f %8d %8d %8d %8d %8d %8.3f %4d %6d %6d %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n',diagnosticdata'); % save statistics
%               continue
%             
%         else
%              %cycle through all snowflake objects in frame
%             good = zeros(1,length(statsall));
%             areafocus = zeros(1,length(statsall));
%             
%             for i = 1:length(statsall);
%                 
%                 flakeloc = statsall(i).PixelList; %indices
%                 flakebox = statsall(i).SubarrayIdx{:}; %box around each flake
%                 flakeareatotal = statsall(i).Area; %total area of each flake
%                 areamask = cat(1,statsall(i).PixelIdxList); %indices for each flake
%                 flakemask = areamask(find(flakebw(areamask) > backthresh)); %flake cross-section indices
%                 flakearea = length(flakemask); %indices that exceed background threshold
%                 intens = mean(flakebw(flakemask))/256; %[0 1] average brightness
%                 rangeintens = mean(rangearry(flakemask))/256; %[0 1] As a complexity measure, the mean interpixel range of normalized intensity within the snowflake bounds
%                 partialarea = length(flakemask)./length(areamask); %fraction of enclosed area that exceeds background
%                 
%                 % On the basis that in focus flakes are both bright and variable, a rough metric that seems to work for estimating degree of focus
%                 focus = intens.*rangeintens;
%               
%                 %length of flake that touches edge of image frame
%                 edgetouch = res*sum([length(find(flakeloc(:,1) == 1 | flakeloc(:,1) == horz)) length(find(flakeloc(:,2) == 1 | flakeloc(:,2) == vert))]);
%                 areafocus(i) = focus*flakearea;
%   
%                 %%%% identify whether a flake in the frame is good (1) or
%                 %%%% bad (0)
%                 if flakearea <= floor((sizemin/res)^2) | ... %must exceed minimum size
%                         mean(mean(rangefilt(flakebw(flakemask)))) < rangefiltthresh | ... %exceed minimum focus
%                         max(flakebw(flakemask)) < minbright*256 | ...%exceed minimum brightness
%                         edgetouch > edgetouchlength | ... % micrometers touching the edge of the frame
%                         (focusreject == 1 & round(focus*100)/100 < focusthresh); %limit bad illumination or focus; %
%                     good(i) = 0;
%                 else       
%                     good(i) = 1;
%                 end
%                 
%             end
%         end
%         
%         idx = find(good == 1);
%         nflakes = length(idx); %Total number of good images
%         
%         % Of the good flakes (idx) find the flake that is most large and in focus
%         maxareafocus = max(areafocus(idx)); 
%         relareafocus = areafocus(idx)./maxareafocus; %relative areafocus
%         relareafocussort = sort(relareafocus,2,'descend'); %Sort from high to low
%         relareafocussort(isnan(relareafocussort) == 1) = []; % Omit NaNs
%         
%         
%         %Send to the image frame to the rejects folder if no good flakes or if more than one good flake
%         %Too many flakes confuse the velocity measurement and might be
%         %blowing snow. However, the largest, most in focus image will be
%         %selected provided AREAFOCUS passes a relative threshold VELTHRESH compared to the
%         %next largest and in focus flake.
%         
%         if  nflakes == 0 | (length(relareafocussort) > 1 ...
%                 & max(relareafocus)/relareafocussort(2) < velthresh); 
%          %   imwrite(flakebw,[rejects '/' filesname],'png',...
%          %       'Title', flakeImInfo.Title, 'Author', flakeImInfo.Author, 'Description', flakeImInfo.Description, ...
%          %       'Copyright', flakeImInfo.Copyright, 'CreationTime', flakeImInfo.CreationTime, 'Source', flakeImInfo.Source, ...
%          %       'BitDepth', 8);
%             imwrite(flakebw,[rejects '/' filesname],'png')
%             acceptid = 0;
%             labelsdata = flakefile;
%             fprintf(fid0,'%s\n',labelsdata'); % save label data 
%             diagnosticdata = [id cam idcam date time acceptid totalflakes nflakes NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN];
%                 fprintf(fid2,'%6d %4d %8.1f %8d %8d %8d %8d %8d %8.3f %4d %6d %6d %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n',diagnosticdata'); % save statistics
%               if displayreject == 1;
%             'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'       
%              %pause(0.3)           
%             hr = figure(1);
%             imshow(flakebw,'InitialMagnification','fit')
%             
%             %truesize;
%             end
%             
%             continue
%             
%         else
%             % if the flake escapes the rejects, save!!!
%             
%             idxgood = idx(find(relareafocus == 1));
%             
%             %%%% The next section calculates the properties of good flake
%             %%%% raw image
%             stats = statsall(idxgood); %select the flake
%             
%             % Select the flake indices
%             areamask = cat(1,stats.PixelIdxList);
%             flakemask = areamask(find(flakebw(areamask) > backthresh)); %flake cross-section indices
%             imagealone = flakebw(flakemask);
%             
%             % Flake physical statistics
%             partialarea = length(flakemask)./length(areamask); %fraction of enclosed area that exceeds background            
%             maxdim = cat(1,stats.MajorAxisLength)*res*1e-3; %Maximum dimension along major axis
%             flakeang = cat(1,stats.Orientation); %Angle from horizontal of major axis
%             asprat = cat(1,stats.MinorAxisLength)/cat(1,stats.MajorAxisLength); %Minor/Major
%             xsec = cat(1,stats.Area)*(res*1e-3)^2; %area in mm^2
%             perim = cat(1,stats.Perimeter)*res*1e-3; %perimeter in mm
%             strucdens = bwarea(BW(flakemask))*res*1e-3/xsec; % internal structure density in edges/mm;            
%             
%             % Flake image statistics
%             intens = mean(imagealone)/255; %[0 1] average brightness
%             minintens = double(min(imagealone))/255; %[0 1] maximum brightness
%             maxintens = double(max(imagealone))/255; %[0 1] minimum brightness
%             medintens = double(median(double(imagealone)))/255; %[0 1] minimum brightness
%             rangeintens = mean(rangearry(flakemask))/255; %[0 1] As a complexity measure, the mean interpixel range of normalized intensity within snowflake bounds            
%             height = cat(1,stats.BoundingBox(4))*res*1e-3; %height in mm
%             width = cat(1,stats.BoundingBox(3))*res*1e-3; %width in mm
%             botloc = cat(1,stats.BoundingBox(2))*res*1e-3; %Distance of bottom of flake from top of frame in mm
%             horzloc = cat(1,stats.BoundingBox(1))*res*1e-3  + width/2;  %Distance of center of flake from left of frame in mm
%             
%             %ADJUST BRIGHTNESS
%             if flakebrighten == 1;
%                 flakebw(flakemask) = brightening(flakebw(flakemask),minintens,maxintens,limitintens);
%             end
%             
%             
%             %CREATE BOUNDING BOX FOR CROPPED FLAKES
%             B = stats.BoundingBox; %top left and width and height
%             cmin = floor(B(1)); cmax = cmin + B(3);
%             rmin = floor(B(2)); rmax = rmin + B(4); %Rows are from top down
%             
%             buf = ceil(100/res); %100 micrometers black space around image
%             rflakemin = max(rmin-buf,1); rflakemax = min(rmax+buf,vert); %omit edge noise
%             cflakemin = max(cmin-buf,1); cflakemax = min(cmax+buf,horz); %omit edge noise
%             flaketrunc = flakebw(rflakemin:rflakemax,cflakemin:cflakemax);
%             if displayaccept == 1; %If the desire is to see the flakes as they are processed
%             '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'       
%             ha = figure(2);
%             
%             %pause(0.3)           
%             imshow(flakebw,'InitialMagnification','fit')
%             %truesize;
%             
%             end
%             
%             %%%%%%%%%%% Identification of candidates for the Triplets
%             %%%%%%%%%%% folder where three good images pass from each camera
%             idtriplet(cam+1) = id;
%             filesnametriplet{cam+1} = filesname;                 
%             triptych{cam+1} = flaketrunc; %for triplet saving
%                                                    
%             
%                 % Write cropped and uncropped accepted images
%                 %imwrite(flaketrunc,[cropcam '/' filesname],'png',...
%                   %  'Title', flakeImInfo.Title, 'Author', flakeImInfo.Author, 'Description', flakeImInfo.Description, ...
%                    % 'Copyright', flakeImInfo.Copyright, 'CreationTime', flakeImInfo.CreationTime, 'Source', flakeImInfo.Source, ...
%                    % 'BitDepth', 8,'XResolution',flakeImInfo.XResolution,'YResolution',flakeImInfo.YResolution);
%                 %imwrite(flakebw,[uncropcam '/' filesname],'png',...
%                  %   'Title', flakeImInfo.Title, 'Author', flakeImInfo.Author, 'Description', flakeImInfo.Description, ...
%                  %   'Copyright', flakeImInfo.Copyright, 'CreationTime', flakeImInfo.CreationTime, 'Source', flakeImInfo.Source, ...
%                  %   'BitDepth', 8,'XResolution',flakeImInfo.XResolution,'YResolution',flakeImInfo.YResolution);
%                 imwrite(flaketrunc,[cropcam '/' filesname],'png')
%                     
%                 imwrite(flakebw,[uncropcam '/' filesname],'png')
%                     
%                 % Write triplet images
%                 % 'triplet'
%                 % [id cam]
%                 % idtriplet
%                 if std(idtriplet) == 0 %all camera ids are the same in idtriplet
%                     
%                     for i = 1:3;
%  %                       imwrite(triptych{i},[triplets '/' filesnametriplet{i}],'png',...
%   %                          'Title', flakeImInfo.Title, 'Author', flakeImInfo.Author, 'Description', flakeImInfo.Description, ...
%    %                         'Copyright', flakeImInfo.Copyright, 'CreationTime', flakeImInfo.CreationTime, 'Source', flakeImInfo.Source, ...
%     %                        'BitDepth', 8,'XResolution',flakeImInfo.XResolution,'YResolution',flakeImInfo.YResolution);
%                         imwrite(triptych{i},[triplets '/' filesnametriplet{i}],'png');
%                         
%                     end
%                 end
%                 
%                 %%%Any adjustments made here MUST be changed correspondingly in
%                 %%%the header labels specified in OUTDIRS.m
%                 
%                 %Write labels parameters to file
%                 labelsdata = flakefile;
% 
%                 %Write physical parameters to file
%                 flakedata = [id cam idcam date time nflakes speed maxdim xsec perim partialarea rangeintens strucdens flakeang asprat];
%                 %Write diagnostic parameters to file
%                 acceptid = 1;
%                 diagnosticdata = [id cam idcam date time acceptid totalflakes nflakes intens minintens maxintens rangeintens strucdens focus maxareafocus height width botloc horzloc];
%                 
%                 if length(flakedata) ~= 19;
%                     error('flakedata is the wrong length for correct printing to file');
%                 end
%                 
%                 if length(diagnosticdata) ~= 23;
%                     error('diagnosticdata is the wrong length for correct printing to file');
%                 end
%                 
%                 fprintf(fid0,'%s\n',labelsdata'); % save label data 
%                 fprintf(fid1,'%6d %4d %8.1f %8d %8d %8d %8d %8d %8.3f %6d %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n',flakedata'); % save statistics
%                 fprintf(fid2,'%6d %4d %8.1f %8d %8d %8d %8d %8d %8.3f %4d %6d %6d %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n',diagnosticdata'); % save statistics
%         end
%                     
% %            fclose(fid0);
% %            fclose(fid1);
% %            fclose(fid2);
% 
%     end
%     
% 
% end
% 
%     delete(poolobj); %Close parallel pool
% 
