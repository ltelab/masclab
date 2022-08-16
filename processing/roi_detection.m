% ============== Compute all the ROI in a BW eroded MASC pic ==============
%
% If the algo is not matching triplets, the function additonnally returns
% the index of the likeliest ROI in the list 
%
% roi = roi_detection_new(data_in,data_in_eroded,process,cam)
%
% inputs
%
%   data_in_eroded : MASC eroded picture (passed through edge_detection.m)
%   process        : struct of process_params
%   cam            : ID of the cam (0,1,2)
%
% output
%
%   all_roi            : struct containing the output ROIs
%   idx_best           : index of the "best" ROI (if triplet=false)
%
%
% Author : Christophe Praz (christophe.praz@epfl.ch)
% Last updates :
% - Jan 2016 : add several infos to the roi to get insight on the
% reliability of the detection (# particules on the frame, ratio to the
% second best particule, ...)
% - Aug 2016 : a lot of changes...
% ========================================================================
function [all_roi, idx_best, area_focus_ratio, flag, status] = roi_detection(data_in,data_in_eroded,process,cam,compute_idx_best)
    
    % Detect all ROI
    all_roi = regionprops(data_in_eroded,'Image','BoundingBox','SubarrayIdx', ...
    'PixelList','PixelIdxList','Perimeter','Area','MajorAxisLength', ...
    'Orientation','MinorAxisLength','Centroid'); 

    % Filter roi smaller than process.min_area
    all_areas = [all_roi.Area];
    all_roi = all_roi(all_areas >process.min_area);
   
    % Filter roi touching the borders of the mask defined by process.discardmat
    all_boundingbox = [all_roi.BoundingBox];
    all_boundingbox = reshape(all_boundingbox,4,length(all_boundingbox)/4); %each bounding box = a column
    all_dist_from_top = all_boundingbox(2,:);
    all_dist_from_bot = size(data_in_eroded,1) - (all_boundingbox(2,:) + all_boundingbox(4,:));
    all_dist_from_left = all_boundingbox(1,:);
    all_dist_from_right = size(data_in_eroded,2) - (all_boundingbox(1,:) + all_boundingbox(3,:));
    if cam == process.camera_order(1) % led on the left
        idx = find(all_dist_from_top > process.discardmat(1)+1 & all_dist_from_bot > process.discardmat(2)+1 & all_dist_from_left > process.discardmat(3)+1 & all_dist_from_right > 1);
    elseif cam == process.camera_order(3) % led on the right
        idx = find(all_dist_from_top > process.discardmat(1)+1 & all_dist_from_bot > process.discardmat(2)+1 & all_dist_from_right > process.discardmat(4)+1 & all_dist_from_left > 1);
    elseif cam == process.camera_order(2) % mid 
        idx = find(all_dist_from_top > process.discardmat(1)+1 & all_dist_from_bot > process.discardmat(2)+1 & all_dist_from_left > 1 & all_dist_from_right > 1);
    end
    if exist('idx','var')
        all_roi = all_roi(idx);
    end

    if ~compute_idx_best || isempty(all_roi) %(process.use_triplet_algo && ismember(cam,process.camera_order(1:3)))
        
        idx_best = NaN;
        area_focus_ratio = NaN;
        flag = NaN;
        status = NaN;
        
    else
        
        % Loop over all ROIs left - flake ROI selection            
        for k=1:length(all_roi)

            % crop initial image around the ROIs
            crop_x = ceil(all_roi(k).BoundingBox(2));
            crop_y = ceil(all_roi(k).BoundingBox(1));
            crop_width = floor(all_roi(k).BoundingBox(3));
            crop_height = floor(all_roi(k).BoundingBox(4));
            
            % compute some basic descriptors on those ROIs 
            all_roi_max_size(k) = max(crop_width,crop_height);         
            all_roi_data = data_in(crop_x:(crop_x+crop_height-1),crop_y:(crop_y+crop_width-1));
            % maximal brightness [0 1]
            all_roi_max_intens(k) = double(max(all_roi_data(all_roi(k).Image)))/255;
            % average brightness [0 1]
            all_roi_mean_intens(k) = mean(all_roi_data(all_roi(k).Image))/255;
            % local variability in the image (using 3x3 box)
            all_roi_range_array = rangefilt(all_roi_data);
            % average local variability [0 1]
            all_roi_range_intens(k) = mean(all_roi_range_array(all_roi(k).Image))/255;        
            % focus parameter
            all_roi_focus(k) = all_roi_mean_intens(k) * all_roi_range_intens(k)^2; % <- I added a ^2 here !!!
            % area X focus parameter (used to select the ROI)
            all_roi_area_focus(k) = all_roi_focus(k) * all_roi(k).Area; % Area is the flake area (~= box size) according to Matlab

        end

        % ROI selection
        [~,idx_best] = max(all_roi_area_focus);
        
        % second maximum
        if length(all_roi) > 1
            sorted_area_focus = sort(all_roi_area_focus,'descend');
            area_focus_ratio = sorted_area_focus(1)/sorted_area_focus(2);
        else
            area_focus_ratio = inf;      
        end 
        
        % some validity checks
        flag = 'GOOD';
        status = '';
        % minimal size of the flake
        if all_roi_max_size(idx_best) < process.sizemin
            flag = 'BAD';
            status = strcat(status,' too small.');
        end
        
        if all_roi_mean_intens(idx_best) < process.minbright
            flag = 'BAD';
            status = strcat(status,' too dark (mean).');
        end

        if all_roi_max_intens(idx_best) < process.max_intensthresh
            flag = 'BAD';
            status = strcat(status,' too dark (max).');
        end
   
    end
        
%         % keep the number of candidates in roi
%         roi.n_roi = length(guess_roi);
%         
%         % find the second maximum
%         if length(guess_roi) > 1
%             sorted_area_focus = sort(guess_roi_area_focus,'descend');
%             roi.area_focus_ratio = sorted_area_focus(1)/sorted_area_focus(2);
%             if length(guess_roi) > 2
%                 roi.area_focus_ratio_2 = sorted_area_focus(1)/sorted_area_focus(3);
%             else
%                 roi.area_focus_ratio_2 = inf;
%             end
%             %idx_second = find(guess_roi_area_focus == sorted_area_focus(end-1));
%             %roi.area_focus_ratio = guess_roi_area_focus(idx_selected)/guess_roi_area_focus(idx_second);
%         else
%             roi.area_focus_ratio = inf;   
%         end
%      
%         % VALIDITY CHECK - rejection criterions
%         roi.flag = 'GOOD';
%         roi.status = '';
%         % minimal size of the flake
%         if guess_roi_max_size(idx_selected) < process.sizemin
%             roi.flag = 'BAD';
%             roi.status = strcat(roi.status,' too small.');
%         end
%         
%         if guess_roi_mean_intens(idx_selected) < process.minbright
%             roi.flag = 'BAD';
%             roi.status = strcat(roi.status,' too dark (mean).');
%         end
% 
%         if guess_roi_max_intens(idx_selected) < process.max_intensthresh
%             roi.flag = 'BAD';
%             roi.status = strcat(roi.status,' too dark (max).');
%         end
%         
%         % I don't use the following tests anymore
%         % 
%         %         if guess_roi_focus(idx_selected) < process.focusthresh
%         % 
%         %         roi.flag = 'BAD';
%         %         roi.status = strcat(roi.status,' low focus.');
%         % 
%         %         end
%         % 
%         %         if guess_roi_area_focus(idx_selected) < process.areafocusthresh
%         % 
%         %         roi.flag = 'BAD';
%         %         roi.status = strcat(roi.status,' low area focus.');
%         % 
%         %         end
%         % 
%         %         if guess_roi_max_intens(idx_selected) < process.max_intensthresh
%         % 
%         %         roi.flag = 'BAD';
%         %         roi.status = strcat(roi.status,' low max intens.');
%         % 
%         %         end
%         % 
%         %         if guess_roi_range_intens(idx_selected) < process.range_intensthresh
%         % 
%         %         roi.flag = 'BAD';
%         %         roi.status = strcat(roi.status,' low range intens.');
%         % 
%         %         end
%         
%         % if the ROI passed all the tests, crop around the snowflake and
%         % save the interesting features
%  
%         if ~strcmp(roi.status,'no ROI detected')
%              
%             crop_x = ceil(guess_roi(idx_selected).BoundingBox(2));
%             crop_y = ceil(guess_roi(idx_selected).BoundingBox(1));
%             crop_width = floor(guess_roi(idx_selected).BoundingBox(3));
%             crop_height = floor(guess_roi(idx_selected).BoundingBox(4));
%             
%             % clean the leftover pixels around the roi in the cropped
%             % rectangle
%             mask = false(size(data_in));
%             mask(guess_roi(idx_selected).PixelIdxList) = true;
%             data_in(~mask) = 0;
%             roi.data = data_in(crop_x:(crop_x+crop_height-1),crop_y:(crop_y+crop_width-1));
%            
%             % retrieve some textural features
%             roi.mean_intens = guess_roi_mean_intens(idx_selected);         
%             roi.max_intens = guess_roi_max_intens(idx_selected);
%             roi.range_intens = guess_roi_range_intens(idx_selected);
%             roi.focus = guess_roi_focus(idx_selected);
%             roi.area_focus = guess_roi_area_focus(idx_selected);
%             roi.area_range = guess_roi_area_range(idx_selected);
%             
%             % retrieve regionprops features
%             roi.x_loc = ceil(guess_roi(idx_selected).BoundingBox(2));
%             roi.y_loc = ceil(guess_roi(idx_selected).BoundingBox(1));
%             
%             % the fitted ellipse
%             roi.E.a = guess_roi(idx_selected).MajorAxisLength/2;
%             roi.E.b = guess_roi(idx_selected).MinorAxisLength/2;
%             roi.E.theta = guess_roi(idx_selected).Orientation;
%             roi.centroid = guess_roi(idx_selected).Centroid;
%             %roi.E.X0 = roi.centroid(1) - roi.y_loc;
%             %roi.E.Y0 = roi.centroid(2) - roi.x_loc;
%             %roi.dist_from_top = guess_roi(idx_selectd).Boundingbox(2);
%             
%         end
%     
%     end
    
    
    
end
  