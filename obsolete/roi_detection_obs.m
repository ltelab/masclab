% === Select the likeliest ROI to be a snowflake in a BW eroded MASC pic =
%
% roi = roi_detection(data_in,data_in_eroded,process,cam)
%
% inputs
%
%   data_in        : MASC picture (2d matrix)
%   data_in_eroded : MASC eroded picture (passed through edge_detection.m)
%   process        : struct of process_params
%   cam            : ID of the cam (0,1,2)
%
% output
%
%   roi            : struct containing the output ROI + some features
%                    roi.good = false means that no satisfying roi has been
%                    found. Check roi.status for to see why.
%
%
% Author : Christophe Praz (christophe.praz@epfl.ch)
% Last updates :
% - Jan 2016 : add several infos to the roi to get insight on the
% reliability of the detection (# particules on the frame, ratio to the
% second best particule, ...)
% ========================================================================
function roi = roi_detection_obs(data_in,data_in_eroded,process,cam)
    
    % Detect all ROI
    all_roi = regionprops(data_in_eroded,'Image','BoundingBox','SubarrayIdx', ...
    'PixelList','PixelIdxList','Perimeter','Area','MajorAxisLength', ...
    'Orientation','MinorAxisLength','Centroid'); 

    % Snowflake size selection criterion
    all_areas = [all_roi.Area];
    % fprintf('Number of particles detected : %u \n',length(all_areas));
    guess_roi = all_roi(all_areas >20);
    % fprintf('Number of particles kept : %u \n',length(guess_roi));
   
    % remove ROI touching borders of the mask (or below the LED array)
    % this is bricolage
    all_boundingbox = [guess_roi.BoundingBox];
    all_boundingbox = reshape(all_boundingbox,4,length(all_boundingbox)/4); %each bounding box = a column
    all_dist_from_top = all_boundingbox(2,:);
    all_dist_from_bot = size(data_in,1) - (all_boundingbox(2,:) + all_boundingbox(4,:));
    all_dist_from_left = all_boundingbox(1,:);
    all_dist_from_right = size(data_in,2) - (all_boundingbox(1,:) + all_boundingbox(3,:));
    if cam == 0
        idx = find(all_dist_from_top > process.discardmat(1)+1 & all_dist_from_bot > process.discardmat(2)+1 & all_dist_from_left > process.discardmat(3)+1);
    elseif cam ==2
        idx = find(all_dist_from_top > process.discardmat(1)+1 & all_dist_from_bot > process.discardmat(2)+1 & all_dist_from_right > process.discardmat(4)+1);
    else
        idx = find(all_dist_from_top > process.discardmat(1)+1 & all_dist_from_bot > process.discardmat(2)+1);
    end
    guess_roi = guess_roi(idx);

    % Check if they are some ROI detected
    if isempty(guess_roi)

        roi.flag = 'BAD';
        roi.status = 'no ROI detected';
        

    else

        roi.status = '';
        guess_roi_max_intens = zeros(length(guess_roi),1);
        guess_roi_mean_intens = zeros(length(guess_roi),1);
        guess_roi_range_intens = zeros(length(guess_roi),1);
        guess_roi_focus = zeros(length(guess_roi),1);
        guess_roi_area_focus = zeros(length(guess_roi),1);
        guess_roi_area_range = zeros(length(guess_roi),1);
        guess_roi_average_size = zeros(length(guess_roi),1);

        % Loop over all ROIs left - flake ROI selection            
        for k=1:length(guess_roi)

            % crop initial image around the ROI
            crop_x = ceil(guess_roi(k).BoundingBox(2));
            crop_y = ceil(guess_roi(k).BoundingBox(1));
            crop_width = floor(guess_roi(k).BoundingBox(3));
            crop_height = floor(guess_roi(k).BoundingBox(4));
            guess_roi_max_size(k) = max(crop_width,crop_height); 
            %guess_roi_average_size(k) = 0.5*double(crop_width + crop_height);        
            guess_roi_data = data_in(crop_x:(crop_x+crop_height-1),crop_y:(crop_y+crop_width-1));
            % maximal brightness [0 1]
            guess_roi_max_intens(k) = double(max(guess_roi_data(guess_roi(k).Image)))/255;
            % average brightness [0 1]
            guess_roi_mean_intens(k) = mean(guess_roi_data(guess_roi(k).Image))/255;
            % local variability in the image (using 3x3 box)
            guess_roi_range_array = rangefilt(guess_roi_data);
            % average local variability [0 1]
            guess_roi_range_intens(k) = mean(guess_roi_range_array(guess_roi(k).Image))/255;        
            % focus parameter
            guess_roi_focus(k) = guess_roi_mean_intens(k) * guess_roi_range_intens(k);
            % area X focus parameter (used to select the ROI)
            guess_roi_area_focus(k) = guess_roi_focus(k) * guess_roi(k).Area; % Area is the flake area (~= box size) according to Matlab
            % area X range_intens
            guess_roi_area_range(k) = guess_roi_range_intens(k) * guess_roi(k).Area;

        end

        % ROI selection
        [~,idx_selected] = max(guess_roi_area_focus);
        
        % keep the number of candidates in roi
        roi.n_roi = length(guess_roi);
        
        % find the second maximum
        if length(guess_roi) > 1
            sorted_area_focus = sort(guess_roi_area_focus,'descend');
            roi.area_focus_ratio = sorted_area_focus(1)/sorted_area_focus(2);
            if length(guess_roi) > 2
                roi.area_focus_ratio_2 = sorted_area_focus(1)/sorted_area_focus(3);
            else
                roi.area_focus_ratio_2 = inf;
            end
            %idx_second = find(guess_roi_area_focus == sorted_area_focus(end-1));
            %roi.area_focus_ratio = guess_roi_area_focus(idx_selected)/guess_roi_area_focus(idx_second);
        else
            roi.area_focus_ratio = inf;   
        end
     
        % VALIDITY CHECK - rejection criterions
        roi.flag = 'GOOD';
        roi.status = '';
        % minimal size of the flake
        if guess_roi_max_size(idx_selected) < process.sizemin
            roi.flag = 'BAD';
            roi.status = strcat(roi.status,' too small.');
        end
        
        if guess_roi_mean_intens(idx_selected) < process.minbright
            roi.flag = 'BAD';
            roi.status = strcat(roi.status,' too dark (mean).');
        end

        if guess_roi_max_intens(idx_selected) < process.max_intensthresh
            roi.flag = 'BAD';
            roi.status = strcat(roi.status,' too dark (max).');
        end
        
        % I don't use the following tests anymore
        % 
        %         if guess_roi_focus(idx_selected) < process.focusthresh
        % 
        %         roi.flag = 'BAD';
        %         roi.status = strcat(roi.status,' low focus.');
        % 
        %         end
        % 
        %         if guess_roi_area_focus(idx_selected) < process.areafocusthresh
        % 
        %         roi.flag = 'BAD';
        %         roi.status = strcat(roi.status,' low area focus.');
        % 
        %         end
        % 
        %         if guess_roi_max_intens(idx_selected) < process.max_intensthresh
        % 
        %         roi.flag = 'BAD';
        %         roi.status = strcat(roi.status,' low max intens.');
        % 
        %         end
        % 
        %         if guess_roi_range_intens(idx_selected) < process.range_intensthresh
        % 
        %         roi.flag = 'BAD';
        %         roi.status = strcat(roi.status,' low range intens.');
        % 
        %         end
        
        % if the ROI passed all the tests, crop around the snowflake and
        % save the interesting features
 
        if ~strcmp(roi.status,'no ROI detected')
             
            crop_x = ceil(guess_roi(idx_selected).BoundingBox(2));
            crop_y = ceil(guess_roi(idx_selected).BoundingBox(1));
            crop_width = floor(guess_roi(idx_selected).BoundingBox(3));
            crop_height = floor(guess_roi(idx_selected).BoundingBox(4));
            
            % clean the leftover pixels around the roi in the cropped
            % rectangle
            mask = false(size(data_in));
            mask(guess_roi(idx_selected).PixelIdxList) = true;
            data_in(~mask) = 0;
            roi.data = data_in(crop_x:(crop_x+crop_height-1),crop_y:(crop_y+crop_width-1));
           
            % retrieve some textural features
            roi.mean_intens = guess_roi_mean_intens(idx_selected);         
            roi.max_intens = guess_roi_max_intens(idx_selected);
            roi.range_intens = guess_roi_range_intens(idx_selected);
            roi.focus = guess_roi_focus(idx_selected);
            roi.area_focus = guess_roi_area_focus(idx_selected);
            roi.area_range = guess_roi_area_range(idx_selected);
            
            % retrieve regionprops features
            roi.x_loc = ceil(guess_roi(idx_selected).BoundingBox(2));
            roi.y_loc = ceil(guess_roi(idx_selected).BoundingBox(1));
            
            % the fitted ellipse
            roi.E.a = guess_roi(idx_selected).MajorAxisLength/2;
            roi.E.b = guess_roi(idx_selected).MinorAxisLength/2;
            roi.E.theta = guess_roi(idx_selected).Orientation;
            roi.centroid = guess_roi(idx_selected).Centroid;
            %roi.E.X0 = roi.centroid(1) - roi.y_loc;
            %roi.E.Y0 = roi.centroid(2) - roi.x_loc;
            %roi.dist_from_top = guess_roi(idx_selectd).Boundingbox(2);
            
        end
    
    end
    
end
    
    
    % old selection system : points attribution
    %
    %             [~,idx] = max(guess_roi_max_intens);
    %             guess_roi_points(idx) = guess_roi_points(idx) + 1;
    %             [~,idx] = max(guess_roi_mean_intens);
    %             guess_roi_points(idx) = guess_roi_points(idx) + 1;
    %             [~,idx] = max(guess_roi_range_intens);
    %             guess_roi_points(idx) = guess_roi_points(idx) + 1;
    %             [~,idx] = max(guess_roi_entropy);
    %             guess_roi_points(idx) = guess_roi_points(idx) + 1;
