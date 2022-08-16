function roi = process_basic_descriptors(data_in,regionprops_roi,process)
    
    crop_x = ceil(regionprops_roi.BoundingBox(2));
    crop_y = ceil(regionprops_roi.BoundingBox(1));
    crop_width = floor(regionprops_roi.BoundingBox(3));
    crop_height = floor(regionprops_roi.BoundingBox(4));

    % clean the leftover pixels around the roi in the cropped
    % rectangle
    mask = false(size(data_in));
    mask(regionprops_roi.PixelIdxList) = true;
    data_in(~mask) = 0;
    roi.data = data_in(crop_x:(crop_x+crop_height-1),crop_y:(crop_y+crop_width-1));
    
    
    % Generate mask of the snowflake + holes
    roi.bw_mask = logical(roi.data>0); % mask + holes
    roi.bw_mask_filled = imfill(roi.bw_mask,'holes'); % mask filled        
    bw_holes_mask = logical(roi.bw_mask_filled - roi.bw_mask);
    holes_roi = regionprops(bw_holes_mask,'PixelIdxList','Perimeter','Area');
    holes_area = [holes_roi.Area];
    idx = find(holes_area > process.min_hole_area);
    roi.holes_mask = false(size(roi.data,1),size(roi.data,2));

    if ~isempty(idx)

        holes_roi = holes_roi(idx);    
        roi.nb_holes = length(holes_roi);


        for m=1:length(holes_roi)

            roi.holes_mask(holes_roi(m).PixelIdxList) = true;

        end

    else 

        roi.nb_holes = 0;

    end

    roi.bw_mask = logical(roi.bw_mask_filled - roi.holes_mask);
    
    % the area
    [roi.y,roi.x] = find(roi.bw_mask_filled > 0);
    roi.area = sum(roi.bw_mask_filled(:));
    roi.area_porous = sum(roi.bw_mask(:));

    % retrieve some textural features
    roi.mean_intens = mean(roi.data(roi.bw_mask_filled))/255;
    roi.max_intens = double(max(roi.data(roi.bw_mask_filled)))/255;
    range_array = rangefilt(roi.data);
    roi.range_intens = mean(range_array(roi.bw_mask_filled))/255;
    roi.focus = roi.mean_intens * roi.range_intens;
    roi.area_focus = roi.area * roi.focus;
    roi.area_range = roi.area * roi.range_intens;

    % retrieve regionprops features
    roi.x_loc = ceil(regionprops_roi.BoundingBox(2));
    roi.y_loc = ceil(regionprops_roi.BoundingBox(1));

    % the fitted ellipse
    roi.E.a = regionprops_roi.MajorAxisLength/2;
    roi.E.b = regionprops_roi.MinorAxisLength/2;
    roi.E.theta = regionprops_roi.Orientation;
    roi.centroid = regionprops_roi.Centroid;
    
    % the perimeter
    B = bwboundaries(roi.bw_mask_filled,'noholes');
    if numel(B) > 1
        disp('Warning: more than 1 particle detected on bw_mask_filled !!!');
    end
    B = B{1};
    roi.x_perim = B(:,2);
    roi.y_perim = B(:,1); 
    roi.perim = length(roi.x_perim);
                       
    % the convex hull
    roi.hull = compute_convex_hull(roi.x,roi.y,roi.perim,0);

    % the dimensions
    roi.width = size(roi.bw_mask,2);
    roi.height = size(roi.bw_mask,1);
    roi.Dmean = 0.5*(roi.width+roi.height);
    [roi.Dmax,roi.Dmax_theta,roi.DmaxA,roi.DmaxB] = compute_Dmax(roi.hull.xh,roi.hull.yh,0);
    roi.eq_radius = sqrt(roi.area/pi);
    roi.D90 = compute_D90(roi.bw_mask_filled,roi.Dmax,roi.Dmax_theta,0,0);
            
    % the complexity
    roi.complex = length(roi.y_perim)/(2*pi*roi.eq_radius);
            
    % the circumscribed circle
    roi.C_out = fit_circle_around(roi.x_perim,roi.y_perim,roi.hull.xh,roi.hull.yh,0);
            
    % the local centroid
    tmp = regionprops(roi.bw_mask_filled,'Centroid');
    roi.centroid_local = tmp.Centroid;
    roi.E.X0 = roi.centroid_local(1);
    roi.E.Y0 = roi.centroid_local(2);
            
    % the inscribed/circumscribed ellipses
    roi.E_in = fit_ellipse_inside(roi.x,roi.y,roi.x_perim,roi.y_perim,roi.E.theta*pi/180,0);
    roi.E_out = fit_ellipse_around(roi.x,roi.y,roi.hull.xh,roi.hull.yh,roi.E.theta*pi/180,0);
            
    % the smallest encompassing rectangle
    roi.Rect = compute_rectangularity(roi.x,roi.y,roi.perim,0);
                
    % descriptors based on those shapes
    roi.compactness = roi.area/(pi*roi.E_out.a*roi.E_out.b);
    roi.roundness = roi.area/roi.C_out.A;

    % skeleton
    roi.skel = skeleton_props(roi.data);
            
    % fractal dimension
    roi.F = fractal_dim(roi.data);
    roi.F_jac = 2*log(roi.perim/4)/log(roi.area); % [1,2] 
    
    % symmetry features
    roi.Sym = compute_symmetry_features(roi.bw_mask_filled,roi.Dmax,roi.eq_radius,0);
            
    % Haralick texture features
    roi.H = haralick_props(roi.data);
            
    % More textural descriptors
    roi.lap = fmeasure(roi.data,'LAPM',[]); roi.area_lap = roi.lap * roi.area;
    roi.hist_entropy = fmeasure(roi.data,'HISE',[]);
    roi.wavs = fmeasure(roi.data,'WAVS',[]);
            
    roi.std = std2(roi.data(roi.bw_mask_filled));
    local_std = stdfilt(roi.data);
    roi.local_std = mean(local_std(roi.bw_mask_filled)); 
    local_std = stdfilt(roi.data,ones(5));
    roi.local_std5 = mean(local_std(roi.bw_mask_filled)); 
    local_std = stdfilt(roi.data,ones(7));
    roi.local_std7 = mean(local_std(roi.bw_mask_filled));

    roi.range_complex = roi.complex * roi.range_intens;
    roi.min_intens = min(roi.data(roi.bw_mask_filled));
    roi.contrast = double(roi.max_intens - roi.min_intens)/roi.mean_intens;
            
    % adjust brightness if desired
    if process.flakebrighten
        %roi.min_intens = double(min(roi.data(roi.bw_mask)))/255;
        roi.new.data = brightening(roi.data);
    end

    % Same textural descriptors but applied on the brightened img
    roi.new.lap = fmeasure(roi.new.data,'LAPM',[]); 
    roi.new.area_lap = roi.new.lap * roi.area;       
    range_array = rangefilt(roi.new.data); %(roi.bw_mask_filled)
    roi.new.range_intens = mean(range_array(roi.bw_mask_filled))/255;
    roi.new.range_complex = roi.complex * roi.new.range_intens;
    roi.new.std = std2(roi.new.data(roi.bw_mask_filled));
    local_std = stdfilt(roi.new.data);
    roi.new.local_std = mean(local_std(roi.bw_mask_filled));
    local_std = stdfilt(roi.new.data,ones(5));
    roi.new.local_std5 = mean(local_std(roi.bw_mask_filled));
    local_std = stdfilt(roi.new.data,ones(7));
    roi.new.local_std7 = mean(local_std(roi.bw_mask_filled));
    roi.new.contrast = double(max(roi.new.data(roi.bw_mask_filled)) - min(roi.new.data(roi.bw_mask_filled))) / mean(roi.new.data(roi.bw_mask_filled));

    % Compute the magig quality parameter
    roi.xhi = log((roi.lap+roi.new.lap)/2 * roi.complex * (roi.local_std + roi.new.local_std)/2 * roi.Dmean); 

end




