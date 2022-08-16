% first script to investigate Massimo images
clear all; close all;

data_dir = '/media/praz/MyData/Massimo_data';
data_time = [2015 01 18 00 15 00];
data_tstr = datestr(data_time,'dd_mm_yyyy_HH_MM');
data_tnum = datenum(data_time);
data_name = sprintf('ICE2_%s.fig',data_tstr);

out_dir = '/media/praz/MyData/Massimo_data';
folder_name = data_name(1:end-4);

% processing parameters
process = process_params_for_Massimo_data;


fprintf('Loading all data found in %s...',fullfile(data_dir,data_name));

% load each particle present on the image
hgload(fullfile(data_dir,data_name));
myhandle = findall(gcf,'type','image');
data = get(myhandle,'cdata');

fprintf(' Done!\n');
fprintf('%u particles found! Processing and saving individual snowflakes...\n',length(data)-1);

for i=1:length(data)-1
    fprintf('%u / %u...',i,length(data)-1);
    cdata = data{i};
       
    % transform into 0-255 uint8 (same as MASC)
    cdata = double(cdata);
    cdata = cdata/65535;
    cdata = uint8(cdata*255);
    
    % "clean" the image background
    cdata_cleaned = cdata;
    cdata_cleaned(cdata_cleaned < process.backthresh) = 0;

    % generate b&w mask
    cdata_bw = false(size(cdata));
    cdata_bw(cdata_cleaned>0) = 1;
    cdata_bw_filled = imfill(cdata_bw,'holes');
    
    % keep the largest, better roi
    roi = roi_detection(cdata,cdata_bw_filled,process,1);

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
    
    % the perimeter
    B = bwboundaries(roi.bw_mask_filled,'noholes');
    if numel(B) > 1
        disp('Warning: more than 1 particle detected on bw_mask_filled !!!');
    end
    B = B{1};
    roi.x_perim = B(:,2);
    roi.y_perim = B(:,1); 
    roi.perim = length(roi.x_perim);

    % the area
    [roi.y,roi.x] = find(roi.bw_mask_filled > 0);
    roi.area = sum(roi.bw_mask_filled(:));
    roi.area_porous = sum(roi.bw_mask(:));

    % the convex hull
    roi.hull = compute_convex_hull(roi.x,roi.y,roi.perim,0);

    % the dimensions
    roi.width = size(roi.bw_mask,2);
    roi.height = size(roi.bw_mask,1);
    roi.Dmean = 0.5*(roi.width+roi.height);
    [roi.Dmax,~,DmaxA,DmaxB] = compute_Dmax(roi.hull.xh,roi.hull.yh,0);
    roi.eq_radius = sqrt(roi.area/pi);

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

    % Haralick texture features
    roi.H = haralick_props(roi.data);

    % More textural descriptors
    roi.lap = fmeasure(roi.data,'LAPM',[]); roi.area_lap = roi.lap * roi.area;
    roi.hist_entropy = fmeasure(roi.data,'HISE',[]);
    roi.wavs = fmeasure(roi.data,'WAVS',[]);

    roi.std = std2(roi.data(roi.bw_mask_filled));
    local_std = stdfilt(roi.data);
    roi.local_std = mean(local_std(roi.bw_mask_filled)); 

    roi.range_complex = roi.complex * roi.range_intens;
    roi.min_intens = min(roi.data(roi.bw_mask_filled));
    roi.contrast = double(roi.max_intens - roi.min_intens)/roi.mean_intens;

    % Compute the magig quality parameter (MODIFIED)
    roi.xhi = log((roi.lap)/2 * roi.complex * (roi.local_std)/2 * roi.Dmean); 

    % name of the particle
    roi.name = sprintf('%s__part%u',data_name(1:end-4),i);  
    roi.tnum = data_tnum;
    
    % order the fields
    roi = orderfields(roi);
    
    
    % ----- illustration -----
    if process.generate_figs
        
        fig_illu=figure('Visible',process.display_figs);
        colormap(gray(255));
        hold on; axis equal; box on;
        image(double(roi.data));
        set(gca,'YDir','reverse');  

        t = linspace(0,2*pi);
        theta = -roi.E.theta*pi/180;
        xt1 = roi.E.X0 + cos(theta).*roi.E.a.*cos(t) - sin(theta).*roi.E.b.*sin(t);
        yt1 = roi.E.Y0 + sin(theta).*roi.E.a.*cos(t) + cos(theta).*roi.E.b.*sin(t);
        xt2 = roi.E_out.X0 + cos(theta).*roi.E_out.a.*cos(t) - sin(theta).*roi.E_out.b.*sin(t);
        yt2 = roi.E_out.Y0 + sin(theta).*roi.E_out.a.*cos(t) + cos(theta).*roi.E_out.b.*sin(t);
        xt3 = roi.E_in.X0 + cos(theta).*roi.E_in.a.*cos(t) - sin(theta).*roi.E_in.b.*sin(t);
        yt3 = roi.E_in.Y0 + sin(theta).*roi.E_in.a.*cos(t) + cos(theta).*roi.E_in.b.*sin(t);
        xt4 = roi.C_out.X0 + roi.C_out.r .*cos(t);
        yt4 = roi.C_out.Y0 + roi.C_out.r .*sin(t);

        plot(roi.x_perim,roi.y_perim,'y-');
        plot(xt1,yt1,'r-','linewidth',2);
        plot(xt2,yt2,'c-','linewidth',2);
        plot(roi.hull.xh,roi.hull.yh,'c--');
        plot(xt3,yt3,'g-','linewidth',2);
        plot(xt4,yt4,'b-','linewidth',2);
        plot(roi.E.X0,roi.E.Y0,'rx');
        plot(roi.E_out.X0,roi.E_out.Y0,'co');
        plot(roi.E_in.X0,roi.E_in.Y0,'gv');
        plot([DmaxA(1) DmaxB(1)],[DmaxA(2) DmaxB(2)],'r--');
        line(roi.Rect.rectx,roi.Rect.recty);
        xlabel('x axis [pixels]');
        ylabel('y axis [pixels]');
        title(sprintf('angle : %2.2f degrees',roi.E.theta));

        drawnow();
    
    end
    
    if process.saveresults
        
        if ~exist(fullfile(out_dir,folder_name),'dir')
            
            mkdir(fullfile(out_dir,folder_name));
            mkdir(fullfile(out_dir,folder_name,'DATA'));
            mkdir(fullfile(out_dir,folder_name,'IMG'));
            mkdir(fullfile(out_dir,folder_name,'FIG'));
            
        end
        
        % save PNG image
        imwrite(roi.data,fullfile(out_dir,folder_name,'IMG',roi.name),'png','BitDepth', 8);
        % save Matlab roi structure
        save_mat(fullfile(out_dir,folder_name,'DATA'),roi.name,roi);
        % save illustration figure (if desired)
        if process.save_figs
            %savefig(fullfile(out_dir,folder_name,'FIG',roi.name),'png');
            saveas(fig_illu,fullfile(out_dir,folder_name,'FIG',roi.name),'png');
        end
      
    end

    fprintf(' Done!\n');
    
end
    
    