% === Core code to process MASC data =====================================
%
% to be completed !!!
%
% ========================================================================
function [tt,flag] = MASC_triplet_process(flake_id,idx_pics,current_dir_list,pic_list,label,process)

    % Load data if possible
    t_loading = tic;
    try
        
        for i=1:length(idx_pics)
            
            if pic_list.cam(idx_pics(i)) == process.camera_order(1)
        
                % left image
                left_filepath = fullfile(current_dir_list,pic_list.files{idx_pics(i)});
                left_filename = pic_list.files{idx_pics(i)};
                left_data = imread(left_filepath);
                left_cam = pic_list.cam(idx_pics(i));
                left_tnum = pic_list.time_num(idx_pics(i));
                
            elseif pic_list.cam(idx_pics(i)) == process.camera_order(2)
        
                % mid image
                mid_filepath = fullfile(current_dir_list,pic_list.files{idx_pics(i)});
                mid_filename = pic_list.files{idx_pics(i)};
                mid_data = imread(mid_filepath);
                mid_cam = pic_list.cam(idx_pics(i));
                mid_tnum = pic_list.time_num(idx_pics(i));
                
            elseif pic_list.cam(idx_pics(i)) == process.camera_order(3)
        
                % right image
                right_filepath = fullfile(current_dir_list,pic_list.files{idx_pics(i)});
                right_filename = pic_list.files{idx_pics(i)};
                right_data = imread(right_filepath);
                right_cam = pic_list.cam(idx_pics(i));
                right_tnum = pic_list.time_num(idx_pics(i));
                
            % additional img for CSU MASCRAD for example    
            elseif length(process.camera_order) > 3 && pic_list.cam(idx_pics(i)) == process.camera_order(4)
                    
                ad1_filepath = fullfile(current_dir_list,pic_list.files{idx_pics(i)});
                ad1_filename = pic_list.files{idx_pics(i)};
                ad1_data = imread(ad1_filepath);
                ad1_cam = pic_list.cam(idx_pics(i));
                ad1_tnum = pic_list.time_num(idx_pics(i));
                
            elseif length(process.camera_order) > 4 && pic_list.cam(idx_pics(i)) == process.camera_order(5)
                
                ad2_filepath = fullfile(current_dir_list,pic_list.files{idx_pics(i)});
                ad2_filename = pic_list.files{idx_pics(i)};
                ad2_data = imread(ad2_filepath);
                ad2_cam = pic_list.cam(idx_pics(i));
                ad2_tnum = pic_list.time_num(idx_pics(i));
                
            end
            
        end
        
        % fallspeed
        flake_fallspeed = pic_list.fallspeed(pic_list.fallid == flake_id);
        

    catch that_ERR
        
        fprintf('Warning : couldn''t load flake %u. Triplet discarded (something went wrong!!) \n',j);
        flag = -1;
        tt = [];
        return;
        
    end
    
    tt.loading = toc(t_loading);
    
    % Remove clutter
    t_clutter = tic;
    left_data = masking(left_data,left_cam,process);
    mid_data = masking(mid_data,mid_cam,process);
    right_data = masking(right_data,right_cam,process);
    if exist('ad1_data','var')
        ad1_data = masking(ad1_data,ad1_cam,process);
    end
    if exist('ad2_data','var')
        ad2_data = masking(ad2_data,ad2_cam,process);
    end
    tt.clutter = toc(t_clutter);
       
    % Create B&W mask
    t_edging = tic;
    left_eroded = false(size(left_data));
    left_eroded(left_data>0) = 1;
    left_eroded = imfill(left_eroded,'holes');    
    mid_eroded = false(size(mid_data));
    mid_eroded(mid_data>0) = 1;
    mid_eroded = imfill(mid_eroded,'holes');  
    right_eroded = false(size(right_data));
    right_eroded(right_data>0) = 1;
    right_eroded = imfill(right_eroded,'holes');
    if exist('ad1_data','var')
        ad1_eroded = false(size(ad1_data));
        ad1_eroded(ad1_data>0) = 1;
        ad1_eroded = imfill(ad1_eroded,'holes'); 
    end
    if exist('ad2_data','var')
        ad2_eroded = false(size(ad2_data));
        ad2_eroded(ad2_data>0) = 1;
        ad2_eroded = imfill(ad2_eroded,'holes'); 
    end
    tt.edging = toc(t_edging);
    
    % Detecting the likeliest ROIs on each image which satisfies the
    % matching criterion (~same vertical location, ~same vertical dimension)
     
    % Detecting the ROIs on each picture of the triplet
    t_roiying = tic;
    left_all_roi = roi_detection(left_data,left_eroded,process,left_cam,false);
    mid_all_roi = roi_detection(mid_data,mid_eroded,process,mid_cam,false);
    right_all_roi = roi_detection(right_data,right_eroded,process,right_cam,false);
    
    if isempty(left_all_roi) || isempty(mid_all_roi) || isempty(right_all_roi)
        
        %fprintf('WARNING : no roi candidate found on (at least) one of the pictures for flake ID %u \n',flake_id);
        flag = 0;
        tt.roiying = toc(t_roiying);
        tt.matching = 0;
        tt.feature = 0;
        tt.plotting = 0;
        tt.saving = 0;
        return;
        
    end
    
    if exist('ad1_data','var')
        [ad1_all_roi,ad1_idx_best] = roi_detection(ad1_data,ad1_eroded,process,ad1_cam,true);
    end
    if exist('ad2_data','var')
        [ad2_all_roi,ad2_idx_best] = roi_detection(ad2_data,ad2_eroded,process,ad2_cam,true);
    end
    
    
    %[all_roi,idx_best,area_focus_ratio,flag,status] = roi_detection(data,data_eroded,process,flake_cam);
    
    
    tt.roiying = toc(t_roiying);
    
    % particles matching (ref=mid camera)
    t_matching = tic;
    
    for k=1:length(mid_all_roi)
        
        mid_height = floor(mid_all_roi(k).BoundingBox(4));   
        mid_vert_pos = ceil(mid_all_roi(k).BoundingBox(2));
        matching_tol_height = max(process.matching_tol_pix,process.matching_tol_percent*mid_height);
        matching_tol_vert_pos = max(process.matching_tol_pix,process.matching_tol_percent*mid_height);

        % check possible match on left camera
        all_boundingbox = [left_all_roi.BoundingBox];
        all_boundingbox = reshape(all_boundingbox,4,length(all_boundingbox)/4); %each bounding box = a column
        all_vert_pos = all_boundingbox(2,:);
        all_height = all_boundingbox(4,:);
        idx_match = find(all_height >= (mid_height-matching_tol_height) & all_height <= (mid_height+matching_tol_height) & all_vert_pos >= (mid_vert_pos-matching_tol_vert_pos) & all_vert_pos <= (mid_vert_pos+matching_tol_vert_pos));
        % keep only the best match
        if length(idx_match) > 1
            height_var = abs(all_height-mid_height) + abs(all_vert_pos-mid_vert_pos);
            [~,idx_left_match(k)] = min(height_var);
        elseif length(idx_match) == 1
            idx_left_match(k) = idx_match;
        else
            idx_left_match(k) = NaN;
        end
        
        % check possible match on right camera 
        all_boundingbox = [right_all_roi.BoundingBox];
        all_boundingbox = reshape(all_boundingbox,4,length(all_boundingbox)/4); %each bounding box = a column
        all_vert_pos = all_boundingbox(2,:);
        all_height = all_boundingbox(4,:);       
        idx_match = find(all_height >= (mid_height-matching_tol_height) & all_height <= (mid_height+matching_tol_height) & all_vert_pos >= (mid_vert_pos-matching_tol_vert_pos) & all_vert_pos <= (mid_vert_pos+matching_tol_vert_pos));
        % keep only the best match
        if length(idx_match) > 1
            height_var = abs(all_height-mid_height) + abs(all_vert_pos-mid_vert_pos);
            [~,idx_right_match(k)] = min(height_var);
        elseif length(idx_match) == 1
            idx_right_match(k) = idx_match;
        else
            idx_right_match(k) = NaN;
        end        
  
    end
    
    % among all matching triplets, we select the sharpest one (on average)
    area_focus_matrix = zeros(length(mid_all_roi),3);
    area_focus_matrix(area_focus_matrix==0) = NaN;
    for k=1:length(mid_all_roi)
   
        if ~isnan(idx_left_match(k)) && ~isnan(idx_right_match(k))
            
            % compute mid area*focus index
            crop_x = ceil(mid_all_roi(k).BoundingBox(2));
            crop_y = ceil(mid_all_roi(k).BoundingBox(1));
            crop_width = floor(mid_all_roi(k).BoundingBox(3));
            crop_height = floor(mid_all_roi(k).BoundingBox(4));
            mid_roi_data = mid_data(crop_x:(crop_x+crop_height-1),crop_y:(crop_y+crop_width-1));
            mid_roi_mean_intens = mean(mid_roi_data(mid_all_roi(k).Image))/255;
            mid_roi_range_array = rangefilt(mid_roi_data);
            mid_roi_range_intens = mean(mid_roi_range_array(mid_all_roi(k).Image))/255;
            mid_roi_focus = mid_roi_mean_intens * mid_roi_range_intens^2;
            area_focus_matrix(k,2) = mid_roi_focus * mid_all_roi(k).Area;
            
            % compute left area*focus index
            kl = idx_left_match(k);
            crop_x = ceil(left_all_roi(kl).BoundingBox(2));
            crop_y = ceil(left_all_roi(kl).BoundingBox(1));
            crop_width = floor(left_all_roi(kl).BoundingBox(3));
            crop_height = floor(left_all_roi(kl).BoundingBox(4));
            left_roi_data = left_data(crop_x:(crop_x+crop_height-1),crop_y:(crop_y+crop_width-1));
            left_roi_mean_intens = mean(left_roi_data(left_all_roi(kl).Image))/255;
            left_roi_range_array = rangefilt(left_roi_data);
            left_roi_range_intens = mean(left_roi_range_array(left_all_roi(kl).Image))/255;
            left_roi_focus = left_roi_mean_intens * left_roi_range_intens^2; %<- I added a range_intens^2 here to give more weight to sharpness vs (size+brightness)
            area_focus_matrix(k,1) = left_roi_focus * left_all_roi(kl).Area;
            
            % compute right area*focus index
            kr = idx_right_match(k);
            crop_x = ceil(right_all_roi(kr).BoundingBox(2));
            crop_y = ceil(right_all_roi(kr).BoundingBox(1));
            crop_width = floor(right_all_roi(kr).BoundingBox(3));
            crop_height = floor(right_all_roi(kr).BoundingBox(4));
            right_roi_data = right_data(crop_x:(crop_x+crop_height-1),crop_y:(crop_y+crop_width-1));
            right_roi_mean_intens = mean(right_roi_data(right_all_roi(kr).Image))/255;
            right_roi_range_array = rangefilt(right_roi_data);
            right_roi_range_intens = mean(right_roi_range_array(right_all_roi(kr).Image))/255;
            right_roi_focus = right_roi_mean_intens * right_roi_range_intens^2;
            area_focus_matrix(k,3) = right_roi_focus * right_all_roi(kr).Area;
            

        end
        
        
        
    end
    
    %disp(area_focus_matrix);   
    tt.matching = toc(t_matching);
    
    % if the area_focus_matrix is full of NaN, we didn't find any potential match
    if sum(isnan(area_focus_matrix(:))) == length(area_focus_matrix(:))
        
        %fprintf('WARNING : no matching triplet found for flake ID %u \n',flake_id);
        flag = 1;
        tt.feature = 0;
        tt.plotting = 0;
        tt.saving = 0;
        return;
        
    end
       
    % disp(area_focus_matrix);
    mean_area_focus = mean(area_focus_matrix,2);
    [~,idx_mid] = max(mean_area_focus);
    idx_left = idx_left_match(idx_mid);
    idx_right = idx_right_match(idx_mid);
    
    t_feature = tic;
    % Now that we selected the ROIs for each image, let's compute some basic descriptors
    % left image
    roi_left = process_basic_descriptors(left_data,left_all_roi(idx_left),process);
    roi_left.n_roi = length(left_all_roi);
    roi_left.name = left_filename;
    roi_left.id = flake_id;
    roi_left.cam = left_cam;
    roi_left.tnum = left_tnum;
    roi_left.fallspeed = flake_fallspeed;
 
    % mid image
    roi_mid = process_basic_descriptors(mid_data,mid_all_roi(idx_mid),process);
    roi_mid.n_roi = length(mid_all_roi);
    roi_mid.name = mid_filename;
    roi_mid.id = flake_id;
    roi_mid.cam = mid_cam;
    roi_mid.tnum = mid_tnum;
    roi_mid.fallspeed = flake_fallspeed;
    
    % right image
    roi_right = process_basic_descriptors(right_data,right_all_roi(idx_right),process);
    roi_right.n_roi = length(right_all_roi);
    roi_right.name = right_filename;
    roi_right.id = flake_id;
    roi_right.cam = right_cam;
    roi_right.tnum = right_tnum;
    roi_right.fallspeed = flake_fallspeed;
    
    % additional images
    if exist('ad1_data','var') && ~isempty(ad1_all_roi)
        roi_ad1 = process_basic_descriptors(ad1_data,ad1_all_roi(ad1_idx_best),process);
        roi_ad1.n_roi = length(ad1_all_roi);
        roi_ad1.name = ad1_filename;
        roi_ad1.id = flake_id;
        roi_ad1.cam = ad1_cam;
        roi_ad1.tnum = ad1_tnum;
        roi_ad1.fallspeed = flake_fallspeed;
    end
    
    if exist('ad2_data','var') && ~isempty(ad2_all_roi)
        roi_ad2 = process_basic_descriptors(ad2_data,ad2_all_roi(ad2_idx_best),process);
        roi_ad2.n_roi = length(ad2_all_roi);
        roi_ad2.name = ad2_filename;
        roi_ad2.id = flake_id;
        roi_ad2.cam = ad2_cam;
        roi_ad2.tnum = ad2_tnum;
        roi_ad2.fallspeed = flake_fallspeed;
    end

    tt.feature = toc(t_feature);
    
  
    % generate plots and display (if desired)
    t_plotting = tic;
    if process.generate_figs
        
        if exist('roi_ad1','var') || exist('roi_ad2','var')
            nline_fig = 2;
        else
            nline_fig = 1;
        end
        
        % triplet
        fig1=figure('Visible',process.display_figs,'Renderer','painters');
        subplot(nline_fig,3,1);
        imshow(roi_left.data);
        xlabel(sprintf('cam %u',roi_left.cam));
        subplot(nline_fig,3,2);
        imshow(roi_mid.data);
        xlabel(sprintf('cam %u',roi_mid.cam));
        title(sprintf('Triplet #%u',flake_id));
        subplot(nline_fig,3,3);
        imshow(roi_right.data);
        xlabel(sprintf('cam %u',roi_right.cam));
        if nline_fig>1 && exist('roi_ad1','var')
            subplot(nline_fig,3,4);
            imshow(roi_ad1.data);
            xlabel(sprintf('add. cam %u',roi_ad1.cam));
        end
%         if nline_fig>1
%             subplot(nline_fig,3,5);
%             title('Additional images');
%         end
        if nline_fig>1 & exist('roi_ad2','var')
            subplot(nline_fig,3,6);
            imshow(roi_ad2.data);
            xlabel(sprintf('add. cam %u',roi_ad2.cam));
        end
        
%         % additional images
%         if exist('roi_ad1','var') || exist('roi_ad2','var')
%             fig9 = figure('Visible',process.display_figs,'Renderer','painters');
%             subplot(121);
%             if exist('roi_ad1','var')
%                 imshow(roi_ad1.data);
%                 xlabel(sprintf('add. cam %u',roi_ad1.cam));
%             end
%             subplot(122);
%             if exist('roi_ad2','var')
%                 imshow(roi_ad2.data);
%                 xlabel(sprintf('add. cam %u',roi_ad2.cam));
%             end
%         end
        
%         % processing left
%         roi = roi_left;
%         fig2=figure('Visible',process.display_figs);
%         colormap(gray(255));
%         hold on; axis equal; box on;
%         image(double(roi.data));
%         set(gca,'YDir','reverse');  
%         t = linspace(0,2*pi);
%         theta = -roi.E.theta*pi/180;
%         xt1 = roi.E.X0 + cos(theta).*roi.E.a.*cos(t) - sin(theta).*roi.E.b.*sin(t);
%         yt1 = roi.E.Y0 + sin(theta).*roi.E.a.*cos(t) + cos(theta).*roi.E.b.*sin(t);
%         xt2 = roi.E_out.X0 + cos(theta).*roi.E_out.a.*cos(t) - sin(theta).*roi.E_out.b.*sin(t);
%         yt2 = roi.E_out.Y0 + sin(theta).*roi.E_out.a.*cos(t) + cos(theta).*roi.E_out.b.*sin(t);
%         xt3 = roi.E_in.X0 + cos(theta).*roi.E_in.a.*cos(t) - sin(theta).*roi.E_in.b.*sin(t);
%         yt3 = roi.E_in.Y0 + sin(theta).*roi.E_in.a.*cos(t) + cos(theta).*roi.E_in.b.*sin(t);
%         xt4 = roi.C_out.X0 + roi.C_out.r .*cos(t);
%         yt4 = roi.C_out.Y0 + roi.C_out.r .*sin(t);
%         plot(roi.x_perim,roi.y_perim,'y-');
%         plot(xt1,yt1,'r-','linewidth',2);
%         plot(xt2,yt2,'c-','linewidth',2);
%         plot(roi.hull.xh,roi.hull.yh,'c--');
%         plot(xt3,yt3,'g-','linewidth',2);
%         plot(xt4,yt4,'b-','linewidth',2);
%         plot(roi.E.X0,roi.E.Y0,'rx');
%         plot(roi.E_out.X0,roi.E_out.Y0,'co');
%         plot(roi.E_in.X0,roi.E_in.Y0,'gv');
%         plot([roi.DmaxA(1) roi.DmaxB(1)],[roi.DmaxA(2) roi.DmaxB(2)],'r--');
%         line(roi.Rect.rectx,roi.Rect.recty);
%         xlabel('x axis [pixels]');
%         ylabel('y axis [pixels]');
%         title(roi.E.theta);
%         
%         % processing mid
%         roi = roi_mid;
%         fig2=figure('Visible',process.display_figs);
%         colormap(gray(255));
%         hold on; axis equal; box on;
%         image(double(roi.data));
%         set(gca,'YDir','reverse');  
%         t = linspace(0,2*pi);
%         theta = -roi.E.theta*pi/180;
%         xt1 = roi.E.X0 + cos(theta).*roi.E.a.*cos(t) - sin(theta).*roi.E.b.*sin(t);
%         yt1 = roi.E.Y0 + sin(theta).*roi.E.a.*cos(t) + cos(theta).*roi.E.b.*sin(t);
%         xt2 = roi.E_out.X0 + cos(theta).*roi.E_out.a.*cos(t) - sin(theta).*roi.E_out.b.*sin(t);
%         yt2 = roi.E_out.Y0 + sin(theta).*roi.E_out.a.*cos(t) + cos(theta).*roi.E_out.b.*sin(t);
%         xt3 = roi.E_in.X0 + cos(theta).*roi.E_in.a.*cos(t) - sin(theta).*roi.E_in.b.*sin(t);
%         yt3 = roi.E_in.Y0 + sin(theta).*roi.E_in.a.*cos(t) + cos(theta).*roi.E_in.b.*sin(t);
%         xt4 = roi.C_out.X0 + roi.C_out.r .*cos(t);
%         yt4 = roi.C_out.Y0 + roi.C_out.r .*sin(t);
%         plot(roi.x_perim,roi.y_perim,'y-');
%         plot(xt1,yt1,'r-','linewidth',2);
%         plot(xt2,yt2,'c-','linewidth',2);
%         plot(roi.hull.xh,roi.hull.yh,'c--');
%         plot(xt3,yt3,'g-','linewidth',2);
%         plot(xt4,yt4,'b-','linewidth',2);
%         plot(roi.E.X0,roi.E.Y0,'rx');
%         plot(roi.E_out.X0,roi.E_out.Y0,'co');
%         plot(roi.E_in.X0,roi.E_in.Y0,'gv');
%         plot([roi.DmaxA(1) roi.DmaxB(1)],[roi.DmaxA(2) roi.DmaxB(2)],'r--');
%         line(roi.Rect.rectx,roi.Rect.recty);
%         xlabel('x axis [pixels]');
%         ylabel('y axis [pixels]');
%         title(roi.E.theta);
%         
%         % processing right
%         roi = roi_right;
%         fig2=figure('Visible',process.display_figs);
%         colormap(gray(255));
%         hold on; axis equal; box on;
%         image(double(roi.data));
%         set(gca,'YDir','reverse');  
%         t = linspace(0,2*pi);
%         theta = -roi.E.theta*pi/180;
%         xt1 = roi.E.X0 + cos(theta).*roi.E.a.*cos(t) - sin(theta).*roi.E.b.*sin(t);
%         yt1 = roi.E.Y0 + sin(theta).*roi.E.a.*cos(t) + cos(theta).*roi.E.b.*sin(t);
%         xt2 = roi.E_out.X0 + cos(theta).*roi.E_out.a.*cos(t) - sin(theta).*roi.E_out.b.*sin(t);
%         yt2 = roi.E_out.Y0 + sin(theta).*roi.E_out.a.*cos(t) + cos(theta).*roi.E_out.b.*sin(t);
%         xt3 = roi.E_in.X0 + cos(theta).*roi.E_in.a.*cos(t) - sin(theta).*roi.E_in.b.*sin(t);
%         yt3 = roi.E_in.Y0 + sin(theta).*roi.E_in.a.*cos(t) + cos(theta).*roi.E_in.b.*sin(t);
%         xt4 = roi.C_out.X0 + roi.C_out.r .*cos(t);
%         yt4 = roi.C_out.Y0 + roi.C_out.r .*sin(t);
%         plot(roi.x_perim,roi.y_perim,'y-');
%         plot(xt1,yt1,'r-','linewidth',2);
%         plot(xt2,yt2,'c-','linewidth',2);
%         plot(roi.hull.xh,roi.hull.yh,'c--');
%         plot(xt3,yt3,'g-','linewidth',2);
%         plot(xt4,yt4,'b-','linewidth',2);
%         plot(roi.E.X0,roi.E.Y0,'rx');
%         plot(roi.E_out.X0,roi.E_out.Y0,'co');
%         plot(roi.E_in.X0,roi.E_in.Y0,'gv');
%         plot([roi.DmaxA(1) roi.DmaxB(1)],[roi.DmaxA(2) roi.DmaxB(2)],'r--');
%         line(roi.Rect.rectx,roi.Rect.recty);
%         xlabel('x axis [pixels]');
%         ylabel('y axis [pixels]');
%         title(roi.E.theta);
%                 
%         % - insert more figures here -         
%         drawnow(); 
        
        
    end
    tt.plotting = toc(t_plotting);
    
    % save processed triplets in label_params.outidr (if desired)
    t_saving = tic;

    if process.saveresults   
                
        % creation of the metadata file (first time)
        if ~exist(label.outdir,'dir')    
            mkdir(label.outdir);
            create_proc_params_file(label.outdir,label,process);
            fprintf('creation of a metadata file proc_params.txt...\n');
        end
                
        date_folder = datestr(roi_mid.tnum,'yyyy.mm.dd');
        hour_folder = datestr(roi_mid.tnum,'HH');  
 
        path2save = fullfile(label.outdir,date_folder,hour_folder);

        if ~exist(path2save,'dir')

            mkdir(path2save);

        end
        
        triplet_mean_intens = mean([roi_left.mean_intens,roi_mid.mean_intens,roi_right.mean_intens]);
        triplet_max_intens = mean([roi_left.max_intens,roi_mid.max_intens,roi_right.max_intens]);
        triplet_max_dim = mean([max(roi_left.width,roi_left.height),max(roi_mid.width,roi_mid.height),max(roi_right.width,roi_right.height)]);
        
        
        if (triplet_mean_intens >= process.minbright) && (triplet_max_intens >= process.max_intensthresh) && (triplet_max_dim >= process.sizemin)

            % save triplet images and data .mat structure
            imwrite(roi_left.data,fullfile(path2save,left_filename),'png','BitDepth', 8);
            imwrite(roi_mid.data,fullfile(path2save,mid_filename),'png','BitDepth', 8);
            imwrite(roi_right.data,fullfile(path2save,right_filename),'png','BitDepth', 8);
            save_mat(fullfile(path2save),strcat(left_filename(1:end-4),'.mat'),roi_left);
            save_mat(fullfile(path2save),strcat(mid_filename(1:end-4),'.mat'),roi_mid);
            save_mat(fullfile(path2save),strcat(right_filename(1:end-4),'.mat'),roi_right);
            if exist('roi_ad1','var')
                imwrite(roi_ad1.data,fullfile(path2save,ad1_filename),'png','BitDepth', 8);
                save_mat(fullfile(path2save),strcat(ad1_filename(1:end-4),'.mat'),roi_ad1);
            end
            if exist('roi_ad2','var')
                imwrite(roi_ad2.data,fullfile(path2save,ad2_filename),'png','BitDepth', 8);
                save_mat(fullfile(path2save),strcat(ad2_filename(1:end-4),'.mat'),roi_ad2);
            end
            
        else
            
            if ~exist(fullfile(path2save,'BAD'),'dir')
                mkdir(fullfile(path2save,'BAD'));
            end
            % save triplet images and data .mat structure
            imwrite(roi_left.data,fullfile(path2save,'BAD',left_filename),'png','BitDepth', 8);
            imwrite(roi_mid.data,fullfile(path2save,'BAD',mid_filename),'png','BitDepth', 8);
            imwrite(roi_right.data,fullfile(path2save,'BAD',right_filename),'png','BitDepth', 8);
            save_mat(fullfile(path2save,'BAD'),strcat(left_filename(1:end-4),'.mat'),roi_left);
            save_mat(fullfile(path2save,'BAD'),strcat(mid_filename(1:end-4),'.mat'),roi_mid);
            save_mat(fullfile(path2save,'BAD'),strcat(right_filename(1:end-4),'.mat'),roi_right);
            if exist('roi_ad1','var')
                imwrite(roi_ad1.data,fullfile(path2save,'BAD',ad1_filename),'png','BitDepth', 8); 
                save_mat(fullfile(path2save,'BAD'),strcat(ad1_filename(1:end-4),'.mat'),roi_ad1);
            end
            if exist('roi_ad2','var')
                imwrite(roi_ad2.data,fullfile(path2save,'BAD',ad2_filename),'png','BitDepth', 8); 
                save_mat(fullfile(path2save,'BAD'),strcat(ad2_filename(1:end-4),'.mat'),roi_ad2);
            end
        end

        % save triplet figures (if desired)
        if process.save_figs
            path2save = fullfile(label.outdir,date_folder,hour_folder,'TRIPLETS');
            if ~exist(path2save,'dir')
                mkdir(path2save);
            end
            saveas(fig1,fullfile(path2save,sprintf('triplet_%u.png',flake_id)));
        end
        
    end
    
    tt.saving = toc(t_saving);
         
    % if everyting went well, return a 2 as flag
    flag = 2;
    
    % manually clear some variables because matlab is not good to
    % do it by itself when parallel processing is enabled
    if strcmp(process.display_figs,'off')
        close all;
    end
    clearvars -except tt flag;

end












