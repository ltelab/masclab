
function process_new_descriptors(dirname,t_str_start,t_str_stop)

    %dirname = '/home/praz/Documents/MASC/sample_snow_20150620/sample_cool/processed_05-Feb-2016/DATA/GOOD';
    %'/media/praz/Masc-Data/APRES3_2015/PROCESSED_20160119/DATA/GOOD';
    
    fprintf('Processing new descriptors (if necessary)... ');

    tmin = datenum(t_str_start,'yyyymmddHHMMSS');
    tmax = datenum(t_str_stop,'yyyymmddHHMMSS');
    
    all_dir = dir(fullfile(dirname,'**','*.mat'));
    all_folders = {all_dir.folder};
    all_files = {all_dir.name};
    n_files = numel(all_files);
    
    fprintf('%u snowflakes found ! \n',n_files); pause(1);

%     file_list = dir(fullfile(dirname,'20*.mat'));
%     % for Massimo data
%     if isempty(file_list)
%         file_list = dir(fullfile(dirname,'ICE*.mat'));
%     end
%     % for Vanderbilt
%     if isempty(file_list)
%         file_list = dir(fullfile(dirname,'FLAKE*.mat'));
%     end
%     file_only_list = {file_list.name};
%     file_list = fullfile(dirname,file_only_list);

    % loop over the snowflakes present in the folder
    for i=1:n_files

        is_modif = false;
        file = fullfile(all_folders{i},all_files{i});
        load(file);
        %disp(all_files{i});
        if mod(i,1000) == 0
            fprintf('%u/%u \n',i,n_files);
        end

        % load only pictures within the time interval
        if roi.tnum >= tmin && roi.tnum <= tmax
           
            %roi.H = haralick_props(roi.data);
            %is_modif = true;
            
            % New descriptors based on convex hull
            if ~isfield(roi,'hull')
                roi.hull = compute_convex_hull(roi.x,roi.y,roi.perim,0);
                is_modif = true;
            end
            
            % New descriptors based on rotational symmetry
            if ~isfield(roi,'Sym')
                roi.Sym = compute_symmetry_features(roi.bw_mask_filled,roi.Dmax,roi.eq_radius,0);
                is_modif = true;
            end
            
            if ~isfield(roi,'wavs')
                roi.wavs = fmeasure(roi.data,'WAVS',[]);
                is_modif = true;
            end
            
            if ~isfield(roi,'hist_entropy')
                roi.hist_entropy = fmeasure(roi.data,'HISE',[]);
                is_modif = true;
            end
            
            if ~isfield(roi,'area_range')
                roi.area_range = roi.area * roi.range_intens;
                is_modif = true;
            end
            
            if ~isfield(roi,'blur_idx')
                blurry2 = imgaussfilt(roi.data,2);
                blurry4 = imgaussfilt(roi.data,4);
                diff2 = roi.data - blurry2;
                diff4 = roi.data - blurry4;
                roi.blur_idx.xhi2 = std2(diff2);
                roi.blur_idx.xhi4 = std2(diff4);
                is_modif = true;
            end
            
            if ~isfield(roi,'D90')
                [roi.Dmax,roi.Dmax_theta,roi.DmaxA,roi.DmaxB] = compute_Dmax(roi.hull.xh,roi.hull.yh,0);
                roi.D90 = compute_D90(roi.bw_mask_filled,roi.Dmax,roi.Dmax_theta,0,0);
                is_modif = true;
            end
                
            
            % Fmeasure descriptors
%             if ~isfield(roi,'fmeas')
%                 roi.fmeas.acmo = fmeasure(roi.data,'ACMO',[]);
%                 roi.fmeas.bren = fmeasure(roi.data,'BREN',[]);
%                 roi.fmeas.cont = fmeasure(roi.data,'CONT',[]);
%                 roi.fmeas.curv = fmeasure(roi.data,'CURV',[]);
%                 roi.fmeas.dctr = fmeasure(roi.data,'DCTR',[]);
%                 roi.fmeas.gder = fmeasure(roi.data,'GDER',[]);
%                 roi.fmeas.glva = fmeasure(roi.data,'GLVA',[]);
%                 roi.fmeas.gllv = fmeasure(roi.data,'GLLV',[]);
%                 roi.fmeas.glvn = fmeasure(roi.data,'GLVN',[]);
%                 roi.fmeas.grae = fmeasure(roi.data,'GRAE',[]);
%                 roi.fmeas.grat = fmeasure(roi.data,'GRAT',[]);
%                 roi.fmeas.gras = fmeasure(roi.data,'GRAS',[]);
%                 roi.fmeas.helm = fmeasure(roi.data,'HELM',[]);
%                 roi.fmeas.hisr = fmeasure(roi.data,'HISR',[]);
%                 roi.fmeas.lape = fmeasure(roi.data,'LAPE',[]);
%                 roi.fmeas.lapv = fmeasure(roi.data,'LAPV',[]);
%                 roi.fmeas.lapd = fmeasure(roi.data,'LAPD',[]);
%                 roi.fmeas.sfil = fmeasure(roi.data,'SFIL',[]);
%                 roi.fmeas.sfrq = fmeasure(roi.data,'SFRQ',[]);
%                 roi.fmeas.teng = fmeasure(roi.data,'TENG',[]);
%                 roi.fmeas.tenv = fmeasure(roi.data,'TENV',[]);
%                 roi.fmeas.vola = fmeasure(roi.data,'VOLA',[]);
%                 roi.fmeas.wavv = fmeasure(roi.data,'WAVV',[]);
%                 roi.fmeas.wavr = fmeasure(roi.data,'WAVR',[]);
%                 is_modif = true;
%             end
            
            
%             if ~isfield(roi,'std_intens')
%                 std_array = stdfilt(roi.data);
%                 roi.std_intens = mean(std_array(roi.bw_mask_filled))/255;
%                 is_modif = true;
%             end
    
            
            
            % New descriptors based on Hu features
%             if ~isfield(roi,'Hu')
%                 roi.Hu = compute_Hu_features(roi.data);
%                 is_modif = true;
%             end
            
            
            % complete set of Haralick features
%             if ~isfield(roi,'H_new')
%                 %roi.glcm_vec32 = compute_GLCM_features(roi.data,32);
%                 roi.H_new = haralick_props_full(roi.data,0);
%                 is_modif = true;
%             end

            % updated skeleton
%             if ~isfield(roi,'skel2')
%                 roi.skel2 = roi.skel;
%                 is_modif = true;
%             end
            %roi.skel2 = skeleton_props(roi.data);
            
            % Dmean, hull (new desc rewritten when updating MASC_process
%             if ~isfield(roi,'Dmean')
%                 roi.Dmean = 0.5*(roi.width+roi.height);
%                 is_modif = true;
%             end
            
%             if ~isfield(roi,'hull')
%                 roi.hull = compute_convex_hull(roi.x,roi.y,roi.perim,0);
%                 is_modif = true;
%             end
            
            % hog (100,100,10,8,0.2)
            %roi.hog1 = compute_hogs(roi.data,roi.orientation,300,300,30,8,1,0);
            %roi.hog2 = compute_hogs(roi.data,roi.orientation,300,300,30,8,1,0);

            % Save new structure
            if is_modif
                save(file,'roi');
            end
      
        end

    end
    
    fprintf('   Done!\n');
    
end












%                 figure;
%                 colormap(gray(255));
%                 hold on; axis equal; box on;
%                 image(double(roi.data));
%                 set(gca,'YDir','reverse');  
%                 
%                 t = linspace(0,2*pi);
%                 theta = -roi.E.theta*pi/180;
%                 xt1 = roi.E.X0 + cos(theta).*roi.E.a.*cos(t) - sin(theta).*roi.E.b.*sin(t);
%                 yt1 = roi.E.Y0 + sin(theta).*roi.E.a.*cos(t) + cos(theta).*roi.E.b.*sin(t);
%                 xt2 = roi.E_out.X0 + cos(theta).*roi.E_out.a.*cos(t) - sin(theta).*roi.E_out.b.*sin(t);
%                 yt2 = roi.E_out.Y0 + sin(theta).*roi.E_out.a.*cos(t) + cos(theta).*roi.E_out.b.*sin(t);
%                 xt3 = roi.E_in.X0 + cos(theta).*roi.E_in.a.*cos(t) - sin(theta).*roi.E_in.b.*sin(t);
%                 yt3 = roi.E_in.Y0 + sin(theta).*roi.E_in.a.*cos(t) + cos(theta).*roi.E_in.b.*sin(t);
%                 xt4 = roi.C_out.X0 + roi.C_out.r .*cos(t);
%                 yt4 = roi.C_out.Y0 + roi.C_out.r .*sin(t);
% 
%                 plot(roi.x_perim,roi.y_perim,'y-');
%                 plot(xt1,yt1,'r-','linewidth',2);
%                 plot(xt2,yt2,'c-','linewidth',2);
%                 plot(roi.hull.xh,roi.hull.yh,'c--');
%                 plot(xt3,yt3,'g-','linewidth',2);
%                 plot(xt4,yt4,'b-','linewidth',2);
%                 plot(roi.E.X0,roi.E.Y0,'rx');
%                 plot(roi.E_out.X0,roi.E_out.Y0,'co');
%                 plot(roi.E_in.X0,roi.E_in.Y0,'gv');
%                 %plot([DmaxA(1) DmaxB(1)],[DmaxA(2) DmaxB(2)],'r--');
%                 %line(roi.Rect.rectx,roi.Rect.recty);
%                 xlabel('x axis [pixels]');
%                 ylabel('y axis [pixels]');
%                 title(roi.E.theta);
















