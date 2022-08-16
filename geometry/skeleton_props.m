% === Compute snowflake skeleton-based features ==========================
%
% Compute image pseudo-skeleton and output a structure SKEL containing a set
% of features derived from this skeleton.
%
% Author : Christophe Praz (christophe.praz@epfl.ch)
% Last update : May 2016
% ========================================================================
function skel = skeleton_props(im,illustration)

    if nargin == 1
        
        illustration = false;
        
    end
  
    % create bw mask
    bw = zeros(size(im)); 
    bw(im>0) = 1;  

    % smooth the silhouette
    bw_close = bwmorph(bw,'close');
    bw_close_open = bwmorph(bw_close,'open');

    % compute simple pseudo-skeleton
    bw_thin = bwmorph(bw_close_open,'thin',Inf);
     
    % remove potential leftover pixels
    all_roi = regionprops(bw_thin,'Area','PixelIdxList');
    
    %remove spur branches (experimental!)
    B = bwmorph(bw_thin, 'branchpoints');
    E = bwmorph(bw_thin, 'endpoints');

    [y,x] = find(E);
    B_loc = find(B);

    Dmask = false(size(bw_thin));
    for k = 1:numel(x)
        D = bwdistgeodesic(bw_thin,x(k),y(k));
        distanceToBranchPt = min(D(B_loc));
        if distanceToBranchPt <=3
            Dmask(D < distanceToBranchPt) =true;
        end
    end
    bw_thin = bw_thin - Dmask;
    
    % can happen after closing/opening if area is very small (2x5 pixels)
    if isempty(all_roi)
        skel.p_ratio = 0;
        skel.A_ratio = 0;
        skel.N_ends = 0;
        skel.N_junctions = 0;
%         skel.F = 0;
%         skel.max_dist = 0;
%         skel.min_dist = 0;
%         skel.med_dist = 0;
%         skel.mean_dist = 0;
%         skel.std_dist = 0;
%         skel.max_interangle = 0;
%         skel.min_interangle = 0;
%         skel.med_interangle = 0;
%         skel.mean_interangle = 0;
%         skel.std_interangle = 0; 
        return;
    end
    
    all_areas = [all_roi.Area];
    [~,idx] = max(all_areas);
    bw_thin(:) = 0;
    bw_thin(all_roi(idx).PixelIdxList) = 1;
    
    % compute associated features
    perim = sum(sum(bwperim(bw))); % bw_close_open before
    area = sum(bw(:));  % bw_close_open before
    eq_radius = sqrt(area/pi);
    perim_skel = sum(bw_thin(:));
    M_end = bwmorph(bw_thin,'endpoints');
    [skel_y,skel_x] = find(bw_thin==1);
    [skel_end_y,skel_end_x] = find(M_end==1);    
   
    skel.p_ratio = perim_skel/perim;
    skel.A_ratio = perim_skel/area;
    skel.N_ends = sum(M_end(:));%sum(sum(bwmorph(bw_thin,'endpoints')));
    skel.N_junctions = sum(sum(bwmorph(bw_thin,'branchpoints')));
    
%     % fractal_dim of the skeleton
%     skel.F = fractal_dim(bw_thin,inf);
%     %skel.F_4box = fractal_dim(bw_thin,4);
%     
%     % computation of the centroid
%     s = regionprops(bw,'centroid');
%     
%     % center the skeleton
%     for i=1:length(skel_end_x)
%         skel_end_x(i) = skel_end_x(i) - s.Centroid(1);
%         skel_end_y(i) = skel_end_y(i) - s.Centroid(2);
%     end
%     
%     % compute distance + angle for each end of the skeleton
%     if length(skel_end_x) > 0
%         for i=1:length(skel_end_x)
%             angle(i) = atand(skel_end_y(i)/skel_end_x(i));
%             dist(i) = sqrt(skel_end_x(i)^2 + skel_end_y(i)^2);
%             % convert angle to 0->360 degrees
%             if skel_end_y(i) >= 0 && skel_end_x(i) >= 0
%             elseif skel_end_y(i) >= 0 && skel_end_x(i) < 0
%                 angle(i) = 180 + angle(i);
%             elseif skel_end_y(i) < 0 && skel_end_x(i) < 0
%                 angle(i) = 180 + angle(i);
%             elseif skel_end_y(i) < 0 && skel_end_x(i) >= 0
%                 angle(i) = 360 + angle(i);
%             end
%         end
% 
%         [angle_sorted,idx_sorted] = sort(angle);
%         dist_sorted = dist(idx_sorted);
%         interdist = diff(angle_sorted);
%         interdist(end+1) = (360-angle_sorted(end)) + angle_sorted(1);
%         skel.max_dist = max(dist_sorted)/eq_radius;
%         skel.min_dist = min(dist_sorted)/eq_radius;
%         skel.med_dist = median(dist_sorted)/eq_radius;
%         skel.mean_dist = mean(dist_sorted)/eq_radius;
%         skel.std_dist = std(dist_sorted)/eq_radius;
%         skel.max_interangle = max(interdist);
%         skel.min_interangle = min(interdist);
%         skel.med_interangle = median(interdist);
%         skel.mean_interangle = mean(interdist);
%         skel.std_interangle = std(interdist);
%     else
%         skel.max_dist = 0;
%         skel.min_dist = 0;
%         skel.med_dist = 0;
%         skel.mean_dist = 0;
%         skel.std_dist = 0;
%         skel.max_interangle = 0;
%         skel.min_interangle = 0;
%         skel.med_interangle = 0;
%         skel.mean_interangle = 0;
%         skel.std_interangle = 0; 
%     end
    
%     tmp_mean = mean(dist_sorted);
%     tmp_std = std(dist_sorted);
%     dist_sorted = (dist_sorted - tmp_mean)/tmp_std;
%     dist_sorted(end+1) = dist_sorted(1);
%     % fft stuff
%     Fs = 1; % sampling frequency
%     L = length(dist_sorted); % length of signal
%     if L>=7
%         f_ = fft(dist_sorted);
%         disp(length(f_));
%         % two-sided spectrum P2
%         P2 = abs(f_/L)
%         P1 = P2(1:(L/2+1))
%         P1(2:end-1) = 2*P1(2:end-1);
%         skel.P6 = P1(6);
%         skel.P6_ratio = P1(6)/max(P1);
%     else
%         skel.P6 = 0;
%         skel.P6_ratio = 0;
%     end
    
    

    if illustration

        figure;
        subplot(221);
        im_plot = im;
        im_plot(bw_thin==1) = 255;
        imshow(im_plot);
        title('thin pseudo-skeleton');
        subplot(222);
        im_plot = bw_close_open;
        im_plot(bw_thin==1) = 0;
        imshow(im_plot);
        title('smoothed particle silhouette');
        h_text = subplot(2,2,3:4);
        xl = xlim(h_text); 
        xPos = xl(1) + diff(xl) / 2; 
        yl = ylim(h_text); 
        yPos = yl(1) + diff(yl) / 2; 
        plot_text = text(xPos, yPos, sprintf('p_{ratio} = %2.2f \n A_{ratio} = %2.2f \n N_{ends} = %u \n N_{junctions} = %u \n',skel.p_ratio,skel.A_ratio,skel.N_ends,skel.N_junctions), 'Parent', h_text);
        set(plot_text, 'HorizontalAlignment', 'center','FontSize',12,'FontWeight','demi');
        set(h_text,'visible','off'); 
        
%         figure;
%         subplot(3,3,1);
%         imshow(im);
%         subplot(3,3,2);
%         imshow(bw_close_open);
%         subplot(3,3,3); hold on; box on;
%         plot(skel_x,skel_y,'ro');
%         %plot(xp2,yp2,'kx');
%         set(gca,'Ydir','reverse');
%         plot(s.Centroid(1),s.Centroid(2),'gx','MarkerSize',9);
%         axis equal;
%         subplot(3,3,4:6); hold on; box on;
%         plot(angle_sorted,dist_sorted,'k-');
%         %     plot(angles,distances_raw,'r-');
%         xlabel('angle');
%         ylabel('distance from centroid');
%         v = axis;
%         axis([0 359 v(3) v(4)]);
%         if L>6
%             subplot(3,3,7:9);
%             plot([0:length(P1)-1],P1);
%             xlabel('frequency domain');
%             ylabel('|P1(f)|');
%         end
        
    end
    
end