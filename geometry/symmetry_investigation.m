    %load example
    clear all; close all;
    imdir = '/home/praz/Documents/MASC/masclab/events/aggregates';
    datafiles = dir(fullfile(imdir,'*.mat'));
    datafiles = {datafiles.name}'; 
    load(fullfile(imdir,datafiles{14}));   
    x = roi.x;
    y = roi.y;
    im = roi.data;
    
    xp2 = roi.x_perim;
    yp2 = roi.y_perim;
    
    % create bw mask
    bw = zeros(size(im)); 
    bw(im>0) = 1;  

    % smooth the silhouette
    bw_close = bwmorph(bw,'close');
    bw_close_open = bwmorph(bw_close,'open');

    % compute simple pseudo-skeleton
    bw_thin = bwmorph(bw_close_open,'thin',Inf);
    %bw_skel = bwmorph(bw,'skel',Inf);
     
    % remove potential leftover pixels
    all_roi = regionprops(bw_thin,'Area','PixelIdxList');
    
    
    % remove spur branches
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
    skelD = bw_thin - Dmask;
    figure;
    imshow(skelD);
    hold all;
    [y,x] = find(B); plot(x,y,'ro')

    
    
    
    % can happen after closing/opening if area is very small (2x5 pixels)
    if isempty(all_roi)
        skel.p_ratio = NaN;
        skel.A_ratio = NaN;
        skel.N_ends = NaN;
        skel.N_junctions = NaN;
        return;
    end
    
    all_areas = [all_roi.Area];
    [~,idx] = max(all_areas);
    bw_thin(:) = 0;
    bw_thin(all_roi(idx).PixelIdxList) = 1;
    [skel_y,skel_x] = find(bw_thin==1);
   
    M_end = bwmorph(bw_thin,'endpoints');
    [skel_end_y,skel_end_x] = find(M_end==1);
    M_ramif = bwmorph(bw_thin,'branchpoints');
    
    % compute associated features
    perim = sum(sum(bwperim(bw_close_open)));
    area = sum(bw_close_open(:)); 
    perim_skel = sum(bw_thin(:));
   
    skel.p_ratio = perim_skel/perim;
    skel.A_ratio = perim_skel/area;
    skel.N_ends = sum(sum(bwmorph(bw_thin,'endpoints')));
    skel.N_junctions = sum(sum(bwmorph(bw_thin,'branchpoints')));
    
    
    %scheme = zeros(size(M_ramif));
    %scheme(bw_close_open) = 255;
    %scheme(bw_thin) = 50;
    %scheme(M_ramif) = 150;


    figure;
    subplot(221);
    im_plot = im;
    im_plot(bw_thin==1) = 255;
    imshow(im_plot);
    title('thin pseudo-skeleton');
    subplot(222);
    im_plot = bw_close_open;
    im_plot(bw_thin==1) = 0;
    %implot(M_ramif==1) = 1;
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
 
    
    % computation of the centroid
    s = regionprops(roi.bw_mask_filled,'centroid');
    
    % center the skeleton
    for i=1:length(skel_end_x)
        skel_end_x(i) = skel_end_x(i) - s.Centroid(1);
        skel_end_y(i) = skel_end_y(i) - s.Centroid(2);
    end
    
    % compute distance + angle for each end of the skeleton
    for i=1:length(skel_end_x)
        angle(i) = atand(skel_end_y(i)/skel_end_x(i));
        dist(i) = sqrt(skel_end_x(i)^2 + skel_end_y(i)^2);
        % convert angle to 0->360 degrees
        if skel_end_y(i) >= 0 && skel_end_x(i) >= 0
            quadrant(i) = 1; 
        elseif skel_end_y(i) >= 0 && skel_end_x(i) < 0
            quadrant(i) = 2;
            angle(i) = 180 + angle(i);
        elseif skel_end_y(i) < 0 && skel_end_x(i) < 0
            quadrant(i) = 3;
            angle(i) = 180 + angle(i);
        elseif skel_end_y(i) < 0 && skel_end_x(i) >= 0
            quadrant(i) = 4;
            angle(i) = 360 + angle(i);
        end
    end
    [angle_sorted,idx_sorted] = sort(angle);
    dist_sorted = dist(idx_sorted);
    interdist = diff(angle_sorted);
    interdist(end+1) = (360-angle_sorted(end)) + angle_sorted(1);
    skel.max_dist = max(dist_sorted)/roi.eq_radius;
    skel.min_dist = min(dist_sorted)/roi.eq_radius;
    skel.med_dist = median(dist_sorted)/roi.eq_radius;
    skel.mean_dist = mean(dist_sorted)/roi.eq_radius;
    skel.max_interdist = max(interdist);
    skel.min_interdist = min(interdist);
    skel.med_interdist = median(interdist);
    skel.mean_interdist = mean(interdist);
    
    tmp_mean = mean(dist_sorted);
    tmp_std = std(dist_sorted);
    dist_sorted = (dist_sorted - tmp_mean)/tmp_std;
    %     % fft stuff
     Fs = 1; % sampling frequency
%     T = 1/Fs; % sampling period
     L = length(dist_sorted); % length of signal
%     t = (0:L-1)*T; % time vector
    f_ = fft(dist_sorted);
%     % two-sided spectrum P2
     P2 = abs(f_/L);
     P1 = P2(1:(L/2+1));
     P1(2:end-1) = 2*P1(2:end-1);
     fvec = 10*Fs*(0:(L/2))/L;
%     %fvec2 = 10*Fs*(1:L)/L;
    
    figure;
    subplot(3,3,1);
    imshow(roi.data);
    subplot(3,3,2);
    imshow(roi.bw_mask_filled);
    subplot(3,3,3); hold on; box on;
    plot(skel_x,skel_y,'ro');
    plot(xp2,yp2,'kx');
    %set(gca,'Ydir','reverse');
    plot(s.Centroid(1),s.Centroid(2),'gx','MarkerSize',9);
    axis equal;
    subplot(3,3,4:6); hold on; box on;
    plot(angle_sorted,dist_sorted,'k-');
%     plot(angles,distances_raw,'r-');
    xlabel('angle');
    ylabel('distance from centroid');
    v = axis;
    axis([0 359 v(3) v(4)]);
     subplot(3,3,7:9); hold on; box on;
     plot([0:length(P1)-1],P1);
     xlabel('frequency domain');
     ylabel('|P1(f)|');
    
%     
%     % computation of the perimeter (sorted)
%     B = bwboundaries(roi.bw_mask_filled);
%     B = B{1};
%     xp = B(:,2);
%     yp = B(:,1);
%     
%     
%     % center the image
%     xp = xp - s.Centroid(1);
%     yp = yp - s.Centroid(2);
%     xp2 = xp2 - s.Centroid(1);
%     yp2 = yp2 - s.Centroid(2);
%     
%     % heuristic shitty approach
%     for i=1:length(xp)
%         angle(i) = atand(yp(i)/xp(i));
%         dist(i) = sqrt(xp(i)^2+yp(i)^2);
%         % convert to 0->360 degrees
%         if yp(i) >= 0 && xp(i) >= 0
%             quadrant(i) = 1; 
%         elseif yp(i) >= 0 && xp(i) < 0
%             quadrant(i) = 2;
%             angle(i) = 180 + angle(i);
%         elseif yp(i) < 0 && xp(i) < 0
%             quadrant(i) = 3;
%             angle(i) = 180 + angle(i);
%         elseif yp(i) < 0 && xp(i) >= 0
%             quadrant(i) = 4;
%             angle(i) = 360 + angle(i);
%         end
%     end
%     [angle_sorted,idx_sorted] = sort(angle);
%     dist_sorted = dist(idx_sorted);
%     % rescales distance values between 0 and 1 
%     tmp_mean = mean(dist_sorted);
%     tmp_std = std(dist_sorted);
%     dist_sorted = (dist_sorted-tmp_mean)/tmp_std;
%     
%     % mega awesome approach
%     angles = 0:1:359;
%     for i=1:length(angles)
%         xl = roi.Dmax * cosd(angles(i));
%         yl = roi.Dmax * sind(angles(i));
%         [allx,ally] = polyxpoly([0 xl],[0 yl],xp,yp);
%         allnorm = sqrt(allx.^2+ally.^2);
%         if ~isempty(allnorm)
%             distances(i) = max(allnorm);
%         else 
%             distances(i) = roi.eq_radius;
%         end
%     end
%     % rescales distance values between 0 and 1
%     distances_raw = distances;
%     tmp_mean = mean(distances);
%     tmp_std = std(distances);
%     distances = (distances-tmp_mean)/tmp_std;
%     
%     % fft stuff
%     Fs = 1; % sampling frequency
%     T = 1/Fs; % sampling period
%     L = 360; % length of signal
%     t = (0:L-1)*T; % time vector
%     f_ = fft(distances);
%     % two-sided spectrum P2
%     P2 = abs(f_/L);
%     P1 = P2(1:L/2+1);
%     P1(2:end-1) = 2*P1(2:end-1);
%     fvec = 10*Fs*(0:(L/2))/L;
%     %fvec2 = 10*Fs*(1:L)/L;
%     
%     figure;
%     subplot(1,3,1);
%     imshow(roi.data);
%     subplot(1,3,2);
%     imshow(roi.bw_mask);
%     subplot(1,3,3);
%     imshow(roi.bw_mask_filled);
%     
%     figure;
%     subplot(3,3,1);
%     imshow(roi.data);
%     subplot(3,3,2);
%     imshow(roi.bw_mask_filled);
%     subplot(3,3,3); hold on; box on;
%     %plot(xp2,yp2,'kx');
%     plot(xp,yp,'r-');
%     plot(0,0,'ro');
%     axis equal;
%     subplot(3,3,4:6); hold on; box on;
%     %plot(angle_sorted,dist_sorted,'k-');
%     plot(angles,distances_raw,'r-');
%     xlabel('angle');
%     ylabel('distance from centroid');
%     v = axis;
%     axis([0 359 v(3) v(4)]);
%     subplot(3,3,7:9); hold on; box on;
%     plot(fvec,P1);
%     xlabel('frequency domain');
%     ylabel('|P1(f)|');
%     

% clc;    % Clear the command window.
% close all;  % Close all figures (except those of imtool.)
% clear;  % Erase all existing variables.
% workspace;  % Make sure the workspace panel is showing.
% format longg;
% format compact;
% fontSize = 20;
% 
% angles = 0 : 10 : 360;
% x = cosd(angles);
% y = sind(angles);
% plot(x, y, 'ro-');
% xlabel('X', 'FontSize', fontSize);
% ylabel('Y', 'FontSize', fontSize);
% % Enlarge figure to full screen.
% set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
% grid on;
% xCenter = 0;
% yCenter = 0;
% for k = 1 : length(angles)
% 	% Draw a line from the center to the edge of the circle.
% 	line([xCenter x(k)], [yCenter y(k)], 'LineWidth', 3);
% 	% Calculate the angle of that line
% 	angle = atand(y(k)/x(k));
% 	% Convert to 0-360
% 	if y(k) >= yCenter && x(k) >= xCenter
% 		quadrant = 1;
% 	elseif y(k) >= yCenter && x(k) <= xCenter
% 		angle = 180 + angle;
% 		quadrant = 2;
% 	elseif y(k) <= yCenter && x(k) < xCenter
% 		angle = 180 + angle;
% 		quadrant = 3;
% 	elseif y(k) <= yCenter && x(k) >= xCenter
% 		angle = 360 + angle;
% 		quadrant = 4;
% 	end
% 	promptMessage = sprintf('This latest line, going from (0,0) to (%.8f, %.8f),\nis at angle %f, and in quadrant %d.', ...
% 		x(k), y(k), angle, quadrant);
% 	fprintf('%s\n', promptMessage);
% 	titleBarCaption = 'Continue?';
% 	button = questdlg(promptMessage, titleBarCaption, 'Continue', 'Cancel', 'Continue');
% 	if strcmp(button, 'Cancel')
% 		break;
% 	end
% end
