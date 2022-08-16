% small script to display a few illustration of descriptors computation,
% for illustration purposes
clear all; close all;

fig_CSU_paper = false;
fig_Dmax_AR = true;

data_dir = '/home/praz/Documents/MASC/charts/CSU_paper';
%'/home/praz/Documents/MASC/sample_snow_20150620/sample_cool/processed_05-Feb-2016/DATA/GOOD';%'/home/praz/Documents/MASC/misc_examples/geometry';
data_filename = '2015.11.26_18.29.18_flake_417_cam_3.mat';%'2015.02.25_05.39.01_flake_52068_cam_2.mat';
load(fullfile(data_dir,data_filename));

if fig_Dmax_AR
    
    figure(99); subplot(131); imshow(roi.data);
    subplot(132); imshow(roi.bw_mask_filled);
    [tmp_Dmax, tmp_angle, tmp_A, tmp_B]  = compute_Dmax(roi.hull.xh,roi.hull.yh,1);
    % computation of Dmax_90
    tmp_im = imrotate(roi.bw_mask_filled,-tmp_angle);
    Dmax_0 = sum(tmp_im,1);
    Dmax_0 = Dmax_0(Dmax_0>0);
    Dmax_0 = length(Dmax_0);
    Dmax_90 = sum(tmp_im,2);
    Dmax_90 = Dmax_90(Dmax_90>0)
    Dmax_90 = length(Dmax_90);
    figure(99); subplot(133); imshow(tmp_im);
    
    
end



if fig_CSU_paper

    figure; hold on;
    %plot([-10 10],[-10 10],'k-');
    %imshow(roi.data); hold on;
    plot(roi.x_perim,roi.y_perim,'k-');
    % ellipse fit
    t = linspace(0,2*pi);
    theta = roi.E.theta * pi / 180;
    xt = roi.E.X0 + cos(-theta).*roi.E.a.*cos(t) - sin(-theta).*roi.E.b.*sin(t);
    yt = roi.E.Y0 + sin(-theta).*roi.E.a.*cos(t) + cos(-theta).*roi.E.b.*sin(t);
    plot(xt,yt,'r-');
    % ellipse out
    xt = roi.E_out.X0 + cos(-theta).*roi.E_out.a.*cos(t) - sin(-theta).*roi.E_out.b.*sin(t);
    yt = roi.E_out.Y0 + sin(-theta).*roi.E_out.a.*cos(t) + cos(-theta).*roi.E_out.b.*sin(t);
    plot(xt,yt,'b-');
    % ellipse in
    xt = roi.E_in.X0 + cos(-theta).*roi.E_in.a.*cos(t) - sin(-theta).*roi.E_in.b.*sin(t);
    yt = roi.E_in.Y0 + sin(-theta).*roi.E_in.a.*cos(t) + cos(-theta).*roi.E_in.b.*sin(t);
    plot(xt,yt,'g-');
    % center
    plot(roi.E.X0,roi.E.Y0,'r^');
    % conv. hull
    k = convhull(roi.x,roi.y);
    xh = roi.x(k);
    yh = roi.y(k);
    plot(xh,yh,'b--');
    % Dmax
    Dmax = 0;
    for i=1:length(xh)-1
        for j=i+1:length(xh)
            tmp_dist = sqrt((xh(i)-xh(j))^2 + (yh(i)-yh(j))^2);
            if tmp_dist > Dmax
                Dmax = tmp_dist;
                idxA = i;
                idxB = j;
            end
        end
    end 
    A = [xh(idxA); yh(idxA)];
    B = [xh(idxB); yh(idxB)];
    plot([A(1) B(1)],[A(2) B(2)],'m-.')

    set(gca,'Ydir','reverse');
    set(gca,'Xtick','');
    set(gca,'Ytick','');
    axis equal;

    % histogram of pixel values
    M = roi.data;
    M = M(:);
    M(M<10) = [];
    figure; hold on; box on;
    histogram(M,'DisplayStyle','stairs');
    plot([87.51 87.51],[0 1300],'r-');
    plot([53.20 53.20],[0 1300],'r--');
    plot([121.82 121.82],[0 1300],'r--');

end





% 
% 
% 
% %skeleton_props(roi.data,1);
% compute_symmetry_features(roi.bw_mask_filled,roi.Dmax,roi.eq_radius,1);
% 
% 
% % computation of the perimeter (sorted by link)
% B = bwboundaries(roi.bw_mask_filled);
% if length(B) > 1
%     disp('Warning : symmetry routine detected more than 1 particle in the filled mask image !');
%     disp('----> results are not reliable !');
% end
% B = B{1};
% xp = B(:,2);
% yp = B(:,1);
% 
% % computation of the centroid
% s = regionprops(roi.bw_mask_filled,'centroid');
% s = s(1,:);
% 
% % center the image (makes the calculations easier)
% %xp = xp - s.Centroid(1);
% %yp = yp - s.Centroid(2);
% 
% [roi.x,roi.y] = find(roi.bw_mask_filled);
% 
% % E = fit_ellipse_matlab(roi.bw_mask_filled,0);
% % E2 = fit_ellipse_around(roi.x,roi.y,0);
% % E3 = fit_ellipse_inside(roi.x,roi.y,roi.x_perim,roi.y_perim,1);
% 
% Dmax = compute_Dmax(roi.x_perim,roi.y_perim,1);
% C = fit_circle_around(roi.x_perim,roi.y_perim,1);
% 
% figure; 
% imshow(roi.data);
% 
% figure; hold on;
% plot(xp,yp,'k-');
% %plot(0,0,'ro');
% % 
% t = linspace(0,2*pi);
% xt = C.x0 + C.r*cos(t);
% yt = C.Y0 + C.r*sin(t);
% plot(xt,yt,'r-');
% 
% t = linspace(0,2*pi);
% xt = C.x0 + roi.eq_radius*cos(t);
% yt = C.Y0 + roi.eq_radius*sin(t);
% plot(xt,yt,'r-');
% 
% axis equal;
% % 
% % xt = s.Centroid(1) + cos(-E.theta).*E2.a.*cos(t) - sin(-E.theta).*E2.b.*sin(t);
% % yt = s.Centroid(2) + sin(-E.theta).*E2.a.*cos(t) + cos(-E.theta).*E2.b.*sin(t);
% % plot(xt,yt,'b-');
% % 
% % xt = s.Centroid(1) + cos(-E.theta).*E3.a.*cos(t) - sin(-E.theta).*E3.b.*sin(t);
% % yt = s.Centroid(2) + sin(-E.theta).*E3.a.*cos(t) + cos(-E.theta).*E3.b.*sin(t);
% % plot(xt,yt,'g-');
% 
% 
% 
% %legend('data','fit');
% set(gca,'Ydir','reverse');



% skel = skeleton_props(roi.data,1);
% 
% [x,y] = find(roi.bw_mask_filled); 
% E_out = fit_ellipse_around(x,y,1);
% E = fit_ellipse_ls(roi.x_perim,roi.y_perim,1);
% E2 = fit_ellipse_pca(x,y,0.9,1);
% E3 = fit_ellipse_matlab(roi.bw_mask_filled,1);