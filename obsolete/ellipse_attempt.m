% I am adding a random line here to see if the commit function of the EPFL git is working correctly
% script sandbox circumbscribed ellipse + ellipse fit + ...
clear all; close all;

% first we load a snowflake example
path = '/home/praz/Documents/MASC/sample_snow_20150620/sample_cool/PROCESSED_old/DATA/GOOD';
path2save = '/home/praz/Documents/MASC/ellipse_fit_examples';
my_dir = dir(fullfile(path,'*.mat'));
name_list = {my_dir.name};
%dataname = '2015.06.20_09.24.32_flake_47067_cam_2.mat';

for i=1:length(name_list)
    
    fprintf('image %u out of %u being processed... \n',i,length(name_list));
    
    load(fullfile(path,name_list{i}));

    [y,x] = find(roi.bw_mask_filled==1);
    xP = roi.x_perim;
    yP = roi.y_perim;

    E1 = fit_ellipse_pca(x,y,0.9,0);
    E2 = fit_ellipse_around(x,y,0);
    E3 = fit_ellipse_inside(x,y,xP,yP,0);


    t = linspace(0,2*pi);
    xt1 = E1.X0 + cos(E1.theta).*E1.a.*cos(t) - sin(E1.theta).*E1.b.*sin(t);
    yt1 = E1.Y0 + sin(E1.theta).*E1.a.*cos(t) + cos(E1.theta).*E1.b.*sin(t);

    xt2 = E2.X0 + cos(E2.theta).*E2.a.*cos(t) - sin(E2.theta).*E2.b.*sin(t);
    yt2 = E2.Y0 + sin(E2.theta).*E2.a.*cos(t) + cos(E2.theta).*E2.b.*sin(t);

    xt3 = E3.X0 + cos(E3.theta).*E3.a.*cos(t) - sin(E3.theta).*E3.b.*sin(t);
    yt3 = E3.Y0 + sin(E3.theta).*E3.a.*cos(t) + cos(E3.theta).*E3.b.*sin(t);

    fig = figure('Visible','off'); hold on; axis equal;
    plot(x,y,'k.');
    plot(xt1,yt1,'r-');
    plot(xt2,yt2,'b-');
    plot(E2.xhull,E2.yhull,'b--');
    plot(xt3,yt3,'g-');
    plot(E1.X0,E1.Y0,'rx');
    plot(E2.X0,E2.Y0,'bo');
    plot(E3.X0,E3.Y0,'gv');

    if isempty(dir(path2save))
        mkdir(path2save);
    end

    print(fig,fullfile(path2save,name_list{i}),'-dpng','-r80');


end


% centroid = [mean(x), mean(y)];
% % centroid(1) = 0.5*(max(x)-min(x));
% % centroid(2) = 0.5*(max(y)-min(y));
% xc = x - centroid(1);
% yc = y - centroid(2);
% [coeff,~,~] = pca([xc,yc]); 
% 
% % retrieve the tilt
% theta = acos(coeff(1,1));
% if coeff(2,1) < 0 
%     theta = -theta;
% end
% % untilt the data
% xc0 = cos(-theta).*xc - sin(-theta).*yc;
% yc0 = sin(-theta).*xc + cos(-theta).*yc;
% % maximal focal distance guess
% maxd = max(xc0) - min(xc0);
% 
% % instead of working on the whole data sample, it is easier and faster to
% % keep only the data points forming a polygon encompassing the data sample
% k = convhull(x,y);
% xh = x(k);
% yh = y(k); %x,y hull start from the leftmost point and then counterclockwise 
% %centering
% xhc = xh - centroid(1);
% yhc = yh - centroid(2);
% % % untilt the data
% xhc0 = cos(-theta).*xhc - sin(-theta).*yhc;
% yhc0 = sin(-theta).*xhc + cos(-theta).*yhc;
% 
% 
% d = 0:1:maxd;
% s = zeros(length(d),1);
% A = zeros(length(d),1);
% for i=1:length(d)
%     s(i) = max(sqrt((xhc0-0.5*d(i)).^2 + yhc0.^2) + sqrt((xhc0+0.5*d(i)).^2 + yhc0.^2));
%     A(i) = pi*s(i)/4 * sqrt(s(i)^2-d(i)^2);
% end
% [A,idx] = min(A);
% s = s(idx);
% d = d(idx);
% 
% % from optimal s, we can derive a and b
% a = s/2;
% b = sqrt((s/2)^2 - (d/2)^2);
% 
% % lets do the same for the biggest inscribed ellipse (thus we have to use
% % perimeter
% xPc = xP - centroid(1);
% yPc = yP - centroid(2);
% xPc0 = cos(-theta).*xPc - sin(-theta).*yPc;
% yPc0 = sin(-theta).*xPc + cos(-theta).*yPc;
% 
% d = 0:1:maxd;
% s2 = zeros(length(d),1);
% A2 = zeros(length(d),1);
% for i=1:length(d)
%     s2(i) = min(sqrt((xPc0-0.5*d(i)).^2 + yPc0.^2) + sqrt((xPc0+0.5*d(i)).^2 + yPc0.^2));
%     A2(i) = pi*s2(i)/4 * sqrt(s2(i)^2-d(i)^2);
% end
% [A2,idx] = max(A2);
% s2 = s2(idx);
% d = d(idx);
% 
% a2 = s2/2;
% b2 = sqrt((s2/2)^2 - (d/2)^2);
% 
% % we can try to shift the ellipse a little bit
% % bottom-top direction
% % idx_top = find(yc0>0);
% % idx_bottom = find(yc0<0);
% % xc0_top = xc0(idx_top);
% % yc0_top = yc0(idx_top);
% % xc0_bottom = xc0(idx_bottom);
% % yc0_bottom = yc0(idx_bottom);
% % 
% % s_top = max(sqrt((xc0_top-0.5*d).^2 + yc0_top.^2) + sqrt((xc0_top+0.5*d).^2 + yc0_top.^2));
% % s_bottom = max(sqrt((xc0_bottom-0.5*d).^2 + yc0_bottom.^2) + sqrt((xc0_bottom+0.5*d).^2 + yc0_bottom.^2));
% % y_offset = 0;
% % 
% % if s_top > s_bottom
% %     while s_top > s_bottom
% %         y_offset = y_offset-1;
% %         s_top = max(sqrt((xc0_top-0.5*d).^2 + (yc0_top+y_offset).^2) + sqrt((xc0_top+0.5*d).^2 + (yc0_top+y_offset).^2));
% %         s_bottom = max(sqrt((xc0_bottom-0.5*d).^2 + (yc0_bottom+y_offset).^2) + sqrt((xc0_bottom+0.5*d).^2 + (yc0_bottom+y_offset).^2));
% %     end
% % end
% 
% % build ellipse for illustration
% t = linspace(0,2*pi);
% xt = centroid(1) + cos(theta).*a.*cos(t) - sin(theta).*b.*sin(t);
% yt = centroid(2) + sin(theta).*a.*cos(t) + cos(theta).*b.*sin(t);
% 
% 
% figure;
% subplot(1,3,1);
% imshow(roi.data);
% subplot(1,3,2:3); hold on; axis equal;
% plot(x,y,'k.');
% plot(xh,yh,'b--');
% plot(xt,yt,'r-');
% set(gca,'Ydir','reverse');
% 
% t = linspace(0,2*pi);
% xt = a.*cos(t);
% yt = b.*sin(t);
% 
% t2 = linspace(0,2*pi);
% xt2 = a2.*cos(t);
% yt2 = b2.*sin(t);
% 
% figure; hold on; axis equal;
% plot(xc0,yc0,'k.');
% plot(xhc0,yhc0,'b--');
% plot(xt,yt,'r-');
% plot(xt2,yt2,'r-');
%set(gca,'Ydir','reverse');




% 
% F1_block = zeros(length(name_list),4);
% accuracy_block = zeros(length(name_list),4);
% 
% tt_pca = 0;
% tt_ls = 0;
% tt_matlab = 0;

% for ii=1:length(name_list)
% 
% % dataname = '2015.06.20_09.12.49_flake_46499_cam_2.mat'; % example with a negative tilt theta
% %dataname = '2015.06.20_09.19.57_flake_46855_cam_1.mat'; % example with a positive tilt theta
% %dataname = '2015.06.20_09.24.32_flake_47067_cam_2.mat';
% 
% dataname = name_list{ii};
% load(fullfile(path,dataname));
% 
% % PCA
% t_pca = tic;
% [y,x] = find(roi.bw_mask_filled==1);
% E1 = fit_ellipse_pca(y,x,0);
% tt_pca = tt_pca + toc(t_pca);
% 
% % Least Squares
% t_ls = tic;
% E2 = fit_ellipse_ls(roi.x_perim,roi.y_perim,0);
% tt_ls = tt_ls + toc(t_ls);
% 
% % Matlab
% t_matlab = tic;
% E3 = fit_ellipse_matlab(roi.bw_mask_filled,0);
% tt_matlab = tt_matlab + toc(t_matlab);
% 
% fprintf('end of iteration %u ! \n',ii);

%end


% figure;
% imshow(roi.data);

% x = roi.x_perim;
% y = roi.y_perim;
% 
% mean_x = mean(x);
% mean_y = mean(y);
% x = x-mean_x;
% y = y-mean_y;
% 
% 
% f = -1;
% %A = sum(X)'\(X'*X);
% X = [x.^2, x.*y, y.^2, x, y];
% A = sum(X)/(X'*X);
% [a,b,c,d,e] = deal(A(1),A(2),A(3),A(4),A(5));
% 
% % ellipse equation :
% % A(1)*x^2 + A(2)*x*y + A(3)*y^2 + A(4)*x + A(5)*y + F = 0
% % let us transform this into a paramteric (plotable) form
% M0 = [f, d/2, e/2; ...
%     d/2, a b/2; ...
%     e/2, b/2, c];
% 
% M = [a, b/2; ...
%     b/2, c];
% 
% 
% 
% %eigenvalues of M
% eigM = eig(M);
% if abs(eigM(1) - a) <= abs(eigM(1) -c)
%     lambda1 = eigM(1);
%     lambda2 = eigM(2); 
% else
%     lambda1 = eigM(2);
%     lambda2 = eigM(1);
% end
% 
% % parameters of the parametric form
% u = sqrt(-det(M0)/(det(M)*lambda1)); % semi grand-axe
% v = sqrt(-det(M0)/(det(M)*lambda2)); % semi petit-axe
% h = (b*e - 2*c*d)/(4*a*c - b^2); % x-position of the center
% k = (b*d - 2*a*e)/(4*a*c - b^2); % y-position of the center
% theta = acot((a-c)/b)/2; % the orientation (tilt from x-axis)
% 
% % reshift the whole ellipse to the original place
% h = h + mean_x;
% k = k + mean_y;
% 
% % illustration
% t = linspace(0,2*pi);
% xt = h + cos(theta).*u.*cos(t) - sin(theta).*v.*sin(t);
% yt = k + sin(theta).*u.*cos(t) + cos(theta).*v.*sin(t);
% 
% % 2 extremities of the semi petit-axe
% P1 = [h - sin(theta)*v, k + cos(theta)*v]; 
% P2 = [h + sin(theta)*v, k - cos(theta)*v];
% 
% % 2 extremitites of the semi grand-axe
% P3 = [h + cos(theta)*u, k + sin(theta)*u];
% P4 = [h - cos(theta)*u, k - sin(theta)*u];
% 
% figure(1);
% hold on;
% plot(xt,yt,'r-');
% plot(roi.x_perim,roi.y_perim,'kx');
% plot(h,k,'ro');
% plot(P1(1),P1(2),'ro');
% plot(P2(1),P2(2),'ro');
% plot(P3(1),P3(2),'ro');
% plot(P4(1),P4(2),'ro');
% plot([P4(1) P3(1)],[P4(2) P3(2)],'r-');
% plot([P2(1) P1(1)],[P2(2) P1(2)],'r-');
% %set(gca,'Ydir','reverse');
% axis equal;
% 
% % compute overlap ellipse - snowflake
% [y_im,x_im] = find(roi.bw_mask_filled > -1);
% y_im = y_im - mean_y;
% x_im = x_im - mean_x;
% x_imr = cos(-theta).*x_im - sin(-theta).*y_im;
% y_imr = sin(-theta).*x_im + cos(-theta).*y_im;
% 
% [y_flake,x_flake] = find(roi.bw_mask_filled == 1);
% y_flake = y_flake - mean_y;
% x_flake = x_flake - mean_x;
% x_flaker = cos(-theta).*x_flake - sin(-theta).*y_flake;
% y_flaker = sin(-theta).*x_flake + cos(-theta).*y_flake;
% 
% figure(2);
% subplot(1,4,1);  hold on; axis equal;
% plot(x_imr,y_imr,'y.');
% plot(x_flaker,y_flaker,'k.');
% t = linspace(0,2*pi);
% xt = u.*cos(t);
% yt = v.*sin(t);
% plot(xt,yt,'b-');
% 
% overlap = zeros(size(x_imr));
% ellipse_only = zeros(size(x_imr));
% flake_only = zeros(size(x_imr));
% nothing = zeros(size(x_imr));
% 
% for i=1:length(x_imr)
%     
%     is_in_flake =  ~isempty(find(x_flaker==x_imr(i),1)) && ~isempty(find(y_flaker==y_imr(i),1));
%     is_in_ellipse = x_imr(i)^2/u^2 + y_imr(i)^2/v^2 <= 1;
%     
%     if is_in_flake && is_in_ellipse
%         
%         overlap(i) = 1;
%         
%     elseif is_in_flake
%         
%         flake_only(i) = 1;
%         
%     elseif is_in_ellipse 
%         
%         ellipse_only(i) = 1;
%         
%     else
%         
%         nothing(i) = 1;
%         
%     end
%       
% end
%     
% % compute accuracy and F1 score based on confusion matrix
% accuracy(1) = (sum(overlap(:)) + sum(nothing(:))) / length(x_imr);
% F1(1) = 2*sum(overlap(:)) / (2*sum(overlap(:)) + sum(ellipse_only(:)) + sum(flake_only(:))); 
% 
% 
% plot(x_imr(overlap==1),y_imr(overlap==1),'g.');
% plot(x_imr(flake_only==1),y_imr(flake_only==1),'m.');
% plot(x_imr(ellipse_only==1),y_imr(ellipse_only==1),'r.');
% plot(x_imr(nothing==1),y_imr(nothing==1),'k.');
% 
% 
% %% Least Squares with centroid applied on all the pixels (not only Perim.)
% clearvars -except roi accuracy F1 ii F1_block accuracy_block name_list path
% 
% x = roi.x_perim;
% y = roi.y_perim;
% 
% [y_all,x_all] = find(roi.bw_mask_filled == 1);
% 
% mean_x = mean(x_all);
% mean_y = mean(y_all);
% x = x-mean_x;
% y = y-mean_y;
% 
% 
% f = -1;
% %A = sum(X)'\(X'*X);
% X = [x.^2, x.*y, y.^2, x, y];
% A = sum(X)/(X'*X);
% [a,b,c,d,e] = deal(A(1),A(2),A(3),A(4),A(5));
% 
% % ellipse equation :
% % A(1)*x^2 + A(2)*x*y + A(3)*y^2 + A(4)*x + A(5)*y + F = 0
% % let us transform this into a paramteric (plotable) form
% M0 = [f, d/2, e/2; ...
%     d/2, a b/2; ...
%     e/2, b/2, c];
% 
% M = [a, b/2; ...
%     b/2, c];
% 
% 
% 
% %eigenvalues of M
% eigM = eig(M);
% if abs(eigM(1) - a) <= abs(eigM(1) -c)
%     lambda1 = eigM(1);
%     lambda2 = eigM(2); 
% else
%     lambda1 = eigM(2);
%     lambda2 = eigM(1);
% end
% 
% % parameters of the parametric form
% u = sqrt(-det(M0)/(det(M)*lambda1)); % semi grand-axe
% v = sqrt(-det(M0)/(det(M)*lambda2)); % semi petit-axe
% h = (b*e - 2*c*d)/(4*a*c - b^2); % x-position of the center
% k = (b*d - 2*a*e)/(4*a*c - b^2); % y-position of the center
% theta = acot((a-c)/b)/2; % the orientation (tilt from x-axis)
% 
% % reshift the whole ellipse to the original place
% h = h + mean_x;
% k = k + mean_y;
% 
% % illustration
% t = linspace(0,2*pi);
% xt = h + cos(theta).*u.*cos(t) - sin(theta).*v.*sin(t);
% yt = k + sin(theta).*u.*cos(t) + cos(theta).*v.*sin(t);
% 
% % 2 extremities of the semi petit-axe
% P1 = [h - sin(theta)*v, k + cos(theta)*v]; 
% P2 = [h + sin(theta)*v, k - cos(theta)*v];
% 
% % 2 extremitites of the semi grand-axe
% P3 = [h + cos(theta)*u, k + sin(theta)*u];
% P4 = [h - cos(theta)*u, k - sin(theta)*u];
% 
% figure(1);
% plot(xt,yt,'r--');
% % plot(roi.x_perim,roi.y_perim,'kx');
% % plot(h,k,'ro');
% % plot(P1(1),P1(2),'ro');
% % plot(P2(1),P2(2),'ro');
% % plot(P3(1),P3(2),'ro');
% % plot(P4(1),P4(2),'ro');
% % plot([P4(1) P3(1)],[P4(2) P3(2)],'r--');
% % plot([P2(1) P1(1)],[P2(2) P1(2)],'r--');
% %set(gca,'Ydir','reverse');
% 
% [y_im,x_im] = find(roi.bw_mask_filled > -1);
% y_im = y_im - mean_y;
% x_im = x_im - mean_x;
% x_imr = cos(-theta).*x_im - sin(-theta).*y_im;
% y_imr = sin(-theta).*x_im + cos(-theta).*y_im;
% 
% [y_flake,x_flake] = find(roi.bw_mask_filled == 1);
% y_flake = y_flake - mean_y;
% x_flake = x_flake - mean_x;
% x_flaker = cos(-theta).*x_flake - sin(-theta).*y_flake;
% y_flaker = sin(-theta).*x_flake + cos(-theta).*y_flake;
% 
% figure(2);
% subplot(1,4,2);  hold on; axis equal;
% plot(x_imr,y_imr,'y.');
% plot(x_flaker,y_flaker,'k.');
% t = linspace(0,2*pi);
% xt = u.*cos(t);
% yt = v.*sin(t);
% plot(xt,yt,'b-');
% 
% overlap = zeros(size(x_imr));
% ellipse_only = zeros(size(x_imr));
% flake_only = zeros(size(x_imr));
% nothing = zeros(size(x_imr));
% 
% for i=1:length(x_imr)
%     
%     is_in_flake =  ~isempty(find(x_flaker==x_imr(i),1)) && ~isempty(find(y_flaker==y_imr(i),1));
%     is_in_ellipse = x_imr(i)^2/u^2 + y_imr(i)^2/v^2 <= 1;
%     
%     if is_in_flake && is_in_ellipse
%         
%         overlap(i) = 1;
%         
%     elseif is_in_flake
%         
%         flake_only(i) = 1;
%         
%     elseif is_in_ellipse 
%         
%         ellipse_only(i) = 1;
%         
%     else
%         
%         nothing(i) = 1;
%         
%     end
%       
% end
%     
% % compute accuracy and F1 score based on confusion matrix
% accuracy(2) = (sum(overlap(:)) + sum(nothing(:))) / length(x_imr);
% F1(2) = 2*sum(overlap(:)) / (2*sum(overlap(:)) + sum(ellipse_only(:)) + sum(flake_only(:))); 
% 
% 
% plot(x_imr(overlap==1),y_imr(overlap==1),'g.');
% plot(x_imr(flake_only==1),y_imr(flake_only==1),'m.');
% plot(x_imr(ellipse_only==1),y_imr(ellipse_only==1),'r.');
% plot(x_imr(nothing==1),y_imr(nothing==1),'k.');
% 
% 
% %% PCA approach
% 
% clearvars -except roi accuracy F1 ii F1_block accuracy_block name_list path
% x_P = roi.x_perim;
% y_P = roi.y_perim;
% [y,x] = find(roi.bw_mask_filled == 1);
% % all points have the same weight, so that the centroid is given by :
% centroid = [mean(x), mean(y)];
% 
% x_center = x - centroid(1);
% y_center = y - centroid(2);
% 
% %[U,S,V] = svd([x,y]);
% [coeff,score,latent] = pca([x_center,y_center]); 
% 
% % retrieve the tilt
% theta = acos(coeff(1,1));
% if coeff(2,1) < 0 
%     theta = -theta;
% end
% angle = theta*180/pi
% 
% figure(1);
% %plot(x_P-centroid(1),y_P-centroid(2),'rx');
% plot([centroid(1),centroid(1)+0.1*latent(1)*coeff(1,1)],[centroid(2),centroid(2)+0.1*latent(1)*coeff(2,1)],'b-');
% plot([centroid(1),centroid(1)+0.1*latent(2)*coeff(1,2)],[centroid(2),centroid(2)+0.1*latent(2)*coeff(2,2)],'b-');
% 
% 
% % I take all my data and I tilt them to match pc1-pc2 with x-y axis
% xr = cos(-theta).*x_center - sin(-theta).*y_center;
% yr = sin(-theta).*x_center + cos(-theta).*y_center;
% 
% figure(3);
% subplot(2,1,1);
% histogram(xr);
% subplot(2,1,2);
% histogram(yr);
% 
% % take a/b ratio as (x-axis maximal distance)/(x-axis maximal distance ratio)
% a_min = quantile(xr,0.90); a_max = quantile(xr,0.10);
% a = 0.5 * (a_max-a_min);
% b_min = quantile(yr,0.90); b_max = quantile(yr,0.10);
% b = 0.5 * (b_max-b_min);
% axis_ratio = a/b;
% 
% % compute final a,b fixing area(ellipse) = area(snowflake)
% a = sqrt(axis_ratio * sum(roi.bw_mask_filled(:)) / pi);
% b = a / axis_ratio;
% 
% % illustration
% t = linspace(0,2*pi);
% xt = centroid(1) + cos(theta).*a.*cos(t) - sin(theta).*b.*sin(t);
% yt = centroid(2) + sin(theta).*a.*cos(t) + cos(theta).*b.*sin(t);
% 
% 
% figure(1);
% plot(xt,yt,'b-');
% 
% [y_im,x_im] = find(roi.bw_mask_filled > -1);
% y_im = y_im - centroid(2);
% x_im = x_im - centroid(1);
% x_imr = cos(-theta).*x_im - sin(-theta).*y_im;
% y_imr = sin(-theta).*x_im + cos(-theta).*y_im;
% 
% [y_flake,x_flake] = find(roi.bw_mask_filled == 1);
% y_flake = y_flake - centroid(2);
% x_flake = x_flake - centroid(1);
% x_flaker = cos(-theta).*x_flake - sin(-theta).*y_flake;
% y_flaker = sin(-theta).*x_flake + cos(-theta).*y_flake;
% 
% figure(2); 
% subplot(1,4,3);  hold on; axis equal;
% plot(x_imr,y_imr,'y.');
% plot(x_flaker,y_flaker,'k.');
% t = linspace(0,2*pi);
% xt = a.*cos(t);
% yt = b.*sin(t);
% plot(xt,yt,'b-');
% 
% overlap = zeros(size(x_imr));
% ellipse_only = zeros(size(x_imr));
% flake_only = zeros(size(x_imr));
% nothing = zeros(size(x_imr));
% 
% for i=1:length(x_imr)
%     
%     is_in_flake =  ~isempty(find(x_flaker==x_imr(i),1)) && ~isempty(find(y_flaker==y_imr(i),1));
%     is_in_ellipse = x_imr(i)^2/a^2 + y_imr(i)^2/b^2 <= 1;
%     
%     if is_in_flake && is_in_ellipse
%         
%         overlap(i) = 1;
%         
%     elseif is_in_flake
%         
%         flake_only(i) = 1;
%         
%     elseif is_in_ellipse 
%         
%         ellipse_only(i) = 1;
%         
%     else
%         
%         nothing(i) = 1;
%         
%     end
%       
% end
%     
% % compute accuracy and F1 score based on confusion matrix
% accuracy(3) = (sum(overlap(:)) + sum(nothing(:))) / length(x_imr);
% F1(3) = 2*sum(overlap(:)) / (2*sum(overlap(:)) + sum(ellipse_only(:)) + sum(flake_only(:))); 
% 
% 
% plot(x_imr(overlap==1),y_imr(overlap==1),'g.');
% plot(x_imr(flake_only==1),y_imr(flake_only==1),'m.');
% plot(x_imr(ellipse_only==1),y_imr(ellipse_only==1),'r.');
% plot(x_imr(nothing==1),y_imr(nothing==1),'k.');
% 
%    
% %% MATLAB approach
% clearvars -except roi accuracy F1 ii F1_block accuracy_block name_list path
% a = 0.5*roi.major_axis_length;
% b = 0.5*roi.minor_axis_length;
% theta = -roi.orientation*pi/180;
% 
% [y,x] = find(roi.bw_mask_filled == 1);
% % all points have the same weight, so that the centroid is given by :
% centroid = [mean(x), mean(y)];
% 
% x_center = x - centroid(1);
% y_center = y - centroid(2);
% 
% [y_im,x_im] = find(roi.bw_mask_filled > -1);
% y_im = y_im - centroid(2);
% x_im = x_im - centroid(1);
% x_imr = cos(-theta).*x_im - sin(-theta).*y_im;
% y_imr = sin(-theta).*x_im + cos(-theta).*y_im;
% 
% [y_flake,x_flake] = find(roi.bw_mask_filled == 1);
% y_flake = y_flake - centroid(2);
% x_flake = x_flake - centroid(1);
% x_flaker = cos(-theta).*x_flake - sin(-theta).*y_flake;
% y_flaker = sin(-theta).*x_flake + cos(-theta).*y_flake;
% 
% figure(2); 
% subplot(1,4,4);  hold on; axis equal;
% plot(x_imr,y_imr,'y.');
% plot(x_flaker,y_flaker,'k.');
% t = linspace(0,2*pi);
% xt = a.*cos(t);
% yt = b.*sin(t);
% plot(xt,yt,'b-');
% 
% overlap = zeros(size(x_imr));
% ellipse_only = zeros(size(x_imr));
% flake_only = zeros(size(x_imr));
% nothing = zeros(size(x_imr));
% 
% for i=1:length(x_imr)
%     
%     is_in_flake =  ~isempty(find(x_flaker==x_imr(i),1)) && ~isempty(find(y_flaker==y_imr(i),1));
%     is_in_ellipse = x_imr(i)^2/a^2 + y_imr(i)^2/b^2 <= 1;
%     
%     if is_in_flake && is_in_ellipse
%         
%         overlap(i) = 1;
%         
%     elseif is_in_flake
%         
%         flake_only(i) = 1;
%         
%     elseif is_in_ellipse 
%         
%         ellipse_only(i) = 1;
%         
%     else
%         
%         nothing(i) = 1;
%         
%     end
%       
% end
%     
% % compute accuracy and F1 score based on confusion matrix
% accuracy(4) = (sum(overlap(:)) + sum(nothing(:))) / length(x_imr);
% F1(4) = 2*sum(overlap(:)) / (2*sum(overlap(:)) + sum(ellipse_only(:)) + sum(flake_only(:))); 
% 
% 
% plot(x_imr(overlap==1),y_imr(overlap==1),'g.');
% plot(x_imr(flake_only==1),y_imr(flake_only==1),'m.');
% plot(x_imr(ellipse_only==1),y_imr(ellipse_only==1),'r.');
% plot(x_imr(nothing==1),y_imr(nothing==1),'k.');
% 
% F1_block(ii,:) = F1;
% accuracy_block(ii,:) = accuracy;
% fprintf('end of iteration %u \n',ii);


%end
