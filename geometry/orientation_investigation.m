clear all; close all;
    
dir_graupels = '/home/praz/Documents/MASC/masclab/events/graupels';
graupels_filenames = dir(fullfile(dir_graupels,'*.mat'));
graupels_filenames = {graupels_filenames.name}';
graupels_picnames = dir(fullfile(dir_graupels,'*.png'));
graupels_picnames = {graupels_picnames.name}';

dir_agg = '/home/praz/Documents/MASC/masclab/events/aggregates';
agg_filenames = dir(fullfile(dir_agg,'*.mat'));
agg_filenames = {agg_filenames.name}';
agg_picnames = dir(fullfile(dir_agg,'*.png'));
agg_picnames = {agg_picnames.name}';

dir_melt = '/home/praz/Documents/MASC/masclab/events/melting_snow';
melt_filenames = dir(fullfile(dir_melt,'*.mat'));
melt_filenames = {melt_filenames.name}';
melt_picnames = dir(fullfile(dir_melt,'*.png'));
melt_picnames = {melt_picnames.name}';

dir_small = '/home/praz/Documents/MASC/masclab/events/small_particles';
small_filenames = dir(fullfile(dir_small,'*.mat'));
small_filenames = {small_filenames.name}';
small_picnames = dir(fullfile(dir_small,'*.png'));
small_picnames = {small_picnames.name}';

dir_col = '/home/praz/Documents/MASC/masclab/events/columns';
col_filenames = dir(fullfile(dir_col,'*.mat'));
col_filenames = {col_filenames.name}';
col_picnames = dir(fullfile(dir_col,'*.png'));
col_picnames = {col_picnames.name}';

dir_dend = '/home/praz/Documents/MASC/masclab/events/dendrites';
dend_filenames = dir(fullfile(dir_dend,'*.mat'));
dend_filenames = {dend_filenames.name}';
dend_picnames = dir(fullfile(dir_dend,'*.png'));
dend_picnames = {dend_picnames.name}';

dir_bulros = '/home/praz/Documents/MASC/masclab/events/rosettes';
bulros_filenames = dir(fullfile(dir_bulros,'*.mat'));
bulros_filenames = {bulros_filenames.name}';
bulros_picnames = dir(fullfile(dir_bulros,'*.png'));
bulros_picnames = {bulros_picnames.name}';

data.theta_pca = [];
data.theta_matlab = [];
data.theta_rect = [];
data.theta_Dmax = [];
data.class_id = [];


for i=1:length(graupels_filenames)  
    file = fullfile(dir_graupels,graupels_filenames{i});
    load(file);
    
    E_pca = fit_ellipse_pca(roi.x,roi.y,0.9,0);
    E_m = fit_ellipse_matlab(roi.bw_mask,0);
    Rect = compute_rectangularity(roi.x,roi.y,roi.perim,0);
    [~,Dmax_theta] = compute_Dmax(roi.x_perim,roi.y_perim,0);
    
    data.theta_pca(end+1,1) = abs(E_pca.theta*180/pi);
    data.theta_matlab(end+1,1) = abs(E_m.theta*180/pi);
    data.theta_rect(end+1,1) = abs(Rect.theta*180/pi);
    data.theta_Dmax(end+1,1) = abs(Dmax_theta*180/pi);
    data.class_id(end+1,1) = 1;
  
end

for i=1:length(agg_filenames)  
    file = fullfile(dir_agg,agg_filenames{i});
    load(file);
    
    E_pca = fit_ellipse_pca(roi.x,roi.y,0.9,0);
    E_m = fit_ellipse_matlab(roi.bw_mask,0);
    Rect = compute_rectangularity(roi.x,roi.y,roi.perim,0);
    [~,Dmax_theta] = compute_Dmax(roi.x_perim,roi.y_perim,0);
    
    data.theta_pca(end+1,1) = abs(E_pca.theta*180/pi);
    data.theta_matlab(end+1,1) = abs(E_m.theta*180/pi);
    data.theta_rect(end+1,1) = abs(Rect.theta*180/pi);
    data.theta_Dmax(end+1,1) = abs(Dmax_theta*180/pi);
    data.class_id(end+1,1) = 2;
  
end

for i=1:length(melt_filenames)  
    file = fullfile(dir_melt,melt_filenames{i});
    load(file);
    
    E_pca = fit_ellipse_pca(roi.x,roi.y,0.9,0);
    E_m = fit_ellipse_matlab(roi.bw_mask,0);
    Rect = compute_rectangularity(roi.x,roi.y,roi.perim,0);
    [~,Dmax_theta] = compute_Dmax(roi.x_perim,roi.y_perim,0);
    
    data.theta_pca(end+1,1) = abs(E_pca.theta*180/pi);
    data.theta_matlab(end+1,1) = abs(E_m.theta*180/pi);
    data.theta_rect(end+1,1) = abs(Rect.theta*180/pi);
    data.theta_Dmax(end+1,1) = abs(Dmax_theta*180/pi);
    data.class_id(end+1,1) = 3;
  
end

for i=1:length(small_filenames)  
    file = fullfile(dir_small,small_filenames{i});
    load(file);
    
    E_pca = fit_ellipse_pca(roi.x,roi.y,0.9,0);
    E_m = fit_ellipse_matlab(roi.bw_mask,0);
    Rect = compute_rectangularity(roi.x,roi.y,roi.perim,0);
    [~,Dmax_theta] = compute_Dmax(roi.x_perim,roi.y_perim,0);
    
    data.theta_pca(end+1,1) = abs(E_pca.theta*180/pi);
    data.theta_matlab(end+1,1) = abs(E_m.theta*180/pi);
    data.theta_rect(end+1,1) = abs(Rect.theta*180/pi);
    data.theta_Dmax(end+1,1) = abs(Dmax_theta*180/pi);
    data.class_id(end+1,1) = 4;
  
end

for i=1:length(col_filenames)  
    file = fullfile(dir_col,col_filenames{i});
    load(file);
    
    E_pca = fit_ellipse_pca(roi.x,roi.y,0.9,0);
    E_m = fit_ellipse_matlab(roi.bw_mask,0);
    Rect = compute_rectangularity(roi.x,roi.y,roi.perim,0);
    [~,Dmax_theta] = compute_Dmax(roi.x_perim,roi.y_perim,0);
    
    data.theta_pca(end+1,1) = abs(E_pca.theta*180/pi);
    data.theta_matlab(end+1,1) = abs(E_m.theta*180/pi);
    data.theta_rect(end+1,1) = abs(Rect.theta*180/pi);
    data.theta_Dmax(end+1,1) = abs(Dmax_theta*180/pi);
    data.class_id(end+1,1) = 5;
  
end

for i=1:length(dend_filenames)  
    file = fullfile(dir_dend,dend_filenames{i});
    load(file);
    
    E_pca = fit_ellipse_pca(roi.x,roi.y,0.9,0);
    E_m = fit_ellipse_matlab(roi.bw_mask,0);
    Rect = compute_rectangularity(roi.x,roi.y,roi.perim,0);
    [~,Dmax_theta] = compute_Dmax(roi.x_perim,roi.y_perim,0);
    
    data.theta_pca(end+1,1) = abs(E_pca.theta*180/pi);
    data.theta_matlab(end+1,1) = abs(E_m.theta*180/pi);
    data.theta_rect(end+1,1) = abs(Rect.theta*180/pi);
    data.theta_Dmax(end+1,1) = abs(Dmax_theta*180/pi);
    data.class_id(end+1,1) = 6;
  
end

for i=1:length(bulros_filenames)  
    file = fullfile(dir_bulros,bulros_filenames{i});
    load(file);
    
    E_pca = fit_ellipse_pca(roi.x,roi.y,0.9,0);
    E_m = fit_ellipse_matlab(roi.bw_mask,0);
    Rect = compute_rectangularity(roi.x,roi.y,roi.perim,0);
    [~,Dmax_theta] = compute_Dmax(roi.x_perim,roi.y_perim,0);
    
    data.theta_pca(end+1,1) = abs(E_pca.theta*180/pi);
    data.theta_matlab(end+1,1) = abs(E_m.theta*180/pi);
    data.theta_rect(end+1,1) = abs(Rect.theta*180/pi);
    data.theta_Dmax(end+1,1) = abs(Dmax_theta*180/pi);
    data.class_id(end+1,1) = 7;
  
end

%%
clearvars -except data; close all;

c = colormap(prism(6)); hold on; box on; grid on;
listplot = [];
gscatter(data.theta_pca,data.theta_matlab,data.class_id);
set(gca,'Fontsize',12)
axis([0 90 0 90]);
xlabel('from ellipse fit [deg]');
ylabel('from Matlab regionprops() [deg]');
legend(listplot,{'graupels','aggregates','melting snow','small particles','columns/needles','dendrites/plates'});
title('Canting angle');

figure; hold on; box on; grid on;
listplot = [];
gscatter(data.theta_pca,data.theta_rect,data.class_id);
set(gca,'Fontsize',12)
axis([0 90 0 90]);
xlabel('from ellipse fit [deg]');
ylabel('from minimal rectangle [deg]');
title('Canting angle');
%legend(listplot,{'graupels','aggregates','melting snow','small particles','columns/needles','dendrites/plates'});

figure; hold on; box on; grid on;
listplot = [];
gscatter(data.theta_pca,data.theta_Dmax,data.class_id);
set(gca,'Fontsize',12)
axis([0 90 0 90]);
xlabel('from ellipse fit [deg]');
ylabel('from Dmax [deg]');
title('Canting angle');
%legend(listplot,{'graupels','aggregates','melting snow','small particles','columns/needles','dendrites/plates'});
;
classes = {'graupels','aggregates','melting snow','small particles','columns/needles','dendrites/plates'};
figure;
for i=1:6
    subplot(3,2,i); hold on; box on; grid on;
    histogram(data.theta_pca(data.class_id==i),14,'Normalization','probability','FaceColor',c(i,:));
    xlabel(classes{i});
    set(gca,'Xlim',[0 90]);
    set(gca,'Fontsize',12);
end

figure;
histogram(data.theta_pca);






%     clear all; close all;
%     imdir = '/home/praz/Documents/MASC/masclab/events/mixed_sample_1';
%     datafiles = dir(fullfile(imdir,'*.mat'));
%     datafiles = {datafiles.name}'; 
%     load(fullfile(imdir,datafiles{5}));   
%     x = roi.x;
%     y = roi.y;
%     im = roi.data;
%     
%     
%     rect = compute_rectangularity(roi.x,roi.y,roi.perim,1);
%     E_m = fit_ellipse_matlab(roi.bw_mask,1);
%     [Dmax,Dmax_theta] = compute_Dmax(roi.x_perim,roi.y_perim,1);
%     figure; imshow(roi.data);
    
    
    
    
    %% part with Dmax
%     xP = roi.x_perim;
%     yP = roi.y_perim;
%     k = convhull(xP,yP);
%     
%     tic;
%     Dmax = 0;
%     for i=1:length(k)-1
%         for j=i+1:length(k)
%             tmp_dist = sqrt((xP(k(i))-xP(k(j)))^2 + (yP(k(i))-yP(k(j)))^2);
%             if tmp_dist > Dmax
%                 Dmax = tmp_dist;
%                 idxA = k(i);
%                 idxB = k(j);
%             end
%         end
%     end
%     toc
%         
%     A = [xP(idxA); yP(idxA)];
%     B = [xP(idxB); yP(idxB)];
%     vec = B-A;
%     vec_0 = [1;0];
%     theta = acos(dot(vec,vec_0)/(norm(vec)*norm(vec_0)));
%     if theta > pi/2
%         theta = pi - theta;
%     end
%    
% 
%     figure; hold on; axis equal; box on;
%     plot(xP,yP,'k.');
%     plot(xP(k),yP(k),'b--');
%     plot([A(1) B(1)],[A(2) B(2)],'r--');
%     xlabel('x [pix]');
%     ylabel('y [pix]');
%     set(gca,'Ydir','reverse','Fontsize',14); 
%     title(sprintf('D_{max} = %2.2f pixels. \theta = %2.1f deg.',Dmax, theta*180/pi));

    
    
    %costheta = dot(a,b)/(norm(a)*norm(b));
    %theta = acos(costheta);
%     tic;
%     Dmax2 = max(pdist([xP(k),yP(k)]));
%     toc
    
%     E1=fit_ellipse_pca(x,y,0.9,1);
%     E2=fit_ellipse_ls(roi.x_perim,roi.y_perim,1);
%     E3=fit_ellipse_matlab(roi.bw_mask,1);