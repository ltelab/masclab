% snowflake selection analysis based on a few features
%close all; clear all;
clear all;  close all;

main_dir = '/media/kiko/Masc-Data/Data_Davos_MayJune_2015/proc_201506_21/DATA/GOOD';
% '/media/praz/Masc-Data/DATA_meteoswiss/PROCESSED_w_boarders_cleaning/DATA/GOOD';  %_w_boarders_cleaning
% '/home/praz/Documents/MASC/sample_snow_20150620/sample_cool/PROCESSED_with_cleaning_v2/DATA/GOOD';

%% GOOD detection
file_list = dir(fullfile(main_dir,'*.mat'));
file_only_list = {file_list.name};
file_list = fullfile(main_dir,file_only_list);

% the images
% good.data = cell(length(file_list),1);

% features based on raw picture
good.area = zeros(length(file_list),1);
good.dim = zeros(length(file_list),1);
good.mean_intens = zeros(length(file_list),1);
good.max_intens = zeros(length(file_list),1);
good.lap = zeros(length(file_list),1);
good.complex = zeros(length(file_list),1);
good.std = zeros(length(file_list),1);
good.local_std = zeros(length(file_list),1);
good.local_std5 = zeros(length(file_list),1);
good.local_std7 = zeros(length(file_list),1);
good.x_perim = cell(length(file_list),1);
good.y_perim = cell(length(file_list),1);
good.P2 = zeros(length(file_list),1);
good.A2 = zeros(length(file_list),1);
good.complex2 = zeros(length(file_list),1);
good.centroid = zeros(length(file_list),2);
good.x = cell(length(file_list),1);
good.y = cell(length(file_list),1);

% features based on enhanced picture
good.new.lap = zeros(length(file_list),1);
good.new.std = zeros(length(file_list),1);
good.new.local_std = zeros(length(file_list),1);
good.new.local_std5 = zeros(length(file_list),1);
good.new.local_std7 = zeros(length(file_list),1);

% Thygan retrieval
good.T = zeros(length(file_list),1);
good.RH = zeros(length(file_list),1);

% other stuff
good.name = cell(length(file_list),1);
good.xhi = zeros(length(file_list),1);
good.tnum = zeros(length(file_list),1);
good.cam = zeros(length(file_list),1);
good.id = zeros(length(file_list),1);

% new studd
good.xhi = zeros(length(file_list),1);

% loop over the snowflakes and save features in a vector
for i=1:length(file_list)
    
    load(file_list{i});
    good.name{i} = file_only_list{i};
    good.area(i) = roi.area;
    good.dim(i) = 0.5*double(roi.width + roi.height);
    good.centroid(i,:) = roi.centroid;
    good.mean_intens(i) = roi.mean_intens;
    good.max_intens(i) = roi.max_intens;
    good.lap(i) = roi.lap;
    good.complex(i) = roi.complex;
    good.std(i) = roi.std;
    good.local_std(i) = roi.local_std;
    good.local_std5(i) = roi.local_std5;
    good.local_std7(i) = roi.local_std7;
    
    good.new.lap(i) = roi.new.lap;
    good.new.std(i) = roi.new.std;
    good.new.local_std(i) = roi.new.local_std;
    good.new.local_std5(i) = roi.new.local_std5;
    good.new.local_std7(i) = roi.new.local_std7;
    
    %good.T(i) = roi.thygan.T;
    %good.RH(i) = roi.thygan.RH;
    good.tnum(i) = roi.tnum;
    
    good.cam(i) = roi.cam;
    good.id(i) = roi.id;
    
    %good.data{i} = roi.data;

    % retrieve all idx (x,y) of pixels within the snowflake
    [good.y{i},good.x{i}] = find(roi.bw_mask_filled ==1);
    good.x_perim{i} = roi.x_perim;
    good.y_perim{i} = roi.y_perim;
    
    
    % computation of new perimeter, area and complexity
    % good.P2(i) = my_perim(roi.x_perim,roi.y_perim);
    good.A2(i) = good.area(i) - 0.5*length(roi.x_perim);% + 0.5*good.P2(i);
    good.complex2(i) = length(roi.x_perim)/(2*sqrt(pi*good.A2(i)));
    
    % the magic parameter
    good.xhi(i) = roi.xhi;
    good.xhi_old(i) = log((good.lap(i)+good.new.lap(i))/2 * good.complex2(i) * (good.local_std(i) + good.new.local_std(i))/2 * good.dim(i)); 
    
end

good.dim_mm = good.dim * 30 / 1000; %[mm]



%% DATA ANALYSIS - HISTOGRAMS
% close all;
clearvars -except good good_old;

figure;
histogram(good.dim_mm);
title('PSD');
xlabel('dim window [mm]');
ylabel('N');

figure;
histogram(good.mean_intens);
title('Mean brightness of snowflakes');
xlabel('Mean brightness [0-1]');
ylabel('N');

figure;
subplot(2,1,1);
histogram(good.lap,linspace(0,30));
title('Lap before CLAHE');
ylabel('N');
subplot(2,1,2);
histogram(good.new.lap,linspace(0,30));
title('Lap after CLAHE');
xlabel('lap intens.');

figure;
subplot(2,1,1);
histogram(good.local_std);
title('local STD before CLAHE');
ylabel('N');
subplot(2,1,2);
histogram(good.new.local_std);
title('local STD after CLAHE');
xlabel('local STD');
ylabel('N');

figure;
subplot(2,1,1);
histogram(good.complex,linspace(0.85,7));
title('Initial complexity parameter');
ylabel('N');
subplot(2,1,2);
histogram(good.complex2,linspace(0.85,7));
title('Modified complexity parameter');
ylabel('N');

figure;
histogram(good.xhi);
xlabel('\xi');
ylabel('N');

figure;
subplot(2,1,1);
histogram(good.centroid(:,1));
title('Flake center location on the X-axis [0 2448]]');
subplot(2,1,2);
histogram(good.centroid(:,2));
title('Flake center location on the Y-axis [0 2048]');



%% DATA ANALYSIS - SCATTERPLOTS
% close all;
clearvars -except good;

figure;
plot(good.local_std,good.dim,'o','MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','r');

%% TEMPERATURE RH DEPENDANCE ANALYSIS
% close all;
clearvars -except good;

figure;
subplot(2,1,1);
histogram(good.T);
title('Local temperature histo');
subplot(2,1,2);
histogram(good.RH);
title('Local RH histo');

% time series
startDate = datenum('02-10-2015');
endDate = datenum('03-13-2015');
xData = linspace(startDate,endDate,10);
%thygan = load_wfj_thygan('20150210000000','20150313000000',0);

figure;
subplot(3,1,1);
histogram(good.tnum);
ylabel('N flakes');
v = axis;
axis([startDate endDate v(3) v(4)]);
ax = gca;
ax.XTick = xData;
datetick('x','keeplimits','keepticks');
% subplot(3,1,2);
% plot(thygan.tnum,thygan.T,'b-');
% v = axis;
% axis([startDate endDate -17 5]);
% ax = gca;
% ax.XTick = xData;
% datetick('x','keeplimits','keepticks');
% ylabel('T [C]');
% subplot(3,1,3);
% plot(thygan.tnum,thygan.RH,'r-');
% v = axis;
% axis([startDate endDate 0 120]);
% ylabel('RH [%]');
% ax = gca;
% ax.XTick = xData;
% datetick('x','keeplimits','keepticks');

% temperature vs geometry
T = -14.5:1:0.5;
for i=1:length(T)
    idx = find(good.T >= T(i)-0.5 & good.T < T(i)+0.5);
    mean_dim(i) = mean(good.dim(idx));
    std_dim(i) = std(good.dim(idx));
    med_dim(i) = median(good.dim(idx));
    mean_complex(i) = mean(good.complex2(idx));
    med_complex(i) = median(good.complex2(idx));
end

figure; hold on;
errorbar(T,mean_dim,std_dim,'ro');
plot(T,med_dim,'bo');

figure; hold on;
plot(T,mean_complex,'ro');
plot(T,med_complex,'bo');

%% SNOWFLAKE - SELECTION
% close all;
clearvars -except good;

%main_dir = '/media/praz/Masc-Data/DATA_meteoswiss/PROCESSED/DATA/GOOD';

figure;
histogram(good.xhi);
xlabel('\xi');
ylabel('N');

idx = find(good.xhi >= 10);% & good.tnum>datenum('20150801','yyyymmdd'));
namelist = cell(length(idx),1);
for i=1:length(idx)
    tmp_name = good.name{idx(i)};
    namelist{i} = tmp_name(1:end-4);
end

im_dir = '/media/kiko/Masc-Data/Data_Davos_MayJune_2015/proc_201506_21/IMAGES/GOOD';
%'/media/praz/Masc-Data/DATA_meteoswiss/PROCESSED_20151109/IMAGES/GOOD';
im_out = '/media/kiko/Masc-Data/Data_Davos_MayJune_2015/proc_201506_21/IMAGES/xi_above_10';
%'/media/praz/Masc-Data/DATA_meteoswiss/PROCESSED_20151109/IMAGES/xi_10-11';
if isempty(dir(im_out))
    mkdir(im_out);
end
data_dir = '/media/kiko/Masc-Data/Data_Davos_MayJune_2015/proc_201506_21/DATA/GOOD';
%'/media/praz/Masc-Data/DATA_meteoswiss/PROCESSED_20151109/DATA/GOOD';
data_out = '/media/kiko/Masc-Data/Data_Davos_MayJune_2015/proc_201506_21/DATA/xi_above_10';
%'/media/praz/Masc-Data/DATA_meteoswiss/PROCESSED_20151109/DATA/xi_10-11';
if isempty(dir(data_out))
    mkdir(data_out);
end
for i=1:length(idx)
    copyfile(fullfile(im_dir,strcat(namelist{i},'.png')),fullfile(im_out,strcat(namelist{i},'.png')));
    copyfile(fullfile(data_dir,strcat(namelist{i},'.mat')),fullfile(data_out,strcat(namelist{i},'.mat')));
end

%% DATA ANALYSIS - stats about ellipses, cam and location of flake in the pic
clearvars -except good;

% add new features to the list
for i=1:length(good.area)
    good.cam(i) = get_cam_id(good.name{i});
    good.id(i) = get_snowflake_id(good.name{i});
    good.E(i) = fit_ellipse_pca(good.x{i},good.y{i},0.9,0);
    good.E_out(i) = fit_ellipse_around(good.x{i},good.y{i},0);
    good.E_in(i) = fit_ellipse_inside(good.x{i},good.y{i},good.x_perim{i},good.y_perim{i},0);
end


%% ILLUSTRATION OF THE ABOVE
clearvars -except good; % close all;

figure;
subplot(2,3,1); hold on; box on;
histogram(good.centroid(good.cam==0,1));
v = axis;
plot([320 320],[v(3) v(4)],'r--');
set(gca,'XLim',[0,2448]);
subplot(2,3,2); hold on; box on;
histogram(good.centroid(good.cam==1,1));
title('distribution of snowflake along the x-axis [0-2448]');
set(gca,'XLim',[0,2448]);
subplot(2,3,3); hold on; box on;
histogram(good.centroid(good.cam==2,1));
v = axis;
plot([2448-320 2448-320],[v(3) v(4)],'r--');
set(gca,'XLim',[0,2448]);
subplot(2,3,4); hold on; box on;
histogram(good.centroid(good.cam==0,2));
ylabel('camera 0');
v = axis;
plot([400 400],[v(3) v(4)],'r--');
plot([2048-400 2048-400],[v(3) v(4)],'r--');
set(gca,'XLim',[0 2048]);
set(gca,'XDir','reverse');
view(90,-90);
subplot(2,3,5); hold on; box on;
histogram(good.centroid(good.cam==1,2));
ylabel('camera 1');
title('distribution of snowflake along the y-axis [0-2048]');
v = axis;
plot([400 400],[v(3) v(4)],'r--');
plot([2048-400 2048-400],[v(3) v(4)],'r--');
set(gca,'XLim',[0 2048]);
set(gca,'XDir','reverse');
view(90,-90);
subplot(2,3,6); hold on; box on;
histogram(good.centroid(good.cam==2,2));
ylabel('camera 2');
v = axis;
plot([400 400],[v(3) v(4)],'r--');
plot([2048-400 2048-400],[v(3) v(4)],'r--');
set(gca,'XLim',[0 2048]);
set(gca,'XDir','reverse');
view(90,-90);


% subplot(2,1,1);
% histogram(good.centroid(:,1));
% title('Flake center location on the X-axis [0 2448]]');
% subplot(2,1,2);
% histogram(good.centroid(:,2));
% title('Flake center location on the Y-axis [0 2048]');






