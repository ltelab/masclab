% small script to make general statistics on the computed features

%% input parameters
%clear all;  close all;

% --- insert path to the dir containing the processed masc particles here
main_dir = '/media/praz/Masc-Data/APRES3_2015/Raw_data/processed_09-Feb-2016/DATA/GOOD';
% --- insert desired time inverval here
tmin_vec = [2015 11 10 00 00 00];
tmax_vec = [2015 11 13 00 00 00];

tmin = datenum(tmin_vec);
tmax = datenum(tmax_vec);
%% loading snowflakes (.mat) in a structure data
% (to be run only once!)
file_list = dir(fullfile(main_dir,'*.mat'));
file_only_list = {file_list.name};
file_list = fullfile(main_dir,file_only_list);

% loop over the snowflakes and save each feature in a vector of the main
% structure data
% as an example, the field data.area will be a vector of size [N_particles]
% containing the area of all the snowflakes found in main_dir
for i=1:length(file_list)
    
    load(file_list{i});
    % load picture saving time
    t = roi.tnum;
    
    % load only pictures within the time interval
    if roi.tnum >= tmin && roi.tnum <= tmax
        
        % general
        data.name{i} = file_only_list{i};
        data.area(i) = roi.area;
        data.dim(i) = 0.5*double(roi.width + roi.height);
        data.Dmax(i) = roi.Dmax;
        data.centroid(i,:) = roi.centroid;
        data.mean_intens(i) = roi.mean_intens;
        data.max_intens(i) = roi.max_intens;
        data.lap(i) = roi.lap;
        data.complex(i) = roi.complex;
        data.std(i) = roi.std;
        data.local_std(i) = roi.local_std;
        data.local_std5(i) = roi.local_std5;
        data.local_std7(i) = roi.local_std7;
        data.bright.lap(i) = roi.new.lap;
        data.bright.std(i) = roi.new.std;
        data.bright.local_std(i) = roi.new.local_std;
        data.bright.local_std5(i) = roi.new.local_std5;
        data.bright.local_std7(i) = roi.new.local_std7;   
        data.tnum(i) = roi.tnum; 
        data.cam(i) = roi.cam;
        data.id(i) = roi.id;
        data.xhi(i) = roi.xhi;
        data.nb_holes(i) = roi.nb_holes;
        % ellipse fit
        data.a(i) = roi.E.a;
        data.b(i) = roi.E.b;
        data.theta(i) = roi.E.theta;
        if roi.E.a >= roi.E.b   
            data.eccentricity(i) = sqrt(1 - roi.E.b^2/roi.E.a^2); 
        else 
            data.eccentricity(i) = sqrt(1 - roi.E.a^2/roi.E.b^2);
        end
        % skeleton
        data.skel_p_ratio(i) = roi.skel.p_ratio;
        data.skel_A_ratio(i) = roi.skel.A_ratio;
        data.skel_N_ends(i) = roi.skel.N_ends;
        data.skel_N_junctions(i) = roi.skel.N_junctions;
        % Haralick
        data.haralick_contrast(i) = roi.H.Contrast;
        data.haralick_corr(i) = roi.H.Correlation;
        data.haralick_energy(i) = roi.H.Energy;
        data.haralick_homog(i) = roi.H.Homogeneity;
        % Relation with other ROIs
        %data.n_roi(i) = roi.n_roi;
        data.ratio(i) = roi.area_focus_ratio;
        % Environmental conditions
        %data.T(i) = roi.thygan.T;
        %data.RH(i) = roi.thygan.RH;
        
        if ~isempty(roi.fallspeed)
            data.fallspeed(i) = roi.fallspeed;
        else
            data.fallspeed(i) = NaN;
        end
    
    end
    
end

% resolution is between ~30 microns/pixels on the x-axis assuming that the
% flakes are in focus
data.dim_mm = data.dim * 30 / 1000; %[mm]
data.Dmax_mm = data.Dmax * 30 / 1000;
data.area_mm = data.dim * (30/1000)^2; % [mm]

% reorder fields of the structure in ABC order
data = orderfields(data);


%% comparison
close all;
clearvars -except data data2
figure; hold on; box on;
histogram(data.dim_mm,linspace(0,10,100),'Normalization','probability');
set(gca,'XGrid','on');
xlabel('D_{max} [mm]');
set(gca,'Fontsize',14);
set(gca,'YTick',[]);
figure; hold on; box on; grid on;
histogram(data2.dim_mm,linspace(0,10,100),'Normalization','probability');
set(gca,'XGrid','on');
xlabel('D_{max} [mm]');
set(gca,'Fontsize',14);
set(gca,'YTick',[]);

figure; hold on; box on;
histogram(data.eccentricity,linspace(0,1,100),'Normalization','probability');
set(gca,'XGrid','on');
xlabel('Aspect ratio');
set(gca,'Fontsize',14);
set(gca,'YTick',[]);
figure; hold on; box on; grid on;
histogram(data2.eccentricity,linspace(0,1,100),'Normalization','probability');
set(gca,'XGrid','on');
xlabel('Aspect ratio');
set(gca,'Fontsize',14);
set(gca,'YTick',[]);

figure; hold on; box on;
histogram(data.complex+0.1,linspace(1,4,100),'Normalization','probability');
set(gca,'XGrid','on');
xlabel('Complexity \chi');
set(gca,'Fontsize',14);
set(gca,'YTick',[]);
figure; hold on; box on; grid on;
histogram(data2.complex+0.1,linspace(1,4,100),'Normalization','probability');
set(gca,'XGrid','on');
xlabel('Complexity \chi');
set(gca,'Fontsize',14);
set(gca,'YTick',[]);

figure; hold on; box on;
histogram(data.mean_intens,linspace(0,1,100),'Normalization','probability');
set(gca,'XGrid','on');
xlabel('Brightness');
set(gca,'Fontsize',14);
set(gca,'YTick',[]);
figure; hold on; box on; grid on;
histogram(data2.mean_intens,linspace(0,1,100),'Normalization','probability');
set(gca,'XGrid','on');
xlabel('Brightness');
set(gca,'Fontsize',14);
set(gca,'YTick',[]);

figure; hold on; box on;
histogram(data.fallspeed,linspace(0,10,100),'Normalization','probability');
set(gca,'XGrid','on');
xlabel('Fallspeed [m/s]');
set(gca,'Fontsize',14);
set(gca,'YTick',[]);
figure; hold on; box on; grid on;
histogram(data2.fallspeed,linspace(0,10,100),'Normalization','probability');
set(gca,'XGrid','on');
xlabel('Fallspeed [m/s]');
set(gca,'Fontsize',14);
set(gca,'YTick',[]);


%% some statistics illustration
close all;
clearvars -except data data2
% you can add more hereafter

% fig 1 : some statistics about size and shape
figure;
subplot(221); hold on; grid on;
histogram(data.dim_mm);
title('PSD');
xlabel('dim. cropped window [mm]');
ylabel('N');
subplot(222); hold on; grid on;
histogram(data.area_mm);
title('Area of flakes');
xlabel('A [mm^2]');
ylabel('N');
subplot(223); hold on; grid on;
histogram(data.eccentricity);
title('Eccentricity');
xlabel('E');
ylabel('N');
subplot(224); hold on; grid on;
histogram(data.complex);
title('Part. complexity');
xlabel('C');
ylabel('N');
grid on;

% fig 2 : average brightness of snoflakes
figure; hold on; grid on;
histogram(data.mean_intens);
title('Mean brightness of snowflakes');
xlabel('Mean brightness [0-1]');
ylabel('N');

% fig 3 : repartition of snowflake location on camera's images
figure;
subplot(2,3,1); hold on; box on;
histogram(data.centroid(data.cam==0,1));
v = axis;
plot([320 320],[v(3) v(4)],'r--');
set(gca,'XLim',[0,2448]);
subplot(2,3,2); hold on; box on;
histogram(data.centroid(data.cam==1,1));
title('distribution of snowflake along the x-axis [0-2448]');
set(gca,'XLim',[0,2448]);
subplot(2,3,3); hold on; box on;
histogram(data.centroid(data.cam==2,1));
v = axis;
plot([2448-320 2448-320],[v(3) v(4)],'r--');
set(gca,'XLim',[0,2448]);
subplot(2,3,4); hold on; box on;
histogram(data.centroid(data.cam==0,2));
ylabel('camera 0');
v = axis;
plot([400 400],[v(3) v(4)],'r--');
plot([2048-400 2048-400],[v(3) v(4)],'r--');
set(gca,'XLim',[0 2048]);
set(gca,'XDir','reverse');
view(90,-90);
subplot(2,3,5); hold on; box on;
histogram(data.centroid(data.cam==1,2));
ylabel('camera 1');
title('distribution of snowflake along the y-axis [0-2048]');
v = axis;
plot([400 400],[v(3) v(4)],'r--');
plot([2048-400 2048-400],[v(3) v(4)],'r--');
set(gca,'XLim',[0 2048]);
set(gca,'XDir','reverse');
view(90,-90);
subplot(2,3,6); hold on; box on;
histogram(data.centroid(data.cam==2,2));
ylabel('camera 2');
v = axis;
plot([400 400],[v(3) v(4)],'r--');
plot([2048-400 2048-400],[v(3) v(4)],'r--');
set(gca,'XLim',[0 2048]);
set(gca,'XDir','reverse');
view(90,-90);

%%
close all;
clearvars -except data data2

startDate = datenum('11-11-2015');
endDate = datenum('12-16-2015');
xData = linspace(startDate,endDate,10);
%thygan = load_wfj_thygan('20150210000000','20150313000000',0);

figure;
histogram(data2.tnum,1000);
ylabel('N flakes');
v = axis;
axis([startDate endDate v(3) v(4)]);
ax = gca;
ax.XTick = xData;
datetick('x','keeplimits','keepticks');

