% script to analyze a snow event more in details
clear all; close all;

% headers
% data_dir = '/media/praz/MyData/MASC/processed_27-Apr-2016/GOOD';
data_dir = '/media/praz/MyData/MASC/Davos_test_LTE/processed_13-May-2016/GOOD';
tstart_vec = [2017 03 06 00 00 00]; 
tstart_str = datestr(tstart_vec,'yyyymmddHHMMSS');
tstart = datenum(tstart_vec);
tstop_vec  = [2017 03 07 00 00 00];
tstop_str  = datestr(tstop_vec,'yyyymmddHHMMSS');
tstop  = datenum(tstop_vec);
res = 35; % 1 pixel = ??? mu

to_load = false;

if to_load

    % load data
    X = [];
    Xname = {};
    Xt = [];
    dir_list = uploaddirs(data_dir,tstart_vec,tstop_vec);
    for i=1:numel(dir_list)

        [X_tmp,Xlab,Xname_tmp,Xt_tmp] = load_processed_data(dir_list{i},tstart_str,tstop_str);
        if ~isempty(X_tmp)
            X = [X; X_tmp];
            Xname = [Xname; Xname_tmp];
            Xt = [Xt; Xt_tmp];
        end

    end

else 
    
    load('prediction/data/newFormat/classif_datastruct_Verbier_20170306.mat');
    Xt = data.Xt;
    Xlab = data.Xlab;
    Xname = data.Xname;
    X = data.X;
    
end

% if a .mat file loaded, filter flakes included in the time window
idx = find(Xt>=tstart & Xt<=tstop);
X = X(idx,:);
Xt = Xt(idx);
Xname = Xname(idx);
fprintf('%u snowflakes found in the desired time interval \n',length(idx));

% keep good snowflakes only
xhi_thresh = 8.5;
idx = find(X(:,6)>=xhi_thresh);
X = X(idx,:);
Xt = Xt(idx);
Xname = Xname(idx);
fprintf('%u remaining after blurry snowflakes filtered out (xhi=%2.1f)\n',length(idx),xhi_thresh);

 
if true
    figure;
    histogram(Xt);
    set(gca,'Xlim',[tstart tstop]);
    set(gca,'XTick',tstart:datenum([0 0 1 0 0 0]):tstop);
    set(gca,'XTicklabel',''); 
    datetick('x','dd.mm-HH:MM','keepticks','keeplimits');
end


%% general illustrations
% shift theta to [0 90]
theta = X(:,10);
theta(theta<0) = -theta(theta<0);

% discard large fallspeed
fs = X(:,7);
fs(fs>10) = 10;

if false 
    figure;
    subplot(2,3,1); hold on; box on; grid on;
    histogram(X(:,3).*res/1000);
    xlabel('Dmax [mm]');
    set(gca,'Xlim',[0 10]);
    subplot(2,3,2); hold on; hold on; box on; grid on;
    title(sprintf('EVENT : %s to %s',datestr(tstart,'yyyy-mm-dd-HHAM'),datestr(tstop,'yyyy-mm-dd-HHAM')));
    histogram(X(:,9));
    xlabel('Aspect ratio');
    subplot(2,3,3); hold on; box on; grid on;
    histogram(X(:,5));
    xlabel('Complexity');
    set(gca,'Xlim',[0.8 3]);
    subplot(2,3,4); hold on; box on; grid on;
    histogram(theta);
    xlabel('Orientation');
    set(gca,'Xlim',[0 90]);
    subplot(2,3,5); hold on; box on; grid on;
    histogram(fs);
    xlabel('fallspeed [m/s]');
    %set(gca,'Xscale','log');
end
%% dynamic time window

% variables of interest
Dmax = X(:,3)*res/1000;
AR = X(:,9);
complex = X(:,5);
snowclass = X(:,1);


% time window
twin_vec = [0 0 0 0 1 0];
twin = datenum(twin_vec);

% how many points minimum we want in the time window to consider it as non-nan
Npts_thresh = 10;

% time shift
tshift_vec = [0 0 0 0 1 0];
tshift = datenum(tshift_vec);

tbot = tstart;
ttop = tbot + twin;

idx_in = find(Xt>=tbot & Xt<ttop);
N(1) = numel(idx_in);
t(1) = tbot + twin/2;

Dmax_t{1} = Dmax(idx_in);
AR_t{1} = AR(idx_in);
complex_t{1} = complex(idx_in);
theta_t{1} = theta(idx_in);
fs_t{1} = fs(idx_in);

for i=1:6
    sclass(i,1) = sum(snowclass(idx_in)==i);
end
N_dry(1) = sum(X(idx_in,15)==0);
N_wet(1) = sum(X(idx_in,15)==1);
Avg_riming(1) = mean(X(idx_in,14));


k = 1;
max_idx_in = length(idx_in);

while tbot < tstop
    
    k = k+1;
    tbot = tbot + tshift;
    ttop = tbot + twin;
    
    idx_in = find(Xt>=tbot & Xt<ttop);
    N(k) = numel(idx_in);
    t(k) = tbot + twin/2;
 
    Dmax_t{k} = Dmax(idx_in);
    AR_t{k} = AR(idx_in);
    complex_t{k} = complex(idx_in);
    theta_t{k} = theta(idx_in);
    fs_t{k} = fs(idx_in);
    
    %snowflake class
    for i=1:6
        sclass(i,k) = sum(snowclass(idx_in)==i);
    end
    
    %melting snow
    N_dry(k) = sum(X(idx_in,15)==0 & X(idx_in,1)==4);
    N_wet(k) = sum(X(idx_in,15)==1 & X(idx_in,1)==4);
    
    %average degree of riming
    Avg_riming(k) = mean(X(idx_in,14));
    
    if length(idx_in) > max_idx_in
        max_idx_in = length(idx_in);
    end
    
end

Dmax_mat = nan(max_idx_in,length(Dmax_t));
AR_mat = nan(max_idx_in,length(AR_t));
complex_mat = nan(max_idx_in,length(complex_t));
fs_mat = nan(max_idx_in,length(fs_t));
theta_mat = nan(max_idx_in,length(theta_t));

for i=1:length(Dmax_t)
    Dmax_vec = Dmax_t{i}';
    if length(Dmax_vec) < max_idx_in
        Dmax_vec(end+1:max_idx_in) = NaN;
    end
    Dmax_mat(:,i) = Dmax_vec;
    
    AR_vec = AR_t{i}';
    if length(AR_vec) < max_idx_in
        AR_vec(end+1:max_idx_in) = NaN;
    end
    AR_mat(:,i) = AR_vec;
    
    complex_vec = complex_t{i}';
    if length(complex_vec) < max_idx_in
        complex_vec(end+1:max_idx_in) = NaN;
    end
    complex_mat(:,i) = complex_vec;
    
    fs_vec = fs_t{i}';
    if length(fs_vec) < max_idx_in
        fs_vec(end+1:max_idx_in) = NaN;
    end
    fs_mat(:,i) = fs_vec;   
    
    theta_vec = theta_t{i}';
    if length(theta_vec) < max_idx_in
        theta_vec(end+1:max_idx_in) = NaN;
    end
    theta_mat(:,i) = theta_vec;
       
end

% remove nan entries
 
% t_nan(idx_nan) = [];
% Dmax_mat(:,idx_nan) = [];
% AR_mat(:,idx_nan) = [];
% complex_mat(:,idx_nan) = [];
% fs_mat(:,idx_nan) = [];
% bright_mat(:,idx_nan) = [];
% theta_mat(:,idx_nan) = [];


% plot every continuous section iteratively
idx_notnan = find(nansum(Dmax_mat)>Npts_thresh);
[~,idxs_stop] = find(diff(idx_notnan)~=1);
idxs_start = [idx_notnan(1)  idx_notnan([idxs_stop+1])];
idxs_stop = [idx_notnan([idxs_stop]) idx_notnan(end)];


% snowflake class : normalization
for i=1:length(N)
    sclass(:,i) = sclass(:,i)/N(i);
end
% riming : normalization between 0 and 1
Avg_riming = (Avg_riming-1)/4;

% creating the suppa-duppa colormap
c = hsv(6);
cmasc = c;
cmasc(4,:) = c(6,:);
cmasc(6,:) = [49,163,84]./255;
cdvd(1,:) = c(6,:);
cdvd(2,:) = c(3,:);
cdvd(3,:) = c(5,:);
cdvd(4,:) = c(1,:);
cdvd(5,:) = c(2,:);
cdvd(6,:) = [241,105,19]./255;
cdvd(7,:) = [140,45,4]./255;
cdvd(8,:) = c(4,:);

% fig1 : MASC classification
figure(1);
subplot(5,1,1:2); box on; hold on;
title(sprintf('EVENT : %s to %s',datestr(tstart,'yyyy-mm-dd-HH:MM'),datestr(tstop,'yyyy-mm-dd-HH:MM')));
h = area(t,sclass'); 
for i=1:6
    h(i).FaceColor = cmasc(i,:);
end
legend('SM','CO','DE','AG','GR','CC');
ylabel('proportions');
set(gca,'Xlim',[tstart tstop]);
set(gca,'XTick',tstart:datenum([0 0 0 1 0 0]):tstop);
set(gca,'XTicklabel',''); 
%datetick('x','HHPM','keepticks','keeplimits');
%ax = gca;
%ax.XTickLabelRotation=55;
set(gca,'Ylim',[0 1]);
% riming degree
subplot(5,1,3); grid on; box on;
plot(t,Avg_riming,'k-');
ylabel('mean riming degree');
set(gca,'Xlim',[tstart tstop]);
set(gca,'XTick',tstart:datenum([0 0 0 1 0 0]):tstop);
set(gca,'XTicklabel',''); 

% melting snow 
subplot(5,1,4); grid on; box on;
plot(t,N_wet./(N_wet+N_dry),'k-');
ylabel('% of wet snow');
set(gca,'Xlim',[tstart tstop]);
set(gca,'XTick',tstart:datenum([0 0 0 1 0 0]):tstop);
set(gca,'XTicklabel',''); 
%datetick('x','HH:MM','keepticks','keeplimits');
%ax = gca;
%ax.XTickLabelRotation=55;
set(gca,'Ylim',[0 1]);

% atmospheric variables
atmdata_filepath = '/home/praz/Documents/MASC/masclab/atm_vars/order_45574_data.txt';
% tstr_start = '20160101000000';
% tstr_stop = '20160131235959';
atm = load_atm_vars(atmdata_filepath,tstart_str,tstop_str,0);

subplot(5,1,5);
plot(atm.tnum,atm.T);
addaxis(atm.tnum,atm.RH);
addaxislabel(2,'RH [%]');
%[hAx,~,~] = plotyy(atm.tnum,atm.T,atm.tnum,atm.RH);
ylabel('T [C]');
%ylabel(hAx(2),'RH [%]');
set(gca,'Xlim',[tstart tstop]);
set(gca,'XTick',tstart:datenum([0 0 0 1 0 0]):tstop);
set(gca,'XTicklabel',''); 
datetick('x','HH:MM','keepticks','keeplimits');
ax = gca;
ax.XTickLabelRotation=55;


% alternative with 3 classes
% sclass2(1,:) = sclass(4,:); % agg
% sclass2(2,:) = sclass(5,:); % grau
% sclass2(3,:) = sclass(2,:)+sclass(3,:)+sclass(6,:); % crystals
% for i=1:length(N)
%     tmp_sum = sum(sclass2(:,i));
%     sclass2(:,i) = sclass2(:,i)/tmp_sum;
% end
% figure; box on;
% h = area(t,sclass2');
% h(1).FaceColor = c(4,:);
% h(2).FaceColor = c(5,:);
% h(3).FaceColor = c(3,:);
% set(gca,'Xlim',[tstart tstop]);
% set(gca,'XTick',tstart:datenum([0 0 0 1 0 0]):tstop);
% set(gca,'XTicklabel',''); 
% datetick('x','HH','keepticks','keeplimits');
% ax = gca;
% ax.XTickLabelRotation=55;
% set(gca,'Ylim',[0 1]);


t_nan = t;

% fig2 : MASC descriptors I
figure(2);
subplot(311); hold on; grid on; box on; 
title(sprintf('EVENT : %s to %s',datestr(tstart,'yyyy-mm-dd-HH:MM'),datestr(tstop,'yyyy-mm-dd-HH:MM')));
plot(t,N,'k.-');
addaxis(atm.tnum,atm.precip);
addaxislabel(2,'precip [mm/h]');
%xlabel('time');
ylabel('number of images');
set(gca,'Xlim',[tstart tstop]);
set(gca,'XTick',tstart:datenum([0 0 0 1 0 0]):tstop);
set(gca,'XTicklabel',''); 
%datetick('x','HH:MM','keepticks','keeplimits');

% ax = gca;
ax.XTickLabelRotation=55;

subplot(312); hold on; box on;
for i=1:length(idxs_start)
plot_area_curve(t(idxs_start(i):idxs_stop(i))',Dmax_mat(:,idxs_start(i):idxs_stop(i))','','Dmax');
end
set(gca,'Xlim',[tstart tstop]);
set(gca,'XTick',tstart:datenum([0 0 0 1 0 0]):tstop);
set(gca,'XTicklabel',''); 
%datetick('x','HHPM','keepticks','keeplimits');
set(gca,'Ylim',[0 4]);
ax = gca;
ax.XTickLabelRotation=55;

subplot(313); hold on; box on;
for i=1:length(idxs_start)
plot_area_curve(t(idxs_start(i):idxs_stop(i))',AR_mat(:,idxs_start(i):idxs_stop(i))','','Aspect ratio',[8,69,148]/255,[66,146,198]/255,[107,174,214]/255);
end
%plot_area_curve(t_nan',AR_mat','','Aspect Ratio',[8,69,148]/255,[66,146,198]/255,[107,174,214]/255);
set(gca,'Xlim',[tstart tstop]);
set(gca,'XTick',tstart:datenum([0 0 0 1 0 0]):tstop);
set(gca,'XTicklabel',''); 
datetick('x','HH:MM','keepticks','keeplimits');
set(gca,'Ylim',[0.3 1]);
ax = gca;
ax.XTickLabelRotation=55;

% fig3: MASC descriptors II
figure(3);
subplot(311); hold on; box on;
title(sprintf('EVENT : %s to %s',datestr(tstart,'yyyy-mm-dd-HH:MM'),datestr(tstop,'yyyy-mm-dd-HH:MM')));
for i=1:length(idxs_start)
plot_area_curve(t(idxs_start(i):idxs_stop(i))',complex_mat(:,idxs_start(i):idxs_stop(i))','','Complexity',[35,139,69]/255,[102,194,164]/255,[178,226,226]/255);
end
%plot_area_curve(t_nan',complex_mat','','Complexity',[35,139,69]/255,[102,194,164]/255,[178,226,226]/255);
set(gca,'Xlim',[tstart tstop]);
set(gca,'XTick',tstart:datenum([0 0 0 1 0 0]):tstop);
set(gca,'XTicklabel',''); 
%datetick('x','HHPM','keepticks','keeplimits');
set(gca,'Ylim',[0.9 2]);
ax = gca;
ax.XTickLabelRotation=55;

subplot(312); hold on; box on;
for i=1:length(idxs_start)
plot_area_curve(t(idxs_start(i):idxs_stop(i))',theta_mat(:,idxs_start(i):idxs_stop(i))','','Orientation',[136,86,167]/255,[158,188,218]/255,[224,236,244]/255);
end
%plot_area_curve(t_nan',theta_mat','','Orientation',[136,86,167]/255,[158,188,218]/255,[224,236,244]/255);
set(gca,'Xlim',[tstart tstop]);
set(gca,'XTick',tstart:datenum([0 0 0 1 0 0]):tstop);
set(gca,'XTicklabel',''); 
%datetick('x','HHPM','keepticks','keeplimits');
set(gca,'Ylim',[0 90]);
ax = gca;
ax.XTickLabelRotation=55;

subplot(313); hold on; box on;
for i=1:length(idxs_start)
plot_area_curve(t(idxs_start(i):idxs_stop(i))',fs_mat(:,idxs_start(i):idxs_stop(i))','','Fallspeed',[99,99,99]/255,[189,189,189]/255,[240,240,240]/255);
end
%plot_area_curve(t_nan',fs_mat','time','Fallspeed',[99,99,99]/255,[189,189,189]/255,[240,240,240]/255);
set(gca,'Xlim',[tstart tstop]);
set(gca,'XTick',tstart:datenum([0 0 0 1 0 0]):tstop);
set(gca,'XTicklabel',''); 
datetick('x','HH:MM','keepticks','keeplimits');
set(gca,'Ylim',[0 2]);
ax = gca;
ax.XTickLabelRotation=55;



%% comparison with 2DVD
close all;
dvd1 = load_2DVD_classif(tstart_str);
dvd2 = load_2DVD_classif(tstop_str);
dvd.t = [dvd1.t; dvd2.t];
dvd.N = [dvd1.N; dvd2.N];
dvd.classif = [dvd1.classif; dvd2.classif];
idx_in = find(dvd.t >= tstart & dvd.t <= tstop);
dvd.t = dvd.t(idx_in);
dvd.N = dvd.N(idx_in);
dvd.classif = dvd.classif(idx_in);
dvd.label = dvd1.label;

% /!\ time shift in 2DVD to match MASC time
tshift_vec = [0 0 0 0 10 0];
tshift_num = datenum(tshift_vec);
dvd.t = dvd.t + tshift_num;

% fig4 : comparison with 2DVD 
figure(4);
% subplot : number of particles
subplot(5,1,1); hold on;
title(sprintf('EVENT : %s to %s',datestr(tstart,'yyyy-mm-dd-HH:MM'),datestr(tstop,'yyyy-mm-dd-HH:MM')));
plot(t,N,'k-');
addaxis(dvd.t,dvd.N);
addaxislabel(2,'# particles 2DVD');
% addaxis(atm.tnum,atm.precip);
% addaxislabel(3,'precip [mm/h]');
%xlabel('time');
ylabel('# images MASC');
set(gca,'Xlim',[tstart tstop]);
set(gca,'XTick',tstart:datenum([0 0 0 1 0 0]):tstop);
set(gca,'XTicklabel',''); 
box on;
% subplot: MASC classification
subplot(5,1,2:3); hold on;
h = area(t,sclass','EdgeColor','none'); 
for i=1:6
    h(i).FaceColor = cmasc(i,:);
end
legend('SM','CO','DE','AG','GR','CC');
ylabel('proportions');
set(gca,'Xlim',[tstart tstop]);
set(gca,'XTick',tstart:datenum([0 0 0 1 0 0]):tstop);
set(gca,'XTicklabel',''); 
set(gca,'Ylim',[0 1]);
set(gca,'YTick',[]);
box on;
% subplot: 2DVD classification
XX = dvd.t;
YY = [0 1];
ZZ = [dvd.classif dvd.classif];
subplot(5,1,4:5); hold on;
plot3(t,Avg_riming,linspace(10,10,length(t)),'k-','linewidth',1.5);
surface(XX,YY,ZZ','EdgeColor','none');
colormap(hsv(8));
colormap(cdvd);
set(gca,'Xlim',[tstart tstop]);
set(gca,'XTick',tstart:datenum([0 0 0 1 0 0]):tstop);
set(gca,'XTicklabel',''); 
datetick('x','HH:MM','keepticks','keeplimits');
set(gca,'Ylim',[0 1]);
set(gca,'YTick',[]);
box on;
legend('Riming degree');
ax = gca;
ax.XTickLabelRotation=55;



%%
% green
% 237,248,251
% 178,226,226
% 102,194,164
% 35,139,69
% purple
% 224,236,244
% 158,188,218
% 136,86,167
% grey
% 240,240,240
% 189,189,189
% 99,99,99



