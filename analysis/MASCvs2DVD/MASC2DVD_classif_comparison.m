clear all; close all;

load_data = 1;

if load_data

    % load MASC data
    xhi_thresh = 0;
    Nmin_shift = 5;
    Nmin_interval = 5;  
    tshift = 6; %in min, for 2dvd (7min for 2016.06.16 15h-04h)  , 6min for 04.23 03h-12h
    
    % load('prediction/data/Xstruct_Davos_winter15-16_extended_last.mat')   >>>>  the one used for the paper!!!
    %masc = load('prediction/data/Xstruct_EastonAirport_single.mat');
    
    load('prediction/data/newFormat/classif_datastruct_Davos2015-16_triplet_withgraupfix_more_entries.mat');
    masc.Xt = data.Xt;
    masc.X = data.X;
    masc.Xlab = data.Xlab;
    masc.Xname = data.Xname;
    clear data;
    
    masc.Xt = datetime(masc.Xt,'ConvertFrom','datenum');
    tstart_vec = [2016 04 23 06 00 00]; % june, 04.23, 11.20
    tday = '20160423';
    %tstart_vec = [2016 04 23 03 00 00];
    tstop_vec  = [2016 04 23 12 00 00];
    %tstop_vec  = [2016 04 23 12 00 00];
    tstart = datetime(tstart_vec);
    tstop = datetime(tstop_vec);
    %tgrid = tstart-minutes(1):minutes(Nmin_interval):tstop; tgrid = tgrid';
    %tgrid2 = tgrid + minutes(1);
    tgrid = tstart:minutes(Nmin_shift):tstop; tgrid = tgrid';
    tgrid2 = tgrid + minutes(Nmin_interval);
    

    % load 2DVD data
    dvd = load_2DVD_classif_campaign(datestr(tstart-minutes(tshift),'yyyymmddHHMMSS'),datestr(tstop,'yyyymmddHHMMSS'));
    dvd.t = datetime(dvd.t,'ConvertFrom','datenum');
    % timeshift between masc and 2dvd
    dvd.t = dvd.t + minutes(tshift);
    
    % load atm variables
    atmdata_filepath = '/home/praz/Documents/MASC/masclab/atm_vars/order_45574_data.txt';
    atm = load_atm_vars(atmdata_filepath,datestr(tstart,'yyyymmddHHMMSS'),datestr(tstop,'yyyymmddHHMMSS'),0);
    
    % merge both on the same grid
    j=1;
    for i=1:length(tgrid)
        idx_masc = find(masc.Xt >= tgrid(i) & masc.Xt <= tgrid2(i) & masc.X(:,6)>=xhi_thresh);
        idx_masc_all = find(masc.Xt >= tgrid(i) & masc.Xt <= tgrid2(i));
        idx_dvd = find(dvd.t >= tgrid(i) & dvd.t <= tgrid2(i)); 
        
        data.t(j,1) = tgrid2(i);
        if ~isempty(idx_masc)
            data.Nmasc(j,1) = numel(idx_masc);
            data.Nmasc_all(j,1) = numel(idx_masc_all); 
            data.class(j,1) = mode(masc.X(idx_masc,1));
            for k=1:6
                data.sclass(j,k) = sum(masc.X(idx_masc,1)==k);
            end         
            data.melting(j,1) = mean(masc.X(idx_masc,16));
            data.riming(j,1) = mean(masc.X(idx_masc,15));
        else
            data.Nmasc(j,1) = 0;
            data.Nmasc_all(j,1) = numel(idx_masc_all);
            data.class(j,1) = NaN;
            for k=1:6
                data.sclass(j,k) = NaN;
            end
            data.melting(j,1) = NaN;
            data.riming(j,1) = NaN;
        end
        if ~isempty(idx_dvd)
            data.Ndvd(j,1) = sum(dvd.N(idx_dvd));
            scores_2dvd = zeros(8,1);
            for k=1:length(idx_dvd)
                that_class = dvd.classif(idx_dvd(k));
                if ~isnan(that_class)
                    scores_2dvd(that_class) = scores_2dvd(that_class)+dvd.N(idx_dvd(k));
                end
            end
            if max(scores_2dvd) > 0
                [~,data.classdvd(j,1)] = max(scores_2dvd);
            else
                data.classdvd(j,1) = NaN;
            end
        else
            data.Ndvd(j,1) = 0;
            data.classdvd(j,1) = NaN;
        end
        j=j+1;
    end
    
    % load pluvio data
    pluvio = load_DFIR_Pluvio2_day(tday);
    %pluvio.t = pluvio.t + minutes(tshift);
       
    % remove NaN in N
    data.Nmasc(isnan(data.Nmasc)) = 0;
    data.Ndvd(isnan(data.Ndvd)) = 0;

    fprintf('%u timesteps of %u minutes in commmon found\n',numel(data.t),Nmin_interval);
    
else  
    load('20151012_20160619__5min_5min_triplet_xhi9.mat');   
    %load('whole_season_MASC_2DVD_5min.mat');
end

% removing timesteps without enough points
N1 = 0;
N2 = 0;
data.sclass(data.Nmasc<=N1 | data.Ndvd<=N2,:) = NaN;
data.classdvd(data.Nmasc<=N1 | data.Ndvd<=N2) = NaN;
data.melting(data.Nmasc<=N1 | data.Ndvd<=N2) = NaN;
data.riming(data.Nmasc<=N1 | data.Ndvd<=N2) = NaN;


%% arranging and plotting

%snowflake class : normalization
for i=1:length(data.Nmasc)
    data.sclass(i,:) = data.sclass(i,:)/data.Nmasc(i);
end

% creating the suppa-duppa colormap
c = hsv(6); % rouge-jaune-vert-turquoise-bleu-rose
cmasc = c; 
cmasc(4,:) = c(6,:);
cmasc(6,:) = [49,163,84]./255;
cdvd(1,:) = c(6,:); %AGG - rose
cdvd(2,:) = c(3,:); %DEN - vert
cdvd(3,:) = c(5,:); %GRA - bleu
cdvd(4,:) = c(1,:); %SP  - rouge
cdvd(5,:) = c(2,:); %COL - jaune
cdvd(6,:) = [241,105,19]./255; %MS - orange
cdvd(7,:) = [140,45,4]./255; %R - brun
cdvd(8,:) = c(4,:); %RIM - turquoise

%cmasc(1,:) = cdvd(7,:);

%%
data.class_agg = data.class;
data.class_agg(data.class == 1) = 7;
data.class_agg(data.class == 2) = 5;
data.class_agg(data.class == 3) = 2;
data.class_agg(data.class == 4) = 1;
data.class_agg(data.class == 5) = 3;
data.class_agg(data.melting >= 0.3) = 6;
 
%% illustrations
close all;




% subplot : number of particles
% subplot(6,1,1); hold on;

% plot(data.t,data.Nmasc,'k-','linewidth',1.25);
% addaxis(data.t,data.Ndvd,'linewidth',1.25);
% addaxislabel(2,'# particles 2DVD');
% % addaxis(atm.tnum,atm.precip);
% % addaxislabel(3,'precip [mm/h]');
% %xlabel('time');
% ylabel('# images MASC');
% set(gca,'Xlim',[datenum(tstart) datenum(tstop)]);
% set(gca,'XTick',datenum(tstart):datenum([0 0 0 1 0 0]):datenum(tstop));
% set(gca,'XTicklabel','');
% datetick('x','HH','keepticks','keeplimits');
% box on;

% figure intensities MASC, 2DVD, Pluvio
if 0
figure;
subplot(211); hold on; grid on; box on;
plot(datenum(data.t),data.Nmasc,'r-','linewidth',1.25);
plot(datenum(data.t),data.Nmasc_all,'k-','linewidth',1.25);
addaxis(datenum(data.t),data.Ndvd,'linewidth',1.25);
addaxislabel(2,'# particles 2DVD');
xlabel('time');
ylabel('# images MASC');
set(gca,'Xlim',[datenum(tstart) datenum(tstop)]);
set(gca,'XTick',datenum(tstart):datenum([0 0 0 1 0 0]):datenum(tstop));
set(gca,'XTicklabel',''); 
datetick('x','HH','keepticks');

subplot(212);
plot(datenum(pluvio.t),pluvio.bucketRT,'k-');
xlabel('time');
ylabel('accuNRT total [mm]');
set(gca,'Xlim',[datenum(tstart) datenum(tstop)]);
set(gca,'XTick',datenum(tstart):datenum([0 0 0 1 0 0]):datenum(tstop));
set(gca,'XTicklabel',''); 
datetick('x','HH','keepticks');

end



% fig1 : comparison with 2DVD

if 0
myfig = figure('Visible','on','units','normalized','outerposition',[0 0 1 1]);
% subplot: MASC classification
subplot(4,1,1:2); hold on;
title(sprintf('%s to %s',datestr(tstart,'yyyy.mm.dd - HH:MM'),datestr(tstop,'yyyy.mm.dd - HH:MM')));
h = area(datenum(data.t),data.sclass,'EdgeColor','none'); 
for i=1:6
        h(i).FaceColor = cmasc(i,:);
end
legend('SP','CC','PC','AG','GR','CPC');
set(gca,'Xlim',[datenum(tstart) datenum(tstop)]);
set(gca,'XTick',datenum(tstart):datenum([0 0 0 1 0 0]):datenum(tstop));
set(gca,'XTicklabel',''); 
set(gca,'Ylim',[0 1]);
set(gca,'YTick',[0 0.5 1]);
ylabel('Proportions');
datetick('x','HH','keepticks','keeplimits')
box on;
set(gca,'Fontsize',30);
set(gca, 'Layer', 'top');

% subplot: T, riming, melting
subplot(4,1,3); hold on;
hline(1) = plot(datenum(data.t),data.riming,'k-','color',c(5,:),'linewidth',3);
hline(2) = plot(datenum(data.t),data.melting,'k-','color',cdvd(6,:),'linewidth',3);
set(gca,'Ylim',[0 1]);
addaxis(atm.tnum,atm.T,'k--','linewidth',3);
%addaxis(atm.tnum,atm.gust3s,'linewidth',2);
addaxislabel(1,'[-]');
addaxislabel(2,'T [C]');
legend('R_i','% wet','T');
% plot(atm.tnum,atm.T,'linewidth',1.25);
% axis([tstart tstop min(atm.T)-1 max(atm.T)+1]);
% %addaxis(atm.tnum,atm.RH);
% addaxis(atm.tnum,atm.winv,'linewidth',1.25);
% addaxislabel(1,'T [C]');
% %addaxislabel(2,'RH [%]');
% addaxislabel(2,'wind [m/s]');
set(gca,'Xlim',[datenum(tstart) datenum(tstop)]);
set(gca,'XTick',datenum(tstart):datenum([0 0 0 1 0 0]):datenum(tstop));
set(gca,'XTicklabel',''); 
datetick('x','HH','keepticks');
box on;
set(gca,'Fontsize',30);



% subplot: MASC aggregated to 2dvd classif.
% subplot(5,1,4);
% XX = [1;2;3;4;5;6;7;8; datenum(data.t)];
% YY = [0 1];
% ZZ = [[1 1; 2 2; 3 3; 4 4; 5 5; 6 6; 7 7; 8 8;]; [data.class_agg data.class_agg]];
% surface(XX,YY,ZZ','EdgeColor','none');
% colormap(cdvd);
% set(gca,'Xlim',[datenum(tstart) datenum(tstop)]);
% set(gca,'XTick',datenum(tstart):datenum([0 0 0 1 0 0]):datenum(tstop));
% set(gca,'XTicklabel',''); 
% datetick('x','HH','keepticks','keeplimits');
% set(gca,'Ylim',[0 1]);
% set(gca,'YTick',[]);
% box on;
% set(gca,'Fontsize',16);

% subplot: 2DVD classification
subplot(4,1,4); hold on;
XX = datenum(data.t);
YY = [0 1];
ZZ = [data.classdvd data.classdvd];
%plot3(datenum(data.t),data.riming,linspace(10,10,length(data.t)),'k-','linewidth',1.25);
% %plot3(t,N_wet./(N_wet+N_dry),linspace(10,10,length(t)),'k--','linewidth',1.25);
surface(XX,YY,ZZ','EdgeColor','none');
%unique_class = unique(data.classdvd);
%unique_class = unique_class(~isnan(unique_class));
colormap(cdvd);
set(gca,'Xlim',[datenum(tstart) datenum(tstop)]);
set(gca,'XTick',datenum(tstart):datenum([0 0 0 1 0 0]):datenum(tstop));
set(gca,'XTicklabel',''); 
datetick('x','HH','keepticks','keeplimits');
set(gca,'Ylim',[0 1]);
set(gca,'YTick',[]);
box on;
xlabel('Time UTC [h]');
set(gca,'Fontsize',30);
set(gca, 'Layer', 'top');
% legend('Riming degree');
% % subplot : melting snow
%subplot(6,1,5); box on;
% plot(t,N_wet./(N_wet+N_dry),'k-','color',[241,105,19]./255,'linewidth',1.25);
% legend('% of wet snow');
% set(gca,'Xlim',[tstart tstop]);
% set(gca,'XTick',tstart:datenum([0 0 0 1 0 0]):tstop);
% set(gca,'XTicklabel',''); 
% %datetick('x','HH:MM','keepticks','keeplimits');
% %ax = gca;
% %ax.XTickLabelRotation=55;
% set(gca,'Ylim',[0 1]);
% % subplot : atm. variables
%subplot(6,1,6); hold on;
% plot([tstart tstop],[0,0],'k--');
% plot(atm.tnum,atm.T,'linewidth',1.25);
% axis([tstart tstop min(atm.T)-1 max(atm.T)+1]);
% %addaxis(atm.tnum,atm.RH);
% addaxis(atm.tnum,atm.winv,'linewidth',1.25);
% addaxislabel(1,'T [C]');
% %addaxislabel(2,'RH [%]');
% addaxislabel(2,'wind [m/s]');
% set(gca,'Xlim',[tstart tstop]);
% set(gca,'XTick',tstart:datenum([0 0 0 1 0 0]):tstop);
% set(gca,'XTicklabel',''); 
% datetick('x','HH:MM','keepticks','keeplimits');
% ax = gca;
% ax.XTickLabelRotation=90;
% set(gcf,'Visible','on');
% box on;

% save figure
%savename = sprintf('MASC2DVD_%s',datestr(tstart,'yyyy-mm-dd'));
%export_fig(fullfile(savepath,savename),'-png');

%close all;

end


% fig2 : classification output + LWC [mm]
myfig = figure('Visible','on','units','normalized','outerposition',[0 0 1 1]);
% subplot: MASC classification
subplot(4,1,1:2); hold on;
title(sprintf('%s to %s',datestr(tstart,'yyyy.mm.dd - HH:MM'),datestr(tstop,'yyyy.mm.dd - HH:MM')));
h = area(datenum(data.t),data.sclass,'EdgeColor','none'); 
for i=1:6
        h(i).FaceColor = cmasc(i,:);
end
hline(1) = plot(datenum(data.t),data.riming,'w-','linewidth',3);
legend('SP','CC','PC','AG','GR','CPC');
set(gca,'Xlim',[datenum(tstart) datenum(tstop)]);
set(gca,'XTick',datenum(tstart):datenum([0 0 0 1 0 0]):datenum(tstop));
set(gca,'XTicklabel',''); 
set(gca,'Ylim',[0 1]);
set(gca,'YTick',[0 0.5 1]);
ylabel('Proportions');
datetick('x','HH','keepticks','keeplimits')
box on;
set(gca,'Fontsize',30);
set(gca, 'Layer', 'top');

subplot(4,1,3:4); hold on; grid on; box on;
idx = find(pluvio.t>=tstart,1,'first');
plot(datenum(pluvio.t),pluvio.bucketRT-pluvio.bucketRT(idx),'k-','linewidth',3);
%addaxis(datenum(pluvio.t),pluvio.accuNRT,'r-','linewidth',3);
plot(datenum(data.t),cumsum(data.Nmasc)./sum(data.Nmasc),'r-','linewidth',3);
%addaxislabel(2,'# particles MASC');
xlabel('Time UTC [h]');
ylabel('[mm]');
set(gca,'Xlim',[datenum(tstart) datenum(tstop)]);
set(gca,'XTick',datenum(tstart):datenum([0 0 0 1 0 0]):datenum(tstop));
set(gca,'XTicklabel',''); 
datetick('x','HH','keepticks');
legend('LW accumulation','cumulative sum of MASC images');
set(gca,'Fontsize',30);
% 
% subplot(4,1,4); hold on;
% plot(datenum(data.t),cumsum(data.Nmasc),'k-','linewidth',3);
% %addaxis(datenum(pluvio.t),pluvio.accuNRT,'r-','linewidth',3);
% %addaxis(datenum(data.t),cumsum(data.Nmasc),'linewidth',3);
% %addaxislabel(2,'# particles MASC');
% xlabel('time');
% ylabel('accuNRT total [mm]');
% set(gca,'Xlim',[datenum(tstart) datenum(tstop)]);
% set(gca,'XTick',datenum(tstart):datenum([0 0 0 1 0 0]):datenum(tstop));
% set(gca,'XTicklabel',''); 
% datetick('x','HH','keepticks');






%
% make MASC vs 2DVD plots for the whole campaign !
% clear all; close all;
% 
% % MASC parameters
% data = load('prediction/data/Xstruct_Davos_winter15-16_extended_last.mat');
% tstart_vec = [2016 02 29 00 00 00]; 
% tstop_vec  = [2016 02 29 23 59 59];
% tend_vec = [2016 03 01 00 00 00];
% xhi_thresh = 8.5;
% twin_vec = [0 0 0 0 5 0]; % time window size
% tshift_vec = [0 0 0 0 1 0]; % time window shift
% Npts_thresh = 10; % Nmin points in interval (MASC) 
% 
% % ATMVAR parameters
% atmdata_filepath = '/home/praz/Documents/MASC/masclab/atm_vars/order_45574_data.txt';
% 
% % 2DVD parameters
% tshift2dvd_vec = [0 0 0 0 10 0]; %time shift in 2DVD to match MASC time
% 
% % path for saving figures
% savepath = '/home/praz/Documents/MASC/MASC2DVD_comparison';
% 
% %%
% 
% while datenum(tstart_vec) ~= datenum(tend_vec)
%    
%     disp(datestr(datenum(tstart_vec),'yyyy-mm-dd'));
%     
%     %% load MASC
%     tstart_str = datestr(tstart_vec,'yyyymmddHHMMSS');
%     tstart = datenum(tstart_vec);
%     tstop_str  = datestr(tstop_vec,'yyyymmddHHMMSS');
%     tstop  = datenum(tstop_vec);
%     res = 35; % 1 pixel = ??? mu
% 
%     % if a .mat file loaded, filter flakes included in the time window
%     idx = find(data.Xt>=tstart & data.Xt<=tstop);
%     
%     if isempty(idx) 
%         tstart_vec = datevec(datenum(tstart_vec) + datenum([0 0 1 0 0 0]));
%         tstop_vec = datevec(datenum(tstop_vec) + datenum([0 0 1 0 0 0]));
%         continue;
%     end
%     
%     X = data.X(idx,:);
%     Xt = data.Xt(idx);
%     Xname = data.Xname(idx);
%     fprintf('%u snowflakes found in the desired time interval \n',length(idx));
% 
%     % keep good snowflakes only
%     idx = find(X(:,6)>=xhi_thresh);
%     X = X(idx,:);
%     Xt = Xt(idx);
%     Xname = Xname(idx);
%     fprintf('%u remaining after blurry snowflakes filtered out (xhi=%2.1f)\n',length(idx),xhi_thresh);
% 
%     %% load atmospheric variables
%     atm = load_atm_vars(atmdata_filepath,tstart_str,tstop_str,0);
% 
%     %% load 2DVD
%     dvd = load_2DVD_classif(tstart_str);
%     % time shift
%     tshift2dvd_num = datenum(tshift2dvd_vec);
%     dvd.t = dvd.t + tshift2dvd_num;
% 
%     %% dynamic time window
%     % variables of interest
%     Dmax = X(:,3)*res/1000;
%     AR = X(:,9);
%     complex = X(:,5);
%     snowclass = X(:,1);
%     theta = X(:,10);
%     theta(theta<0) = -theta(theta<0);
%     fs = X(:,7);
%     fs(fs>10) = 10;
% 
%     % time window
%     twin = datenum(twin_vec);
%     tshift = datenum(tshift_vec);
% 
%     tbot = tstart;
%     ttop = tbot + twin;
% 
%     % first bin
%     idx_in = find(Xt>=tbot & Xt<ttop);
%     N(1) = numel(idx_in);
%     t(1) = tbot + twin/2;
%     Dmax_t{1} = Dmax(idx_in);
%     AR_t{1} = AR(idx_in);
%     complex_t{1} = complex(idx_in);
%     theta_t{1} = theta(idx_in);
%     fs_t{1} = fs(idx_in);
%     for i=1:6
%         sclass(i,1) = sum(snowclass(idx_in)==i);
%     end
%     N_dry(1) = sum(X(idx_in,15)==0);
%     N_wet(1) = sum(X(idx_in,15)==1);
%     Avg_riming(1) = mean(X(idx_in,15));
% 
%     % while loop for the following bins
%     k = 1;
%     max_idx_in = length(idx_in);
%     while tbot < tstop
% 
%         k = k+1;
%         tbot = tbot + tshift;
%         ttop = tbot + twin;
% 
%         idx_in = find(Xt>=tbot & Xt<ttop);
%         N(k) = numel(idx_in);
%         t(k) = tbot + twin/2;
% 
%         Dmax_t{k} = Dmax(idx_in);
%         AR_t{k} = AR(idx_in);
%         complex_t{k} = complex(idx_in);
%         theta_t{k} = theta(idx_in);
%         fs_t{k} = fs(idx_in);
% 
%         %snowflake class
%         for i=1:6
%             sclass(i,k) = sum(snowclass(idx_in)==i);
%         end
% 
%         %melting snow
%         N_dry(k) = sum(X(idx_in,16)==0);% & X(idx_in,1)==4);
%         N_wet(k) = sum(X(idx_in,16)==1);% & X(idx_in,1)==4);
% 
%         %average degree of riming
%         Avg_riming(k) = mean(X(idx_in,15));
% 
%         if length(idx_in) > max_idx_in
%             max_idx_in = length(idx_in);
%         end
% 
%     end
% 
%     % tweak to use nice area plots
%     Dmax_mat = nan(max_idx_in,length(Dmax_t));
%     AR_mat = nan(max_idx_in,length(AR_t));
%     complex_mat = nan(max_idx_in,length(complex_t));
%     fs_mat = nan(max_idx_in,length(fs_t));
%     theta_mat = nan(max_idx_in,length(theta_t));
% 
%     for i=1:length(Dmax_t)
%         Dmax_vec = Dmax_t{i}';
%         if length(Dmax_vec) < max_idx_in
%             Dmax_vec(end+1:max_idx_in) = NaN;
%         end
%         Dmax_mat(:,i) = Dmax_vec;
% 
%         AR_vec = AR_t{i}';
%         if length(AR_vec) < max_idx_in
%             AR_vec(end+1:max_idx_in) = NaN;
%         end
%         AR_mat(:,i) = AR_vec;
% 
%         complex_vec = complex_t{i}';
%         if length(complex_vec) < max_idx_in
%             complex_vec(end+1:max_idx_in) = NaN;
%         end
%         complex_mat(:,i) = complex_vec;
% 
%         fs_vec = fs_t{i}';
%         if length(fs_vec) < max_idx_in
%             fs_vec(end+1:max_idx_in) = NaN;
%         end
%         fs_mat(:,i) = fs_vec;   
% 
%         theta_vec = theta_t{i}';
%         if length(theta_vec) < max_idx_in
%             theta_vec(end+1:max_idx_in) = NaN;
%         end
%         theta_mat(:,i) = theta_vec;
% 
%     end
% 
%     % plot every continuous section iteratively
%     %idx_notnan = find(nansum(Dmax_mat)>Npts_thresh);
%     %[~,idxs_stop] = find(diff(idx_notnan)~=1);
%     %idxs_start = [idx_notnan(1)  idx_notnan([idxs_stop+1])];
%     %idxs_stop = [idx_notnan([idxs_stop]) idx_notnan(end)];
% 
% 
%     % snowflake class : normalization
%     for i=1:length(N)
%         sclass(:,i) = sclass(:,i)/N(i);
%     end
%     % riming : normalization between 0 and 1
%     % Avg_riming = (Avg_riming-1)/4;
% 
%     % creating the suppa-duppa colormap
%     c = hsv(6);
%     cmasc = c;
%     cmasc(4,:) = c(6,:);
%     cmasc(6,:) = [49,163,84]./255;
%     cdvd(1,:) = c(6,:);
%     cdvd(2,:) = c(3,:);
%     cdvd(3,:) = c(5,:);
%     cdvd(4,:) = c(1,:);
%     cdvd(5,:) = c(2,:);
%     cdvd(6,:) = [241,105,19]./255;
%     cdvd(7,:) = [140,45,4]./255;
%     cdvd(8,:) = c(4,:);
% 
% 
%     %% illustrations
% 
%     % fig1 : comparison with 2DVD
%     myfig = figure('Visible','on','units','normalized','outerposition',[0 0 1 1]);
%     % subplot : number of particles
%     subplot(6,1,1); hold on;
%     title(sprintf('EVENT : %s to %s',datestr(tstart,'yyyy-mm-dd-HH:MM'),datestr(tstop,'yyyy-mm-dd-HH:MM')));
%     plot(t,N,'k-','linewidth',1.25);
%     addaxis(dvd.t,dvd.N,'linewidth',1.25);
%     addaxislabel(2,'# particles 2DVD');
%     % addaxis(atm.tnum,atm.precip);
%     % addaxislabel(3,'precip [mm/h]');
%     %xlabel('time');
%     ylabel('# images MASC');
%     set(gca,'Xlim',[tstart tstop]);
%     set(gca,'XTick',tstart:datenum([0 0 0 1 0 0]):tstop);
%     set(gca,'XTicklabel',''); 
%     box on;
%     % subplot: MASC classification
%     subplot(6,1,2:3); hold on;
%     h = area(t,sclass','EdgeColor','none'); 
%     for i=1:6
%         h(i).FaceColor = cmasc(i,:);
%     end
%     legend('SM','CO','DE','AG','GR','CC');
%     set(gca,'Xlim',[tstart tstop]);
%     set(gca,'XTick',tstart:datenum([0 0 0 1 0 0]):tstop);
%     set(gca,'XTicklabel',''); 
%     set(gca,'Ylim',[0 1]);
%     set(gca,'YTick',[]);
%     box on;
%     % subplot: 2DVD classification
%     XX = dvd.t;
%     YY = [0 1];
%     ZZ = [dvd.classif dvd.classif];
%     subplot(6,1,4); hold on;
%     plot3(t,Avg_riming,linspace(10,10,length(t)),'k-','linewidth',1.25);
%     %plot3(t,N_wet./(N_wet+N_dry),linspace(10,10,length(t)),'k--','linewidth',1.25);
%     surface(XX,YY,ZZ','EdgeColor','none');
%     unique_class = unique(ZZ);
%     unique_class = unique_class(~isnan(unique_class));
%     colormap(cdvd(unique_class,:));
%     set(gca,'Xlim',[tstart tstop]);
%     set(gca,'XTick',tstart:datenum([0 0 0 1 0 0]):tstop);
%     set(gca,'XTicklabel',''); 
%     %datetick('x','HH:MM','keepticks','keeplimits');
%     set(gca,'Ylim',[0 1]);
%     set(gca,'YTick',[]);
%     box on;
%     legend('Riming degree');
%     % subplot : melting snow
%     subplot(6,1,5); box on;
%     plot(t,N_wet./(N_wet+N_dry),'k-','color',[241,105,19]./255,'linewidth',1.25);
%     legend('% of wet snow');
%     set(gca,'Xlim',[tstart tstop]);
%     set(gca,'XTick',tstart:datenum([0 0 0 1 0 0]):tstop);
%     set(gca,'XTicklabel',''); 
%     %datetick('x','HH:MM','keepticks','keeplimits');
%     %ax = gca;
%     %ax.XTickLabelRotation=55;
%     set(gca,'Ylim',[0 1]);
%     % subplot : atm. variables
%     subplot(6,1,6); hold on;
%     plot([tstart tstop],[0,0],'k--');
%     plot(atm.tnum,atm.T,'linewidth',1.25);
%     axis([tstart tstop min(atm.T)-1 max(atm.T)+1]);
%     %addaxis(atm.tnum,atm.RH);
%     addaxis(atm.tnum,atm.winv,'linewidth',1.25);
%     addaxislabel(1,'T [C]');
%     %addaxislabel(2,'RH [%]');
%     addaxislabel(2,'wind [m/s]');
%     set(gca,'Xlim',[tstart tstop]);
%     set(gca,'XTick',tstart:datenum([0 0 0 1 0 0]):tstop);
%     set(gca,'XTicklabel',''); 
%     datetick('x','HH:MM','keepticks','keeplimits');
%     ax = gca;
%     ax.XTickLabelRotation=90;
%     set(gcf,'Visible','on');
%     box on;
% 
%     % save figure
%     %savename = sprintf('MASC2DVD_%s',datestr(tstart,'yyyy-mm-dd'));
%     %export_fig(fullfile(savepath,savename),'-png');
% 
%     %close all;
% 
%     if false
%     % fig2 : MASC descriptors I
%     figure(2);
%     subplot(311); hold on; grid on; box on; 
%     title(sprintf('EVENT : %s to %s',datestr(tstart,'yyyy-mm-dd-HH:MM'),datestr(tstop,'yyyy-mm-dd-HH:MM')));
%     plot(t,N,'k.-');
%     addaxis(atm.tnum,atm.precip);
%     addaxislabel(2,'precip [mm/h]');
%     %xlabel('time');
%     ylabel('number of images');
%     set(gca,'Xlim',[tstart tstop]);
%     set(gca,'XTick',tstart:datenum([0 0 0 1 0 0]):tstop);
%     set(gca,'XTicklabel',''); 
%     %datetick('x','HH:MM','keepticks','keeplimits');
% 
%     % ax = gca;
%     ax.XTickLabelRotation=55;
% 
%     subplot(312); hold on; box on;
%     for i=1:length(idxs_start)
%     plot_area_curve(t(idxs_start(i):idxs_stop(i))',Dmax_mat(:,idxs_start(i):idxs_stop(i))','','Dmax');
%     end
%     set(gca,'Xlim',[tstart tstop]);
%     set(gca,'XTick',tstart:datenum([0 0 0 1 0 0]):tstop);
%     set(gca,'XTicklabel',''); 
%     %datetick('x','HHPM','keepticks','keeplimits');
%     set(gca,'Ylim',[0 4]);
%     ax = gca;
%     ax.XTickLabelRotation=55;
% 
%     subplot(313); hold on; box on;
%     for i=1:length(idxs_start)
%     plot_area_curve(t(idxs_start(i):idxs_stop(i))',AR_mat(:,idxs_start(i):idxs_stop(i))','','Aspect ratio',[8,69,148]/255,[66,146,198]/255,[107,174,214]/255);
%     end
%     %plot_area_curve(t_nan',AR_mat','','Aspect Ratio',[8,69,148]/255,[66,146,198]/255,[107,174,214]/255);
%     set(gca,'Xlim',[tstart tstop]);
%     set(gca,'XTick',tstart:datenum([0 0 0 1 0 0]):tstop);
%     set(gca,'XTicklabel',''); 
%     datetick('x','HH:MM','keepticks','keeplimits');
%     set(gca,'Ylim',[0.3 1]);
%     ax = gca;
%     ax.XTickLabelRotation=55;
% 
%     % fig3: MASC descriptors II
%     figure(3);
%     subplot(311); hold on; box on;
%     title(sprintf('EVENT : %s to %s',datestr(tstart,'yyyy-mm-dd-HH:MM'),datestr(tstop,'yyyy-mm-dd-HH:MM')));
%     for i=1:length(idxs_start)
%     plot_area_curve(t(idxs_start(i):idxs_stop(i))',complex_mat(:,idxs_start(i):idxs_stop(i))','','Complexity',[35,139,69]/255,[102,194,164]/255,[178,226,226]/255);
%     end
%     %plot_area_curve(t_nan',complex_mat','','Complexity',[35,139,69]/255,[102,194,164]/255,[178,226,226]/255);
%     set(gca,'Xlim',[tstart tstop]);
%     set(gca,'XTick',tstart:datenum([0 0 0 1 0 0]):tstop);
%     set(gca,'XTicklabel',''); 
%     %datetick('x','HHPM','keepticks','keeplimits');
%     set(gca,'Ylim',[0.9 2]);
%     ax = gca;
%     ax.XTickLabelRotation=55;
% 
%     subplot(312); hold on; box on;
%     for i=1:length(idxs_start)
%     plot_area_curve(t(idxs_start(i):idxs_stop(i))',theta_mat(:,idxs_start(i):idxs_stop(i))','','Orientation',[136,86,167]/255,[158,188,218]/255,[224,236,244]/255);
%     end
%     %plot_area_curve(t_nan',theta_mat','','Orientation',[136,86,167]/255,[158,188,218]/255,[224,236,244]/255);
%     set(gca,'Xlim',[tstart tstop]);
%     set(gca,'XTick',tstart:datenum([0 0 0 1 0 0]):tstop);
%     set(gca,'XTicklabel',''); 
%     %datetick('x','HHPM','keepticks','keeplimits');
%     set(gca,'Ylim',[0 90]);
%     ax = gca;
%     ax.XTickLabelRotation=55;
% 
%     subplot(313); hold on; box on;
%     for i=1:length(idxs_start)
%     plot_area_curve(t(idxs_start(i):idxs_stop(i))',fs_mat(:,idxs_start(i):idxs_stop(i))','','Fallspeed',[99,99,99]/255,[189,189,189]/255,[240,240,240]/255);
%     end
%     %plot_area_curve(t_nan',fs_mat','time','Fallspeed',[99,99,99]/255,[189,189,189]/255,[240,240,240]/255);
%     set(gca,'Xlim',[tstart tstop]);
%     set(gca,'XTick',tstart:datenum([0 0 0 1 0 0]):tstop);
%     set(gca,'XTicklabel',''); 
%     datetick('x','HH:MM','keepticks','keeplimits');
%     set(gca,'Ylim',[0 2]);
%     ax = gca;
%     ax.XTickLabelRotation=55;
%     end
% 
%     tstart_vec = datevec(datenum(tstart_vec) + datenum([0 0 1 0 0 0]));
%     tstop_vec = datevec(datenum(tstop_vec) + datenum([0 0 1 0 0 0]));
%     
%     
% end
% 
% 
% 
% 
% 
% 




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



