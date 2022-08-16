% small script to analyse EastonAirport data
clear all; close all;

load('prediction/data/newFormat/classif_datastruct_CSU_triplet_withgraupfix_oldHCS.mat');
idx = find(data.X(:,5)<1);
data.X(idx,5) =1;

%% time series

% create a structure masc in the required shape for load_MASC_classif_2
masc.Xt = data.Yt;
masc.xhi = data.Y(:,6);
masc.label_ID = data.Y(:,1);
masc.riming_idx = data.Y(:,15);
masc.melting = data.Y(:,16);
% time serie
TMP = load_MASC_classif_2(masc,[2015 11 27 02 00 00],[2015 11 27 05 00 00], 0, 10, 8, 0, true); % window size, window shift
 
%% some histograms
close all;
c1 = [153 102 255]./255;
c2 = [255 153  51]./255;
fs = 14;

% riming histograms
t_thresh = datenum([2015 11 27 00 00 00]);
idx_before = find(data.Yt < t_thresh & data.Y(:,1)~=1);
idx_after = find(data.Yt > t_thresh & data.Y(:,1)~=1);
Ridx_before = data.Y(idx_before,15);
Ridx_after = data.Y(idx_after,15);
hspace = 0.1;
figure; 
hold on; grid on; box on;
histogram(Ridx_before,[0:hspace:1],'DisplayStyle','stairs','linewidth',2);
histogram(Ridx_after,[0:hspace:1],'DisplayStyle','stairs','linewidth',2);
xlabel('Riming index \in [0,1]');
ylabel('# images');
set(gca,'Fontsize',16);
legend('R_i 26 Nov','R_i 27 Nov');

% size (pixel) for cameras (2 3 4) histograms
idx_before = find(data.Xt < t_thresh & data.X(:,18)>1);
idx_after = find(data.Xt > t_thresh & data.X(:,18)>1);

figure; 
subplot(211); hold on; grid on; box on;
title('PSD [pix] for side cameras');
histogram(data.X(idx_before,3),[0:25:400],'FaceColor',c1);
legend('26 Nov');
subplot(212); hold on; grid on; box on;
histogram(data.X(idx_after,3),[0:25:400],'FaceColor',c2);
legend('27 Nov');

figure; hold on; grid on; box on;
[N_Dmax,X_Dmax] = hist(data.X(idx_before,3),[0:25:400]);
N_Dmax = N_Dmax./(sum(N_Dmax)*25);
plot(X_Dmax,N_Dmax,'Color',c1,'linewidth',2);
[N_Dmax,X_Dmax] = hist(data.X(idx_after,3),[0:25:400]);
N_Dmax = N_Dmax./(sum(N_Dmax)*25);
plot(X_Dmax,N_Dmax,'Color',c2,'linewidth',2);
xlabel('PSD [pix] for side cameras');
ylabel('Normalized density');
legend('26 Nov','27 Nov');
set(gca,'Fontsize',fs);



% size (pixel) for cameras (0 1) histograms
idx_before = find(data.Xt < t_thresh & data.X(:,18)<2);
idx_after = find(data.Xt > t_thresh & data.X(:,18)<2);

figure; 
subplot(211); hold on; grid on; box on;
title('PSD [pix] for downwards pointing cameras');
histogram(data.X(idx_before,3),[0:20:200],'FaceColor',c1);
set(gca,'Ylim',[0 80]);
legend('26 Nov');
subplot(212); hold on; grid on; box on;
histogram(data.X(idx_after,3),[0:20:200],'FaceColor',c2);
set(gca,'Ylim',[0 80]);
legend('27 Nov');

figure; hold on; grid on; box on;
[N_Dmax,X_Dmax] = hist(data.X(idx_before,3),[0:20:200]);
N_Dmax = N_Dmax./(sum(N_Dmax)*20);
plot(X_Dmax,N_Dmax,'Color',c1,'linewidth',2);
[N_Dmax,X_Dmax] = hist(data.X(idx_after,3),[0:20:200]);
N_Dmax = N_Dmax./(sum(N_Dmax)*20);
plot(X_Dmax,N_Dmax,'Color',c2,'linewidth',2);
xlabel('PSD [pix] for downwards pointing cameras');
ylabel('Normalized density');
legend('26 Nov','27 Nov');
set(gca,'Fontsize',fs);

% orientation on cameras (2 3 4) histograms
idx_before = find(data.Xt < t_thresh & data.X(:,18)>1);
idx_after = find(data.Xt > t_thresh & data.X(:,18)>1);

figure; 
subplot(211); hold on; grid on; box on;
title('Particle orientation [degrees from horizontal]');
histogram(abs(data.X(idx_before,10)),'FaceColor',c1);
%set(gca,'Ylim',[0 80]);
legend('26 Nov');
subplot(212); hold on; grid on; box on;
histogram(abs(data.X(idx_after,10)),'FaceColor',c2);
%set(gca,'Ylim',[0 80]);
legend('27 Nov');

figure; hold on; grid on; box on;
[N,X] = hist(abs(data.X(idx_before,10)),[0:10:90]);
N = N./(sum(N)*10);
plot(X,N,'Color',c1,'linewidth',2);
[N,X] = hist(abs(data.X(idx_after,10)),[0:10:90]);
N = N./(sum(N)*10);
plot(X,N,'Color',c2,'linewidth',2);
xlabel('Orientation [degrees from horizontal]');
ylabel('Normalized density');
legend('26 Nov','27 Nov');
set(gca,'Fontsize',fs,'Xlim',[0 90]);

% aspect ratio (all cameras)
idx_before = find(data.Xt < t_thresh & data.X(:,18)>-1);
idx_after = find(data.Xt > t_thresh & data.X(:,18)>-1);

figure; 
subplot(211); hold on; grid on; box on;
title('Particle Aspect Ratio');
histogram(data.X(idx_before,9),'FaceColor',c1);
%set(gca,'Ylim',[0 80]);
legend('26 Nov');
subplot(212); hold on; grid on; box on;
histogram(data.X(idx_after,9),'FaceColor',c2);
%set(gca,'Ylim',[0 80]);
legend('27 Nov');

figure; hold on; grid on; box on;
[N,X] = hist(data.X(idx_before,9),[0:0.1:1]);
N = N./(sum(N)*0.1);
plot(X,N,'Color',c1,'linewidth',2);
[N,X] = hist(data.X(idx_after,9),[0:0.1:1]);
N = N./(sum(N)*0.1);
plot(X,N,'Color',c2,'linewidth',2);
xlabel('Aspect ratio');
ylabel('Normalized density');
legend('26 Nov','27 Nov');
set(gca,'Fontsize',fs);

% complexity 
idx_before = find(data.Xt < t_thresh & data.X(:,18)>-1);
idx_after = find(data.Xt > t_thresh & data.X(:,18)>-1);

figure; 
subplot(211); hold on; grid on; box on;
title('Particle complexity');
histogram(data.X(idx_before,5),[1:0.25:5],'FaceColor',c1);
%set(gca,'Ylim',[0 80]);
legend('26 Nov');
subplot(212); hold on; grid on; box on;
histogram(data.X(idx_after,5),[1:0.25:5],'FaceColor',c2);
%set(gca,'Ylim',[0 80]);
legend('27 Nov');

figure; hold on; grid on; box on;
[N,X] = hist(data.X(idx_before,5),[1:0.25:5]);
N = N./(sum(N)*0.25);
plot(X,N,'Color',c1,'linewidth',2);
[N,X] = hist(data.X(idx_after,5),[1:0.25:5]);
N = N./(sum(N)*0.25);
plot(X,N,'Color',c2,'linewidth',2);
xlabel('Particles complexity');
ylabel('Normalized density');
legend('26 Nov','27 Nov');
set(gca,'Fontsize',fs);

% fallspeed
% idx_before = find(data.Yt < t_thresh);
% idx_after = find(data.Yt > t_thresh);
% 
% figure; hold on; grid on; box on;
% [N,X] = hist(data.Y(idx_before,7));
% N = N./(sum(N)*0.25);
% plot(X,N,'Color',c1,'linewidth',2);
% [N,X] = hist(data.Y(idx_after,7));
% N = N./(sum(N)*0.25);
% plot(X,N,'Color',c2,'linewidth',2);
% xlabel('Fallspeed');
% ylabel('Normalized density');
% legend('26 Nov','27 Nov');
% set(gca,'Fontsize',fs)




%% riming time series (different ways)

c1 = [153 102 255]./255;
c2 = [255 153  51]./255;

t_thresh = datenum([2015 11 27 00 00 00]);
idx_before = find(data.Yt < t_thresh & data.Y(:,1)~=1);
idx_after = find(data.Yt > t_thresh & data.Y(:,1)~=1);
xx1 = data.Yt(idx_before);
yy1 = data.Y(idx_before,15);
xx2 = data.Yt(idx_after);
yy2 = data.Y(idx_after,15);

figure; fs=12;

% lines with markers
subplot(121); hold on; grid on; box on;
plot(xx1,yy1,'bo','MarkerFaceColor','b','MarkerEdgeColor',c1,'MarkerFaceColor',c1,'Markersize',4);
set(gca,'Xlim',[datenum([2015 11 26 17 00 00]) datenum([2015 11 26 19 00 00])]);
set(gca,'XTick',datenum([2015 11 26 17 00 00]):datenum([0 0 0 0 30 0]):datenum([2015 11 26 19 00 00]));
%set(gca,'XTicklabel','');
set(gca,'Ylim',[0.1 1]);
datetick('x','HH:MM','keepticks','keeplimits');
ylabel('all pts')
set(gca,'Fontsize',fs);
title('26 Nov. 17h-19h');

subplot(122); hold on; grid on; box on;
plot(xx2,yy2,'ro','MarkerFaceColor','r','MarkerEdgeColor',c2,'MarkerFaceColor',c2,'Markersize',4);
set(gca,'Xlim',[datenum([2015 11 27 02 00 00]) datenum([2015 11 27 05 00 00])]);
set(gca,'XTick',datenum([2015 11 27 02 00 00]):datenum([0 0 0 0 30 0]):datenum([2015 11 27 05 00 00]));
datetick('x','HH:MM','keepticks','keeplimits');
set(gca,'Fontsize',fs);
set(gca,'Ylim',[0.1 1]);
title('27 Nov. 02h-05h');

% stems
% subplot(423); hold on; grid on; box on;
% stem(xx1,yy1,'filled','Color','b','Markersize',4);
% set(gca,'Xlim',[datenum([2015 11 26 17 00 00]) datenum([2015 11 26 19 00 00])]);
% set(gca,'XTick',datenum([2015 11 26 17 00 00]):datenum([0 0 0 0 30 0]):datenum([2015 11 26 19 00 00]));
% set(gca,'Ylim',[0.1 1]);
% datetick('x','HH:MM','keepticks','keeplimits');
% set(gca,'Fontsize',fs);
% 
% subplot(424); hold on; grid on; box on;
% stem(xx2,yy2,'filled','Color','r','Markersize',4);
% set(gca,'Xlim',[datenum([2015 11 27 02 00 00]) datenum([2015 11 27 05 00 00])]);
% set(gca,'XTick',datenum([2015 11 27 02 00 00]):datenum([0 0 0 0 30 0]):datenum([2015 11 27 05 00 00]));
% datetick('x','HH:MM','keepticks','keeplimits');
% set(gca,'Fontsize',fs);
% set(gca,'Ylim',[0.1 1]);

% moving window
Nmin_shift = 5;
Nmin_interval = 10;

tstart = datetime([2015 11 26 16 50 00]);
tstop = datetime([2015 11 26 19 00 00]);
tgrid = tstart:minutes(Nmin_shift):tstop; tgrid = tgrid';
tgrid2 = tgrid + minutes(Nmin_interval);
xt1 = datetime(xx1,'ConvertFrom','datenum');
for i=1:length(tgrid)
    idx_in = find(xt1>tgrid(i) & xt1<=tgrid2(i));
    yr1(i) = mean(yy1(idx_in));
end

figure;
subplot(121); hold on; grid on; box on;
plot(datenum(mean([tgrid tgrid2],2)),yr1,'b-','Color',c1,'linewidth',1.5);
plot(xx1,yy1,'bx','MarkerEdgeColor',c1,'MarkerFaceColor',c1,'markersize',3);
set(gca,'Xlim',[datenum([2015 11 26 17 00 00]) datenum([2015 11 26 19 00 00])]);
set(gca,'XTick',datenum([2015 11 26 17 00 00]):datenum([0 0 0 0 30 0]):datenum([2015 11 26 19 00 00]));
set(gca,'Ylim',[0.1 1]);
datetick('x','HH:MM','keepticks','keeplimits');
ylabel('Riming index');
title('26 Nov. 17h-19h');
set(gca,'Fontsize',fs);

tstart = datetime([2015 11 27 01 50 00]);
tstop = datetime([2015 11 27 05 00 00]);
tgrid = []; tgrid2 = [];
tgrid = tstart:minutes(Nmin_shift):tstop; tgrid = tgrid';
tgrid2 = tgrid + minutes(Nmin_interval);
xt2 = datetime(xx2,'ConvertFrom','datenum');
for i=1:length(tgrid)
    idx_in = find(xt2>tgrid(i) & xt2<=tgrid2(i));
    yr2(i) = mean(yy2(idx_in));
end

subplot(122); hold on; grid on; box on;
plot(datenum(mean([tgrid tgrid2],2)),yr2,'r-','Color',c2,'linewidth',1.5);
plot(xx2,yy2,'rx','MarkerEdgeColor',c2,'MarkerFaceColor',c2,'markersize',3);
set(gca,'Xlim',[datenum([2015 11 27 02 00 00]) datenum([2015 11 27 05 00 00])]);
set(gca,'XTick',datenum([2015 11 27 02 00 00]):datenum([0 0 0 0 30 0]):datenum([2015 11 27 05 00 00]));
set(gca,'Ylim',[0.1 1]);
datetick('x','HH:MM','keepticks','keeplimits');
title('27 Nov. 02h-05h');
set(gca,'Fontsize',fs);

% average over N points
n = 5;
for i=ceil(n/2):1:length(xx1)-floor(n/2)
    
    avg_x1(i) = xx1(i);
    avg_y1(i) = mean(yy1(i-floor(n/2):i+floor(n/2)));
    
end
figure;
subplot(121); hold on; grid on; box on;
plot(avg_x1,avg_y1,'bo','MarkerFaceColor',c1,'MarkerEdgeColor',c1,'Markersize',4);
set(gca,'Xlim',[datenum([2015 11 26 17 00 00]) datenum([2015 11 26 19 00 00])]);
set(gca,'XTick',datenum([2015 11 26 17 00 00]):datenum([0 0 0 0 30 0]):datenum([2015 11 26 19 00 00]));
datetick('x','HH:MM','keepticks','keeplimits');
ylabel('Riming degree');
title('26 Nov. 17h-19h');
set(gca,'Fontsize',fs);
set(gca,'Ylim',[0.1 1]);


for i=ceil(n/2):1:length(xx2)-floor(n/2)
    
    avg_x2(i) = xx2(i);
    avg_y2(i) = mean(yy2(i-floor(n/2):i+floor(n/2)));
    
end
subplot(122); hold on; grid on; box on;
plot(avg_x2,avg_y2,'ro','MarkerFaceColor',c2,'MarkerEdgeColor',c2,'Markersize',4);
set(gca,'Xlim',[datenum([2015 11 27 02 00 00]) datenum([2015 11 27 05 00 00])]);
set(gca,'XTick',datenum([2015 11 27 02 00 00]):datenum([0 0 0 0 30 0]):datenum([2015 11 27 05 00 00]));
datetick('x','HH:MM','keepticks','keeplimits');
set(gca,'Fontsize',fs);
title('27 Nov. 02h-05h');
set(gca,'Ylim',[0.1 1]);


if 0

% saving results in ASCII file
filename = 'EastonAirport_masc_classification_all_particles_newclassifier.dat';
fid = fopen(filename,'w');

fprintf(fid,'# particle-by-particle MASC classification output \n');
fprintf(fid,'# columns header: 1) particle timestamp [UTC], 2) particle ID, 3) particle degree of riming 4) particle type \n');
fprintf(fid,'# columns delimiter = tabulation \n');
fprintf(fid,'# particle types : 1 = small particle (SP), 2 = columnar crystal (CC), 3 = planar crystal (PC), 4 = aggregate (AG), 5 = graupel (GR), 6 = combination of columnar and planar crystals (CPC) \n');
for i=1:length(data.Yt)
    fprintf(fid,'%s \t %u \t %2.2f \t %u \n',datestr(data.Yt(i),'yyyy.mm.dd HH:MM:SS'),data.Y(i,28),data.Y(i,15),data.Y(i,1));
end
fclose(fid);

% results with aggregation in time
% filename = 'EastonAirport_masc_classification_aggregated.dat';
% fid = fopen(filename,'w');
% fprintf(fid,'# MASC classification output aggregated over 30min intervals \n');
% fprintf(fid,'# columns header: 1) beginning of the time interval [UTC], 2) end of the time interval [UTC], 3) number of particles classified within the time interval 4) averaged degree of riming \n');
% fprintf(fid,'5) %% of small particles (SP) 6) %% of columnar crystals (CC), 7) %% of planar crystals (PC), 8) %% of aggregates (AG), 9) %% of graupels (GR), 10) %% of combination of columnar and planar crystals (CPC) \n');
% fprintf(fid,'# columns delimiter = tabulation \n');
% for i=1:length(TMP.t)-1
%     fprintf(fid,'%s \t %s \t %u \t %2.2f \t %2.3f \t %2.3f \t %2.3f \t %2.3f \t %2.3f \t %2.3f \n',datestr(TMP.t(i)-minutes(30),'yyyy.mm.dd HH:MM:SS'),datestr(TMP.t(i),'yyyy.mm.dd HH:MM:SS'),TMP.Nmasc(i), TMP.riming(i), ...
%         TMP.sclass(i,1), TMP.sclass(i,2), TMP.sclass(i,3), TMP.sclass(i,4), TMP.sclass(i,5), TMP.sclass(i,6));
% end
%fclose(fid);

end






% t_start_vec = [2015 11 26 17 00 00];
% t_stop_vec = [2015 11 27 05 00 00];
% xhi_thresh = 9;
% Nmin_interval = 5;
% Nmin_shift = 5;
% Nimg_min = 0;
% illu = 1;
% %data = load_MASC_classif(datastruct,t_start_vec,t_stop_vec,xhi_thresh,Nmin_interval,Nmin_shift,Nimg_min,illu);
% 
%  
% load(datastruct);
% t_thresh = datenum([2015 11 27 00 00 00]);
% idx_before = find(Xt < t_thresh);
% idx_after = find(Xt > t_thresh);
% Ridx_before = X(idx_before,15);
% Ridx_after = X(idx_after,15);
% 
% 
% % for each flake ID, we find all views and average stuff over it
% triplet.flake_ID = unique(X(:,22));
% for i=1:length(triplet.flake_ID)
%     idx_views = find(X(:,22)==triplet.flake_ID(i));
%     triplet.avg_label_prob(i,:) = nansum(Xfullprob_label(idx_views,:));
%     triplet.avg_label_prob(i,:) = triplet.avg_label_prob(i,:)./sum(triplet.avg_label_prob(i,:));
%     [~,triplet.label_ID(i,1)] = nanmax(triplet.avg_label_prob(i,:));
%     triplet.avg_riming_index(i,1) = nanmean(X(idx_views,15));
%     if  (triplet.avg_riming_index(i,1) < 0.7) && (triplet.label_ID(i,1) == 5)
%         triplet.avg_riming_index(i,1) = 0.7 + 0.2*rand;
%     end
%     triplet.avg_melting_index(i,1) = nanmean(X(idx_views,17));
%     if triplet.avg_melting_index(i,1) > 0.5
%         triplet.is_melting(i,1) = true;
%     else
%         triplet.is_melting(i,1) = false;
%     end
%     triplet.t(i,1) = Xt(idx_views(1));
%     triplet.xhi(i,1) = nanmean(X(idx_views,6));
% end
%     
% idx_before = find(triplet.t < t_thresh);
% idx_after = find(triplet.t > t_thresh);
% triplet_Ridx_before =  triplet.avg_riming_index(idx_before);
% triplet_Ridx_after = triplet.avg_riming_index(idx_after);
% 
% hspace = 0.1;
% figure; 
% subplot(211); hold on; grid on; box on;
% histogram(Ridx_before,[0:hspace:1],'DisplayStyle','stairs','linewidth',2);
% histogram(Ridx_after,[0:hspace:1],'DisplayStyle','stairs','linewidth',2);
% xlabel('Riming index \in [0,1]');
% ylabel('# images');
% set(gca,'Fontsize',16);
% legend('R_i 26 Nov','R_i 27 Nov');
% subplot(212); hold on; grid on; box on;
% histogram(triplet_Ridx_before,[0:hspace:1],'DisplayStyle','stairs','linewidth',2);
% histogram(triplet_Ridx_after,[0:hspace:1],'DisplayStyle','stairs','linewidth',2);
% xlabel('Riming index \in [0,1]');
% ylabel('# images');
% set(gca,'Fontsize',16);
% legend('R_i 26 Nov','R_i 27 Nov');
% 
% 
% %% temporal analysis of the event
% 
% % create a structure masc in the required shape for load_MASC_classif_2
% masc.Xt = triplet.t;
% masc.xhi = triplet.xhi;
% masc.label_ID = triplet.label_ID;
% masc.riming_idx = triplet.avg_riming_index;
% masc.melting = triplet.is_melting;
% 
% TMP = load_MASC_classif_2(masc,[2015 11 26 17 00 00],[2015 11 27 05 00 00], 0, 30, 30, 0, true);
% 
% 
% if 0
% 
% % saving results in ASCII file
% filename = 'EastonAirport_masc_classification_all_particles.dat';
% fid = fopen(filename,'w');
% fprintf(fid,'# particle-by-particle MASC classification output \n');
% fprintf(fid,'# columns header: 1) particle timestamp [UTC], 2) particle ID, 3) particle degree of riming 4) particle type \n');
% fprintf(fid,'# columns delimiter = tabulation \n');
% fprintf(fid,'# particle types : 1 = small particle (SP), 2 = columnar crystal (CC), 3 = planar crystal (PC), 4 = aggregate (AG), 5 = graupel (GR), 6 = combination of columnar and planar crystals (CPC) \n');
% for i=1:length(triplet.t)
%     fprintf(fid,'%s \t %u \t %2.2f \t %u \n',datestr(triplet.t(i),'yyyy.mm.dd HH:MM:SS'),triplet.flake_ID(i),triplet.avg_riming_index(i),triplet.label_ID(i));
% end
% fclose(fid);
% 
% 
% % results with aggregation in time
% filename = 'EastonAirport_masc_classification_aggregated.dat';
% fid = fopen(filename,'w');
% fprintf(fid,'# MASC classification output aggregated over 30min intervals \n');
% fprintf(fid,'# columns header: 1) beginning of the time interval [UTC], 2) end of the time interval [UTC], 3) number of particles classified within the time interval 4) averaged degree of riming \n');
% fprintf(fid,'5) %% of small particles (SP) 6) %% of columnar crystals (CC), 7) %% of planar crystals (PC), 8) %% of aggregates (AG), 9) %% of graupels (GR), 10) %% of combination of columnar and planar crystals (CPC) \n');
% fprintf(fid,'# columns delimiter = tabulation \n');
% for i=1:length(TMP.t)-1
%     fprintf(fid,'%s \t %s \t %u \t %2.2f \t %2.3f \t %2.3f \t %2.3f \t %2.3f \t %2.3f \t %2.3f \n',datestr(TMP.t(i)-minutes(30),'yyyy.mm.dd HH:MM:SS'),datestr(TMP.t(i),'yyyy.mm.dd HH:MM:SS'),TMP.Nmasc(i), TMP.riming(i), ...
%         TMP.sclass(i,1), TMP.sclass(i,2), TMP.sclass(i,3), TMP.sclass(i,4), TMP.sclass(i,5), TMP.sclass(i,6));
% end
% fclose(fid);
% 
% end






