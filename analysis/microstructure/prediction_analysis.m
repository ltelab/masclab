% script to analyze snowflakes classification output over a period of time
clear all; close all;

% headers
% data_dir = '/media/praz/MyData/MASC/processed_27-Apr-2016/GOOD';
%data_dir = '/media/praz/MyData/MASC/processed_Davos_winter_2015-16/processed_13-Jun-2016/GOOD';
%data_dir = '/media/praz/MyData/MASC/processed_Davos_winter_2015-16/processed_single_full/GOOD';
%data_dir = '/media/praz/MASC_Vanderbilt/part2_processed/nocroptopbot/processed_2016-11-01_11:20/GOOD';
%data_dir = '/media/praz/MASC_Vanderbilt/data_cleaned_processed/single/GOOD';

tstart_vec = [2010 11 21 00 00 00]; 
tstart_str = datestr(tstart_vec,'yyyymmddHHMMSS');
tstart = datenum(tstart_vec);
tstop_vec  = [2020 08 05 22 00 00];
tstop_str  = datestr(tstop_vec,'yyyymmddHHMMSS');
tstop  = datenum(tstop_vec);
res = 33.5; % 1 pixel = ??? mu

to_load = false;

if to_load

    % load data
    X = [];
    Xname = {};
    Xt = [];
    Xfullprob_label = [];
    Xfullprob_riming = [];
    dir_list = uploaddirs(data_dir,tstart_vec,tstop_vec);
    perc_old = 0;
    for i=1:numel(dir_list)

        fprintf('%u / %u\n',i,numel(dir_list));
        [X_tmp,Xlab,Xname_tmp,Xt_tmp,Xfullprob_label_tmp,Xfullprob_riming_tmp] = load_processed_labels(dir_list{i},tstart_str,tstop_str,1);
        if ~isempty(X_tmp)
            X = [X; X_tmp];
            Xname = [Xname; Xname_tmp];
            Xt = [Xt; Xt_tmp];
            Xfullprob_label = [Xfullprob_label; Xfullprob_label_tmp];
            Xfullprob_riming = [Xfullprob_riming; Xfullprob_riming_tmp];
        end
        
        perc_current = floor(100*i/numel(dir_list));
        if perc_current > perc_old
            fprintf('%u %% done \n',perc_current);
            perc_old = perc_current;
        end

    end

else 
    
    load('APRES3_winter_2015-16.mat'); 
    
    % classif_datastruct_APRES3_triplet_withgraupfix.mat <- paper Cryosphere
    
end

% newFormat to oldFormat conversion
X = data.X;
Xt = data.Xt;
Xlab = data.Xlab;
Xname = data.Xname;

% select time interval
idx2keep = find(Xt>=tstart & Xt<=tstop);
X = X(idx2keep,:);
Xt = Xt(idx2keep);
Xname = Xname(idx2keep);

%% small analysis graupels
close all;

idx_gr = find(X(:,1) == 5 & X(:,28) >=0);
gr_xhi = X(idx_gr,6);
gr_riming_idx = X(idx_gr,15);
gr_area = X(idx_gr,2);
gr_dmax = X(idx_gr,3);
gr_dmean = X(idx_gr,4);
gr_xhi2 = X(idx_gr,28);
gr_xhi4 = X(idx_gr,29);

% figure;
% subplot(311);
% histogram(X(idx_gr,6));
% subplot(312);
% histogram(X(idx_gr,28));
% subplot(313);
% histogram(X(idx_gr,29));
% 
% figure;
% plot(X(:,2),X(:,28),'k.');



%figure;
%plot(gr_area,gr_riming_idx,'k.');
% 
% 
ncolor = 10000;
c = jet(ncolor);
d = floor(logspace(1,4,100));
d = ncolor+1 - flipud(d');

c2 = c(d,:);

figure;
n = histogram2(gr_xhi4,gr_riming_idx,'DisplayStyle','tile','ShowEmptyBins','on');
colormap(c2);
colorbar;
title('graupels density map');
xlabel('some parameter');
ylabel('riming index');
set(gca,'Fontsize',14);
%caxis(my_clim)
% 
% 
% % Now change the pseudocolor axis to a log scale.
% %caxis(log10(my_clim));
% %scatter(gr_xhi,gr_riming_idx):
figure;
histogram(gr_riming_idx);
xlabel('riming index');
ylabel('count');
set(gca,'Fontsize',14);
fprintf('Pourcentage of GR with RI below 0.5 : %2.2f%%\n',sum(gr_riming_idx<0.5)/length(gr_riming_idx)*100);


%% selecting only "real" falling snow events (APRES3)
% t2_vec = [2015 11 22 00 00 00];
% t2_str = datenum(t2_vec);
% idx_blow = find(Xt < t2_str);
% X(idx_blow,:) = [];
% Xname(idx_blow) = [];
% Xt(idx_blow) = [];

%select events
% t1_start = [2015 12 26 00 00 00];
% t1_stop = [2015 12 28 00 00 00];
% idx2keep = find(Xt >= datenum(t1_start) & Xt <= datenum(t1_stop));
% X = X(idx2keep,:);
% Xname = Xname(idx2keep);


%% riming study
% riming = X(:,12);
% close all;
% c = jet(5);
% labels = {'1','2','3','4','5'};
% figure(50); hold on; grid on; box on;
% for i=1:5
%     if i<5
%         idx = find(X(:,6) >= 9 & X(:,13) > 0 & X(:,12) == i & X(:,1)~=5);
%     elseif i==5
%         idx = find(X(:,6) >= 9 & X(:,11) > 0.9 & X(:,1)==5);%(X(:,12) == i | X(:,1)~=5));
%     end
%     Xfilt = X(idx,:);
%     [N_fs,X_fs] = hist(Xfilt(:,7),logspace(-1.2,1.2,30));
%     N_fs = N_fs./nansum(N_fs);
%     figure(50);
%     plot(X_fs,N_fs,'Color',c(i,:),'linewidth',2);
% end
% 
% figure(50);
% xlabel('[m/s]');
% title('Fallspeed');
% set(gca,'Xscale','log');
% legend(labels);
% set(gca,'Fontsize',14);
% set(gca,'Xlim',[0.1 10]);
% legend(labels);


%% pie plots as a function of xhi

if 1

xhi_lim = [0 9 9.25 9.5];
max_intens_lim = 0.1;
mean_intens_lim = 0.1;
size_min = 15;

for i=1:4
    idx_xhi = find(X(:,6) >= xhi_lim(i) & X(:,26) >= mean_intens_lim & X(:,27) >= max_intens_lim & X(:,4) >= size_min);
    X1 = X(idx_xhi,:);
    slices_labels = {'small','col','dend','agg','grau','capcol'}; %{'grau','agg','melt','small','dend','col'};
    spaces = {' ',' ',' ',' ',' ',' '};
    perc_symbols = {'%','%','%','%','%','%'};
    slices(1) = sum(X1(:,1)==1);
    slices(2) = sum(X1(:,1)==2);
    slices(3) = sum(X1(:,1)==3);
    slices(4) = sum(X1(:,1)==4);
    slices(5) = sum(X1(:,1)==5);
    slices(6) = sum(X1(:,1)==6);

    for j=1:6
        perc(j) = round(100*sum(slices(j))/length(X1(:,1)));
    end
    slices_labels = strcat(slices_labels,spaces,num2strs(perc),perc_symbols);
    figure(104);
    subplot(2,2,i);
    pie(slices,slices_labels);
    title(sprintf('Davos xhi > %2.1f',xhi_lim(i)));
    
    % plot particles location on screen distribution
%     figure;
%     subplot(2,3,1); hold on; box on;
%     idx = find(X(:,18)==0 & X(:,6)>= xhi_lim(i));
%     histogram(X(idx,19));
%     v = axis;
%     plot([320 320],[v(3) v(4)],'r--');
%     set(gca,'XLim',[0,2448]);
%     subplot(2,3,2); hold on; box on;
%     idx = find(X(:,18)==1 & X(:,6)>= xhi_lim(i));
%     histogram(X(idx,19));
%     title('distribution of snowflake along the x-axis [0-2448]');
%     set(gca,'XLim',[0,2448]);
%     subplot(2,3,3); hold on; box on;
%     idx = find(X(:,18)==2 & X(:,6)>= xhi_lim(i));
%     histogram(X(idx,19));
%     v = axis;
%     plot([2448-320 2448-320],[v(3) v(4)],'r--');
%     set(gca,'XLim',[0,2448]);
%     subplot(2,3,4); hold on; box on;
%     idx = find(X(:,18)==0 & X(:,6)>= xhi_lim(i));
%     histogram(X(idx,20));
%     ylabel('camera 0');
%     v = axis;
%     plot([400 400],[v(3) v(4)],'r--');
%     plot([2048-400 2048-400],[v(3) v(4)],'r--');
%     set(gca,'XLim',[0 2048]);
%     set(gca,'XDir','reverse');
%     view(90,-90);
%     subplot(2,3,5); hold on; box on;
%     idx = find(X(:,18)==1 & X(:,6)>= xhi_lim(i));
%     histogram(X(:,20));
%     ylabel('camera 1');
%     title('distribution of snowflake along the y-axis [0-2048]');
%     v = axis;
%     plot([400 400],[v(3) v(4)],'r--');
%     plot([2048-400 2048-400],[v(3) v(4)],'r--');
%     set(gca,'XLim',[0 2048]);
%     set(gca,'XDir','reverse');
%     view(90,-90);
%     subplot(2,3,6); hold on; box on;
%     idx = find(X(:,18)==2 & X(:,6)>= xhi_lim(i));
%     histogram(X(:,20));
%     ylabel('camera 2');
%     v = axis;
%     plot([400 400],[v(3) v(4)],'r--');
%     plot([2048-400 2048-400],[v(3) v(4)],'r--');
%     set(gca,'XLim',[0 2048]);
%     set(gca,'XDir','reverse');
%     view(90,-90);
    
end

end

%% analysis riming for Christophe Genton paper
%close all;

xhi_lim = 9.5;
max_intens_lim = 0.1;
mean_intens_lim = 0.1;
size_min = 0;

% no filter
Ri_0 = X(:,15);

% filter small particles
idx_1 = find(X(:,1) > 1);
Ri_1 = X(idx_1,15);
filt_1 = (numel(Ri_0)-numel(Ri_1))/numel(Ri_0) * 100;

% filter small particles and low quality particles
idx_2 = find(X(:,1) > 1 & X(:,6) >= xhi_lim & X(:,26) >= mean_intens_lim & X(:,27) >= max_intens_lim & X(:,4) >= size_min);
Ri_2 = X(idx_2,15);
filt_2 = (numel(Ri_0)-numel(Ri_2))/numel(Ri_0) * 100;

figure;
subplot(311); hold on; box on; grid on;
title('All particles');
histogram(Ri_0,[0:0.025:1],'Normalization','probability');
subplot(312); hold on; box on; grid on;
title('Small particles discarded');
histogram(Ri_1,[0:0.025:1],'Normalization','probability');
subplot(313); hold on; box on; grid on;
title('Small discarded, xi>9.25, mean brightness > 0.1');
histogram(Ri_2,[0:0.025:1],'Normalization','probability');
xlabel('Riming index');





%% analysis of distribution of riming per class
close all;

xhi_lim = 9.25;
max_intens_lim = 0;
mean_intens_lim = 0;
size_min = 15;

idx_xhi = find(X(:,6) >= xhi_lim & X(:,26) >= mean_intens_lim & X(:,27) >= max_intens_lim & X(:,4) >= size_min);
X1 = X(idx_xhi,:);

figure;
slices_labels = {'small','col','plan','agg','grau','col+plan'}; %{'grau','agg','melt','small','dend','col'};
spaces = {' ',' ',' ',' ',' ',' '};
perc_symbols = {'%','%','%','%','%','%'};
slices(1) = sum(X1(:,1)==1);
slices(2) = sum(X1(:,1)==2);
slices(3) = sum(X1(:,1)==3);
slices(4) = sum(X1(:,1)==4);
slices(5) = sum(X1(:,1)==5);
slices(6) = sum(X1(:,1)==6);

for j=1:6
    perc(j) = round(100*sum(slices(j))/length(X1(:,1)));
end
slices_labels = strcat(slices_labels,spaces,num2strs(perc),perc_symbols);
pie(slices,slices_labels);
title('Hydrometeor class proportions');
set(gca,'fontsize',12);
set(gca,'linewidth',2);

figure;
for j=1:6
    idx = find(X1(:,1)==j);
    riming_index = X1(idx,15);
    xlim([-0.05 1.05]);
    subplot(2,3,j); hold on; box on;
    title(slices_labels{j});
    histogram(riming_index,30,'Normalization','probability');
end

figure;
for j=1:6
    idx = find(X1(:,1)==j);
    riming_degree = X1(idx,12);
    subplot(2,3,j); hold on; box on;
    title(slices_labels{j});
    histogram(riming_degree,'Normalization','probability');
    xlim([-1 6]);
end


% plot particles location on screen distribution
discardmat = [400 550 200 250]; % [t,b,l,r]
figure;
subplot(2,3,1); hold on; box on;
idx = find(X1(:,18)==0 & X1(:,6)>= xhi_lim);
histogram(X1(idx,19));
v = axis;
plot([discardmat(3) discardmat(3)],[v(3) v(4)],'r--');
set(gca,'XLim',[0,2448]);
subplot(2,3,2); hold on; box on;
idx = find(X1(:,18)==1 & X1(:,6)>= xhi_lim);
histogram(X1(idx,19));
title('distribution of snowflake along the x-axis [0-2448]');
set(gca,'XLim',[0,2448]);
subplot(2,3,3); hold on; box on;
idx = find(X1(:,18)==2 & X1(:,6)>= xhi_lim);
histogram(X1(idx,19));
v = axis;
plot([2448-discardmat(4) 2448-discardmat(4)],[v(3) v(4)],'r--');
set(gca,'XLim',[0,2448]);
subplot(2,3,4); hold on; box on;
idx = find(X1(:,18)==0 & X1(:,6)>= xhi_lim);
histogram(X1(idx,20));
ylabel('camera 0');
v = axis;
plot([discardmat(1) discardmat(1)],[v(3) v(4)],'r--');
plot([2048-discardmat(2) 2048-discardmat(2)],[v(3) v(4)],'r--');
set(gca,'XLim',[0 2048]);
set(gca,'XDir','reverse');
view(90,-90);
subplot(2,3,5); hold on; box on;
idx = find(X1(:,18)==1 & X1(:,6)>= xhi_lim);
histogram(X1(idx,20));
ylabel('camera 1');
title('distribution of snowflake along the y-axis [0-2048]');
v = axis;
plot([discardmat(1) discardmat(1)],[v(3) v(4)],'r--');
plot([2048-discardmat(2) 2048-discardmat(2)],[v(3) v(4)],'r--');
set(gca,'XLim',[0 2048]);
set(gca,'XDir','reverse');
view(90,-90);
subplot(2,3,6); hold on; box on;
idx = find(X1(:,18)==2 & X1(:,6)>= xhi_lim);
histogram(X1(idx,20));
ylabel('camera 2');
v = axis;
plot([discardmat(1) discardmat(1)],[v(3) v(4)],'r--');
plot([2048-discardmat(2) 2048-discardmat(2)],[v(3) v(4)],'r--');
set(gca,'XLim',[0 2048]);
set(gca,'XDir','reverse');
view(90,-90);



%% global distributions : Dmax, AR, canting angle, FS
%close all;
c = hsv(5);
% old scheme
% labels = {'grau','agg','melt','small','dend','col'};
% new scheme
labels = {'small','col','plan','agg','grau','col+plan'};

% pdf
figure(199); hold on; grid on; box on;
figure(1100); hold on; grid on; box on;
figure(1101); hold on; grid on; box on;
figure(1102); hold on; grid on; box on;
for i=1:5
    idx = find(X(:,6) >= 9 & X(:,11) > 0.0 & X(:,1) == i);
    Xfilt = X(idx,:);
    Xfilt(:,3) = Xfilt(:,3)*res/1000;
    Xfilt(:,4) = Xfilt(:,4)*res/1000;
    Xfilt(:,8) = Xfilt(:,8)*res/1000;
    [N_Dmax,X_Dmax] = hist(Xfilt(:,3),[0:0.2:19.75]);
    [N_AR,X_AR] = hist(Xfilt(:,9),[0.025:0.05:0.975]); % 0.05
    [N_angle,X_angle] = hist(abs(Xfilt(:,10)),[2.5:5:87.5]); % [2.5:5:87.5] if abs(..) [-87.5:5:87.5] otherwise
    [N_fs,X_fs] = hist(Xfilt(:,7),logspace(-1.2,1.2,30));
    N_Dmax = N_Dmax./(sum(N_Dmax)*0.2);
    N_AR = N_AR./(sum(N_AR)*0.05);
    N_angle = N_angle./(sum(N_angle)*5);
    N_fs = N_fs./(nansum(N_fs));

    figure(199);
    plot(X_Dmax,N_Dmax,'Color',c(i,:),'linewidth',2);
    figure(1100);
    plot(X_AR,N_AR,'Color',c(i,:),'linewidth',2);
    figure(1101);
    plot(X_angle,N_angle,'Color',c(i,:),'linewidth',2);
    figure(1102);
    plot(X_fs,N_fs,'Color',c(i,:),'linewidth',2);
end
figure(199);
xlabel('Dmax [mm]');
set(gca,'Xlim',[0 12]);
title('PSDs');
legend(labels);
set(gca,'Fontsize',14);
figure(1100);
title('Aspect Ratio');
legend(labels);
set(gca,'Fontsize',14);
figure(1101);
xlabel('angle from horizontal [degree]');
title('Canting Angle');
legend(labels);
set(gca,'Fontsize',14);
figure(1102);
xlabel('[m/s]');
title('Fallspeed');
set(gca,'Xscale','log');
legend(labels);
set(gca,'Fontsize',14);
set(gca,'Xlim',[0.1 10]);


% % axis ratio
% figure(100); hold on; grid on; box on;
% for i=1:6
% find(X(:,6) >= 9 & X(:,11) > 0.7 & X(:,1) == i);
% Xfilt = X(idx,:);





%% distributions as a function of the diameter

fill_between_lines = @(X,Y1,Y2,C) fill( [X fliplr(X)],  [Y1 fliplr(Y2)], C, 'edgecolor','none','facealpha',.2);

close all;
c = hsv(5);
%c(3,:) = [];
labels = {'small','col','dend','agg','grau','capcol'};

n_tresh = 10; % Davos 90; DDU 40
Dmax_range = 0:0.1:20;
Dmax_vec = Dmax_range+0.1;%double(unique(diff(Dmax_range)))/2;
Dmax_vec = Dmax_vec(1:end-1);

figure(103); hold on; grid on; box on;
figure(104); hold on; grid on; box on;
figure(105); hold on; grid on; box on;

yaplot = [];
yoplot = [];

for i=1:5

    if i~=10 % discard dendrites for APRES3
    
    % Davos none, DDU xi9 and prob0.9
    idx = find(X(:,6) >= 0 & X(:,11)>0.0 & X(:,1) == i);
    Xfilt = X(idx,:);
    Xfilt(:,3) = Xfilt(:,3)*0.035;

    for j=1:numel(Dmax_range)-1
        lim_inf = Dmax_range(j);
        lim_sup = Dmax_range(j+1);
        idx = find(Xfilt(:,3) >= lim_inf & Xfilt(:,3) < lim_sup);
        if ~isempty(idx)
            tmp.mean_ar(j) = nanmean(Xfilt(idx,9));
            tmp.med_ar(j) = nanmedian(Xfilt(idx,9));
            tmp.std_ar(j) = nanstd(Xfilt(idx,9));
            tmp.q25_ar(j) = quantile(Xfilt(idx,9),0.25);
            tmp.q75_ar(j) = quantile(Xfilt(idx,9),0.75);
            
            tmp.mean_or(j) = nanmean(abs(Xfilt(idx,10)));
            tmp.med_or(j) = nanmedian(abs(Xfilt(idx,10)));
            
            tmp.mean_fs(j) = nanmean(Xfilt(idx,7));
            tmp.med_fs(j) = nanmedian(Xfilt(idx,7));
            tmp.q25_fs(j) = quantile(Xfilt(idx,7),0.25);
            tmp.q75_fs(j) = quantile(Xfilt(idx,7),0.75);
            
            tmp.n_part(j) = length(idx);
        else
            tmp.mean_ar(j) = NaN;
            tmp.med_ar(j) = NaN;
            
            tmp.mean_or(j) = NaN;
            tmp.med_or(j) = NaN;
            
            tmp.mean_fs(j) = NaN;
            tmp.med_fs(j) = NaN;
            
            tmp.n_part(j) = 0;
        end
    end
                  
    figure(103);
    fill_between_lines(Dmax_vec(tmp.n_part>n_tresh),tmp.q25_ar(tmp.n_part>n_tresh),tmp.q75_ar(tmp.n_part>n_tresh),c(i,:));
    yaplot(end+1) = plot(Dmax_vec(tmp.n_part>n_tresh),tmp.med_ar(tmp.n_part>n_tresh),'Color',c(i,:),'linewidth',2);
    figure(104);
    plot(Dmax_vec(tmp.n_part>n_tresh),tmp.med_or(tmp.n_part>n_tresh),'Color',c(i,:),'linewidth',2);
    figure(105);
    fill_between_lines(Dmax_vec(tmp.n_part>n_tresh),tmp.q25_fs(tmp.n_part>n_tresh),tmp.q75_fs(tmp.n_part>n_tresh),c(i,:));
    yoplot(end+1) = plot(Dmax_vec(tmp.n_part>n_tresh),tmp.med_fs(tmp.n_part>n_tresh),'Color',c(i,:),'linewidth',2);
    
    end
    
end

figure(103);
xlabel('Dmax [mm]');
ylabel('Aspect ratio');
axis([0 12 0.1 1]);
legend(yaplot,labels);
%legend(labels);
set(gca,'Fontsize',14);
figure(104);
xlabel('Dmax [mm]');
ylabel('Orientation [deg]');
legend(labels);
set(gca,'Fontsize',14);
figure(105);
xlabel('Dmax [mm]');
ylabel('Fallspeed [m/s]');
set(gca,'Fontsize',14);
legend(yoplot,labels);



% errorbar(Dmax_vec,med_ar,med_ar-q25_ar,q75_ar-med_ar,'k.-');
% title('Aspect Ratio');
% subplot(3,1,2); box on; hold on;
% errorbar(Dmax_vec,mean_orientation,std_orientation,'k.-');
% title('Canting angle');
% subplot(3,1,3); box on; hold on;
% plot(Dmax_vec,n_part,'k-.');
% title('# particles');
% set(gca,'yscale','log');
% 

%% Structure for Daniel including D_eq, canting angle, AR for graupels and aggregates
idx = find(X(:,6) >= 9.5 & X(:,11)>0.9 & X(:,1) == 2);
Xfilt = X(idx,:);
vec_id = 1:1:size(Xfilt,1);
vec_id = vec_id';
output_mat = [vec_id Xfilt(:,3)*0.035 Xfilt(:,9) Xfilt(:,10)];
idx_weird = find(output_mat(:,2) < 0.6);
Xname_weird = Xname(idx);
Xname_weird = Xname_weird(idx_weird);










