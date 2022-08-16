% script to analyze snowflake microstructural properties as a function of
% the hydrometeor type and degree of riming
clearvars -except fs_data; close all;

% user params
% classif_datastruct_Verbier_triplet.mat
% classif_datastruct_Davos2015-16_triplet_withgraupfix_more_entries.mat
% classif_datastruct_APRES3_triplet_withgraupfix.mat

campaign = 'APRES3_2016-17';

switch campaign
    case 'Verbier2016'
        datafile = '/home/praz/Documents/MASC/masc_classification_datfiles/old/Verbier_winter_2016-17.mat';
    case 'Davos_SPICE'
        datafile = '/home/praz/Documents/MASC/masc_classification_datfiles/old/classif_datastruct_Davos2015-16_triplet_smartMerge.mat';
    case 'ICEPOP2018'
        datafile = '/ltedata/ICEPOP_2018/MASC/Proc_data/ICEPOP2018_MASC_classification_full.mat'; % moved from BKC to MH on 2018.02.21
    case 'APRES3_2016-17'
        datafile = '/ltedata/APRES3/2016-17/MASC/Proc_data/APRES3_2016-17_MASC_classification_full.mat'; % summer = 2017.01.05 -> 2017.02.21 // winter = 2017.03.23 -> 2017.07.18
    otherwise
        fprintf('Error: campaign %s not found.\n',campaign);
end

tstart_vec = [2017 02 21 00 00 00]; 
tstop_vec  = [2019 02 21 00 00 00];
res = 33.5; % 1 pixel = ??? mu

% loading
tstart_str = datestr(tstart_vec,'yyyymmddHHMMSS');
tstart = datenum(tstart_vec);
tstop_str  = datestr(tstop_vec,'yyyymmddHHMMSS');
tstop  = datenum(tstop_vec);
load(datafile); 

X = data.X;
Xt = data.Xt;
Xlab = data.Xlab;
Xname = data.Xname;

% select time interval
idx2keep = find(Xt>=tstart & Xt<=tstop);
X = X(idx2keep,:);
Xt = Xt(idx2keep);
Xname = Xname(idx2keep);

%% filtering
close all;

% (re)loading
ht = X(:,1);
fs = X(:,7);
Dmax = X(:,3) * 33.5/1000;
ypos = X(:,20);
xhi = X(:,6);
rc = X(:,12);
ri = X(:,15);
melt = X(:,16);
melt(isnan(melt)) = 0;
n_roi = X(:,28);

% what to display
fig.disp_filt_roi_effect = false;
fig.disp_fs_ht = false;
fig.disp_fs_ri = false;
fig.disp_overview = true;

% FS filters
switch campaign
    case 'Verbier2016'
        fs_max = 2;
        fs_max_fit = 2;
        fs_min = 0.5;
        ymin = 600;
        ymax = 800;
        xhi_thresh = 8;
        Dmax_min = 0.5;
        Dmax_max = 5;
        n_roi_max = 3;
    case 'Davos_SPICE'
        fs_max = 2;
        fs_max_fit = 2;
        fs_min = 0.5;
        ymin = 1000;
        ymax = 1200;
        xhi_thresh = 8;
        Dmax_min = 0.5;
        Dmax_max = 5;
        n_roi_max = Inf;
    case 'ICEPOP2018'
        fs_max = 2;
        fs_max_fit = 2;
        fs_min = 0.5;
        ymin = 1000;
        ymax = 1200;
        xhi_thresh = 8;
        Dmax_min = 0.5;
        Dmax_max = 5;
        n_roi_max = 10;
    case 'APRES3_2016-17'
        fs_max = 2;
        fs_max_fit = 2;
        fs_min = 0.5;
        ymin = 1000;
        ymax = 1200;
        xhi_thresh = 8;
        Dmax_min = 0.5;
        Dmax_max = 5;
        n_roi_max = 10;
    otherwise
        fprintf('Error: campaign %s not found.\n',campaign);
end


% data filtering
idx2keep = find(fs < fs_max & fs > fs_min & ~isnan(fs) & ~isnan(Dmax) & ~isnan(ri) & ypos >= ymin & ypos <= ymax & xhi >= xhi_thresh & Dmax >= Dmax_min & Dmax <= Dmax_max & n_roi <= n_roi_max);
Dmax = Dmax(idx2keep);
ht = ht(idx2keep);
fs = fs(idx2keep);
ypos = ypos(idx2keep);
xhi = xhi(idx2keep);
rc = rc(idx2keep);
ri = ri(idx2keep);
melt = melt(idx2keep);
n_roi = n_roi(idx2keep);

% keep fs_data in a structure (for campaign intercomparison purposes)
switch campaign
    case 'Verbier2016'
        fs_data.verbier.fs = fs;
    case 'Davos_SPICE'
        fs_data.davos.fs = fs;
    case 'ICEPOP2018'
        fs_data.icepop.fs = fs;
    otherwise
        fprintf('Error: campaign %s not found.\n',campaign);
end

%% illustrations

if fig.disp_filt_roi_effect
    
    fig_axis = [Dmax_min Dmax_max fs_min fs_max];
    % raw plot
    figure; hold on; grid on; box on;
    dscatter(Dmax,fs);
    xlabel('Dmax [mm]');
    ylabel('Fallspeed [m/s]');
    set(gca,'Fontsize',12);
    title('Davos - raw data');
    axis(fig_axis);

    figure;
    subplot(221);
    dscatter(Dmax,fs);
    xlabel('Dmax [mm]');
    ylabel('Fallspeed [m/s]');
    set(gca,'Fontsize',12);
    title('Verbier - no filter');
    axis(fig_axis);

    subplot(222);
    dscatter(Dmax(n_roi <= 5),fs(n_roi <= 5));
    xlabel('Dmax [mm]');
    ylabel('Fallspeed [m/s]');
    set(gca,'Fontsize',12);
    title('Verbier - N_{roi} <= 5');
    axis(fig_axis);

    subplot(223);
    dscatter(Dmax(n_roi <= 2),fs(n_roi <= 2));
    xlabel('Dmax [mm]');
    ylabel('Fallspeed [m/s]');
    set(gca,'Fontsize',12);
    title('Verbier - N_{roi} <= 2');
    axis(fig_axis);

    subplot(224);
    dscatter(Dmax(n_roi <= 1),fs(n_roi <= 1));
    xlabel('Dmax [mm]');
    ylabel('Fallspeed [m/s]');
    set(gca,'Fontsize',12);
    title('Verbier - N_{roi} <= 1');
    axis(fig_axis);
    
end

%fallspeed as a function of ht
if fig.disp_fs_ht

    idx_c2 = find(ht==2);
    idx_c4 = find(ht==4);
    idx_c5 = find(ht==5 & xhi>8.5);
    idx_c236 = find(ht==3 | ht==6);

    bin_size = 0.1;
    bins_center = Dmax_min:bin_size:Dmax_max;
    fs_mean = [];
    fs_med = [];
    fs_25 = [];

    xx = linspace(Dmax_min,Dmax_max,200);
    yy = 0.8.*xx.^0.16; % LH dry agg
    yy2 = 1.1.*xx.^0.57; % LH graupel

    fig_axis = [Dmax_min Dmax_max fs_min fs_max];
    figure('Visible','on','units','normalized','outerposition',[0.25 .25 0.35 1]); hold on; grid on; box on;
    subplot(4,1,1); hold on; grid on; box on;
    dscatter(Dmax(idx_c2),fs(idx_c2));
    plot(xx,yy,'r--');
    plot(xx,yy2,'r--');
    title('Columns');
    ylabel('Fallspeed [m/s]');
    axis(fig_axis);
    for i=1:length(bins_center)
        tmp_Dmax = Dmax(idx_c2);
        tmp_fs = fs(idx_c2);
        tmp_idx = find(tmp_Dmax>=bins_center(i)-bin_size/2 & tmp_Dmax<bins_center(i)+bin_size/2);
        fs_mean(i) = mean(tmp_fs(tmp_idx));
        fs_med(i) = median(tmp_fs(tmp_idx));   
        fs_25(i) = quantile(tmp_fs(tmp_idx),0.25);
        fs_75(i) = quantile(tmp_fs(tmp_idx),0.75);
    end
    plot(bins_center,fs_med,'ko','MarkerFaceColor','k');
    plot(bins_center,fs_25,'k--','linewidth',2);
    plot(bins_center,fs_75,'k--','linewidth',2);
    b = robustfit(log(bins_center),log(fs_med),'bisquare');
    fit_param_ht(1,1) = exp(b(1));
    fit_param_ht(1,2) = b(2);
    plot(xx,yy,'r--');
    plot(xx,yy2,'r--');
    % bidule pour plotter les traitillés ou ya peu de données
    qlow = 0.01;
    qhigh = 0.99;
    qlow_idx_c2 = quantile(Dmax(idx_c2),qlow);
    xx_c2_low = find(xx>=qlow_idx_c2,1,'first');
    qhigh_idx_c2 = quantile(Dmax(idx_c2),qhigh);
    xx_c2_high = find(xx>qhigh_idx_c2,1,'first');
    plot(xx(1:xx_c2_low),exp(b(1)).*xx(1:xx_c2_low).^b(2),'g--','linewidth',2);
    plot(xx(xx_c2_low:xx_c2_high),exp(b(1)).*xx(xx_c2_low:xx_c2_high).^b(2),'g-','linewidth',2);
    plot(xx(xx_c2_high:end),exp(b(1)).*xx(xx_c2_high:end).^b(2),'g--','linewidth',2);
    set(gca,'Fontsize',12);

    subplot(4,1,2); hold on; grid on; box on;
    dscatter(Dmax(idx_c236),fs(idx_c236));
    title('Planar Crystals');
    ylabel('Fallspeed [m/s]');
    axis(fig_axis);
    for i=1:length(bins_center)
        tmp_Dmax = Dmax(idx_c236);
        tmp_fs = fs(idx_c236);
        tmp_idx = find(tmp_Dmax>=bins_center(i)-bin_size/2 & tmp_Dmax<bins_center(i)+bin_size/2);
        fs_mean(i) = mean(tmp_fs(tmp_idx));
        fs_med(i) = median(tmp_fs(tmp_idx));   
        fs_25(i) = quantile(tmp_fs(tmp_idx),0.25);
        fs_75(i) = quantile(tmp_fs(tmp_idx),0.75);
    end
    plot(bins_center,fs_med,'ko','MarkerFaceColor','k');
    plot(bins_center,fs_25,'k--','linewidth',2);
    plot(bins_center,fs_75,'k--','linewidth',2);
    b = robustfit(log(bins_center),log(fs_med),'bisquare');
    fit_param_ht(2,1) = exp(b(1));
    fit_param_ht(2,2) = b(2);
    plot(xx,yy,'r--');
    plot(xx,yy2,'r--');
    qlow_idx_c236 = quantile(Dmax(idx_c236),qlow);
    xx_c236_low = find(xx>=qlow_idx_c236,1,'first');
    qhigh_idx_c236 = quantile(Dmax(idx_c236),qhigh);
    xx_c236_high = find(xx>qhigh_idx_c236,1,'first');
    plot(xx,exp(b(1)).*xx.^b(2),'g-','linewidth',2);
    set(gca,'Fontsize',12);



    subplot(4,1,3); hold on; grid on; box on;
    dscatter(Dmax(idx_c4),fs(idx_c4));
    title('Aggregates');
    ylabel('Fallspeed [m/s]');
    axis(fig_axis);
    for i=1:length(bins_center)
        tmp_Dmax = Dmax(idx_c4);
        tmp_fs = fs(idx_c4);
        tmp_idx = find(tmp_Dmax>=bins_center(i)-bin_size/2 & tmp_Dmax<bins_center(i)+bin_size/2);
        fs_mean(i) = mean(tmp_fs(tmp_idx));
        fs_med(i) = median(tmp_fs(tmp_idx));   
        fs_25(i) = quantile(tmp_fs(tmp_idx),0.25);
        fs_75(i) = quantile(tmp_fs(tmp_idx),0.75);
    end
    plot(bins_center,fs_med,'ko','MarkerFaceColor','k');
    plot(bins_center,fs_25,'k--','linewidth',2);
    plot(bins_center,fs_75,'k--','linewidth',2);
    b = robustfit(log(bins_center),log(fs_med),'bisquare');
    fit_param_ht(3,1) = exp(b(1));
    fit_param_ht(3,2) = b(2);
    plot(xx,yy,'r--');
    plot(xx,yy2,'r--');
    qlow_idx_c4 = quantile(Dmax(idx_c4),qlow);
    xx_c4_low = find(xx>=qlow_idx_c4,1,'first');
    qhigh_idx_c4 = quantile(Dmax(idx_c4),qhigh);
    xx_c4_high = find(xx>qhigh_idx_c4,1,'first');
    plot(xx,exp(b(1)).*xx.^b(2),'g-','linewidth',2);
    set(gca,'Fontsize',12);



    subplot(4,1,4); hold on; grid on; box on;
    dscatter(Dmax(idx_c5),fs(idx_c5));
    xlabel('Dmax [mm]');
    title('Graupel');
    ylabel('Fallspeed [m/s]');
    axis(fig_axis);
    for i=1:length(bins_center)
        tmp_Dmax = Dmax(idx_c5);
        tmp_fs = fs(idx_c5);
        tmp_idx = find(tmp_Dmax>=bins_center(i)-bin_size/2 & tmp_Dmax<bins_center(i)+bin_size/2);
        fs_mean(i) = mean(tmp_fs(tmp_idx));
        fs_med(i) = median(tmp_fs(tmp_idx));   
        fs_25(i) = quantile(tmp_fs(tmp_idx),0.25);
        fs_75(i) = quantile(tmp_fs(tmp_idx),0.75);
    end
    plot(bins_center,fs_med,'ko','MarkerFaceColor','k');
    plot(bins_center,fs_25,'k--','linewidth',2);
    plot(bins_center,fs_75,'k--','linewidth',2);
    b = robustfit(log(bins_center),log(fs_med),'bisquare');
    fit_param_ht(4,1) = exp(b(1));
    fit_param_ht(4,2) = b(2);
    plot(xx,yy,'r--');
    plot(xx,yy2,'r--');
    qlow_idx_c5 = quantile(Dmax(idx_c5),qlow);
    xx_c5_low = find(xx>=qlow_idx_c5,1,'first');
    qhigh_idx_c5 = quantile(Dmax(idx_c5),qhigh);
    xx_c5_high = find(xx>qhigh_idx_c5,1,'first');
    plot(xx,exp(b(1)).*xx.^b(2),'g-','linewidth',2);
    set(gca,'Fontsize',12);

    disp('Hydrometeor type : a and b params in v = a*Dmax^b');
    disp(fit_param_ht);


    % Show fs-Dmax fits for the hydrometeor type
    c = hsv(5);
    %labels = {'small','col','plan','agg','grau','col+plan'};

    figure; hold on; grid on; box on;
    bibi = [];
    bibi(end+1) = plot(xx,yy,'r--','linewidth',2);
    bibi(end+1) = plot(xx,yy2,'r--','linewidth',2);
    plot(xx(1:xx_c2_low),fit_param_ht(1,1).*xx(1:xx_c2_low).^fit_param_ht(1,2),'k--','color',c(2,:),'linewidth',2);
    bibi(end+1) = plot(xx(xx_c2_low:xx_c2_high),fit_param_ht(1,1).*xx(xx_c2_low:xx_c2_high).^fit_param_ht(1,2),'k-','color',c(2,:),'linewidth',2);
    plot(xx(xx_c2_high:end),fit_param_ht(1,1).*xx(xx_c2_high:end).^fit_param_ht(1,2),'k--','color',c(2,:),'linewidth',2);
    plot(xx(1:xx_c236_low),fit_param_ht(2,1).*xx(1:xx_c236_low).^fit_param_ht(2,2),'k--','color',c(3,:),'linewidth',2);
    bibi(end+1) = plot(xx(xx_c236_low:xx_c236_high),fit_param_ht(2,1).*xx(xx_c236_low:xx_c236_high).^fit_param_ht(2,2),'k-','color',c(3,:),'linewidth',2);
    plot(xx(xx_c236_high:end),fit_param_ht(2,1).*xx(xx_c236_high:end).^fit_param_ht(2,2),'k--','color',c(3,:),'linewidth',2);
    plot(xx(1:xx_c4_low),fit_param_ht(3,1).*xx(1:xx_c4_low).^fit_param_ht(3,2),'k--','color',c(4,:),'linewidth',2);
    bibi(end+1) = plot(xx(xx_c4_low:xx_c4_high),fit_param_ht(3,1).*xx(xx_c4_low:xx_c4_high).^fit_param_ht(3,2),'k-','color',c(4,:),'linewidth',2);
    plot(xx(xx_c4_high:end),fit_param_ht(3,1).*xx(xx_c4_high:end).^fit_param_ht(3,2),'k--','color',c(4,:),'linewidth',2);
    plot(xx(1:xx_c5_low),fit_param_ht(4,1).*xx(1:xx_c5_low).^fit_param_ht(4,2),'k--','color',c(5,:),'linewidth',2);
    bibi(end+1) = plot(xx(xx_c5_low:xx_c5_high),fit_param_ht(4,1).*xx(xx_c5_low:xx_c5_high).^fit_param_ht(4,2),'k-','color',c(5,:),'linewidth',2);
    plot(xx(xx_c5_high:end),fit_param_ht(4,1).*xx(xx_c5_high:end).^fit_param_ht(4,2),'k--','color',c(5,:),'linewidth',2);
    axis(fig_axis);
    bobo = {'LH dry snow','LH graupel','Columns','Planar crystals','Aggregates','Graupel'};
    legend(bibi,bobo);
    xlabel('Dmax [mm]');
    ylabel('Fallspeed [m/s]');
    set(gca,'Fontsize',12);

end

% fallspeed as a function of ri
if fig.disp_fs_ri
 
    ri_div = [0.25 0.5 0.75];
    idx_r1 = find(ri <= ri_div(1) & ht>1);
    idx_r2 = find(ri> ri_div(1) & ri <= ri_div(2) & ht>1);
    idx_r3 = find(ri > ri_div(2) & ri <= ri_div(3) & ht>1);
    idx_r4 = find(ri > ri_div(3) & ht>1 & xhi > 9);

    %Dmax_min = 0.5;
    %Dmax_max = 5;
    xx = linspace(Dmax_min,Dmax_max,50);
    yy = 0.8.*xx.^0.16;
    yy2 = 1.1.*xx.^0.57;
    fig_axis = [Dmax_min Dmax_max fs_min fs_max];

    %modelfun = @(b,x)b(1).*x.^b(2);


    % median fallspeed as a function of RI
    bin_size = 0.1;
    bins_center = Dmax_min:bin_size:Dmax_max;
    fs_mean = [];
    fs_med = [];
    fs_25 = [];



    % RI1
    figure('Visible','on','units','normalized','outerposition',[0.25 .25 0.25 0.35]); hold on; grid on; box on;
    dscatter(Dmax(idx_r1),fs(idx_r1));
    plot(xx,yy,'r--','linewidth',2);
    plot(xx,yy2,'r--','linewidth',2);
    b = robustfit(log(Dmax(intersect(idx_r1,find(fs<fs_max_fit)))),log(fs(intersect(idx_r1,find(fs<fs_max_fit)))),'bisquare');
    %plot(xx,exp(b(1)).*xx.^b(2),'g-','linewidth',2);
    axis(fig_axis);

    for i=1:length(bins_center)
        tmp_Dmax = Dmax(idx_r1);
        tmp_fs = fs(idx_r1);
        tmp_idx = find(tmp_Dmax>=bins_center(i)-bin_size/2 & tmp_Dmax<bins_center(i)+bin_size/2);
        fs_mean(i) = mean(tmp_fs(tmp_idx));
        fs_med(i) = median(tmp_fs(tmp_idx));   
        fs_25(i) = quantile(tmp_fs(tmp_idx),0.25);
        fs_75(i) = quantile(tmp_fs(tmp_idx),0.75);
    end
    plot(bins_center,fs_med,'ko','MarkerFaceColor','k');
    plot(bins_center,fs_25,'k--','linewidth',2);
    plot(bins_center,fs_75,'k--','linewidth',2);
    b = robustfit(log(bins_center),log(fs_med),'bisquare');
    fit_param(1,1) = exp(b(1));
    fit_param(1,2) = b(2);
    plot(xx,exp(b(1)).*xx.^b(2),'g-','linewidth',2);
    set(gca,'Fontsize',12);
    xlabel('Dmax [mm]');
    ylabel('Fallspeed [m/s]');
    title('R_i < 0.25');




    % RI2
    figure('Visible','on','units','normalized','outerposition',[0.25 .25 0.25 0.35]); hold on; grid on; box on;
    dscatter(Dmax(idx_r2),fs(idx_r2),'SMOOTHING',20);
    plot(xx,yy,'r--','linewidth',2);
    plot(xx,yy2,'r--','linewidth',2);
    b = robustfit(log(Dmax(intersect(idx_r2,find(fs<fs_max_fit)))),log(fs(intersect(idx_r2,find(fs<fs_max_fit)))),'bisquare');
    %plot(xx,exp(b(1)).*xx.^b(2),'g-','linewidth',2);
    axis(fig_axis);

    for i=1:length(bins_center)
        tmp_Dmax = Dmax(idx_r2);
        tmp_fs = fs(idx_r2);
        tmp_idx = find(tmp_Dmax>=bins_center(i)-bin_size/2 & tmp_Dmax<bins_center(i)+bin_size/2);
        fs_mean(i) = mean(tmp_fs(tmp_idx));
        fs_med(i) = median(tmp_fs(tmp_idx)); 
        fs_25(i) = quantile(tmp_fs(tmp_idx),0.25);
        fs_75(i) = quantile(tmp_fs(tmp_idx),0.75);
    end
    plot(bins_center,fs_med,'ko','MarkerFaceColor','k');
    plot(bins_center,fs_25,'k--','linewidth',2);
    plot(bins_center,fs_75,'k--','linewidth',2);
    b = robustfit(log(bins_center),log(fs_med),'bisquare');
    fit_param(2,1) = exp(b(1));
    fit_param(2,2) = b(2);
    plot(xx,exp(b(1)).*xx.^b(2),'g-','linewidth',2);
    set(gca,'Fontsize',12);
    xlabel('Dmax [mm]');
    ylabel('Fallspeed [m/s]');
    title('R_i \in [0.25; 0.5]');

    % RI3
    figure('Visible','on','units','normalized','outerposition',[0.25 .25 0.25 0.35]); hold on; grid on; box on;
    dscatter(Dmax(idx_r3),fs(idx_r3),'SMOOTHING',20);
    plot(xx,yy,'r--','linewidth',2);
    plot(xx,yy2,'r--','linewidth',2);
    b = robustfit(log(Dmax(intersect(idx_r3,find(fs<fs_max_fit)))),log(fs(intersect(idx_r3,find(fs<fs_max_fit)))),'bisquare');
    %plot(xx,exp(b(1)).*xx.^b(2),'k-','linewidth',2);
    b = robustfit(log(fs(intersect(idx_r3,find(fs<fs_max_fit)))),log(Dmax(intersect(idx_r3,find(fs<fs_max_fit)))),'bisquare');
    A = exp(-b(1)/b(2));
    B = 1/b(2);
    %plot(xx,A.*xx.^B,'k--','linewidth',2);
    axis(fig_axis);

    for i=1:length(bins_center)
        tmp_Dmax = Dmax(idx_r3);
        tmp_fs = fs(idx_r3);
        tmp_idx = find(tmp_Dmax>=bins_center(i)-bin_size/2 & tmp_Dmax<bins_center(i)+bin_size/2);
        fs_mean(i) = mean(tmp_fs(tmp_idx));
        fs_med(i) = median(tmp_fs(tmp_idx)); 
        fs_25(i) = quantile(tmp_fs(tmp_idx),0.25);
        fs_75(i) = quantile(tmp_fs(tmp_idx),0.75);
    end
    plot(bins_center,fs_med,'ko','MarkerFaceColor','k');
    plot(bins_center,fs_25,'k--','linewidth',2);
    plot(bins_center,fs_75,'k--','linewidth',2);
    b = robustfit(log(bins_center),log(fs_med),'bisquare');
    fit_param(3,1) = exp(b(1));
    fit_param(3,2) = b(2);
    plot(xx,exp(b(1)).*xx.^b(2),'g-','linewidth',2);
    set(gca,'Fontsize',12);
    xlabel('Dmax [mm]');
    ylabel('Fallspeed [m/s]');
    title('R_i \in [0.5; 0.75]');

    % R4
    figure('Visible','on','units','normalized','outerposition',[0.25 .25 0.25 0.35]); hold on; grid on; box on;
    dscatter(Dmax(idx_r4),fs(idx_r4),'SMOOTHING',20);
    plot(xx,yy,'r--','linewidth',2);
    plot(xx,yy2,'r--','linewidth',2);
    b = robustfit(log(Dmax(intersect(idx_r4,find(fs<fs_max_fit)))),log(fs(intersect(idx_r4,find(fs<fs_max_fit)))),'bisquare');
    %plot(xx,exp(b(1)).*xx.^b(2),'k-','linewidth',2);
    b = robustfit(log(fs(intersect(idx_r4,find(fs<fs_max_fit)))),log(Dmax(intersect(idx_r4,find(fs<fs_max_fit)))),'bisquare');
    A = exp(-b(1)/b(2));
    B = 1/b(2);
    %plot(xx,A.*xx.^B,'k--','linewidth',2);
    axis(fig_axis);

    for i=1:length(bins_center)
        tmp_Dmax = Dmax(idx_r4);
        tmp_fs = fs(idx_r4);
        tmp_idx = find(tmp_Dmax>=bins_center(i)-bin_size/2 & tmp_Dmax<bins_center(i)+bin_size/2);
        fs_mean(i) = mean(tmp_fs(tmp_idx));
        fs_med(i) = median(tmp_fs(tmp_idx));  
        fs_25(i) = quantile(tmp_fs(tmp_idx),0.25);
        fs_75(i) = quantile(tmp_fs(tmp_idx),0.75);
    end
    plot(bins_center,fs_med,'ko','MarkerFaceColor','k');
    plot(bins_center,fs_25,'k--','linewidth',2);
    plot(bins_center,fs_75,'k--','linewidth',2);
    b = robustfit(log(bins_center),log(fs_med),'bisquare');
    fit_param(4,1) = exp(b(1));
    fit_param(4,2) = b(2);
    plot(xx,exp(b(1)).*xx.^b(2),'g-','linewidth',2);
    set(gca,'Fontsize',12);
    xlabel('Dmax [mm]');
    ylabel('Fallspeed [m/s]');
    title('R_i > 0.75');

    disp('Rming degree : a and b params in v = a*Dmax^b');
    disp(fit_param);


    % RI4 but Dmax(fs) fitted
    if 0
    xx_fs = linspace(fs_min,fs_max,50);
    figure('Visible','on','units','normalized','outerposition',[0.25 .25 0.25 0.35]); hold on; grid on; box on;
    dscatter(Dmax(idx_r4),fs(idx_r4),'SMOOTHING',20);
    plot(xx,yy,'r-','linewidth',2);

    fs_bin_size = 0.05;
    fs_bins_center = fs_min:fs_bin_size:fs_max;
    Dmax_mean = [];
    Dmax_med = [];
    Dmax_25 = [];

    for i=1:length(fs_bins_center)
        tmp_Dmax = Dmax(idx_r4);
        tmp_fs = fs(idx_r4);
        tmp_idx = find(tmp_fs>=fs_bins_center(i)-fs_bin_size/2 & tmp_fs<fs_bins_center(i)+fs_bin_size/2);
        Dmax_mean(i) = mean(tmp_Dmax(tmp_idx));
        Dmax_med(i) = median(tmp_Dmax(tmp_idx));  
        Dmax_25(i) = quantile(tmp_Dmax(tmp_idx),0.25);
        Dmax_75(i) = quantile(tmp_Dmax(tmp_idx),0.75);
    end
    plot(Dmax_med,fs_bins_center,'ko','MarkerFaceColor','k');
    plot(Dmax_25,fs_bins_center,'k--','linewidth',2);
    plot(Dmax_75,fs_bins_center,'k--','linewidth',2);
    b = robustfit(log(fs_bins_center),log(Dmax_med),'bisquare');
    plot(exp(b(1)).*xx_fs.^b(2),xx_fs,'g-','linewidth',2)

    axis(fig_axis);

    end

    % show the different Dmax - FS relationships for the riming degree
    c1 = [189,201,225]./255;
    c2 = [116,169,207]./255;
    c3 = [43,140,190]./255;
    c4 = [4,90,141]./255;


    figure; hold on; grid on; box on;
    plot(xx,yy,'r--','linewidth',2);
    plot(xx,yy2,'r--','linewidth',2);
    plot(xx,fit_param(1,1).*xx.^fit_param(1,2),'k-','color',c1,'linewidth',2);
    plot(xx,fit_param(2,1).*xx.^fit_param(2,2),'k-','color',c2,'linewidth',2);
    plot(xx,fit_param(3,1).*xx.^fit_param(3,2),'k-','color',c3,'linewidth',2);
    plot(xx,fit_param(4,1).*xx.^fit_param(4,2),'k-','color',c4,'linewidth',2);
    axis(fig_axis);
    legend('LH dry snow','LH graupel','R_i < 0.25','R_i \in [0.25,0.5]','R_i \in [0.5,0.75]','R_i > 0.75');
    xlabel('Dmax [mm]');
    ylabel('Fallspeed [m/s]');
    set(gca,'Fontsize',12);

end


% general plots
if fig.disp_overview
    
    xhi_lim = [8.5 9 9.5 10];
    for i=1:4
        idx_xhi = find(X(:,6) >= xhi_lim(i));
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

        %plot particles location on screen distribution
        figure;
        subplot(2,3,1); hold on; box on;
        idx = find(X(:,18)==0 & X(:,6)>= xhi_lim(i));
        histogram(X(idx,19));
        v = axis;
        plot([320 320],[v(3) v(4)],'r--');
        set(gca,'XLim',[0,2448]);
        subplot(2,3,2); hold on; box on;
        idx = find(X(:,18)==1 & X(:,6)>= xhi_lim(i));
        histogram(X(idx,19));
        title('distribution of snowflake along the x-axis [0-2448]');
        set(gca,'XLim',[0,2448]);
        subplot(2,3,3); hold on; box on;
        idx = find(X(:,18)==2 & X(:,6)>= xhi_lim(i));
        histogram(X(idx,19));
        v = axis;
        plot([2448-320 2448-320],[v(3) v(4)],'r--');
        set(gca,'XLim',[0,2448]);
        subplot(2,3,4); hold on; box on;
        idx = find(X(:,18)==0 & X(:,6)>= xhi_lim(i));
        histogram(X(idx,20));
        ylabel('camera 0');
        v = axis;
        plot([400 400],[v(3) v(4)],'r--');
        plot([2048-400 2048-400],[v(3) v(4)],'r--');
        set(gca,'XLim',[0 2048]);
        set(gca,'XDir','reverse');
        view(90,-90);
        subplot(2,3,5); hold on; box on;
        idx = find(X(:,18)==1 & X(:,6)>= xhi_lim(i));
        histogram(X(:,20));
        ylabel('camera 1');
        title('distribution of snowflake along the y-axis [0-2048]');
        v = axis;
        plot([400 400],[v(3) v(4)],'r--');
        plot([2048-400 2048-400],[v(3) v(4)],'r--');
        set(gca,'XLim',[0 2048]);
        set(gca,'XDir','reverse');
        view(90,-90);
        subplot(2,3,6); hold on; box on;
        idx = find(X(:,18)==2 & X(:,6)>= xhi_lim(i));
        histogram(X(:,20));
        ylabel('camera 2');
        v = axis;
        plot([400 400],[v(3) v(4)],'r--');
        plot([2048-400 2048-400],[v(3) v(4)],'r--');
        set(gca,'XLim',[0 2048]);
        set(gca,'XDir','reverse');
        view(90,-90);

    end

end

%% small fs intercomparison (MH vs SPICE, both DFIR but different altitude)

if 0 
fprintf('*********************************************************\n');
fprintf('Davos snowflake fallspeed    : mean = %2.2f, median = %2.2f, mode = %2.2f \n',mean(fs_data.davos.fs),median(fs_data.davos.fs),mode(fs_data.davos.fs));
fprintf('MayHills snowflake fallspeed : mean = %2.2f, median = %2.2f, mode = %2.2f \n',mean(fs_data.icepop.fs),median(fs_data.icepop.fs),mode(fs_data.icepop.fs));
%fprintf('Verbier snowflake fallspeed  : mean = %2.2f, median = %2.2f, mode = %2.2f \n',mean(fs_data.verbier.fs),median(fs_data.verbier.fs),mode(fs_data.verbier.fs));
fprintf('*********************************************************\n');

xbins = 0.05:0.1:5;
nbins_davos = hist(fs_data.davos.fs,xbins);
nbins_mh = hist(fs_data.icepop.fs,xbins);

%normalization
nbins_davos = nbins_davos./sum(nbins_davos);
nbins_mh = nbins_mh./sum(nbins_mh);

figure; hold on; grid on; box on;
plot(xbins,nbins_davos,'b-','linewidth',2);
plot(xbins,nbins_mh,'r-','linewidth',2);
xlabel('Fallspeed [m/s]');
ylabel('Normalized occurence');
set(gca,'xlim',[0.5 3]);
legend('Davos DFIR (2520m asl)','MayHills DFIR (780m asl');
set(gca,'Fontsize',14);

end







%% analysis of distribution of riming per class

xhi_lim = 9.25;
idx_xhi = find(X(:,6) >= xhi_lim);
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
% 
% 
% % plot particles location on screen distribution
% discardmat = [400 550 200 250]; % [t,b,l,r]
% figure;
% subplot(2,3,1); hold on; box on;
% idx = find(X(:,18)==0 & X(:,6)>= xhi_lim);
% histogram(X(idx,19));
% v = axis;
% plot([discardmat(3) discardmat(3)],[v(3) v(4)],'r--');
% set(gca,'XLim',[0,2448]);
% subplot(2,3,2); hold on; box on;
% idx = find(X(:,18)==1 & X(:,6)>= xhi_lim);
% histogram(X(idx,19));
% title('distribution of snowflake along the x-axis [0-2448]');
% set(gca,'XLim',[0,2448]);
% subplot(2,3,3); hold on; box on;
% idx = find(X(:,18)==2 & X(:,6)>= xhi_lim);
% histogram(X(idx,19));
% v = axis;
% plot([2448-discardmat(4) 2448-discardmat(4)],[v(3) v(4)],'r--');
% set(gca,'XLim',[0,2448]);
% subplot(2,3,4); hold on; box on;
% idx = find(X(:,18)==0 & X(:,6)>= xhi_lim);
% histogram(X(idx,20));
% ylabel('camera 0');
% v = axis;
% plot([discardmat(1) discardmat(1)],[v(3) v(4)],'r--');
% plot([2048-discardmat(2) 2048-discardmat(2)],[v(3) v(4)],'r--');
% set(gca,'XLim',[0 2048]);
% set(gca,'XDir','reverse');
% view(90,-90);
% subplot(2,3,5); hold on; box on;
% idx = find(X(:,18)==1 & X(:,6)>= xhi_lim);
% histogram(X(idx,20));
% ylabel('camera 1');
% title('distribution of snowflake along the y-axis [0-2048]');
% v = axis;
% plot([discardmat(1) discardmat(1)],[v(3) v(4)],'r--');
% plot([2048-discardmat(2) 2048-discardmat(2)],[v(3) v(4)],'r--');
% set(gca,'XLim',[0 2048]);
% set(gca,'XDir','reverse');
% view(90,-90);
% subplot(2,3,6); hold on; box on;
% idx = find(X(:,18)==2 & X(:,6)>= xhi_lim);
% histogram(X(idx,20));
% ylabel('camera 2');
% v = axis;
% plot([discardmat(1) discardmat(1)],[v(3) v(4)],'r--');
% plot([2048-discardmat(2) 2048-discardmat(2)],[v(3) v(4)],'r--');
% set(gca,'XLim',[0 2048]);
% set(gca,'XDir','reverse');
% view(90,-90);
% 
% 
% 
% %% global distributions : Dmax, AR, canting angle, FS
% %close all;
% c = hsv(5);
% % old scheme
% % labels = {'grau','agg','melt','small','dend','col'};
% % new scheme
% labels = {'small','col','plan','agg','grau','col+plan'};
% 
% % pdf
% figure(199); hold on; grid on; box on;
% figure(1100); hold on; grid on; box on;
% figure(1101); hold on; grid on; box on;
% figure(1102); hold on; grid on; box on;
% for i=1:5
%     idx = find(X(:,6) >= 9 & X(:,11) > 0.0 & X(:,1) == i);
%     Xfilt = X(idx,:);
%     Xfilt(:,3) = Xfilt(:,3)*res/1000;
%     Xfilt(:,4) = Xfilt(:,4)*res/1000;
%     Xfilt(:,8) = Xfilt(:,8)*res/1000;
%     [N_Dmax,X_Dmax] = hist(Xfilt(:,3),[0:0.2:19.75]);
%     [N_AR,X_AR] = hist(Xfilt(:,9),[0.025:0.05:0.975]); % 0.05
%     [N_angle,X_angle] = hist(abs(Xfilt(:,10)),[2.5:5:87.5]); % [2.5:5:87.5] if abs(..) [-87.5:5:87.5] otherwise
%     [N_fs,X_fs] = hist(Xfilt(:,7),logspace(-1.2,1.2,30));
%     N_Dmax = N_Dmax./(sum(N_Dmax)*0.2);
%     N_AR = N_AR./(sum(N_AR)*0.05);
%     N_angle = N_angle./(sum(N_angle)*5);
%     N_fs = N_fs./(nansum(N_fs));
% 
%     figure(199);
%     plot(X_Dmax,N_Dmax,'Color',c(i,:),'linewidth',2);
%     figure(1100);
%     plot(X_AR,N_AR,'Color',c(i,:),'linewidth',2);
%     figure(1101);
%     plot(X_angle,N_angle,'Color',c(i,:),'linewidth',2);
%     figure(1102);
%     plot(X_fs,N_fs,'Color',c(i,:),'linewidth',2);
% end
% figure(199);
% xlabel('Dmax [mm]');
% set(gca,'Xlim',[0 12]);
% title('PSDs');
% legend(labels);
% set(gca,'Fontsize',14);
% figure(1100);
% title('Aspect Ratio');
% legend(labels);
% set(gca,'Fontsize',14);
% figure(1101);
% xlabel('angle from horizontal [degree]');
% title('Canting Angle');
% legend(labels);
% set(gca,'Fontsize',14);
% figure(1102);
% xlabel('[m/s]');
% title('Fallspeed');
% set(gca,'Xscale','log');
% legend(labels);
% set(gca,'Fontsize',14);
% set(gca,'Xlim',[0.1 10]);
% 
% 
% % % axis ratio
% % figure(100); hold on; grid on; box on;
% % for i=1:6
% % find(X(:,6) >= 9 & X(:,11) > 0.7 & X(:,1) == i);
% % Xfilt = X(idx,:);
% 
% 
% 
% 
% 
% %% distributions as a function of the diameter
% 
% fill_between_lines = @(X,Y1,Y2,C) fill( [X fliplr(X)],  [Y1 fliplr(Y2)], C, 'edgecolor','none','facealpha',.2);
% 
% close all;
% c = hsv(5);
% %c(3,:) = [];
% labels = {'small','col','dend','agg','grau','capcol'};
% 
% n_tresh = 10; % Davos 90; DDU 40
% Dmax_range = 0:0.1:20;
% Dmax_vec = Dmax_range+0.1;%double(unique(diff(Dmax_range)))/2;
% Dmax_vec = Dmax_vec(1:end-1);
% 
% figure(103); hold on; grid on; box on;
% figure(104); hold on; grid on; box on;
% figure(105); hold on; grid on; box on;
% 
% yaplot = [];
% yoplot = [];
% 
% for i=1:5
% 
%     if i~=10 % discard dendrites for APRES3
%     
%     % Davos none, DDU xi9 and prob0.9
%     idx = find(X(:,6) >= 0 & X(:,11)>0.0 & X(:,1) == i);
%     Xfilt = X(idx,:);
%     Xfilt(:,3) = Xfilt(:,3)*0.035;
% 
%     for j=1:numel(Dmax_range)-1
%         lim_inf = Dmax_range(j);
%         lim_sup = Dmax_range(j+1);
%         idx = find(Xfilt(:,3) >= lim_inf & Xfilt(:,3) < lim_sup);
%         if ~isempty(idx)
%             tmp.mean_ar(j) = nanmean(Xfilt(idx,9));
%             tmp.med_ar(j) = nanmedian(Xfilt(idx,9));
%             tmp.std_ar(j) = nanstd(Xfilt(idx,9));
%             tmp.q25_ar(j) = quantile(Xfilt(idx,9),0.25);
%             tmp.q75_ar(j) = quantile(Xfilt(idx,9),0.75);
%             
%             tmp.mean_or(j) = nanmean(abs(Xfilt(idx,10)));
%             tmp.med_or(j) = nanmedian(abs(Xfilt(idx,10)));
%             
%             tmp.mean_fs(j) = nanmean(Xfilt(idx,7));
%             tmp.med_fs(j) = nanmedian(Xfilt(idx,7));
%             tmp.q25_fs(j) = quantile(Xfilt(idx,7),0.25);
%             tmp.q75_fs(j) = quantile(Xfilt(idx,7),0.75);
%             
%             tmp.n_part(j) = length(idx);
%         else
%             tmp.mean_ar(j) = NaN;
%             tmp.med_ar(j) = NaN;
%             
%             tmp.mean_or(j) = NaN;
%             tmp.med_or(j) = NaN;
%             
%             tmp.mean_fs(j) = NaN;
%             tmp.med_fs(j) = NaN;
%             
%             tmp.n_part(j) = 0;
%         end
%     end
%                   
%     figure(103);
%     fill_between_lines(Dmax_vec(tmp.n_part>n_tresh),tmp.q25_ar(tmp.n_part>n_tresh),tmp.q75_ar(tmp.n_part>n_tresh),c(i,:));
%     yaplot(end+1) = plot(Dmax_vec(tmp.n_part>n_tresh),tmp.med_ar(tmp.n_part>n_tresh),'Color',c(i,:),'linewidth',2);
%     figure(104);
%     plot(Dmax_vec(tmp.n_part>n_tresh),tmp.med_or(tmp.n_part>n_tresh),'Color',c(i,:),'linewidth',2);
%     figure(105);
%     fill_between_lines(Dmax_vec(tmp.n_part>n_tresh),tmp.q25_fs(tmp.n_part>n_tresh),tmp.q75_fs(tmp.n_part>n_tresh),c(i,:));
%     yoplot(end+1) = plot(Dmax_vec(tmp.n_part>n_tresh),tmp.med_fs(tmp.n_part>n_tresh),'Color',c(i,:),'linewidth',2);
%     
%     end
%     
% end
% 
% figure(103);
% xlabel('Dmax [mm]');
% ylabel('Aspect ratio');
% axis([0 12 0.1 1]);
% legend(yaplot,labels);
% %legend(labels);
% set(gca,'Fontsize',14);
% figure(104);
% xlabel('Dmax [mm]');
% ylabel('Orientation [deg]');
% legend(labels);
% set(gca,'Fontsize',14);
% figure(105);
% xlabel('Dmax [mm]');
% ylabel('Fallspeed [m/s]');
% set(gca,'Fontsize',14);
% legend(yoplot,labels);
% 
% 
% 
% % errorbar(Dmax_vec,med_ar,med_ar-q25_ar,q75_ar-med_ar,'k.-');
% % title('Aspect Ratio');
% % subplot(3,1,2); box on; hold on;
% % errorbar(Dmax_vec,mean_orientation,std_orientation,'k.-');
% % title('Canting angle');
% % subplot(3,1,3); box on; hold on;
% % plot(Dmax_vec,n_part,'k-.');
% % title('# particles');
% % set(gca,'yscale','log');
% % 
% 
% %% Structure for Daniel including D_eq, canting angle, AR for graupels and aggregates
% idx = find(X(:,6) >= 9.5 & X(:,11)>0.9 & X(:,1) == 2);
% Xfilt = X(idx,:);
% vec_id = 1:1:size(Xfilt,1);
% vec_id = vec_id';
% output_mat = [vec_id Xfilt(:,3)*0.035 Xfilt(:,9) Xfilt(:,10)];
% idx_weird = find(output_mat(:,2) < 0.6);
% Xname_weird = Xname(idx);
% Xname_weird = Xname_weird(idx_weird);
% 









