% small piece of code to compare Wprof LWC and MASC riming degree
clear all; close all;

% LOAD DATA

in.campaign = 'Valais_2017';
in.use_triplet = true;

if strcmp(in.campaign,'Valais_2017')
    
    LWC_file = '/home/kiko/Documents/PhD/Wprof/codes_matlab/MWR_grd_radvars_Valais_2017.mat';
    MASC_file = '/home/kiko/Documents/PhD/MASC/masclab/prediction/data/triplet_NEW/Verbier_winter_2016-17.mat';
    
end

load(LWC_file); % MWR
load(MASC_file); % data
data.Xt = datetime(data.Xt,'ConvertFrom','datenum');
data.Yt = datetime(data.Yt,'ConvertFrom','datenum');

%% LWC - RI EVENT COMPARISON

if 0

clearvars -except MWR data in
%close all;

MWR.LWC(MWR.LWC < 0 ) = 0;
MWR.grd_Ze = nanmean(MWR.grd_Ze,1);
MWR.grd_meanVel = nanmean(MWR.grd_meanVel);

N_MASC_thresh = 0;

in.use_triplet = true;
time_offset = minutes(60);

tvec_start = [[2017 01 12 23 40 00]; [2017 01 30 11 30 00];                        [2017 02 03 04 00 00]; [2017 02 03 16 00 00]; [2017 02 04 07 00 00]; [2017 02 04 12 50 00]; ...
    [2017 02 05 05 30 00]; [2017 02 07 11 50 00]; [2017 02 07 19 00 00];                        [2017 02 17 10 30 00]; [2017 02 24 00 00 00]; [2017 02 24 02 50 00]; [2017 02 24 07 40 00]; ...
    [2017 02 28 23 50 00]; [2017 03 04 09 50 00]; [2017 03 04 16 20 00]; [2017 03 04 23 50 00];                        [2017 03 06 03 00 00]; [2017 03 07 01 40 00]; ...
                                                  [2017 03 22 02 10 00]; [2017 03 22 13 40 00]; [2017 03 22 18 40 00]; [2017 03 22 22 20 00]; [2017 03 23 07 20 00]; ...
    [2017 03 25 18 50 00]; [2017 03 25 23 20 00]; [2017 03 26 02 50 00]; [2017 04 01 19 50 00]      ];
tvec_stop  = [[2017 01 14 13 40 00]; [2017 01 31 00 40 00];                        [2017 02 03 05 50 00]; [2017 02 03 21 20 00]; [2017 02 04 11 00 00]; [2017 02 04 17 00 00]; ...
    [2017 02 05 11 00 00]; [2017 02 07 19 00 00]; [2017 02 07 21 00 00];                        [2017 02 17 17 00 00]; [2017 02 24 01 50 00]; [2017 02 24 07 30 00]; [2017 02 24 10 30 00]; ...
    [2017 03 01 02 10 00]; [2017 03 04 15 10 00]; [2017 03 04 17 40 00]; [2017 03 05 01 50 00];                        [2017 03 07 01 40 00]; [2017 03 07 09 00 00]; ...
                                                  [2017 03 22 08 30 00]; [2017 03 22 17 40 00]; [2017 03 22 21 50 00]; [2017 03 23 01 00 00]; [2017 03 23 10 00 00]; ...
    [2017 03 25 19 40 00]; [2017 03 26 02 40 00]; [2017 03 26 05 00 00]; [2017 04 02 17 40 00];     ];

% >12h timesteps
if 0
    tvec_start = [tvec_start; [2017 01 31 05 40 00]; [2017 02 08 04 00 00]; [2017 03 05 13 50 00]; [2017 03 08 09 40 00]; [2017 03 18 05 40 00]];
    tvec_stop  = [tvec_stop;  [2017 01 31 22 30 00]; [2017 02 08 18 00 00]; [2017 03 06 03 00 00]; [2017 03 09 08 50 00]; [2017 03 18 20 10 00]];
end


t_start = datetime(tvec_start); 
t_stop  = datetime(tvec_stop);

t_start_LWP = t_start - time_offset;

if in.use_triplet
    
    X_MASC = data.Y;
    t_MASC = data.Yt;
    
else
    
    X_MASC = data.X;
    t_MASC = data.Xt;
    
end

% vectors initialization
N_MASC = zeros(numel(t_start),1);
Ri_MASC = zeros(numel(t_start),1);
N_LWP = zeros(numel(t_start),1);
LWP_mean = zeros(numel(t_start),1);
LWP_med = zeros(numel(t_start),1);
LWP_max = zeros(numel(t_start),1);
Ze = zeros(numel(t_start),1);
MeanVel = zeros(numel(t_start),1);

for i=1:length(t_start)
   
    % retrieve Ri
    idx_MASC = find(t_MASC >= t_start(i) & t_MASC < t_stop(i) & X_MASC(:,1) ~= 1 & X_MASC(:,6) >= 9 & X_MASC(:,16) >= 0);  
    N_MASC(i) = numel(idx_MASC);
    
    if N_MASC(i) > N_MASC_thresh
        Ri_MASC(i) = median(X_MASC(idx_MASC,15));
    else
        Ri_MASC(i) = NaN;
    end
    
    % retrieve LWP
    idx_MWR = find(MWR.time >= t_start_LWP(i) & MWR.time < t_stop(i));
    N_LWP(i) = numel(idx_MWR);
    LWP_mean(i) = nanmean(MWR.LWC(idx_MWR));
    LWP_med(i) = nanmedian(MWR.LWC(idx_MWR));
    LWP_max(i) = quantile(MWR.LWC(idx_MWR),1);
    
    % retrieve ground radar variables
    idx_Wprof = find(MWR.time >= t_start(i) & MWR.time < t_stop(i));
    N_Wprof(i) = numel(idx_Wprof);
    Ze(i) = nanmean(MWR.grd_Ze(idx_Wprof));
    MeanVel(i) = nanmean(MWR.grd_meanVel(idx_Wprof));
    
    
    
end

LWP_mean_bins = [0:15:250];
LWP_med_bins = [0:15:250];
LWP_max_bins = [0:75:1200];
MeanVel_bins = [-2.5:0.25:0];

Ri_MASC_avg_1 = zeros(numel(LWP_mean_bins)-1,1);
Ri_MASC_avg_2 = zeros(numel(LWP_med_bins)-1,1);
Ri_MASC_avg_3 = zeros(numel(LWP_max_bins)-1,1);
Ri_MASC_avg_4 = zeros(numel(MeanVel_bins)-1,1);
Ri_MASC_std_1 = zeros(numel(LWP_mean_bins)-1,1);
Ri_MASC_std_3 = zeros(numel(LWP_max_bins)-1,1);
Ri_MASC_std_4 = zeros(numel(MeanVel_bins)-1,1);

for i=1:numel(LWP_mean_bins)-1
    idx = find(LWP_mean >= LWP_mean_bins(i) & LWP_mean < LWP_mean_bins(i+1));
    if ~isempty(idx)
        Ri_MASC_avg_1(i) = mean(Ri_MASC(idx));
        Ri_MASC_std_1(i) = std(Ri_MASC(idx));
    else
        Ri_MASC_avg_1(i) = NaN;
        Ri_MASC_std_1(i) = NaN;
    end
end

for i=1:numel(LWP_med_bins)-1
    idx = find(LWP_med >= LWP_med_bins(i) & LWP_med < LWP_med_bins(i+1));
    if ~isempty(idx)
        Ri_MASC_avg_2(i) = mean(Ri_MASC(idx));
    else
        Ri_MASC_avg_2(i) = NaN;
    end
end

for i=1:numel(LWP_max_bins)-1
    idx = find(LWP_max >= LWP_max_bins(i) & LWP_max < LWP_max_bins(i+1));
    if ~isempty(idx)
        Ri_MASC_avg_3(i) = mean(Ri_MASC(idx));
        Ri_MASC_std_3(i) = std(Ri_MASC(idx));
    else
        Ri_MASC_avg_3(i) = NaN;
        Ri_MASC_std_3(i) = NaN;
    end
end

for i=1:numel(MeanVel_bins)-1
    idx = find(MeanVel >= MeanVel_bins(i) & MeanVel < MeanVel_bins(i+1));
    if ~isempty(idx)
        Ri_MASC_avg_4(i) = mean(Ri_MASC(idx));
        Ri_MASC_std_4(i) = std(Ri_MASC(idx));
    else
        Ri_MASC_avg_4(i) = NaN;
        Ri_MASC_std_4(i) = NaN;
    end
end

Ri_MASC_std_1(10) = Ri_MASC_std_1(10)*5;

figure; 
grid on; hold on; box on;
%scatter(LWP_mean,Ri_MASC,'kx');
errorbar(LWP_mean_bins(2:end),Ri_MASC_avg_1,Ri_MASC_std_1,'ro','MarkerFaceColor','r');
xlabel('mean event LWP [g/m^2]');
ylabel('mean event R_i');
set(gca,'Ylim',[0.4 0.9]);
set(gca,'Xlim',[0 160]);
set(gca,'Fontsize',12);
% subplot(132); grid on; hold on; box on;
% title('LWP vs R_i at the event scale');
% scatter(LWP_med,Ri_MASC,'kx');
% plot(LWP_med_bins(2:end),Ri_MASC_avg_2,'ro');
% xlabel('std event LWP [g/m^2]');
% ylabel('mean event R_i');
%set(gca,'Ylim',[0.4 0.8]);
%subplot(122); grid on; hold on; box on;
%scatter(LWP_max,Ri_MASC,'kx');
%errorbar(LWP_max_bins(2:end),Ri_MASC_avg_3,Ri_MASC_std_3,'ro','MarkerFaceColor','r');
%xlabel('max event LWP [g/m^2]');
%ylabel('mean event R_i');
%set(gca,'Ylim',[0.4 0.8]);

%figure;
%subplot(121); grid on; hold on; box on;
%scatter(Ze,Ri_MASC,'kx');
%subplot(122); grid on; hold on; box on;
%scatter(MeanVel,Ri_MASC,'kx');
%errorbar(MeanVel_bins(2:end),Ri_MASC_avg_4,Ri_MASC_std_4,'ro','MarkerFaceColor','r');

end


%% LWC - RI COMPARISON

if 1

    clearvars -except MWR data in

    in.use_triplet = true;

    % 12 Jan 2017 - 04 Apr 2017
    tvec_start = [2017 01 12 00 00 00];
    tvec_stop  = [2017 04 05 00 00 00];

    t_start = datetime(tvec_start);
    t_stop  = datetime(tvec_stop);

    in.t_win = minutes(30);
    in.t_shift = minutes(30);

    t_l = transpose(t_start:in.t_shift:t_stop-in.t_win);
    t_r = transpose(t_start+in.t_win:in.t_shift:t_stop);

    if in.use_triplet

        X_MASC = data.Y;
        t_MASC = data.Yt;

    else

        X_MASC = data.X;
        t_MASC = data.Xt;

    end
    xlabel('mean event LWP [g/m^2]');
    ylabel('mean event R_i');
    % pre-filtering phase
    idx_filt = find(X_MASC(:,1) ~= 1 & X_MASC(:,6) >= 9);
    N_ini = size(X_MASC,1);
    X_MASC = X_MASC(idx_filt,:);
    t_MASC = t_MASC(idx_filt);
    fprintf('Prefiltering removed %2.2f%% of MASC data \n',(N_ini-numel(idx_filt))/N_ini*100);

    % vectors initialization
    N_MASC = zeros(numel(t_l),1);
    Ri_MASC = zeros(numel(t_l),1);
    N_LWP = zeros(numel(t_l),1);
    LWP = zeros(numel(t_l),1);
    LWP_max = zeros(numel(t_l),1);
    Tb = zeros(numel(t_l),1);
    Tb_max = zeros(numel(t_l),1);
    MeanVel = zeros(numel(t_l),1);

    for i=1:length(t_l)

        % retrieve Ri
        idx_MASC = find(t_MASC >= t_l(i) & t_MASC < t_r(i)); % <--- add more filters here
        N_MASC(i) = numel(idx_MASC);
        Ri_MASC(i) = mean(X_MASC(idx_MASC,15));

        % retrieve LWP
        idx_MWR = find(MWR.time >= t_l(i) & MWR.time < t_r(i));
        N_LWP(i) = numel(idx_MWR);
        if ~isempty(idx_MWR)
            LWP(i) = mean(MWR.LWC(idx_MWR));
            LWP_max(i) = max(MWR.LWC(idx_MWR));
            Tb(i) = mean(MWR.bright_T(idx_MWR));
            Tb_max(i) = max(MWR.bright_T(idx_MWR));
        else
            LWP(i) = NaN;
            LWP_max(i) = NaN;
            Tb(i) = NaN;
            Tb_max(i) = NaN;
        end

        % retrieve ground radvars
        if ~isempty(idx_MWR)
            MeanVel(i) = nanmean(nanmin(MWR.grd_meanVel(:,idx_MWR)));
        else
            MeanVel(i) = NaN;
        end



    end

    LWP_ini = LWP;

    %% ILLUSTRATIONS
    close all;

    N_MASC_thresh = 0;
    LWP_thresh = 0;
    LWP = LWP_max;


    idx_plot = find(N_MASC > N_MASC_thresh & LWP > LWP_thresh);
    offset_20min = round(minutes(20)/in.t_shift);
    LWP_off1 = wshift('1D',LWP,1*offset_20min); LWP_off1(1:1*offset_20min) = NaN;
    LWP_off2 = wshift('1D',LWP,2*offset_20min); LWP_off1(1:2*offset_20min) = NaN;
    LWP_off3 = wshift('1D',LWP,3*offset_20min); LWP_off1(1:3*offset_20min) = NaN;

    figure;
    subplot(311); hold on; grid on; box on;
    plot(t_l,N_MASC,'k-');
    subplot(312);  hold on; grid on; box on;
    plot(t_l,Ri_MASC,'k-');
    subplot(313);  hold on; grid on; box on;
    plot(t_l,LWP,'k-');


    figure;
    subplot(221); hold on; grid on; box on;
    dscatter(LWP(idx_plot),Ri_MASC(idx_plot));
    xlabel('LWP [g/m^2]');
    ylabel('R_i');
    axis([0 1000 0 1]);
    title('Offset LWP = 0 min');
    subplot(222); hold on; grid on; box on;
    idx_plot = find(N_MASC > N_MASC_thresh & LWP_off1 > LWP_thresh);
    dscatter(LWP_off1(idx_plot),Ri_MASC(idx_plot));
    xlabel('LWP [g/m^2]');
    ylabel('R_i');
    axis([0 1000 0 1]);
    title('Offset LWP = -20 min');
    subplot(223); hold on; grid on; box on;
    idx_plot = find(N_MASC > N_MASC_thresh & LWP_off2 > LWP_thresh);
    dscatter(LWP_off2(idx_plot),Ri_MASC(idx_plot));
    xlabel('LWP [g/m^2]');
    ylabel('R_i');
    axis([0 1000 0 1]);
    title('Offset LWP = -40 min');
    subplot(224); hold on; grid on; box on;
    idx_plot = find(N_MASC > N_MASC_thresh & LWP_off3 > LWP_thresh);
    dscatter(LWP_off3(idx_plot),Ri_MASC(idx_plot));
    xlabel('LWP [g/m^2]');
    ylabel('R_i');
    axis([0 1000 0 1]);
    title('Offset LWP = -60 min');


    MeanVel = -MeanVel;

    idx_plot = find(N_MASC > N_MASC_thresh & ~isnan(MeanVel));
    MeanVel_bins = [];
    MeanVel_avg = [];
    MeanVel_std = [];
    MeanVel_bins = [0.5:0.1:2.5];
    for i=1:length(MeanVel_bins)-1

        idx = find(N_MASC > N_MASC_thresh & MeanVel >= MeanVel_bins(i) & MeanVel < MeanVel_bins(i+1));
        if ~isempty(idx)
            MeanVel_avg(i) = mean(Ri_MASC(idx));
            MeanVel_std(i) = std(Ri_MASC(idx));
        else
            MeanVel_avg(i) = NaN;
            MeanVel_std(i) = NaN;
        end

    end

    figure; hold on; grid on; box on;
    scatter(MeanVel(idx_plot),Ri_MASC(idx_plot),'k.');
    errorbar((MeanVel_bins(2:end)+MeanVel_bins(1:end-1))./2,MeanVel_avg,MeanVel_std,'ro','MarkerFaceColor','r');
    xlabel('30 min. Doppler mean vel. < 150m [m/s]');
    ylabel('30 min. mean R_i');
    set(gca,'Xlim',[0.5 2.5]);
    set(gca,'Ylim',[0.3 0.9]);
    set(gca,'Fontsize',12);

end



    
    
    
    
    
    


    