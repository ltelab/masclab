% script to analyze snowflake microstructural properties as a function of
% the hydrometeor type and degree of riming
clear all; close all;

% user params
% classif_datastruct_Verbier_triplet.mat
% classif_datastruct_Davos2015-16_triplet_withgraupfix_more_entries.mat
% classif_datastruct_APRES3_triplet_withgraupfix.mat


figure_title = '';

datafile = 'Davos_winter_2015-16';
tstart_vec = [2010 06 22 00 00 00]; 
tstop_vec  = [2020 08 05 22 00 00];
res = 33.5; % 1 pixel = ??? mu

% loading
tstart_str = datestr(tstart_vec,'yyyymmddHHMMSS');
tstart = datenum(tstart_vec);
tstop_str  = datestr(tstop_vec,'yyyymmddHHMMSS');
tstop  = datenum(tstop_vec);
load(datafile); 

for k=1:2

    if k == 1
        X = data.Y;
        Xt = data.Yt;
        Xlab = data.Ylab;
        Xname = data.Yname;
    elseif k == 2
        X = data.X;
        Xt = data.Xt;
        Xlabl = data.Xlab;
        Xname = data.Xname;
    end

    % select time interval !! avoid melting/rain events in Davos !
    idx2keep = find(Xt>=tstart & Xt<=tstop);
    X = X(idx2keep,:);
    Xt = Xt(idx2keep);
    Xname = Xname(idx2keep);

    %% Dmax,AR,angle,complexity analysis
    %close all;

    ht = X(:,1);
    fs = X(:,7);
    Dmax = X(:,3) * 33.5/1000;
    ypos = X(:,20);
    xhi = X(:,6);
    rc = X(:,12);
    ri = X(:,15);
    melt = X(:,16);
    melt(isnan(melt)) = 0;
    AR = X(:,9);
    angle = abs(X(:,10));
    cplx = X(:,5);
    %cplx(cplx<1) = 1;

    % filters
    xhi_thresh = 8.5;

    % data filtering
    idx2keep = find(xhi > xhi_thresh);
    Dmax = Dmax(idx2keep);
    ht = ht(idx2keep);
    fs = fs(idx2keep);
    ypos = ypos(idx2keep);
    xhi = xhi(idx2keep);
    rc = rc(idx2keep);
    ri = ri(idx2keep);
    melt = melt(idx2keep);
    AR = AR(idx2keep);
    angle = angle(idx2keep);
    cplx = cplx(idx2keep);


    % general distributions (pdf)
    c = hsv(5);
    labels = {'small','col','plan','agg','grau'};
    Dmax_bins_size = 0.2;
    Dmax_bins_center = [0:Dmax_bins_size:19.75];
    AR_bins_size = 0.05;
    AR_bins_center = [0.025:AR_bins_size:0.975];
    angle_bins_size = 10;
    angle_bins_center = [5:angle_bins_size:85]; %  [2.5:5:87.5] if abs(..) [-87.5:5:87.5] otherwise
    cplx_bins_size = 0.2;
    cplx_bins_center = [1:cplx_bins_size:5];

    % pdf
    figure(1); hold on; grid on; box on;
    figure(2); hold on; grid on; box on;
    figure(3); hold on; grid on; box on;
    figure(4); hold on; grid on; box on;
    for i=1:6

        % filter data
        if i < 6
            idx = find(ht == i);
        elseif i == 6
            idx = find(melt == 1);
        end

        fprintf('i = %u, Nb flakes = %u \n',i,length(idx));
        if i==1
            disp(mode(Dmax(idx)));
            disp(mean(Dmax(idx)));
            disp(quantile(Dmax(idx),0.5));
            disp(quantile(Dmax(idx),0.75));
        end
        [N_Dmax,X_Dmax] = hist(Dmax(idx),Dmax_bins_center);
        [N_AR,X_AR] = hist(AR(idx),AR_bins_center); % 0.05
        [N_angle,X_angle] = hist(angle(idx),angle_bins_center); 
        [N_cplx,X_cplx] = hist(cplx(idx),cplx_bins_center);
        
        fprintf('k = %u : mean orientation for class %u : %2.2f \n',k,i,mean(angle(idx)));

        % normalization
        N_Dmax = N_Dmax./(sum(N_Dmax)*Dmax_bins_size);
        N_AR = N_AR./(sum(N_AR)*AR_bins_size);
        N_angle = N_angle./(sum(N_angle)*angle_bins_size);
        N_cplx = N_cplx./(sum(N_cplx)*cplx_bins_size);

        % plot
        if  i~= 6 && k == 1
            
            figure(1); hold on;
            plot(X_Dmax,N_Dmax,'Color',c(i,:),'linewidth',2);
            figure(2);
            plot(X_AR,N_AR,'Color',c(i,:),'linewidth',2);
            figure(3);
            plot(X_angle,N_angle,'Color',c(i,:),'linewidth',2);
            figure(4);
            plot(X_cplx,N_cplx,'Color',c(i,:),'linewidth',2);
            %figure(1102);
            %plot(X_fs,N_fs,'Color',c(i,:),'linewidth',2);
        
        elseif i~=6 && k == 2
            
            figure(1); hold on;
            plot(X_Dmax,N_Dmax,':','Color',c(i,:),'linewidth',1.5);
            figure(2);
            plot(X_AR,N_AR,':','Color',c(i,:),'linewidth',1.5);
            figure(3);
            plot(X_angle,N_angle,':','Color',c(i,:),'linewidth',1.5);
            figure(4);
            plot(X_cplx,N_cplx,':','Color',c(i,:),'linewidth',1.5);
            %figure(1102);
            %plot(X_fs,N_fs,'Color',c(i,:),'linewidth',2);
            
        elseif i==-99 % insert 6 to display melting snow
            
            figure(1); hold on;
            plot(X_Dmax,N_Dmax,'k--','linewidth',1.5);
            figure(2);
            plot(X_AR,N_AR,'k--','linewidth',1.5);
            figure(3);
            plot(X_angle,N_angle,'k--','linewidth',1.5);
            figure(4);
            plot(X_cplx,N_cplx,'k--','linewidth',1.5);
            
        end
    end
    
end

figure(1);
xlabel('Dmax [mm]');
ylabel('dN/d(Dmax)');
set(gca,'Xlim',[0 8]);
title(figure_title);
legend(labels);
set(gca,'Fontsize',12);
figure(2);
xlabel('AR');
ylabel('dN/d(AR)');
title(figure_title);
legend(labels);
set(gca,'Fontsize',12);
figure(3);
xlabel('\theta [°]');
ylabel('dN/d\theta');
title(figure_title);
set(gca,'Xlim',[5 85]);
legend(labels);
set(gca,'Fontsize',12);
figure(4);
xlabel('\chi');
ylabel('dN/d\chi');
set(gca,'Xlim',[1 4]);
legend(labels);
title(figure_title);
set(gca,'Fontsize',12);


%% distributions as a function of the particle type and diameter
if 1

    close all;

    fill_between_lines = @(X,Y1,Y2,C) fill( [X fliplr(X)],  [Y1 fliplr(Y2)], C, 'edgecolor','none','facealpha',.2);

    n_tresh = 100; % Davos 90; DDU 40 || PROB & x : Davos none, DDU xi9 and prob0.9
    Dmax_range = 0:0.2:20;
    Dmax_vec = Dmax_range+0.1;%double(unique(diff(Dmax_range)))/2;
    Dmax_vec = Dmax_vec(1:end-1);

    figure(11); hold on; grid on; box on;
    figure(12); hold on; grid on; box on;
    figure(13); hold on; grid on; box on;

    yaplot = [];
    yiplot = [];
    yoplot = [];

    for i=1:6

        % maybe add proba

        for j=1:numel(Dmax_range)-1
            lim_inf = Dmax_range(j);
            lim_sup = Dmax_range(j+1);
            
            if i<6
                idx_in = find(Dmax >= lim_inf & Dmax < lim_sup & ht == i);
            elseif i==6
                idx_in = find(Dmax >= lim_inf & Dmax < lim_sup & melt == 1);
            end
            
            if ~isempty(idx_in)
                % AR
                tmp.mean_ar(j) = nanmean(AR(idx_in));
                tmp.med_ar(j) = nanmedian(AR(idx_in));
                tmp.std_ar(j) = nanstd(AR(idx_in));
                tmp.q25_ar(j) = quantile(AR(idx_in),0.25);
                tmp.q75_ar(j) = quantile(AR(idx_in),0.75);

                % angle
                tmp.mean_angle(j) = nanmean(angle(idx_in));
                tmp.med_angle(j) = nanmedian(angle(idx_in));
                tmp.std_angle(j) = nanstd(angle(idx_in));
                tmp.q25_angle(j) = quantile(angle(idx_in),0.25);
                tmp.q75_angle(j) = quantile(angle(idx_in),0.75);

                % complexity
                tmp.mean_cplx(j) = nanmean(cplx(idx_in));
                tmp.med_cplx(j) = nanmedian(cplx(idx_in));
                tmp.std_cplx(j) = nanstd(cplx(idx_in));
                tmp.q25_cplx(j) = quantile(cplx(idx_in),0.25);
                tmp.q75_cplx(j) = quantile(cplx(idx_in),0.75);

                tmp.n_part(j) = length(idx_in);
            else
                % AR
                tmp.mean_ar(j) = NaN;
                tmp.med_ar(j) = NaN;
                tmp.std_ar(j) = NaN;
                tmp.q25_ar(j) = NaN;
                tmp.q75_ar(j) = NaN;

                % angle
                tmp.mean_angle(j) = NaN;
                tmp.med_angle(j) = NaN;
                tmp.std_angle(j) = NaN;
                tmp.q25_angle(j) = NaN;
                tmp.q75_angle(j) = NaN;

                % complexity
                tmp.mean_cplx(j) = NaN;
                tmp.med_cplx(j) = NaN;
                tmp.std_cplx(j) = NaN;
                tmp.q25_cplx(j) = NaN;
                tmp.q75_cplx(j) = NaN;

                tmp.n_part(j) = 0;
            end
        end

        if sum(find(tmp.n_part>n_tresh)) > 0
            if i < 6
                figure(11);
                fill_between_lines(Dmax_vec(tmp.n_part>n_tresh),tmp.q25_ar(tmp.n_part>n_tresh),tmp.q75_ar(tmp.n_part>n_tresh),c(i,:));
                yaplot(end+1) = plot(Dmax_vec(tmp.n_part>n_tresh),tmp.med_ar(tmp.n_part>n_tresh),'Color',c(i,:),'linewidth',2);
                figure(12);
                fill_between_lines(Dmax_vec(tmp.n_part>n_tresh),tmp.q25_angle(tmp.n_part>n_tresh),tmp.q75_angle(tmp.n_part>n_tresh),c(i,:));
                yiplot(end+1) = plot(Dmax_vec(tmp.n_part>n_tresh),tmp.med_angle(tmp.n_part>n_tresh),'Color',c(i,:),'linewidth',2);
                figure(13);
                fill_between_lines(Dmax_vec(tmp.n_part>n_tresh),tmp.q25_cplx(tmp.n_part>n_tresh),tmp.q75_cplx(tmp.n_part>n_tresh),c(i,:));
                yoplot(end+1) = plot(Dmax_vec(tmp.n_part>n_tresh),tmp.med_cplx(tmp.n_part>n_tresh),'Color',c(i,:),'linewidth',2);
            elseif i == 6
                figure(11);
                fill_between_lines(Dmax_vec(tmp.n_part>n_tresh),tmp.q25_ar(tmp.n_part>n_tresh),tmp.q75_ar(tmp.n_part>n_tresh),'k');
                yaplot(end+1) = plot(Dmax_vec(tmp.n_part>n_tresh),tmp.med_ar(tmp.n_part>n_tresh),'k--','linewidth',2);
                figure(12);
                fill_between_lines(Dmax_vec(tmp.n_part>n_tresh),tmp.q25_angle(tmp.n_part>n_tresh),tmp.q75_angle(tmp.n_part>n_tresh),'k');
                yiplot(end+1) = plot(Dmax_vec(tmp.n_part>n_tresh),tmp.med_angle(tmp.n_part>n_tresh),'k--','linewidth',2);
                figure(13);
                fill_between_lines(Dmax_vec(tmp.n_part>n_tresh),tmp.q25_cplx(tmp.n_part>n_tresh),tmp.q75_cplx(tmp.n_part>n_tresh),'k');
                yoplot(end+1) = plot(Dmax_vec(tmp.n_part>n_tresh),tmp.med_cplx(tmp.n_part>n_tresh),'k--','linewidth',2);
            end
        end

    end

    figure(11);
    xlabel('Dmax [mm]');
    ylabel('AR');
    title(figure_title);
    axis([0 14 0.1 1]);
    legend(yaplot,labels);
    set(gca,'Fontsize',14);
    figure(12);
    xlabel('Dmax [mm]');
    ylabel('\theta [°]');
    axis([0 14 0 90]);
    title(figure_title);
    legend(yiplot,labels);
    set(gca,'Fontsize',14);
    figure(13);
    xlabel('Dmax [mm]');
    ylabel('\chi');
    axis([0 14 0 4]);
    title(figure_title);
    set(gca,'Fontsize',14);
    legend(yoplot,labels);

end








