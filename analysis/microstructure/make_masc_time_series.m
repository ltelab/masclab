% attempt to create a general script to load MASC processed data and
% generate time series illustration, add other instrument output if
% desired, etc.
function make_masc_time_series(data,tstart_vec,tstop_vec,options,disp)

    MASC_pixres = options.pixres; % MASC pixel resolution in mm/pixel
    xi_thresh = options.xi_thresh; % quality parameter threshold
    Nmin_interval = options.Nmin_interval; % time window (min)
    Nmin_shift = options.Nmin_shift;
    Nclasses_masc = options.Nclasses_masc;
    MASC_classes = options.MASC_classes;
    MASC_classes_desired = options.MASC_classes_desired;
    N_MascSamples_min = options.N_MascSamples_min;
    use_triplet = options.use_triplet;
    verbose = options.verbose;
    OR_180 = options.OR_180;
    
    fs_quicklook1 = 18;
    fs_quicklook2 = 16;

    % 2dvd : not implemented
    
    % pluvio2 (Davos) : not implemented

    %% loading & gridding
    % load masc data
    % load(MASC_datafile);
    if use_triplet
        masc.Xt = datetime(data.Yt,'ConvertFrom','datenum');
        masc.X = data.Y;
        masc.Xlab = data.Ylab;
        masc.Xname = data.Yname;
    else
        masc.Xt = datetime(data.Xt,'ConvertFrom','datenum');
        masc.X = data.X;
        masc.Xlab = data.Xlab;
        masc.Xname = data.Xname;
    end

    clear data;

    tstart = datetime(tstart_vec);
    tstop = datetime(tstop_vec);
    tgrid = tstart:minutes(Nmin_shift):tstop; tgrid = tgrid';
    tgrid2 = tgrid + minutes(Nmin_interval);

    % load 2DVD data 

    % load atmospheric variables

    % load pluvio2

    % grid data (MASC done, need to add 2DVD is available)
    for i=1:length(tgrid)
        idx_masc_all = find(masc.Xt >= tgrid(i) & masc.Xt <= tgrid2(i)); % all images, used for geometric information
        if ~isempty(idx_masc_all)
            idx_masc_filtered = find(masc.Xt >= tgrid(i) & masc.Xt <= tgrid2(i) & masc.X(:,6) >= xi_thresh); 
            idx_masc_filtered_wo_SP = find(masc.Xt >= tgrid(i) & masc.Xt <= tgrid2(i) & masc.X(:,6) >= xi_thresh & masc.X(:,1) ~= 1);
            % filter fallspeed > 15 m/s or < 0.1 m/s
            idx_fallspeed_outlier = find(masc.X(:,7) >= 15 | masc.X(:,7) < 0.1);
            masc.X(idx_fallspeed_outlier,7) = NaN;
        else
            idx_masc_filtered = [];
            idx_masc_filtered_wo_SP = [];
        end
        if isempty(idx_masc_filtered) && verbose
            fprintf('No images found between %s and %s. \n',tgrid(i),tgrid2(i));
        end
        
        data.Nmasc_all(i,1) = numel(idx_masc_all);
        data.Nmasc_filtered(i,1) = numel(idx_masc_filtered);
        % classification related info
        if ~isempty(idx_masc_filtered)
            data.dom_class(i,1) = mode(masc.X(idx_masc_filtered,1));
            for k=1:Nclasses_masc
                data.sclass(i,k) = sum(masc.X(idx_masc_filtered,1)==k);
            end
            data.melting(i,1) = mean(masc.X(idx_masc_filtered,16));
            data.riming(i,1) = mean(masc.X(idx_masc_filtered_wo_SP,15));
        else
            data.dom_class(i,1) = NaN;
            data.sclass(i,1:Nclasses_masc) = NaN;
            data.melting(i,1) = NaN;
            data.riming(i,1) = NaN;
        end
        % geometry related info
        if ~isempty(idx_masc_all)
            % Dmax
            data.Dmax.mean(i,1) = nanmean(masc.X(idx_masc_all,3))*MASC_pixres;
            data.Dmax.med(i,1) = nanmedian(masc.X(idx_masc_all,3))*MASC_pixres;
            data.Dmax.std(i,1) = nanstd(masc.X(idx_masc_all,3))*MASC_pixres;
            % Nroi
            data.Nroi.mean(i,1) = nanmean(masc.X(idx_masc_all,28));
            data.Nroi.med(i,1) = nanmedian(masc.X(idx_masc_all,28));
            data.Nroi.std(i,1) = nanstd(masc.X(idx_masc_all,28));
            % aspect ratio
            data.AR.mean(i,1)  = nanmean(masc.X(idx_masc_all,9));
            data.AR.med(i,1) = nanmedian(masc.X(idx_masc_all,9));
            data.AR.std(i,1) = nanstd(masc.X(idx_masc_all,9));
            % area ratio
            data.ArR.mean(i,1) = nanmean(masc.X(idx_masc_all,29));
            data.ArR.med(i,1) = nanmedian(masc.X(idx_masc_all,29));
            data.ArR.std(i,1) = nanstd(masc.X(idx_masc_all,29));
            % orientation
            if OR_180
                OR_vec = masc.X(idx_masc_all,10);
            else
                OR_vec = abs(masc.X(idx_masc_all,10));
            end
            data.OR.mean(i,1) = nanmean(OR_vec);
            data.OR.med(i,1) = nanmedian(OR_vec);
            data.OR.std(i,1) = nanstd(OR_vec);
            % complexity
            data.cplx.mean(i,1) = nanmean(masc.X(idx_masc_all,5));
            data.cplx.med(i,1) = nanmedian(masc.X(idx_masc_all,5));
            data.cplx.std(i,1) = nanstd(masc.X(idx_masc_all,5));
            % fallspeed
            data.fs.mean(i,1) = nanmean(masc.X(idx_masc_all,7));
            data.fs.med(i,1) = nanmedian(masc.X(idx_masc_all,7));
            data.fs.std(i,1) = nanstd(masc.X(idx_masc_all,7));       


        else
            data.Dmax.mean(i,1) = NaN;
            data.Dmax.med(i,1) = NaN;
            data.Dmax.std(i,1) = NaN;
            data.Nroi.mean(i,1) = NaN;
            data.Nroi.med(i,1) = NaN;
            data.Nroi.std(i,1) = NaN;
            data.AR.mean(i,1)  = NaN;
            data.AR.med(i,1) = NaN;
            data.AR.std(i,1) = NaN;
            data.ArR.mean(i,1) = NaN;
            data.ArR.med(i,1) = NaN;
            data.ArR.std(i,1) = NaN;
            data.OR.mean(i,1) = NaN;
            data.OR.med(i,1) = NaN;
            data.OR.std(i,1) = NaN;
            data.cplx.mean(i,1) = NaN;
            data.cplx.med(i,1) = NaN;
            data.cplx.std(i,1) = NaN;
            data.fs.mean(i,1) = NaN;
            data.fs.med(i,1) = NaN;
            data.fs.std(i,1) = NaN;
        end

    end

    % normalization of the class proportion in every time interval
    data.sclass_norm = zeros(size(data.sclass,1),numel(MASC_classes_desired));
    for i=1:length(data.Nmasc_filtered)
       data.sclass_norm(i,:) = data.sclass(i,MASC_classes_desired)/sum(data.sclass(i,MASC_classes_desired));
    end

    % remove timesteps without enough samples
    data.dom_class(data.Nmasc_filtered < N_MascSamples_min) = NaN;
    data.sclass(data.Nmasc_filtered < N_MascSamples_min,:) = NaN;
    data.melting(data.Nmasc_filtered < N_MascSamples_min) = NaN;
    data.riming(data.Nmasc_filtered < N_MascSamples_min) = NaN;
    
    % check if there are some data within the interval to plot
    if sum(data.Nmasc_all) == 0
        options.savefigs = false;
    end
           
    %% illustration
    close all;

    % definition of the colormap used for hydrometeor classes 
    c = hsv(6); % rouge-jaune-vert-turquoise-bleu-rose
    cmasc = c; 
    cmasc(4,:) = c(6,:);
    cmasc(6,:) = [49,163,84]./255;
    cmelting = [241,105,19]./255;
    criming = 'b';

    % definition of a date string for figures title
    date_title = sprintf('%s - %s',datestr(tstart_vec,'yyyy.mm.dd HH:MM'),datestr(tstop_vec,'yyyy.mm.dd HH:MM'));
    
    if disp.now
        tnow = datetime('now','TimeZone','UTC');
        tgrid.TimeZone = 'UTC';
    end

    % Fig1 : the hydrometeor classification proportions
    if disp.ht_props

        fig1 = figure('Outerposition',[100 1000 1000 400]); hold on; box on;
        h = area(tgrid,data.sclass_norm,'EdgeColor','none'); 
        for i=1:numel(MASC_classes_desired)
            h(i).FaceColor = cmasc(MASC_classes_desired(i),:);
        end
        if disp.riming
            hline = []; hlabel = {};
            hline(end+1) = plot(tgrid,data.riming,'k-.','color','w','linewidth',3);
            hlabel{end+1} = 'R_i';
        end
        legend(MASC_classes(MASC_classes_desired));
        set(gca,'Ylim',[0 1]);
        set(gca,'YTick',[0 0.5 1]);
        ylabel('Proportions');
        title(date_title);
        set(gca,'Fontsize',14)

    end

    % Fig2 : the average degree of riming and proportions of melting snowflakes
    if disp.riming || disp.melting

        hline = []; hlabel = {};
        fig2 = figure('Outerposition',[150 650 1000 400]); hold on; grid on; box on;
        if disp.riming
            hline(end+1) = plot(tgrid,data.riming,'color',criming,'linewidth',2);
            hlabel{end+1} = 'R_i';
        end
        if disp.melting
            hline(end+1) = plot(tgrid,data.melting,'color',cmelting,'linewidth',2);
            hlabel{end+1} = '% wet';
        end
        ylabel('[-]');
        set(gca,'Ylim',[0 1]);
        set(gca,'YTick',[0 0.5 1]);
        title(date_title);
        legend(hline,hlabel);

    end

    % Fig3 : the total number of MASC images (all & filtered) in each time bin
    if disp.n_MASC

        fig3 = figure('Outerposition',[200 600 1000 400]); hold on; grid on; box on;
        plot(tgrid,data.Nmasc_all,'k-','linewidth',2);
        plot(tgrid,data.Nmasc_filtered,'k--','linewidth',2);
        if use_triplet
            ylabel('# MASC triplets');
        else
            ylabel('# MASC images');
        end
        legend('total',sprintf('filtered xi>%1.1f',xi_thresh));
        title(date_title);
    end

    % Fig4 : 1,2,3 on the same graph
    if disp.overview_classif

        fig4 = figure('Visible',disp.mode,'units','normalized','outerposition',[0 0 1 1]);
        subplot(4,1,1); hold on; grid on; box on;
        plot(tgrid,data.Nmasc_all,'k-','linewidth',2); 
        if xi_thresh > 0
            plot(tgrid,data.Nmasc_filtered,'k--','linewidth',2)
        end
        if use_triplet
            ylabel('# MASC triplets');
        else
            ylabel('# MASC images');
        end
        if disp.now
            yaxis = ylim;
           plot([tnow tnow],[yaxis(1) yaxis(2)],'r--','linewidth',2); 
        end
        if xi_thresh > 0
            legend('total',sprintf('filtered xi>%1.1f',xi_thresh));
        else
            legend('total');
        end
        title(date_title);
        set(gca,'Fontsize',fs_quicklook1);
        set(gca,'Xlim',[tgrid(1) tgrid(end)]);
        

        subplot(4,1,2:3); hold on; grid on; box on;
        h = area(tgrid,data.sclass_norm,'EdgeColor','none'); 
        for i=1:numel(MASC_classes_desired)
            h(i).FaceColor = cmasc(MASC_classes_desired(i),:);
        end
        if disp.riming
            hline = []; hlabel = {};
            hline(end+1) = plot(tgrid,data.riming,'k-.','color','w','linewidth',3);
            hlabel{end+1} = 'R_i';
        end
        set(gca,'Xlim',[tgrid(1) tgrid(end)]);
        set(gca,'Ylim',[0 1]);
        set(gca,'YTick',[0 0.5 1]);
        ylabel('Proportions');
        set(gca,'Fontsize',fs_quicklook1);
        if disp.now
            yaxis = ylim;
           plot([tnow tnow],[yaxis(1) yaxis(2)],'r--','linewidth',2); 
        end
        legend(MASC_classes(MASC_classes_desired));

        subplot(4,1,4); hold on; grid on; box on;
        hline = []; hlabel = {};
        hline(end+1) = plot(tgrid,data.riming,'color',criming,'linewidth',2);
        hlabel{end+1} = 'R_i';
        hline(end+1) = plot(tgrid,data.melting,'color',cmelting,'linewidth',2);
        hlabel{end+1} = '% wet';
        ylabel('[-]');
        set(gca,'Xlim',[tgrid(1) tgrid(end)]);
        set(gca,'Ylim',[0 1]);
        set(gca,'YTick',[0 0.5 1]);
        set(gca,'Fontsize',fs_quicklook1);
        if disp.now
            yaxis = ylim;
           plot([tnow tnow],[yaxis(1) yaxis(2)],'r--','linewidth',2); 
        end
        legend(hline,hlabel);

        if options.savefigs               
            saveas(fig4,fullfile(options.savedir,strcat(datestr(tgrid(1),'yyyymmdd'),'_quicklook#1')),'png');
        end
            
            
            
            
            %options.savepath options.save

    end

    % Fig5 : microstructural properties
    if disp.overview_microstruct

        CMAP = jet(20);

        fig5 = figure('Visible',disp.mode,'units','normalized','outerposition',[0 0 1 1]);

        % Dmax
        subplot(321); hold on; grid on; box on;
        hline = []; hlabel = {};
        plot(tgrid,nan(numel(tgrid)));
        fill_btw_curves(tgrid, data.Dmax.mean - data.Dmax.std, data.Dmax.mean + data.Dmax.std, data.Dmax.mean, 'r', 0.35);
        set(gca,'Xlim',[tgrid(1) tgrid(end)]);
        title('Dmax');
        ylabel('[mm]');
        if disp.now
           yaxis = ylim;
           hline(end+1) = plot([tnow tnow],[yaxis(1) yaxis(2)],'r--','linewidth',2); 
           hlabel{end+1} = 'Now';
           legend(hline,hlabel);
        end
        set(gca,'Fontsize',fs_quicklook2);

        % Fallspeed
        subplot(322); hold on; grid on; box on;
        plot(tgrid,nan(numel(tgrid)));
        fill_btw_curves(tgrid, data.fs.mean - data.fs.std, data.fs.mean + data.fs.std, data.fs.mean, 'b', 0.35);
        set(gca,'Xlim',[tgrid(1) tgrid(end)]);
        title('Fallspeed');
        ylabel('[m/s]');
        if disp.now
            yaxis = ylim;
           plot([tnow tnow],[yaxis(1) yaxis(2)],'r--','linewidth',2); 
        end
        set(gca,'Fontsize',fs_quicklook2);

        % Aspect Ratio
        subplot(323); hold on; grid on; box on;
        plot(tgrid,nan(numel(tgrid)));
        fill_btw_curves(tgrid, data.AR.mean - data.AR.std, data.AR.mean + data.AR.std, data.AR.mean, CMAP(15 ,:), 0.35);
        set(gca,'Xlim',[tgrid(1) tgrid(end)]);
        title('Aspect Ratio');
        ylabel('[-]');
        if disp.now
            yaxis = ylim;
           plot([tnow tnow],[yaxis(1) yaxis(2)],'r--','linewidth',2); 
        end
        set(gca,'Fontsize',fs_quicklook2);

        %fig6 = figure('Visible','on','units','normalized','outerposition',[0.5 0 0.5 1]);

        % Orientation 
        subplot(324); hold on; grid on; box on;
        plot(tgrid,nan(numel(tgrid)));
        fill_btw_curves(tgrid, data.OR.mean - data.OR.std, data.OR.mean + data.OR.std, data.OR.mean, 'r', 0.35);
        set(gca,'Xlim',[tgrid(1) tgrid(end)]);
        if OR_180
            title('Orientation [-90; 90]');
        else
            title('Orientation [0; 90]');
        end
        ylabel('[°]');
        if disp.now
            yaxis = ylim;
           plot([tnow tnow],[yaxis(1) yaxis(2)],'r--','linewidth',2); 
        end
        set(gca,'Fontsize',fs_quicklook2);

        % Complexity
        subplot(326); hold on; grid on; box on;
        plot(tgrid,nan(numel(tgrid)));
        fill_btw_curves(tgrid, data.cplx.mean - data.cplx.std, data.cplx.mean + data.cplx.std, data.cplx.mean, 'b', 0.35);
        set(gca,'Xlim',[tgrid(1) tgrid(end)]);
        title('Complexity');
        ylabel('[-]');
        if disp.now
            yaxis = ylim;
           plot([tnow tnow],[yaxis(1) yaxis(2)],'r--','linewidth',2); 
        end
        set(gca,'Fontsize',fs_quicklook2);


        % Area Ratio
        subplot(325); hold on; grid on; box on;
        plot(tgrid,nan(numel(tgrid)));
        fill_btw_curves(tgrid, data.ArR.mean - data.ArR.std, data.ArR.mean + data.ArR.std, data.ArR.mean, CMAP(15 ,:), 0.35);
        set(gca,'Xlim',[tgrid(1) tgrid(end)]);
        title('Area Ratio');
        ylabel('[-]');
        if disp.now
            yaxis = ylim;
           plot([tnow tnow],[yaxis(1) yaxis(2)],'r--','linewidth',2); 
        end
        set(gca,'Fontsize',fs_quicklook2);

        %plot_area_curve(t(idxs_start(i):idxs_stop(i))',Dmax_mat(:,idxs_start(i):idxs_stop(i))','','Dmax');
        
        if options.savefigs               
            saveas(fig5,fullfile(options.savedir,strcat(datestr(tgrid(1),'yyyymmdd'),'_quicklook#2')),'png');
        end

    end

end




% Fig2 : a pie chart of the HT proportions
% if disp.ht_piechart
% 
%     figure('Outerposition',[250 450 500 500]);
%     slices_labels = MASC_classes; 
%     spaces = {' ',' ',' ',' ',' ',' '};
%     perc_symbols = {'%','%','%','%','%','%'};
%     slices(1) = sum(data.sclass(:,1));
%     slices(2) = sum(data.sclass(:,2));
%     slices(3) = sum(data.sclass(:,3));
%     slices(4) = sum(data.sclass(:,4));
%     slices(5) = sum(data.sclass(:,5));
%     slices(6) = sum(data.sclass(:,6));
%     for j=1:6
%         perc(j) = round(100*slices(j)/sum(slices));
%     end
%     slices_labels = strcat(num2strs(perc),perc_symbols);
%     h = pie(slices,slices_labels);
%     hpie = findobj(h,'Type','patch');
%     for i=1:Nclasses_masc    
%         set(hpie(i),'FaceColor',cmasc(i,:));
%     end
%     title(date_title);
%     set(gca,'fontsize',12);
%     legend(MASC_classes);
% 
% end



%% specific illustrations


    

%title(sprintf('%s to %s',datestr(tstart,'yyyy.mm.dd - HH:MM'),datestr(tstop,'yyyy.mm.dd - HH:MM')));

% for Jacopo : Dmax, Nmasc and Nroi for 2 events (blowing vs snowfall)
% figure(756);
% subplot(311); hold on; grid on; box on;
% plot(tgrid,data.Nmasc_all,'k-');
% ylabel('# MASC img');
% set(gca,'Ylim',[0 1500]);
% subplot(312); hold on; grid on; box on;
% plot(tgrid, data.Nroi.med,'k-');
% ylabel('mean # ROIs/img');
% set(gca,'Ylim',[0 30]);
% subplot(313); hold on; grid on; box on;
% plot(tgrid,data.Dmax.mean,'k-');
% plot(tgrid,data.Dmax.mean+data.Dmax.std,'k--');
% plot(tgrid,data.Dmax.mean-data.Dmax.std,'k--');
% xlabel('Time [UTC]');
% ylabel('Dmax [mm]');
% set(gca,'Ylim',[0 3]);


% % subplot: MASC classification


