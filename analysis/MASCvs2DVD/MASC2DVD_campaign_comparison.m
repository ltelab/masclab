clear all; close all;

load_data = 1;

if load_data

    % load MASC data
    xhi_thresh = 9;
    Nmin_interval = 5;  
    tshift = 11; % T2dvd + shift = tmasc
    
    %masc = load('prediction/data/Xstruct_Davos_winter15-16_extended_last.mat');
    %masc.Xt = datetime(masc.Xt,'ConvertFrom','datenum');
    
    load('prediction/data/newFormat/classif_datastruct_Davos2015-16_triplet_withgraupfix.mat');
    masc.X = data.X;
    masc.Xt = data.Xt;
    masc.Xname = data.Xname;
    masc.Xlab = data.Xlab;
    masc.Xt = datetime(masc.Xt,'ConvertFrom','datenum');
    clear data;
    
    
    tstart_vec = [2015 10 12 00 00 00]; % june, 04.23, 11.20
    tstop_vec  = [2016 06 20 00 00 00];
    tstart = datetime(tstart_vec);
    tstop = datetime(tstop_vec);
    tgrid = tstart:minutes(Nmin_interval):tstop; tgrid = tgrid';
    tgrid2 = tgrid + minutes(5);

    % load 2DVD data
    dvd = load_2DVD_classif_campaign(datestr(tstart-minutes(tshift),'yyyymmddHHMMSS'),datestr(tstop,'yyyymmddHHMMSS'));
    dvd.t = datetime(dvd.t,'ConvertFrom','datenum');
    % timeshift between masc and 2dvd
    dvd.t = dvd.t + minutes(tshift);
    

    % merge both on the same grid
    j=1;
    for i=1:length(tgrid)
        idx_masc = find(masc.Xt >= tgrid(i) & masc.Xt <= tgrid2(i) & masc.X(:,6)>=xhi_thresh);
        idx_dvd = find(dvd.t >= tgrid(i) & dvd.t <= tgrid2(i));      
        if ~isempty(idx_masc) && ~isempty(idx_dvd)
            data.t(j,1) = tgrid2(i); 
            data.Nmasc(j,1) = numel(idx_masc);
            data.class(j,1) = mode(masc.X(idx_masc,1));
            data.melting(j,1) = mean(masc.X(idx_masc,16));
            data.riming(j,1) = mean(masc.X(idx_masc,15));
            data.Ndvd(j,1) = sum(dvd.N(idx_dvd));
            scores_2dvd = zeros(8,1);
            for k=1:length(idx_dvd)
                that_class = dvd.classif(idx_dvd(k));
                if ~isnan(that_class)
                    scores_2dvd(that_class) = scores_2dvd(that_class)+dvd.N(idx_dvd(k));
                end
            end
            [~,data.classdvd(j,1)] = max(scores_2dvd);
            %data.classdvd(j,1) = mode(dvd.classif(idx_dvd));
            j=j+1;
        end
    end   

    fprintf('%u timesteps of %u minutes in commmon found\n',numel(data.t),Nmin_interval);
    
else  
    %load('20151012_20160619__5min_1min_triplet_xhi9.mat');   
    %load('whole_season_MASC_2DVD_5min.mat');
    load('20151012_20160619__5min_5min_triplet_xhi9.mat');
end

%% fine tune the time shift
% shift = -2;
% data.t2 = data.t;
% data.t2 = data.t2 + minutes(shift);
% 
% [lia,locb] = ismember(data.t,data.t2);
% locb(locb==0) = [];
% data.Nmasc = data.Nmasc(lia);
% data.class = data.class(lia);
% data.melting = data.melting(lia);
% data.riming = data.riming(lia);
% data.Ndvd = data.Ndvd(locb);
% data.classdvd = data.classdvd(locb);



%% 
close all;
label = {'AG','DE','GR','SP+rain','CO','MS','R','RIM'};
Nmin_pts_masc = 30;
Nmin_pts_dvd = 400;

%1min: 12/200 -> OA:70.7    HSS:57  1340 datapoints
%5min: 30/500 -> OA:74.4    HSS:60  831 datapoints
%5/1min: 40/600 -> OA:75.4  HSS:63  2388 datapoints

% 30/200 pour jan-mars

% copy original data
data.cnew = data.class;
data.cnew2 = data.classdvd;

%  merge hydrometeor type
data.cnew(data.class==1) = 4;  %SP->rain+SP
data.cnew(data.class==2) = 5;  %CC->CO
data.cnew(data.class==3) = 2;  %PC->DE
data.cnew(data.class==4) = 1;  %AG->AG
data.cnew(data.class==5) = 3;  %GR->GR
data.cnew(data.class==6) = NaN;%CPC->NAN

% SP 2dvd -> NAN
data.cnew2(data.classdvd==7) = 4;

% if avg melting degree > thresh => MS
melting_thresh = 0.5; %0.25ini
data.cnew(data.melting >= melting_thresh) = 6;

% if avg riming degree > thresh => RIM
riming_thresh = 0.5;
%data.cnew(data.riming >= riming_thresh & data.cnew ~= 3) = 8;

% additionally merging rimed particles and aggregates
%data.cnew(data.cnew==8) = 1;
data.cnew2(data.classdvd==8) = 1;

% remove if not enough pts
data.cnew(data.Nmasc<Nmin_pts_masc) = NaN;
data.cnew2(data.Ndvd<Nmin_pts_dvd) = NaN;
   
% remove fucking nan-values before computing scores (GRRR)
idx = find(~isnan(data.cnew) & ~isnan(data.cnew2));
data.cnew = data.cnew(idx);
data.cnew2 = data.cnew2(idx);

% compute scores
[M,val] = confusionmat(data.cnew,data.cnew2);
OA = computeOA(data.cnew,data.cnew2);
HSS = computeKAPPA(data.cnew,data.cnew2);

figure('units','pixels','Position',[100 100 762 638]);
heatmap(M,label(val),label(val),1,'Colormap','red','ShowAllTicks',1,'UseLogColorMap',true,'Colorbar',true,'Fontsize',14);
set(gca,'Fontsize',14);
title('left:MASC bot:2DVD');

fprintf('OA : %2.1f    HSS : %2.1f     DIAG : %u \n',100*OA,100*HSS,sum(diag(M)));









