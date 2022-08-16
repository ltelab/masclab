% script to investigate different aspect ratio definition :
% - globally first
% - as a function of the hydrometeor size then
% - as a function of the hydrometeor type finally
% based on Davos winter 2015-2016

%% load data
clear all; close all;
datafile = 'data/triplet_NEW/Davos_winter_2015-16_36desc.mat';
tstart_vec = [2010 06 22 00 00 00]; 
tstop_vec  = [2020 08 05 22 00 00];
tstart_str = datestr(tstart_vec,'yyyymmddHHMMSS');
tstart = datenum(tstart_vec);
tstop_str  = datestr(tstop_vec,'yyyymmddHHMMSS');
tstop  = datenum(tstop_vec);
load(datafile); 

%% select data and filter

% user parameters
res = 33.5; % 1 pixel = "res" micrometers
use_triplet = false;
avoid_AR_lt_1 = true;

if use_triplet    
    X = data.Y;
    Xt = data.Yt;
    Xlab = data.Ylab;
    Xname = data.Yname;
    
else   
    X = data.X;
    Xt = data.Xt;
    Xlab = data.Xlab;
    Xname = data.Xname;
    
end

% select time interval !! avoid melting/rain events in Davos !
idx2keep = find(Xt>=tstart & Xt<=tstop);
X = X(idx2keep,:);
Xt = Xt(idx2keep);
Xname = Xname(idx2keep);

AR1 = data.X(:,9); % ellipse fit
AR2 = data.X(:,30); % D90
AR3 = data.X(:,31); % rectangle fit
AR4 = data.X(:,32); % ellipse out

if avoid_AR_lt_1
    AR1(AR1>1) = 1./AR1(AR1>1);
    AR2(AR2>1) = 1./AR2(AR2>1);
    AR3(AR3>1) = 1./AR3(AR3>1);
    AR4(AR4>1) = 1./AR4(AR4>1);
end

%% illustrations
close all;

% ARs raw scatterplots
MS_size = 3;
MS_color = 'b';

figure;
subplot(221); hold on; box on; grid on;
dscatter(AR1,AR2);
plot([0 1],[0 1],'r--');
[B1,~,~,~,STATS1] = regress(AR2,[ones(numel(AR1),1) AR1]);
plot([0 1],[B1(1) B1(1)+B1(2)*1],'g-');
axis([0 1 0 1]);
xlabel('ellipse fit AR');
ylabel('Dmax_{90}/Dmax AR');
%plot(AR1,AR2,'bo','MarkerSize',MS_size,'MarkerEdgeColor',MS_color,'MarkerFaceColor',MS_color);

subplot(222); hold on; box on; grid on;
dscatter(AR1,AR3);
plot([0 1],[0 1],'r--');
[B2,~,~,~,STATS2] = regress(AR3,[ones(numel(AR1),1) AR1]);
plot([0 1],[B2(1) B2(1)+B2(2)*1],'g-');
axis([0 1 0 1]);
xlabel('ellipse fit AR');
ylabel('rectangle out AR');
%plot(AR1,AR3,'bo','MarkerSize',MS_size,'MarkerEdgeColor',MS_color,'MarkerFaceColor',MS_color);

subplot(223); hold on; box on; grid on;
dscatter(AR1,AR4);
plot([0 1],[0 1],'r--');
[B3,~,~,~,STATS3] = regress(AR4,[ones(numel(AR1),1) AR1]);
plot([0 1],[B3(1) B3(1)+B3(2)*1],'g-');
axis([0 1 0 1]);
xlabel('ellipse fit AR');
ylabel('ellipse out AR');
%plot(AR1,AR4,'bo','MarkerSize',MS_size,'MarkerEdgeColor',MS_color,'MarkerFaceColor',MS_color);





