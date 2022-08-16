% small script to compare Davos and Vanderbilt images 
clear all, close all;
load('classif_datastruct_Vanderbilt_single_withgraupfix.mat');
data1 = data;
clear data;
load('classif_datastruct_Davos2015-16_triplet_withgraupfix.mat');
data2 = data;
clear data;

%%

bright1 = data1.X(:,26);
bright2 = data2.X(:,26);
bright2(bright2 <0.1) = NaN;
xhi1 = data1.X(:,6);
xhi2 = data2.X(:,6);

figure; hold on; box on; grid on;
histogram(bright1,'Normalization','probability');
histogram(bright2,'Normalization','probability');
legend('Greenland','Davos');
set(gca,'Yticklabel',[]);
set(gca,'Xlim',[0 1]);
xlabel('Mean snowflake brightness');
title('Davos (Swiss Alps) - Greenland comparison');
set(gca,'Fontsize',14);

% plot particles location on screen distribution
xhi_lim = 8;
X = data1.X;

discardmat = [400 550 320 410]; % [t,b,l,r]
figure;
subplot(2,3,1); hold on; box on;
idx = find(X(:,18)==0 & X(:,6)>= xhi_lim);
histogram(X(idx,19));
v = axis;
plot([discardmat(3) discardmat(3)],[v(3) v(4)],'r--');
set(gca,'XLim',[0,2448]);
subplot(2,3,2); hold on; box on;
idx = find(X(:,18)==1 & X(:,6)>= xhi_lim);
histogram(X(idx,19));
title('distribution of snowflake along the x-axis [0-2448]');
set(gca,'XLim',[0,2448]);
subplot(2,3,3); hold on; box on;
idx = find(X(:,18)==2 & X(:,6)>= xhi_lim);
histogram(X(idx,19));
v = axis;
plot([2448-discardmat(4) 2448-discardmat(4)],[v(3) v(4)],'r--');
set(gca,'XLim',[0,2448]);
subplot(2,3,4); hold on; box on;
idx = find(X(:,18)==0 & X(:,6)>= xhi_lim);
histogram(X(idx,20));
ylabel('camera 0');
v = axis;
plot([discardmat(1) discardmat(1)],[v(3) v(4)],'r--');
plot([2048-discardmat(2) 2048-discardmat(2)],[v(3) v(4)],'r--');
set(gca,'XLim',[0 2048]);
set(gca,'XDir','reverse');
view(90,-90);
subplot(2,3,5); hold on; box on;
idx = find(X(:,18)==1 & X(:,6)>= xhi_lim);
histogram(X(idx,20));
ylabel('camera 1');
title('distribution of snowflake along the y-axis [0-2048]');
v = axis;
plot([discardmat(1) discardmat(1)],[v(3) v(4)],'r--');
plot([2048-discardmat(2) 2048-discardmat(2)],[v(3) v(4)],'r--');
set(gca,'XLim',[0 2048]);
set(gca,'XDir','reverse');
view(90,-90);
subplot(2,3,6); hold on; box on;
idx = find(X(:,18)==2 & X(:,6)>= xhi_lim);
histogram(X(idx,20));
ylabel('camera 2');
v = axis;
plot([discardmat(1) discardmat(1)],[v(3) v(4)],'r--');
plot([2048-discardmat(2) 2048-discardmat(2)],[v(3) v(4)],'r--');
set(gca,'XLim',[0 2048]);
set(gca,'XDir','reverse');
view(90,-90);