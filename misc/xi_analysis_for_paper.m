%% load data
clear all; %close all;

training_dir = '/home/praz/Documents/MASC/MASC_labelling/mixed_sample_merged_12345';
files = dir(fullfile(training_dir,'*.mat'));
files = {files.name};
files = files';

train = zeros(numel(files),2);
for i=1:numel(files)
    load(fullfile(training_dir,files{i}));
        
    roi.Dmean = 0.5*(roi.width+roi.height);
    roi.lap = fmeasure(roi.data,'LAPM',[]);   
    local_std = stdfilt(roi.data);
    roi.local_std = mean(local_std(roi.bw_mask_filled)); 
    
    roi.new.data = brightening(roi.data);
    roi.new.lap = fmeasure(roi.new.data,'LAPM',[]); 
    local_std = stdfilt(roi.new.data);
    roi.new.local_std = mean(local_std(roi.bw_mask_filled));
    
    xhi = log((roi.lap+roi.new.lap)/2 * roi.complex * (roi.local_std + roi.new.local_std)/2 * roi.Dmean);
    
    train(i,1) = xhi;
    train(i,2) = (size(roi.data,1)*size(roi.data,2))/2;
end

load('Xstruct_Davos_winter15-16_extended_last.mat');

full = zeros(length(Xt),2);
full(:,1) = X(:,6);
full(:,2) = X(:,4);

%% compute histograms
%nbins = 20;

%train(:,1) = train(:,1) -0.5;
%idx = find(train(:,1) < 9);
%mini = min(train(:,1));
%train(idx,1) = 9+rand/2.1;

idx_full = find(full(:,1)>=7);
[N1,X1] = hist(full(idx_full,1),[7:0.2:13]);

idx_train = find(train(:,1)>=7);
[N2,X2] = hist(train(idx_train,1),[7:0.2:13]);

figure; hold on; box on; grid on;
c1 = [43 140 190]./255;
area(X1,N1/(sum(N1)*0.2),'FaceAlpha',0.4,'FaceColor',c1);
idx_start = find(X2==9);
area(X2(idx_start:end),N2(idx_start:end)/(sum(N2)*0.2),'FaceColor','r','FaceAlpha',0.4);
plot([9 9], [0 1],'k--');
axis([7 13 0 1]);
xlabel('\xi');
legend('Davos 2015-16','training set');
set(gca,'Fontsize',14);
%histogram(full(:,1),[7:0.1:13],'Normalization','probability','FaceColor','b','FaceAlpha',0.5);
%histogram(train(:,1),[7:0.1:13],'Normalization','probability','FaceColor','r','FaceAlpha',0.5);