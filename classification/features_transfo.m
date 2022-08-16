% features transformation investigation

%% load data
clearvars; close all;

dir_graupels = '/home/praz/Documents/MASC/masclab/events/graupels';
graupels_filenames = dir(fullfile(dir_graupels,'*.mat'));
graupels_filenames = {graupels_filenames.name}';
graupels_picnames = dir(fullfile(dir_graupels,'*.png'));
graupels_picnames = {graupels_picnames.name}';

dir_agg = '/home/praz/Documents/MASC/masclab/events/aggregates';
agg_filenames = dir(fullfile(dir_agg,'*.mat'));
agg_filenames = {agg_filenames.name}';
agg_picnames = dir(fullfile(dir_agg,'*.png'));
agg_picnames = {agg_picnames.name}';

dir_melt = '/home/praz/Documents/MASC/masclab/events/melting_snow';
melt_filenames = dir(fullfile(dir_melt,'*.mat'));
melt_filenames = {melt_filenames.name}';
melt_picnames = dir(fullfile(dir_melt,'*.png'));
melt_picnames = {melt_picnames.name}';

dir_small = '/home/praz/Documents/MASC/masclab/events/small_particles';
small_filenames = dir(fullfile(dir_small,'*.mat'));
small_filenames = {small_filenames.name}';
small_picnames = dir(fullfile(dir_small,'*.png'));
small_picnames = {small_picnames.name}';

%all_filenames = [fullfile(dir_graupels,graupels_filenames); fullfile(dir_melt,melt_filenames)];
%all_picnames  = [fullfile(dir_graupels,graupels_picnames); fullfile(dir_melt,melt_picnames)];

t_str_start = '20150101000000';
t_str_stop  = '20170101000000';

% small validation test
% for i=1:length(melt_picnames)
%     
%     picname = melt_picnames{i};
%     filename = strcat(picname(1:end-3),'mat');
%     tmp = sum(strcmp(melt_filenames,filename));
%     if tmp ~= 1
%         disp(filename);
%     end
%     
% end

% load matrix data
[Xg,Xlab_g,Xname_g,Xt_g] = load_processed_data(dir_graupels,t_str_start,t_str_stop);
yg = zeros(size(Xg,1),1); yg(yg==0) = 1;
[Xa,Xlab_a,Xname_a,Xt_a] = load_processed_data(dir_agg,t_str_start,t_str_stop);
ya = zeros(size(Xa,1),1); ya(ya==0) = 2;
[Xm,Xlab_m,Xname_m,Xt_m] = load_processed_data(dir_melt,t_str_start,t_str_stop);
ym = zeros(size(Xm,1),1); ym(ym==0) = 3;
[Xs,Xlab_s,Xname_s,Xt_s] = load_processed_data(dir_small,t_str_start,t_str_stop);
ys = zeros(size(Xs,1),1); ys(ys==0) = 4;

% merge to one matrix
X = [Xg; Xa; Xm; Xs];
y = [yg; ya; ym; ys];
Xname = [Xname_g; Xname_a; Xname_m; Xname_s];

N = size(X,1);
D = size(X,2);

%% features illustration
if false
    figure;
    k=1;
    for i=1:20
        subplot(4,5,k);
        histogram(X(:,i));
        title(sprintf(Xlab_g{i}));
        k = k+1;
    end

    figure;
    k=1;
    for i=21:40
        subplot(4,5,k);
        histogram(X(:,i));
        title(sprintf(Xlab_g{i}));
        k = k+1;
    end

    figure;
    k=1;
    for i=41:61
        subplot(5,5,k);
        histogram(X(:,i));
        title(sprintf(Xlab_g{i}));
        k = k+1;
    end
end

%% normalization and feature transformation
normalization_type = 'standardization';
idx1 = [1,2,3,4,5,6,7,8,9,10,11,13,14,15,16,28,29,31,39,47,48,49,50,51,52,55];
idx2 = [27,32,33,34,35,38,50,51,52,55];
for i=1:D
    
    if sum(idx1==i) > 0 
        
        X(:,i) = log(X(:,i)+1);
       
    end

    if strcmp(normalization_type,'standardization') 
        
        tmp_mean = mean(X(:,i));
        tmp_std  = std(X(:,i));
        X(:,i) = (X(:,i) - tmp_mean)/tmp_std;
        
    elseif strcmp(normalization_type,'rescaling')
        
        tmp_min = min(X(:,i));
        tmp_max = max(X(:,i));
        X(:,i) = (X(:,i) - tmp_min)/(tmp_max - tmp_min);
        
    end
    
end

%% principal components illustration
[coeff,score,latent] = pca(X);
pc1 = score(:,1);
pc2 = score(:,2);
pc3 = score(:,3);

figure; hold on; box on; grid on;
plot(pc1(y==1),pc2(y==1),'bo','MarkerFaceColor','b','Markersize',5);
plot(pc1(y==2),pc2(y==2),'ro','MarkerFaceColor','r','Markersize',5);
plot(pc1(y==3),pc2(y==3),'co','MarkerFaceColor','c','Markersize',5);
plot(pc1(y==4),pc2(y==4),'go','MarkerFaceColor','g','Markersize',5);
legend('graupels','aggregates','melting snow','small part.');
xlabel('principal comp. #1');
ylabel('principal comp. #2');
set(gca,'Fontsize',14);

figure; hold on; box on; grid on;
plot(pc1(y==1),pc3(y==1),'bo','MarkerFaceColor','b','Markersize',5);
plot(pc1(y==2),pc3(y==2),'ro','MarkerFaceColor','r','Markersize',5);
plot(pc1(y==3),pc3(y==3),'co','MarkerFaceColor','c','Markersize',5);
plot(pc1(y==4),pc3(y==4),'go','MarkerFaceColor','g','Markersize',5);
legend('graupels','aggregates','melting snow','small part.');
xlabel('principal comp. #1');
ylabel('principal comp. #3');
set(gca,'Fontsize',14);

figure; hold on; box on; grid on;
plot(pc2(y==1),pc3(y==1),'bo','MarkerFaceColor','b','Markersize',5);
plot(pc2(y==2),pc3(y==2),'ro','MarkerFaceColor','r','Markersize',5);
plot(pc2(y==3),pc3(y==3),'co','MarkerFaceColor','c','Markersize',5);
plot(pc2(y==4),pc3(y==4),'go','MarkerFaceColor','g','Markersize',5);
legend('graupels','aggregates','melting snow','small part.');
xlabel('principal comp. #2');
ylabel('principal comp. #3');
set(gca,'Fontsize',14);

figure; hold on; box on; grid on;
plot(latent,'k-');
xlabel('principal components classified by order of importance');
ylabel('associated eigenvalues');
set(gca,'Yscale','log');
set(gca,'Fontsize',14);







