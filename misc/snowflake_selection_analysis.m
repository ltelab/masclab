% snowflake selection analysis based on a few features
%close all; clear all;
clear all; close all;

main_dir = '/home/praz/Documents/MASC/sample_snow_20150620/sample_cool/PROCESSED_TO_CLASSIFY_NEW/DATA';

%% GOOD detection
file_list = dir(fullfile(main_dir,'NICE','*.mat'));
file_list = {file_list.name};
file_list = fullfile(main_dir,'NICE',file_list);

% features based on raw picture
good.area = zeros(length(file_list),1);
good.dim = zeros(length(file_list),1);
good.mean_intens = zeros(length(file_list),1);
good.max_intens = zeros(length(file_list),1);
good.lap = zeros(length(file_list),1);
good.complex = zeros(length(file_list),1);
good.std = zeros(length(file_list),1);
good.local_std = zeros(length(file_list),1);
good.local_std5 = zeros(length(file_list),1);
good.local_std7 = zeros(length(file_list),1);


% % features based on enhanced picture
good.new.lap = zeros(length(file_list),1);
good.new.std = zeros(length(file_list),1);
good.new.local_std = zeros(length(file_list),1);
good.new.local_std5 = zeros(length(file_list),1);
good.new.local_std7 = zeros(length(file_list),1);

good.flag = zeros(length(file_list),1);

% loop over the snowflakes and save features in a vector
for i=1:length(file_list)
    
    load(file_list{i});
    good.area(i) = roi.area;
    good.dim(i) = 0.5*double(roi.width + roi.height);
    good.mean_intens(i) = roi.mean_intens;
    good.max_intens(i) = roi.max_intens;
    good.lap(i) = roi.lap;
    good.complex(i) = roi.complex;
    good.std(i) = roi.std;
    good.local_std(i) = roi.local_std;
    good.local_std5(i) = roi.local_std5;
    good.local_std7(i) = roi.local_std7;
    
    good.new.lap(i) = roi.new.lap;
    good.new.std(i) = roi.new.std;
    good.new.local_std(i) = roi.new.local_std;
    good.new.local_std5(i) = roi.new.local_std5;
    good.new.local_std7(i) = roi.new.local_std7;
    
    good.flag(i) = 4;
    
end

%% BAD detection
file_list = dir(fullfile(main_dir,'BAD','*.mat'));
file_list = {file_list.name};
file_list = fullfile(main_dir,'BAD',file_list);

% features based on raw picture
bad.area = zeros(length(file_list),1);
bad.dim = zeros(length(file_list),1);
bad.mean_intens = zeros(length(file_list),1);
bad.max_intens = zeros(length(file_list),1);
bad.lap = zeros(length(file_list),1);
bad.complex = zeros(length(file_list),1);
bad.std = zeros(length(file_list),1);
bad.local_std = zeros(length(file_list),1);
bad.local_std5 = zeros(length(file_list),1);
bad.local_std7 = zeros(length(file_list),1);


% % features based on enhanced picture
bad.new.lap = zeros(length(file_list),1);
bad.new.std = zeros(length(file_list),1);
bad.new.local_std = zeros(length(file_list),1);
bad.new.local_std5 = zeros(length(file_list),1);
bad.new.local_std7 = zeros(length(file_list),1);

bad.flag = zeros(length(file_list),1);

% loop over the snowflakes and save features in a vector
for i=1:length(file_list)
    
    load(file_list{i});
    bad.area(i) = roi.area;
    bad.dim(i) = 0.5*double(roi.width + roi.height);
    bad.mean_intens(i) = roi.mean_intens;
    bad.max_intens(i) = roi.max_intens;
    bad.lap(i) = roi.lap;
    bad.complex(i) = roi.complex;
    bad.std(i) = roi.std;
    bad.local_std(i) = roi.local_std;
    bad.local_std5(i) = roi.local_std5;
    bad.local_std7(i) = roi.local_std7;
    
    bad.new.lap(i) = roi.new.lap;
    bad.new.std(i) = roi.new.std;
    bad.new.local_std(i) = roi.new.local_std;
    bad.new.local_std5(i) = roi.new.local_std5;
    bad.new.local_std7(i) = roi.new.local_std7;
    
    bad.flag(i) = 0;
    
end

%% GREAT detection
file_list = dir(fullfile(main_dir,'AWESOME','*.mat'));
file_list = {file_list.name};
file_list = fullfile(main_dir,'AWESOME',file_list);

% features based on raw picture
great.area = zeros(length(file_list),1);
great.dim = zeros(length(file_list),1);
great.mean_intens = zeros(length(file_list),1);
great.max_intens = zeros(length(file_list),1);
great.lap = zeros(length(file_list),1);
great.complex = zeros(length(file_list),1);
great.std = zeros(length(file_list),1);
great.local_std = zeros(length(file_list),1);
great.local_std5 = zeros(length(file_list),1);
great.local_std7 = zeros(length(file_list),1);

great.flag = zeros(length(file_list),1);


% % features based on enhanced picture
great.new.lap = zeros(length(file_list),1);
great.new.std = zeros(length(file_list),1);
great.new.local_std = zeros(length(file_list),1);
great.new.local_std5 = zeros(length(file_list),1);
great.new.local_std7 = zeros(length(file_list),1);

% loop over the snowflakes and save features in a vector
for i=1:length(file_list)
    
    load(file_list{i});
    great.area(i) = roi.area;
    great.dim(i) = 0.5*double(roi.width + roi.height);
    great.mean_intens(i) = roi.mean_intens;
    great.max_intens(i) = roi.max_intens;
    great.lap(i) = roi.lap;
    great.complex(i) = roi.complex;
    great.std(i) = roi.std;
    great.local_std(i) = roi.local_std;
    great.local_std5(i) = roi.local_std5;
    great.local_std7(i) = roi.local_std7; 
    
    great.new.lap(i) = roi.new.lap;
    great.new.std(i) = roi.new.std;
    great.new.local_std(i) = roi.new.local_std;
    great.new.local_std5(i) = roi.new.local_std5;
    great.new.local_std7(i) = roi.new.local_std7;
    
    great.flag(i) = 5;
    
end

%% BLURRY detection
file_list = dir(fullfile(main_dir,'BLURRY','*.mat'));
file_list = {file_list.name};
file_list = fullfile(main_dir,'BLURRY',file_list);

% features based on raw picture
blurry.area = zeros(length(file_list),1);
blurry.dim = zeros(length(file_list),1);
blurry.mean_intens = zeros(length(file_list),1);
blurry.max_intens = zeros(length(file_list),1);
blurry.lap = zeros(length(file_list),1);
blurry.complex = zeros(length(file_list),1);
blurry.std = zeros(length(file_list),1);
blurry.local_std = zeros(length(file_list),1);
blurry.local_std5 = zeros(length(file_list),1);
blurry.local_std7 = zeros(length(file_list),1);


% % features based on enhanced picture
blurry.new.lap = zeros(length(file_list),1);
blurry.new.std = zeros(length(file_list),1);
blurry.new.local_std = zeros(length(file_list),1);
blurry.new.local_std5 = zeros(length(file_list),1);
blurry.new.local_std7 = zeros(length(file_list),1);

blurry.flag = zeros(length(file_list),1);

% loop over the snowflakes and save features in a vector
for i=1:length(file_list)
    
    load(file_list{i});
    blurry.area(i) = roi.area;
    blurry.dim(i) = 0.5*double(roi.width + roi.height);
    blurry.mean_intens(i) = roi.mean_intens;
    blurry.max_intens(i) = roi.max_intens;
    blurry.lap(i) = roi.lap;
    blurry.complex(i) = roi.complex;
    blurry.std(i) = roi.std;
    blurry.local_std(i) = roi.local_std;
    blurry.local_std5(i) = roi.local_std5;
    blurry.local_std7(i) = roi.local_std7; 
    
    blurry.new.lap(i) = roi.new.lap;
    blurry.new.std(i) = roi.new.std;
    blurry.new.local_std(i) = roi.new.local_std;
    blurry.new.local_std5(i) = roi.new.local_std5;
    blurry.new.local_std7(i) = roi.new.local_std7;
    
    blurry.flag(i) = 3;
    
end

%% VERY_BLURRY detection
file_list = dir(fullfile(main_dir,'VERY_BLURRY','*.mat'));
file_list = {file_list.name};
file_list = fullfile(main_dir,'VERY_BLURRY',file_list);

% features based on raw picture
v_blurry.area = zeros(length(file_list),1);
v_blurry.dim = zeros(length(file_list),1);
v_blurry.mean_intens = zeros(length(file_list),1);
v_blurry.max_intens = zeros(length(file_list),1);
v_blurry.lap = zeros(length(file_list),1);
v_blurry.complex = zeros(length(file_list),1);
v_blurry.std = zeros(length(file_list),1);
v_blurry.local_std = zeros(length(file_list),1);
v_blurry.local_std5 = zeros(length(file_list),1);
v_blurry.local_std7 = zeros(length(file_list),1);


% % features based on enhanced picture
v_blurry.new.lap = zeros(length(file_list),1);
v_blurry.new.std = zeros(length(file_list),1);
v_blurry.new.local_std = zeros(length(file_list),1);
v_blurry.new.local_std5 = zeros(length(file_list),1);
v_blurry.new.local_std7 = zeros(length(file_list),1);

v_blurry.flag = zeros(length(file_list),1);

% loop over the snowflakes and save features in a vector
for i=1:length(file_list)
    
    load(file_list{i});
    v_blurry.area(i) = roi.area;
    v_blurry.dim(i) = 0.5*double(roi.width + roi.height);
    v_blurry.mean_intens(i) = roi.mean_intens;
    v_blurry.max_intens(i) = roi.max_intens;
    v_blurry.lap(i) = roi.lap;
    v_blurry.complex(i) = roi.complex;
    v_blurry.std(i) = roi.std;
    v_blurry.local_std(i) = roi.local_std;
    v_blurry.local_std5(i) = roi.local_std5;
    v_blurry.local_std7(i) = roi.local_std7; 
    
    v_blurry.new.lap(i) = roi.new.lap;
    v_blurry.new.std(i) = roi.new.std;
    v_blurry.new.local_std(i) = roi.new.local_std;
    v_blurry.new.local_std5(i) = roi.new.local_std5;
    v_blurry.new.local_std7(i) = roi.new.local_std7;
    
    v_blurry.flag(i) = 2;
    
end

%% SHIT detection
file_list = dir(fullfile(main_dir,'SHIT','*.mat'));
file_list = {file_list.name};
file_list = fullfile(main_dir,'SHIT',file_list);

% features based on raw picture
shit.area = zeros(length(file_list),1);
shit.dim = zeros(length(file_list),1);
shit.mean_intens = zeros(length(file_list),1);
shit.max_intens = zeros(length(file_list),1);
shit.lap = zeros(length(file_list),1);
shit.complex = zeros(length(file_list),1);
shit.std = zeros(length(file_list),1);
shit.local_std = zeros(length(file_list),1);
shit.local_std5 = zeros(length(file_list),1);
shit.local_std7 = zeros(length(file_list),1);


% % features based on enhanced picture
shit.new.lap = zeros(length(file_list),1);
shit.new.std = zeros(length(file_list),1);
shit.new.local_std = zeros(length(file_list),1);
shit.new.local_std5 = zeros(length(file_list),1);
shit.new.local_std7 = zeros(length(file_list),1);

shit.flag = zeros(length(file_list),1);

% loop over the snowflakes and save features in a vector
for i=1:length(file_list)
    
    load(file_list{i});
    shit.area(i) = roi.area;
    shit.dim(i) = 0.5*double(roi.width + roi.height);
    shit.mean_intens(i) = roi.mean_intens;
    shit.max_intens(i) = roi.max_intens;
    shit.lap(i) = roi.lap;
    shit.complex(i) = roi.complex;
    shit.std(i) = roi.std;
    shit.local_std(i) = roi.local_std;
    shit.local_std5(i) = roi.local_std5;
    shit.local_std7(i) = roi.local_std7;
    
    shit.new.lap(i) = roi.new.lap;
    shit.new.std(i) = roi.new.std;
    shit.new.local_std(i) = roi.new.local_std;
    shit.new.local_std5(i) = roi.new.local_std5;
    shit.new.local_std7(i) = roi.new.local_std7;
    
    shit.flag(i) = 1;
    
end

%% ALL in the same vector
all.area = [great.area; good.area; blurry.area; v_blurry.area; shit.area];
all.dim = [great.dim; good.dim; blurry.dim; v_blurry.dim; shit.dim];
all.mean_intens = [great.mean_intens; good.mean_intens; blurry.mean_intens; v_blurry.mean_intens; shit.mean_intens];
all.max_intens = [great.max_intens; good.max_intens; blurry.max_intens; v_blurry.max_intens; shit.max_intens];
all.lap = [great.lap; good.lap; blurry.lap; v_blurry.lap; shit.lap];
all.complex = [great.complex; good.complex; blurry.complex; v_blurry.complex; shit.complex];
all.std = [great.std; good.std; blurry.std; v_blurry.std; shit.std];
all.local_std = [great.local_std; good.local_std; blurry.local_std; v_blurry.local_std; shit.local_std];
all.local_std5 = [great.local_std5; good.local_std5; blurry.local_std5; v_blurry.local_std5; shit.local_std5];
all.local_std7 = [great.local_std7; good.local_std7; blurry.local_std7; v_blurry.local_std7; shit.local_std7];
all.newlap = [great.new.lap; good.new.lap; blurry.new.lap; v_blurry.new.lap; shit.new.lap];
all.newstd = [great.new.std; good.new.std; blurry.new.std; v_blurry.new.std; shit.new.std];
all.newlocalstd = [great.new.local_std; good.new.local_std; blurry.new.local_std; v_blurry.new.local_std; shit.new.local_std];
all.newlocalstd5 = [great.new.local_std5; good.new.local_std5; blurry.new.local_std5; v_blurry.new.local_std5; shit.new.local_std5];
all.newlocalstd7 = [great.new.local_std7; good.new.local_std7; blurry.new.local_std7; v_blurry.new.local_std7; shit.new.local_std7];
all.flag = [great.flag; good.flag; blurry.flag; v_blurry.flag; shit.flag];

%% Clustering - Illustration
c = winter(5);

x = all.newlocalstd5;
y = all.dim;

figure(1);
hold on; box on;
for i=1:5
    plot(x(all.flag==i),y(all.flag==i),'o','MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor',c(i,:));
end
xlabel('mean of local std');
ylabel('dimension');
legend('bad','very blurry','blurry','good','great');

%% 
x = log(all.newlocalstd5);
y = log(all.dim);

figure;
hold on; box on;
for i=1:5
    plot(x(all.flag==i),y(all.flag==i),'o','MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor',c(i,:));
end
xlabel('log(mean of local stdev)');
ylabel('log(dimension)');
legend('bad','very blurry','blurry','good','great');

%% 
x = (all.lap+all.newlap)./2.*all.complex.*(all.local_std+all.newlocalstd)./2.*all.dim;
y = rand(length(x),1);%all.dim.*all.mean_intens.*all.local_std;

figure;
hold on; box on;
for i=1:5
    plot(x(all.flag==i),y(all.flag==i),'o','MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor',c(i,:));
end
%xlabel('mean of local stdev');
%ylabel('dimension');
legend('bad','very blurry','blurry','good','great');

% hold on; box on;
% plot(shit.local_std,shit.new.local_std,'r.');
% plot(v_blurry.local_std,v_blurry.new.local_std,'m.');
% plot(blurry.local_std,blurry.new.local_std,'c.');
% plot(good.local_std,good.new.local_std,'b.');
% plot(great.local_std,great.new.local_std,'g.');



%% PCA
pca_input = [all.dim,all.mean_intens,all.complex,all.lap,all.newlap,all.local_std,all.newlocalstd,all.local_std5,all.local_std7,all.newlocalstd5,all.newlocalstd7];
[coeff,score,latent] = pca(pca_input);

pc_1 = zeros(length(all.area),1);
pc_2 = zeros(length(all.area),1);
pc_3 = zeros(length(all.area),1);

for i=1:size(coeff,1)
    
    pc_1 = pc_1 + coeff(i,1) .* pca_input(:,i);
    pc_2 = pc_2 + coeff(i,2) .* pca_input(:,i);
    pc_3 = pc_3 + coeff(i,3) .* pca_input(:,i);
    
    %coeff(1,1) .* pca_input(:,1) + coeff(2,1) .* pca_input(:,2) + coeff(3,1) .* pca_input(:,3) + coeff(4,1) .* pca_input(:,4) + coeff(5,1) .* pca_input(:,5) + ...
    %    coeff(6,1) .* pca_input(:,6) + coeff(7,1) .* pca_input(:,7) + coeff(8,1) .* pca_input(:,8) + coeff(9,1) .* pca_input(:,9);

    %pc_2 = coeff(1,2) .* pca_input(:,1) + coeff(2,2) .* pca_input(:,2) + coeff(3,2) .* pca_input(:,3) + coeff(4,2) .* pca_input(:,4) + coeff(5,2) .* pca_input(:,5) + ...
    %    coeff(6,2) .* pca_input(:,6) + coeff(7,2) .* pca_input(:,7) + coeff(8,2) .* pca_input(:,8) + coeff(9,2) .* pca_input(:,9);
end


figure;
hold on; box on;
for i=1:5
    plot(pc_1(all.flag==i),pc_2(all.flag==i),'o','MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor',c(i,:));
end
xlabel('pc_1');
ylabel('pc_2');
legend('bad','very blurry','blurry','good','great');

%% HISTOGRAMS
% trial and error conclusions :
% - all.complex does not really improve the results
% - working on the PC is not giving good results
% - brigthness increases gap between bad and rest but can lead to some
% blurry detected as good

var1 = log((all.lap+all.newlap)./2.*all.complex.*(all.local_std+all.newlocalstd)./2.*all.dim);
bins_center = linspace(7.5,12.5,25) ;
M_count = []; M_center = [];
w = [1 0.9 0.8 0.7 0.6];
my_legend  = {'bad','much blurry','blurry','good','great'};

figure('Renderer','opengl');
hold on;
for i=1:5
   
    subplot(5,1,i);
    histogram(var1(all.flag==i),bins_center,'FaceColor',c(i,:),'FaceAlpha',1,'Normalization','probability'); 
    title(my_legend{i});
    if i < 5
        set(gca,'XTickLabel','');
    else
        xlabel('Descriptor \xi');
    end   
end


var1 = log((all.lap+all.newlap)./2.*all.complex.*(all.local_std+all.newlocalstd)./2.*all.dim.*all.mean_intens);
bins_center = linspace(5,12.5,25) ;
M_count = []; M_center = [];
w = [1 0.9 0.8 0.7 0.6];
my_legend  = {'bad','much blurry','blurry','good','great'};

figure('Renderer','opengl');
hold on;
for i=1:5
    
    subplot(5,1,i);
    histogram(var1(all.flag==i),bins_center,'FaceColor',c(i,:),'FaceAlpha',1,'Normalization','probability'); 
    title(my_legend{i});
    if i < 5
        set(gca,'XTickLabel','');
    else
        xlabel('Magic descriptor');
    end
end


%% EFFECT OF MASK SIZE ON LOCAL STD

% local STD
figure;
subplot(1,2,1);
hold on; box on;
for i=1:5
    plot(all.local_std(all.flag==i),all.local_std5(all.flag==i),'o','MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor',c(i,:));
end
maxtick = max(max(all.local_std(:)),max(all.local_std5(:)));
axis equal;
axis([0 maxtick 0 maxtick]);
xlabel('mask 3x3');
ylabel('mask5x5');

subplot(1,2,2);
hold on; box on;
for i=1:5
    plot(all.local_std(all.flag==i),all.local_std7(all.flag==i),'o','MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor',c(i,:));
end
maxtick = max(max(all.local_std(:)),max(all.local_std7(:)));
axis equal;
axis([0 maxtick 0 maxtick]);
xlabel('mask 3x3');
ylabel('mask 7x7');
%legend('bad','very blurry','blurry','good','great');

% new local STD
figure;
subplot(1,2,1);
hold on; box on;
for i=1:5
    plot(all.newlocalstd(all.flag==i),all.newlocalstd5(all.flag==i),'o','MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor',c(i,:));
end
maxtick = max(max(all.newlocalstd(:)),max(all.newlocalstd5(:)));
axis equal;
axis([0 maxtick 0 maxtick]);
xlabel('mask 3x3');
ylabel('mask5x5');

subplot(1,2,2);
hold on; box on;
for i=1:5
    plot(all.newlocalstd(all.flag==i),all.newlocalstd7(all.flag==i),'o','MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor',c(i,:));
end
maxtick = max(max(all.newlocalstd(:)),max(all.newlocalstd7(:)));
axis equal;
axis([0 maxtick 0 maxtick]);
xlabel('mask 3x3');
ylabel('mask 7x7');


%% PIECE OF CODE TO SAVE DATA FOR NNSTART

M_input = [all.dim, all.area, all.mean_intens, all.complex, all.lap, all.local_std, all.local_std5, all.local_std7, all.newlap, all.newlocalstd, all.newlocalstd5, all.newlocalstd7];
M_targets = zeros(length(all.flag),1);
M_targets(all.flag==1,1) = 1;
% M_targets(all.flag==2,2) = 1;
% M_targets(all.flag==3,2) = 1;
% M_targets(all.flag==4,2) = 1;
% M_targets(all.flag==5,2) = 1;


X = [all.dim, all.area, all.mean_intens, all.complex, all.lap, all.local_std, all.local_std5, all.local_std7, all.newlap, all.newlocalstd, all.newlocalstd5, all.newlocalstd7];
Y = cell(length(all.flag),1);
for i=1:length(Y)
    if all.flag(i) == 5
        Y{i} = 'great';
    else
        Y{i} = 'others';
    end
end




% % [g_count,g_center] = hist(good.focus,bins_center); g_countN = g_count./sum(g_count);
% % [ba_count,ba_center] = hist(bad.focus,bins_center); ba_countN = ba_count./sum(ba_count);
% % % [bl_count,bl_center] = hist(blurry.focus,bins_center); bl_countN = bl_count./sum(bl_count);
% % % [a_count,a_center] = hist(ambi.focus,bins_center); a_countN = a_count./sum(a_count);
% % % 
% % figure;
% % subplot(2,1,1); hold on; box on;
% bar(g_center,g_count,1,'g');


% %% Principal Component Analysis
% all.area_focus_ratio(all.area_focus_ratio == Inf) = 100;
% pca_input = [all.max_intens all.range_intens.*all.area./all.mean_intens all.area all.lap.*all.area./all.mean_intens all.new.wav.*all.area./all.mean_intens];
% all_variables = [all.max_intens all.range_intens.*all.area./all.mean_intens all.area all.max_intens.*all.lap.*all.area./all.mean_intens all.new.wav.*all.area./all.mean_intens];
% [coeff,score,latent] = pca(pca_input);
% 
% pc_1 = zeros(length(all.area),1);
% pc_2 = zeros(length(all.area),1);
% pc_3 = zeros(length(all.area),1);
% 
% for i=1:size(coeff,1)
%     
%     pc_1 = pc_1 + coeff(i,1) .* all_variables(:,i);
%     pc_2 = pc_2 + coeff(i,2) .* all_variables(:,i);
%     pc_3 = pc_3 + coeff(i,3) .* all_variables(:,i);
%     
%     %coeff(1,1) .* pca_input(:,1) + coeff(2,1) .* pca_input(:,2) + coeff(3,1) .* pca_input(:,3) + coeff(4,1) .* pca_input(:,4) + coeff(5,1) .* pca_input(:,5) + ...
%     %    coeff(6,1) .* pca_input(:,6) + coeff(7,1) .* pca_input(:,7) + coeff(8,1) .* pca_input(:,8) + coeff(9,1) .* pca_input(:,9);
% 
%     %pc_2 = coeff(1,2) .* pca_input(:,1) + coeff(2,2) .* pca_input(:,2) + coeff(3,2) .* pca_input(:,3) + coeff(4,2) .* pca_input(:,4) + coeff(5,2) .* pca_input(:,5) + ...
%     %    coeff(6,2) .* pca_input(:,6) + coeff(7,2) .* pca_input(:,7) + coeff(8,2) .* pca_input(:,8) + coeff(9,2) .* pca_input(:,9);
% end
% 
% pc_1_good = pc_1(all.flag == 1);
% pc_1_bad = pc_1(all.flag == 0);
% pc_2_good = pc_2(all.flag == 1);
% pc_2_bad = pc_2(all.flag == 0);
% pc_3_good = pc_3(all.flag ==1);















% % %% BAD detection
% clear file_list;
% file_list = dir(fullfile(main_dir,'bad','*.mat'));
% file_list = {file_list.name};
% file_list = fullfile(main_dir,'bad',file_list);
% 
% % features based on raw picture
% bad.area_focus_ratio = zeros(length(file_list),1);
% bad.area = zeros(length(file_list),1);
% bad.dim = zeros(length(file_list),1);
% bad.mean_intens = zeros(length(file_list),1);
% bad.max_intens = zeros(length(file_list),1);
% bad.range_intens = zeros(length(file_list),1);
% bad.focus = zeros(length(file_list),1);
% bad.area_focus = zeros(length(file_list),1);
% bad.area_range = zeros(length(file_list),1);
% bad.fallspeed = zeros(length(file_list),1);
% bad.lap = zeros(length(file_list),1);
% bad.area_lap = zeros(length(file_list),1);
% bad.wav = zeros(length(file_list),1);
% bad.area_wav = zeros(length(file_list),1);
% bad.complex = zeros(length(file_list),1);
% bad.range_complex = zeros(length(file_list),1);
% 
% % features based on enhanced picture
% bad.new.lap = zeros(length(file_list),1);
% bad.new.area_lap = zeros(length(file_list),1);
% bad.new.wav = zeros(length(file_list),1);
% bad.new.area_wav = zeros(length(file_list),1);
% bad.new.range_intens = zeros(length(file_list),1);
% bad.new.range_complex = zeros(length(file_list),1);
% 
% % loop over the snowflakes and save features in a vector
% for i=1:length(file_list)
%     
%     load(file_list{i});
%     bad.area_focus_ratio(i) = roi.area_focus_ratio;
%     bad.area(i) = roi.area;
%     bad.dim(i) = 0.5*(roi.width + roi.height);
%     bad.mean_intens(i) = roi.mean_intens;
%     bad.max_intens(i) = roi.max_intens;
%     bad.range_intens(i) = roi.range_intens;
%     bad.focus(i) = roi.focus;
%     bad.area_focus(i) = roi.area_focus; 
%     bad.area_range(i) = roi.area_range;
%     if ~isempty(roi.fallspeed)
%         bad.fallspeed(i) = roi.fallspeed;
%     else
%         bad.fallspeed(i) = NaN;
%     end
%     bad.lap(i) = roi.lap;
%     bad.area_lap(i) = roi.area_lap;
%     bad.wav(i) = roi.wav;
%     bad.area_wav(i) = roi.area_wav;
%     bad.complex(i) = roi.complex;
%     bad.range_complex(i) = roi.range_complex;
%     
%     bad.new.lap(i) = roi.new.lap;
%     bad.new.area_lap(i) = roi.new.area_lap;
%     bad.new.wav(i) = roi.new.wav;
%     bad.new.area_wav(i) = roi.new.area_wav;
%     bad.new.range_intens(i) = roi.new.range_intens;
%     bad.new.range_complex(i) = roi.new.range_complex;
%     
% end
% 
% %% struct containing all data
% all.area_focus_ratio = [good.area_focus_ratio; bad.area_focus_ratio];
% all.area = [good.area; bad.area];
% all.dim = [good.dim; bad.dim];
% all.mean_intens = [good.mean_intens; bad.mean_intens];
% all.max_intens = [good.max_intens; bad.max_intens];
% all.range_intens = [good.range_intens; bad.range_intens];
% all.focus = [good.focus; bad.focus];
% all.area_focus = [good.area_focus; bad.area_focus];
% all.area_range = [good.area_range; bad.area_range];
% all.fallspeed = [good.fallspeed; bad.fallspeed];
% all.lap = [good.lap; bad.lap];
% all.area_lap = [good.area_lap; bad.area_lap];
% all.wav = [good.wav; bad.wav];
% all.area_wav = [good.area_wav; bad.area_wav];
% all.complex = [good.complex; bad.complex];
% all.range_complex = [good.range_complex; bad.range_complex];
% all.new.lap = [good.new.lap; bad.new.lap];
% all.new.area_lap = [good.new.area_lap; bad.new.area_lap];
% all.new.wav = [good.new.wav; bad.new.wav];
% all.new.area_wav = [good.new.area_wav; bad.new.area_wav];
% all.new.range_intens = [good.new.range_intens; bad.new.range_intens];
% all.new.range_complex = [good.new.range_complex; bad.new.range_complex];
% all.flag = zeros(length(all.area),1);
% all.flag(1:length(good.area)) = 1;
% 
% %% Principal Component Analysis
% all.area_focus_ratio(all.area_focus_ratio == Inf) = 100;
% pca_input = [all.max_intens all.range_intens.*all.area./all.mean_intens all.area all.lap.*all.area./all.mean_intens all.new.wav.*all.area./all.mean_intens];
% all_variables = [all.max_intens all.range_intens.*all.area./all.mean_intens all.area all.max_intens.*all.lap.*all.area./all.mean_intens all.new.wav.*all.area./all.mean_intens];
% [coeff,score,latent] = pca(pca_input);
% 
% pc_1 = zeros(length(all.area),1);
% pc_2 = zeros(length(all.area),1);
% pc_3 = zeros(length(all.area),1);
% 
% for i=1:size(coeff,1)
%     
%     pc_1 = pc_1 + coeff(i,1) .* all_variables(:,i);
%     pc_2 = pc_2 + coeff(i,2) .* all_variables(:,i);
%     pc_3 = pc_3 + coeff(i,3) .* all_variables(:,i);
%     
%     %coeff(1,1) .* pca_input(:,1) + coeff(2,1) .* pca_input(:,2) + coeff(3,1) .* pca_input(:,3) + coeff(4,1) .* pca_input(:,4) + coeff(5,1) .* pca_input(:,5) + ...
%     %    coeff(6,1) .* pca_input(:,6) + coeff(7,1) .* pca_input(:,7) + coeff(8,1) .* pca_input(:,8) + coeff(9,1) .* pca_input(:,9);
% 
%     %pc_2 = coeff(1,2) .* pca_input(:,1) + coeff(2,2) .* pca_input(:,2) + coeff(3,2) .* pca_input(:,3) + coeff(4,2) .* pca_input(:,4) + coeff(5,2) .* pca_input(:,5) + ...
%     %    coeff(6,2) .* pca_input(:,6) + coeff(7,2) .* pca_input(:,7) + coeff(8,2) .* pca_input(:,8) + coeff(9,2) .* pca_input(:,9);
% end
% 
% pc_1_good = pc_1(all.flag == 1);
% pc_1_bad = pc_1(all.flag == 0);
% pc_2_good = pc_2(all.flag == 1);
% pc_2_bad = pc_2(all.flag == 0);
% pc_3_good = pc_3(all.flag ==1);
% pc_3_bad = pc_3(all.flag == 0);
% 
% figure; hold on;
% plot3(pc_1_good,pc_2_good,pc_3_good,'g.');
% plot3(pc_1_bad,pc_2_bad,pc_3_bad,'r.');
% 
% % nice separation
% % [all.mean_intens all.max_intens all.range_intens all.focus all.area_focus all.area_range all.range_complex all.area all.lap all.wav all.complex all.new.lap all.new.wav all.new.range_intens];
% 
% 
% % %% FOCUS
% % bins_center = (0:0.01:0.35);
% % % 
% % [g_count,g_center] = hist(good.focus,bins_center); g_countN = g_count./sum(g_count);
% % [ba_count,ba_center] = hist(bad.focus,bins_center); ba_countN = ba_count./sum(ba_count);
% % % [bl_count,bl_center] = hist(blurry.focus,bins_center); bl_countN = bl_count./sum(bl_count);
% % % [a_count,a_center] = hist(ambi.focus,bins_center); a_countN = a_count./sum(a_count);
% % % 
% % figure;
% % subplot(2,1,1); hold on; box on;
% bar(g_center,g_count,1,'g');
% bar(ba_center,ba_count,0.8,'r');
% % bar(bl_center,bl_count,0.6,'b');
% % bar(a_center,a_count,0.4,'y');
% v = axis;
% axis([bins_center(1) bins_center(end) v(3) v(4)]);
% ylabel('N');
% title('FOCUS');
% subplot(2,1,2); hold on; box on;
% plot(g_center,g_countN,'g-');
% plot(ba_center,ba_countN,'r-');
% % plot(bl_center,bl_countN,'b-');
% % plot(a_center,a_countN,'y');
% legend('good','bad');
% ylabel('N_{norm}');
% xlabel('bins');
% 
% %% LAP
% 
% [g_count,g_center] = hist(good.new.area_lap); g_countN = g_count./sum(g_count);
% [ba_count,ba_center] = hist(bad.new.area_lap); ba_countN = ba_count./sum(ba_count);
% % [bl_count,bl_center] = hist(blurry.area_focus,bins_center); bl_countN = bl_count./sum(bl_count);
% % [a_count,a_center] = hist(ambi.area_focus,bins_center); a_countN = a_count./sum(a_count);
% % 
% figure;
% subplot(2,1,1); hold on; box on;
% bar(g_center,g_count,1,'g');
% bar(ba_center,ba_count,0.8,'r');
% % bar(bl_center,bl_count,0.6,'b');
% % bar(a_center,a_count,0.4,'y');
% % v = axis;
% % axis([bins_center(1) 2000 v(3) v(4)]);
% % ylabel('N');
% % title('AREA FOCUS');
% subplot(2,1,2); hold on; box on;
% plot(g_center,g_countN,'g-');
% plot(ba_center,ba_countN,'r-');
% % plot(bl_center,bl_countN,'b-');
% % plot(a_center,a_countN,'y');
% % legend('good','bad','blurry','ambig.');
% % ylabel('N_{norm}');
% % xlabel('bins');
% % v = axis;
% % axis([bins_center(1) bins_center(end) v(3) v(4)]);
% 
% %% MEAN_INTENS
% bins_center = linspace(0,1,35);
% 
% [g_count,g_center] = hist(good.mean_intens,bins_center); g_countN = g_count./sum(g_count);
% [ba_count,ba_center] = hist(bad.mean_intens,bins_center); ba_countN = ba_count./sum(ba_count);
% [bl_count,bl_center] = hist(blurry.mean_intens,bins_center); bl_countN = bl_count./sum(bl_count);
% [a_count,a_center] = hist(ambi.mean_intens,bins_center); a_countN = a_count./sum(a_count);
% 
% figure;
% subplot(2,1,1); hold on; box on;
% bar(g_center,g_count,1,'g');
% bar(ba_center,ba_count,0.8,'r');
% bar(bl_center,bl_count,0.6,'b');
% bar(a_center,a_count,0.4,'y');
% v = axis;
% axis([bins_center(1) bins_center(end) v(3) v(4)]);
% ylabel('N');
% title('MEAN BRIGHTNESS');
% subplot(2,1,2); hold on; box on;
% plot(g_center,g_countN,'g-');
% plot(ba_center,ba_countN,'r-');
% plot(bl_center,bl_countN,'b-');
% plot(a_center,a_countN,'y');
% legend('good','bad','blurry','ambig.');
% ylabel('N_{norm}');
% xlabel('bins');
% v = axis;
% axis([bins_center(1) bins_center(end) v(3) v(4)]);
% 
% %% MAX_INTENS 
% bins_center = linspace(0,1,35);
% 
% [g_count,g_center] = hist(good.max_intens,bins_center); g_countN = g_count./sum(g_count);
% [ba_count,ba_center] = hist(bad.max_intens,bins_center); ba_countN = ba_count./sum(ba_count);
% [bl_count,bl_center] = hist(blurry.max_intens,bins_center); bl_countN = bl_count./sum(bl_count);
% [a_count,a_center] = hist(ambi.max_intens,bins_center); a_countN = a_count./sum(a_count);
% 
% figure;
% subplot(2,1,1); hold on; box on;
% bar(g_center,g_count,1,'g');
% bar(ba_center,ba_count,0.8,'r');
% bar(bl_center,bl_count,0.6,'b');
% bar(a_center,a_count,0.4,'y');
% v = axis;
% axis([bins_center(1) bins_center(end) v(3) v(4)]);
% ylabel('N');
% title('MAX BRIGHTNESS');
% subplot(2,1,2); hold on; box on;
% plot(g_center,g_countN,'g-');
% plot(ba_center,ba_countN,'r-');
% plot(bl_center,bl_countN,'b-');
% plot(a_center,a_countN,'y');
% legend('good','bad','blurry','ambig.');
% ylabel('N_{norm}');
% xlabel('bins');
% v = axis;
% axis([bins_center(1) bins_center(end) v(3) v(4)]);
% 
% %% RANGE_INTENS
% bins_center = linspace(0,1,35);
% 
% [g_count,g_center] = hist(good.range_intens,bins_center); g_countN = g_count./sum(g_count);
% [ba_count,ba_center] = hist(bad.range_intens,bins_center); ba_countN = ba_count./sum(ba_count);
% [bl_count,bl_center] = hist(blurry.range_intens,bins_center); bl_countN = bl_count./sum(bl_count);
% [a_count,a_center] = hist(ambi.range_intens,bins_center); a_countN = a_count./sum(a_count);
% 
% figure;
% subplot(2,1,1); hold on; box on;
% bar(g_center,g_count,1,'g');
% bar(ba_center,ba_count,0.8,'r');
% bar(bl_center,bl_count,0.6,'b');
% bar(a_center,a_count,0.4,'y');
% v = axis;
% axis([bins_center(1) bins_center(end) v(3) v(4)]);
% ylabel('N');
% title('RANGE INTENS');
% subplot(2,1,2); hold on; box on;
% plot(g_center,g_countN,'g-');
% plot(ba_center,ba_countN,'r-');
% plot(bl_center,bl_countN,'b-');
% plot(a_center,a_countN,'y');
% legend('good','bad','blurry','ambig.');
% ylabel('N_{norm}');
% xlabel('bins');
% v = axis;
% axis([bins_center(1) bins_center(end) v(3) v(4)]);
% 
% %% AREA
% bins_center = linspace(0,max(good.area),50);
% 
% [g_count,g_center] = hist(good.area,bins_center); g_countN = g_count./sum(g_count);
% [ba_count,ba_center] = hist(bad.area,bins_center); ba_countN = ba_count./sum(ba_count);
% [bl_count,bl_center] = hist(blurry.area,bins_center); bl_countN = bl_count./sum(bl_count);
% [a_count,a_center] = hist(ambi.area,bins_center); a_countN = a_count./sum(a_count);
% 
% figure;
% subplot(2,1,1); hold on; box on;
% bar(g_center,g_count,1,'g');
% bar(ba_center,ba_count,0.8,'r');
% bar(bl_center,bl_count,0.6,'b');
% bar(a_center,a_count,0.4,'y');
% v = axis;
% axis([bins_center(1) bins_center(end) v(3) v(4)]);
% ylabel('N');
% title('AREA');
% subplot(2,1,2); hold on; box on;
% plot(g_center,g_countN,'g-');
% plot(ba_center,ba_countN,'r-');
% plot(bl_center,bl_countN,'b-');
% plot(a_center,a_countN,'y');
% legend('good','bad','blurry','ambig.');
% ylabel('N_{norm}');
% xlabel('bins');
% v = axis;
% axis([bins_center(1) bins_center(end) v(3) v(4)]);
% 
% %% AREA_FOCUS_RATIO
% bins_center = linspace(0,30,50);
% 
% [g_count,g_center] = hist(good.area_focus_ratio,bins_center); g_countN = g_count./sum(g_count);
% [ba_count,ba_center] = hist(bad.area_focus_ratio,bins_center); ba_countN = ba_count./sum(ba_count);
% [bl_count,bl_center] = hist(blurry.area_focus_ratio,bins_center); bl_countN = bl_count./sum(bl_count);
% [a_count,a_center] = hist(ambi.area_focus_ratio,bins_center); a_countN = a_count./sum(a_count);
% 
% figure;
% subplot(2,1,1); hold on; box on;
% bar(g_center,g_count,1,'g');
% bar(ba_center,ba_count,0.8,'r');
% bar(bl_center,bl_count,0.6,'b');
% bar(a_center,a_count,0.4,'y');
% v = axis;
% axis([bins_center(1) bins_center(end) v(3) 50]);
% ylabel('N');
% title('AREA FOCUS RATIO');
% subplot(2,1,2); hold on; box on;
% plot(g_center,g_countN,'g-');
% plot(ba_center,ba_countN,'r-');
% plot(bl_center,bl_countN,'b-');
% plot(a_center,a_countN,'y');
% legend('good','bad','blurry','ambig.');
% ylabel('N_{norm}');
% xlabel('bins');
% v = axis;
% axis([bins_center(1) bins_center(end) v(3) 0.20]);
% 
% %% ENTROPY
% bins_center = linspace(0,10,50);
% 
% [g_count,g_center] = hist(good.entropy,bins_center); g_countN = g_count./sum(g_count);
% [ba_count,ba_center] = hist(bad.entropy,bins_center); ba_countN = ba_count./sum(ba_count);
% [bl_count,bl_center] = hist(blurry.entropy,bins_center); bl_countN = bl_count./sum(bl_count);
% [a_count,a_center] = hist(ambi.entropy,bins_center); a_countN = a_count./sum(a_count);
% 
% figure;
% subplot(2,1,1); hold on; box on;
% bar(g_center,g_count,1,'g');
% bar(ba_center,ba_count,0.8,'r');
% bar(bl_center,bl_count,0.6,'b');
% bar(a_center,a_count,0.4,'y');
% v = axis;
% axis([bins_center(1) bins_center(end) v(3) 50]);
% ylabel('N');
% title('ENTROPY');
% subplot(2,1,2); hold on; box on;
% plot(g_center,g_countN,'g-');
% plot(ba_center,ba_countN,'r-');
% plot(bl_center,bl_countN,'b-');
% plot(a_center,a_countN,'y');
% legend('good','bad','blurry','ambig.');
% ylabel('N_{norm}');
% xlabel('bins');
% v = axis;
% axis([bins_center(1) bins_center(end) v(3) 0.20]);
% 
% 
