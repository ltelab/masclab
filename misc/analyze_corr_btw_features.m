%% Load data
% script to analyze the correlation between features
clear all; %close all;

dir_data = '/home/praz/Documents/MASC/MASC_labelling/mixed_sample_merged_1234';
data_filenames = dir(fullfile(dir_data,'20*.mat'));
data_filenames = {data_filenames.name}';
data_picnames = dir(fullfile(dir_data,'*.png'));
data_picnames = {data_picnames.name}';

% time interval to take into consideration
t_str_start = '20150101000000';
t_str_stop  = '20170101000000';

% path to labels
dir_labels = 'home/praz/Documents/MASC/masclab/labelling_GUI/labels_v.1.1';
labels_filename = 'FINAL2_corrected.mat';

% list of descriptors wanted
feat_vec = [1:1:96]';
%feat_vec(101) = [];
%feat_vec(41:48) = [];

% compute new descriptors (if necessary)
process_new_descriptors(dir_data,t_str_start,t_str_stop)

% load the training matrix X
[X,Xlab,Xname,Xt] = load_processed_data(dir_data,t_str_start,t_str_stop,feat_vec);


%% Features transformation

D = size(X,2);
normalization_type = 'standardization';

for i=1:D
      
    skew = skewness(X(:,i));
    % save the skewness in model
    classif_params.skew(i,1) = skew;
    
    if skew > 1
        
        X(:,i) = log(abs(X(:,i)+1));
        
    elseif skew > 0.75
        
        X(:,i) = sqrt(abs(X(:,i)));
        
    elseif skew < -1
        
        X(:,i) = exp(X(:,i));
        
    elseif skew < -0.75
        
        X(:,i) = X(:,i).^2;
        
    end

    if strcmp(normalization_type,'standardization') 
        
        tmp_mean = mean(X(:,i));
        tmp_std  = std(X(:,i));
        if i<1000000000
            X(:,i) = (X(:,i) - tmp_mean)/tmp_std;
        else
            X(:,i) = (X(:,i) - tmp_mean)/(5*tmp_std);
        end   
        
        classif_params.mean(i,1) = tmp_mean;
        classif_params.std(i,1) = tmp_std;
        
    elseif strcmp(normalization_type,'rescaling')
        
        tmp_min = min(X(:,i));
        tmp_max = max(X(:,i));
        X(:,i) = (X(:,i) - tmp_min)/(tmp_max - tmp_min);
        
    end
    
end

%% stuff

% remove cols 41-48
X(:,41:48) = NaN;
%X(:,[7 17 18 92 2 14 34 37 38 39 53 80 51 26 70 77 15 33 79 6]) = NaN;
C = corrcoef(X);
C = triu(C,1);
corr_cap = 0.96;

for i=1:size(C,1)
    for j=1:size(C,2)
        if C(i,j) >= corr_cap
            fprintf('Feature %u (%s) is very correlated with feature %u (%s) \n',i,Xlab{i},j,Xlab{j});
        end
    end
end



















% %% Features discretization
% Nbins = 200;
% N = size(X,1);
% for i=1:D
%    
%     mini = min(X(:,i));
%     maxi = max(X(:,i));
%     grid = linspace(mini,maxi,Nbins);
%     
%     for j=1:N
%         idx = find(grid>=X(j,i),1,'first');
%         X(j,i) = grid(idx)-(maxi-mini)/(Nbins*2);      
%     end
%     
% end
% 
% 
% 
% %% Do stuff
% 
% MDL = [];
% 
% for k=1:D
% 
% feat = X(:,k);
% class = X(:,1);
% 
% 
% N = size(X,1);
% D = size(X,2);
% C = length(unique(class));
% classes = unique(class);
% nc = [];
% 
% for i=1:C
%     nc(i) = sum(class==classes(i));
% end
%     
% % Prior MDL
% p1 = log_facto(N);
% p2 = 0;
% for i=1:C
%     p2 = p2 + log_facto(nc(i));
% end
% p3 = log_facto(N+C-1);
% p4 = log_facto(N);
% p5 = log_facto(C-1);
% prior_MDL = p1-p2+p3-p4-p5;
% 
% % Post MDL
% p1 = 0;
% p2 = 0;
% p3 = 0;
% p4 = 0;
% p5 = 0;
% 
% ufeat = unique(feat);
% 
% for j=1:length(ufeat)
% 
%     nj = sum(feat==ufeat(j));
%     p1 = p1 + log_facto(nj);
%     for i=1:C
%         ncj(i) = sum(feat==ufeat(j) & class==classes(i));
%         p2 = p2 + log_facto(ncj(i));
%     end
%     p3 = p3 + log_facto(nj+C-1);
%     p4 = p4 + log_facto(nj);
%     p5 = p5 + log_facto(C-1);
%        
% end
% post_MDL = p1-p2+p3-p4-p5;
% 
% MDL(k) = (prior_MDL - post_MDL)/N;
% %MDL(k) = MDL(k)/(prior_MDL/N);
% 
% end








%[RANKED,WEIGHT] = relieff(X(:,1:10),y,1,'method','classification');



% MDL analysis








% 
% add normalization (or not)
% 
% linear correlation
% C = corrcoef(X);
% figure;
% heatmap(C,1:1:D,1:1:D,1,'Colormap','red');
% 
% 
% 
% function a=compute_corr_with_class(X,y,'method')
% 
%     N = size(X,1);
%     D = size(X,2);
%     classes = unique(y);
%     N_classes = length(classes);
%     centroids = [];
%     all_dist = [];
%     score = [];
%     
%     for i=1:D
%         
%         feat = X(:,i);
%         
%         for j=1:N_classes
%             centroids(j,1) = median(feat(y==classes(j)));
%         end
%         
%         all_dist = pdist2(feat,centroids);
%         [~,pred_class] = min(all_dist,[],2);
%         pred_class = classes(pred_class);
%         score(i) = sum(pred_class==y)/N;      
%         
%     end
%         
%     compute average features correlation
%     C = corrcoef(X);
%     C = triu(C,1);
%     C = C(:);
%     C = C(C~=0);
% 
%     merit = D*mean(score)/ sqrt(D + D*(D-1)*mean(C))
% 
% end
% 











