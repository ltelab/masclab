%% User parameters & initialisation
% User parameters
K=4; 
N_it=10; 
method='logistic';
type_classif='multiclass'; 
use_cost_weights = true;
target = 'riming';
normalization_type = 'standardization';
dynamic_feat_transfo = false;

if strcmp(target,'geometry')
    parameters_method =  {0.0001,1,1000,0,10000};
elseif strcmp(target,'riming')
    parameters_method = {0.0001,0.01,1000,0,10000};
elseif strcmp(target,'melting')
    parameters_method = {0.001,0.01,1000,0,10000};
end

% Path to training data
dir_data = '/home/praz/Documents/MASC/MASC_labelling/mixed_sample_merged_12345';
data_filenames = dir(fullfile(dir_data,'20*.mat'));
data_filenames = {data_filenames.name}';
data_picnames = dir(fullfile(dir_data,'*.png'));
data_picnames = {data_picnames.name}';

% Path to labels
label_version = '1.1';
dir_labels = 'home/praz/Documents/MASC/masclab/labelling_GUI/labels_v.1.1';
labels_filename = 'FINAL3.mat';

% Time interval 
t_str_start = '20150101000000';
t_str_stop  = '20180101000000';

% Selection of descriptors
Dwanted = 25;
if strcmp(target,'geometry')
    load('out_geo_10it_rand_FS_CV4_HSS_D80.mat');
elseif strcmp(target,'riming')
    load('out_riming_10it_rand_FS_CV4_RMSE_D80.mat');
elseif strcmp(target,'melting')
    load('out_melting_10it_rand_FS_CV4_HSS_80D.mat');
end
    
feat_vec = feat_mat(:,Dwanted+1);
feat_vec = feat_vec(feat_vec~=0);
feat_vec = sort(feat_vec,'ascend');

% simple code to remove textural descriptors (HM type)
if 1
    feat_text = [28 29 30 36 39 40 51 53 54 56 75 78];
    for i=1:length(feat_vec)
        if ~isempty(find(feat_text==feat_vec(i)))
            feat_vec(i) = -1;
        end
    end
    feat_vec(feat_vec==-1) = [];
end

%% Data loading
% Compute new descriptors (if necessary)
process_new_descriptors(dir_data,t_str_start,t_str_stop)

% Load the training matrix X
[X,Xlab,Xname,Xt] = load_processed_data(dir_data,t_str_start,t_str_stop,feat_vec);

% Load the labels vector y
[y,yR,yM] = load_labels(dir_data,dir_labels,labels_filename,label_version);

if strcmp(target,'riming')
    y = yR;
elseif strcmp(target,'melting')
    y = yM;
end

% Remove "undetermined", "not labelled yet" or "ambiguous" labels
if strcmp(target,'melting')
    idx_unknown = find(y<0);
else
    idx_unknown = find(y<=0);
end
y(idx_unknown) = [];
yR(idx_unknown) = [];
X(idx_unknown,:) = [];
data_picnames(idx_unknown) = [];
data_filenames(idx_unknown) = [];
fprintf('%u samples discarded because labelled as "undetermined" or "ambiguous" \n',length(idx_unknown));

%% merging the classes into a simpler scheme
if strcmp(target,'geometry')
    fprintf('Merging the classes for "geometry" classification. \n');
    % 1: small particle
    % 2: column & needle
    y(y==3) = 2;
    % 3: sectored plates & dendrites
    y(y==5 | y==6) = 3;
    % 4: remove plates, add aggregates
    idx_plates = find(y==4);
    y(idx_plates) = [];
    yR(idx_plates) = [];
    X(idx_plates,:) = [];
    data_picnames(idx_plates) = [];
    data_filenames(idx_plates) = [];
    y(y==7) = 4;
    % 5: graupels
    y(y==11) = 5;
    % 6: capped columns
    y(y==10) = 6;
    % for now : discard combin of crystals
    idx_combcrystals = find(y>6);
    y(idx_combcrystals) = [];
    yR(idx_combcrystals) = [];
    X(idx_combcrystals,:) = [];
    data_picnames(idx_combcrystals) = [];
    data_filenames(idx_combcrystals) = [];
    % add graupel-like particles into graupel
    %y(y==4 & (yR==4 | yR==5)) = 5;
    labels = {'SP','CC','PC','AG','GR','CPC'};
end

%% Features transformation
if ~dynamic_feat_transfo

D = size(X,2);

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

end

%% Classification 
N = size(X,1);
D = size(X,2);
N_classes = length(unique(y));

Cout = CrossValidation(method,type_classif,parameters_method,K,N_it,X,y,0,1,1,use_cost_weights,dynamic_feat_transfo,feat_vec);

    