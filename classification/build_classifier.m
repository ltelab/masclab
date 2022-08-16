clearvars ; close all;

                % alpha = 0.0001, lambda = 0.01, mini-batch size = 100
                % momentum = 0.9, Max no of iterations = 5000
                %parameters={0.00001,0.01,100,0.9,5000}; 

% USER DEFINED PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

normalization_type = 'standardization';
method='logistic'; % can be either 'logistic', 'naive_bayes,'nn','svm','random_forest' or 'baseline'
parameters_method =  {0.0001,1,1000,0,10000};
%parameters_method = {1,1,'linear'};
% best for svm before Hu : parameters_method =  {10000,0.0001,'rbf'};
% default for logistic : {0.001,0,779,0.9,20000}; 
% default for svm-binary (6 classes) : {100,0.1,'rbf'} or {1,1,'linear'}
type_classif='multiclass'; % can be binary or multiclass
save_log=false; % Can be true of false

features_ranking = false;
display_confmat = true;
label_version = '1.1';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% path to training data
dir_data = '/home/praz/Documents/MASC/MASC_labelling/mixed_sample_merged_12345';
%dir_data = '/media/praz/MyData/Massimo_data/ALL';

data_filenames = dir(fullfile(dir_data,'20*.mat'));
% for Massimo data
if isempty(data_filenames)
    data_filenames = dir(fullfile(dir_data,'ICE*.mat'));
end
data_filenames = {data_filenames.name}';
data_picnames = dir(fullfile(dir_data,'*.png'));
data_picnames = {data_picnames.name}';

% time interval to take into consideration
t_str_start = '20150101000000';
t_str_stop  = '20180101000000';

% path to labels
dir_labels = 'home/praz/Documents/MASC/masclab/labelling_GUI/labels_v.1.1';
labels_filename = 'FINAL3.mat';
%dir_labels = '/media/praz/MyData/Massimo_data';
%labels_filename = 'CP_labels_all_sorted.mat';

%% load flakes descriptors

% list of descriptors wanted
% initial trial & error
%feat_vec = [1 2 3 4 5 9 10 16 19 25 28 29 34 37 49 50 53 54 55 63 64 65 70 75 83 84 87];
% 20 features resulting from forward selection - 1it, norand
%feat_vec = [5 93 87 59 54 35 25 31 28 72 75 84 83 9 82 39 49 92 1 32];
%feat_vec = [1 2 3 4 5 9 10 16 19 25 28 29 34 37 49 50 53 54 55 59 63 64 65 70 75 83 84 87];
%feat_vec = [1 2 3 4 5 9 10 16 19 25 49 50 59 63 64 65 70 83 84 87]; % only geometry
%feat_vec = [28,29,34,37,53,54,55,75]; % only texture
%feat_vec = 1:1:120;
%feat_vec(:,101) = [];

Dwanted = 25;
F = load('out_geo_10it_rand_FS_CV4_HSS_D80.mat');
feat_vec = F.feat_mat(:,Dwanted+1);
feat_vec = feat_vec(feat_vec~=0);
feat_vec = sort(feat_vec,'ascend');

% compute new descriptors (if necessary)
process_new_descriptors(dir_data,t_str_start,t_str_stop)

% load the training matrix X
[X,Xlab,Xname,Xt] = load_processed_data(dir_data,t_str_start,t_str_stop,feat_vec);

% load the labels vector y
[y,yR,yM] = load_labels(dir_data,dir_labels,labels_filename,label_version);

% remove "undetermined", "not labelled yet" or "ambiguous" labels
idx_unknown = find(y<=0);
y(idx_unknown) = [];
yR(idx_unknown) = [];
X(idx_unknown,:) = [];
data_picnames(idx_unknown) = [];
data_filenames(idx_unknown) = [];
fprintf('%u samples discarded because labelled as "undetermined" or "ambiguous" \n',length(idx_unknown));


%% merging the classes into a simpler scheme
if 1
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
    labels = {'smal','col','dend','agg','graup','capcol'};
end

% MASC
%labels = {'graup','agg','melt','smal','dend','col'};
% Massimo
%labels = {'agg','smal','col','bulros'};
% y(y==2) = 1;
% y(y==4) = 2;
% y(y==7) = 3;
% y(y==8) = 4;

N = size(X,1);
D = size(X,2);
N_classes = length(unique(y));

% computing the weights, if we want to penalize the cost
for i=1:N_classes
    Nc(i,1) = sum(y==i);
end
[Ncmax,Ncmaxidx] = max(Nc);
for i=1:N_classes
    weights(i,1) = Ncmax/Nc(i);
end
parameters_method{end+1} = weights;


%% Features transformation
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
  
% Initialize global confusion matrix
global_confmat_Te = zeros(N_classes);
global_confmat_Tr = zeros(N_classes);

% Initialize beta ranking vectors
beta_sum_mat = [];

% Building classifier: XTr = Xtot
XTr = X;
yTr = y;

% Fit model using the wrapping function
model = trainClassifier(method, type_classif, yTr, XTr, parameters_method);
        
% Features ranking
if features_ranking
    betas = abs(model);
    beta_sum = sum(betas,2);
    beta_sum = beta_sum(2:end,:);
    beta_sum_mat(:,end+1) = beta_sum;
end
               
% Predict
[predTr,scoresTr] = predictClassifier(method,type_classif,model,XTr);  
        
% Compute accuracy
BER_Tr = computeBER(yTr,predTr);
OA_Tr = computeOA(yTr,predTr);
kappa_Tr = computeKAPPA(yTr,predTr);

fprintf('Training BER: %.2f%%\n\n', BER_Tr * 100 );
fprintf('Training OA: %.2f%%\n\n', OA_Tr * 100 );
fprintf('Training kappa: %.2f%%\n\n', kappa_Tr * 100 );

% Confusion matrix
if display_confmat
    global_confmat = confusionmat(yTr,predTr);   
    global_confmat = round(global_confmat./sum(global_confmat(:)) * 100,2);
    figure;
    heatmap(global_confmat,labels,labels,1,'Colormap','red','ShowAllTicks',1,'UseLogColorMap',true,'Colorbar',true);
    title('Train confusion matrix');
end



classifier.type = method;
classifier.type_classif = type_classif;
classifier.parameters_method = parameters_method;
classifier.normalization = 'standardization';
classifier.model = model;
classifier.N = N;
classifier.N_classes = N_classes;
classifier.N_labels = labels;
classifier.label_type = labels_filename;
classifier.D = D;
classifier.Xlab = Xlab;
classifier.normalization_params = classif_params;
classifier.BER = BER_Tr;
classifier.OA = OA_Tr;
classifier.kappa = kappa_Tr;
classifier.feat_vec = feat_vec;



