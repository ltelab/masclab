clearvars; close all;

                % alpha = 0.0001, lambda = 0.01, mini-batch size = 100
                % momentum = 0.9, Max no of iterations = 5000
                %parameters={0.00001,0.01,100,0.9,5000}; 

% USER DEFINED PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K=4; % Number of subsamples in the K-fold cross-validation
N_it=20; % Number of iterations of K-fold cross-validation
normalization_type = 'standardization';
method='logistic'; % can be either 'logistic', 'naive_bayes,'nn','svm','random_forest' or 'baseline'
parameters_method =  {0.001,0.1,1000,0,10000};
% best for svm before Hu : parameters_method =  {10000,0.0001,'rbf'};
% default for logistic : {0.001,0,779,0.9,20000}; 
% default for svm-binary (6 classes) : {100,0.1,'rbf'} or {1,1,'linear'}
type_classif='binary'; % can be binary or multiclass
save_log=false; % Can be true of false

features_ranking = false;
display_confmat = true;
label_version = '1.1';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% path to training data
dir_data = '/home/praz/Documents/MASC/MASC_labelling/mixed_sample_merged_1234';
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
t_str_stop  = '20170101000000';

% path to labels
dir_labels = 'home/praz/Documents/MASC/masclab/labelling_GUI/labels_v.1.1';
labels_filename = 'FINAL2_corrected.mat';
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
load('out_geo_10it_rand_FS_CV4_HSS_D80.mat');
feat_vec = feat_mat(:,Dwanted+1);
feat_vec = feat_vec(feat_vec~=0);
feat_vec = sort(feat_vec,'ascend');

% compute new descriptors (if necessary)
process_new_descriptors(dir_data,t_str_start,t_str_stop)

% load the training matrix X
[X,Xlab,Xname,Xt] = load_processed_data(dir_data,t_str_start,t_str_stop,feat_vec);

% load the labels vector y
[y,yR,yM] = load_labels(dir_data,dir_labels,labels_filename,label_version);
%y = yM;

% remove "undetermined", "not labelled yet" or "ambiguous" labels
idx_unknown = find(y<=0);
y(idx_unknown) = [];
yR(idx_unknown) = [];
X(idx_unknown,:) = [];
data_picnames(idx_unknown) = [];
data_filenames(idx_unknown) = [];
fprintf('%u samples discarded because labelled as "undetermined" or "ambiguous" \n',length(idx_unknown));

y = real(y);

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
    Nc(i,1) = sum(y==(i-1));
end
[Ncmax,Ncmaxidx] = max(Nc);
for i=1:N_classes
    weights(i,1) = Ncmax/Nc(i);
end
parameters_method{end+1} = weights;

%% Features transformation

D = size(X,2);

% for i=1:D
%       
%     skew = skewness(X(:,i));
%     % save the skewness in model
%     classif_params.skew(i,1) = skew;
%     
%     if skew > 1
%         
%         X(:,i) = log(abs(X(:,i)+1));
%         
%     elseif skew > 0.75
%         
%         X(:,i) = sqrt(abs(X(:,i)));
%         
%     elseif skew < -1
%         
%         X(:,i) = exp(X(:,i));
%         
%     elseif skew < -0.75
%         
%         X(:,i) = X(:,i).^2;
%         
%     end
% 
%     if strcmp(normalization_type,'standardization') 
%         
%         tmp_mean = mean(X(:,i));
%         tmp_std  = std(X(:,i));
%         if i<1000000000
%             X(:,i) = (X(:,i) - tmp_mean)/tmp_std;
%         else
%             X(:,i) = (X(:,i) - tmp_mean)/(5*tmp_std);
%         end   
%         
%         classif_params.mean(i,1) = tmp_mean;
%         classif_params.std(i,1) = tmp_std;
%         
%     elseif strcmp(normalization_type,'rescaling')
%         
%         tmp_min = min(X(:,i));
%         tmp_max = max(X(:,i));
%         X(:,i) = (X(:,i) - tmp_min)/(tmp_max - tmp_min);
%         
%     end
%     
% end


%% Algorithm training & testing

% split between train and test
setSeed(1);
idx_perm = randperm(N);
Nk = floor(N/K);
idxTe = idx_perm(1:Nk);
idxTr = idx_perm(Nk+1:end);
yTe = y(idxTe);
XTe_ini = X(idxTe,:);
yTr_tot = y(idxTr);
XTr_tot = X(idxTr,:);
%train_sampling = [0.05 0.10 0.15 0.20 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95];
train_sampling = [0.02:0.02:0.98];
%[0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.10 0.15 0.20 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];

for i=1:N_it
    
    fprintf('\n Iteration %u/%u : \n\n',i,N_it);

    idx_perm = randperm(length(yTr_tot));
    for k=1:length(train_sampling)

        idx2keep =idx_perm(1:floor(train_sampling(k)*length(yTr_tot)));
        yTr = yTr_tot(idx2keep);
        XTr = XTr_tot(idx2keep,:);
        
        XTe = XTe_ini;
        
        % Normalization based on training set only
        for d=1:D

            skew = skewness(XTr(:,d));

            if skew > 1

                XTr(:,d) = log(abs(XTr(:,d)+1));
                XTe(:,d) = log(abs(XTe(:,d)+1));

            elseif skew > 0.75

                XTr(:,d) = sqrt(abs(XTr(:,d)));
                XTe(:,d) = sqrt(abs(XTe(:,d)));

            elseif skew < -1

                XTr(:,d) = exp(XTr(:,d));
                XTe(:,d) = exp(XTe(:,d));

            elseif skew < -0.75

                XTr(:,d) = XTr(:,d).^2;
                XTe(:,d) = XTe(:,d).^2;

            end

            if strcmp(normalization_type,'standardization') 

                tmp_mean = mean(XTr(:,d));
                tmp_std  = std(XTr(:,d));
                
                XTr(:,d) = (XTr(:,d) - tmp_mean)/tmp_std;
                XTe(:,d) = (XTe(:,d) - tmp_mean)/tmp_std;
                
            end
 
        end

        % Fit model using the wrapping function
        model = trainClassifier(method, type_classif, yTr, XTr, parameters_method);
        
%         % Features ranking
%         betas = abs(model);
%         beta_sum = sum(betas,2);
%         beta_mean = mean(betas,2);
%         [~,beta_rank_sum] = sort(beta_sum,'descend');
%         [~,beta_rank_mean] = sort(beta_mean,'descend');
%         beta_rankvec_sum = [beta_rankvec_sum beta_rank_sum];
%         beta_rankvec_mean = [beta_rankvec_mean beta_rank_mean];
        
        
        % Predict
        [predTe,scoresTe] = predictClassifier(method,type_classif,model,XTe);
        [predTr,scoresTr] = predictClassifier(method,type_classif,model,XTr);  

        % Compute accuracy
         BER_Te(i,k) = computeBER(yTe,predTe);
         BER_Tr(i,k) = computeBER(yTr,predTr);
         OA_Te(i,k) = computeOA(yTe,predTe);
         OA_Tr(i,k) = computeOA(yTr,predTr);
         kappa_Te(i,k) = computeKAPPA(yTe,predTe);
         kappa_Tr(i,k) = computeKAPPA(yTr,predTr);

         fprintf('size of training sample : %2.2f%%',train_sampling(k)*100);
         fprintf('Testing BER: %.2f%%\n', BER_Te(i,k) * 100 );
         fprintf('Training BER: %.2f%%\n\n', BER_Tr(i,k) * 100 );

    end
    
end

        
%% diagnostic
figure;
plot_learning_curves(train_sampling',kappa_Te','Test',kappa_Tr','Train','% of training sample used','Kappa');
set(gca,'Fontsize',14);
figure;
plot_learning_curves(train_sampling',BER_Te','Test',BER_Tr','Train','% of training sample used','BER');
set(gca,'Fontsize',14);

