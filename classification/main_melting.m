clearvars; close all;

                % alpha = 0.0001, lambda = 0.01, mini-batch size = 100
                % momentum = 0.9, Max no of iterations = 5000
                %parameters={0.00001,0.01,100,0.9,5000}; 

% USER DEFINED PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K=4; % Number of subsamples in the K-fold cross-validation
N_it=1; % Number of iterations of K-fold cross-validation
normalization_type = 'standardization';
method='logistic'; % can be either 'logistic', 'naive_bayes,'nn','svm','random_forest' or 'baseline'
parameters_method =  {0.001,0.001,1000,0,10000};
%parameters_method = {1,1,'linear'};
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
%all except absurd
%feat_vec = [1:1:96]';
%feat_vec([32 33 37 38 41:48 67 77 79 80]) = [];

Dwanted = 25;
load('out_melting_10it_rand_FS_CV4_HSS_80D.mat');
feat_vec = feat_mat(:,Dwanted+1);
feat_vec = feat_vec(feat_vec~=0);
%feat_vec = sort(feat_vec,'ascend');


% compute new descriptors (if necessary)
process_new_descriptors(dir_data,t_str_start,t_str_stop)

% load the training matrix X
[X,Xlab,Xname,Xt] = load_processed_data(dir_data,t_str_start,t_str_stop,feat_vec);

% load the labels vector y
[~,~,y] = load_labels(dir_data,dir_labels,labels_filename,label_version);

% convert in double
y = double(y);

% remove "undetermined", "not labelled yet" or "ambiguous" labels
idx_unknown = find(y<0);
y(idx_unknown) = [];
X(idx_unknown,:) = [];
data_picnames(idx_unknown) = [];
data_filenames(idx_unknown) = [];
fprintf('%u samples discarded because labelled as "undetermined" or "ambiguous" \n',length(idx_unknown));

N = size(X,1);
D = size(X,2);
N_classes = length(unique(y));
labels = {'dry','melting'};

% computing the weights, if we want to penalize the cost
for i=1:N_classes
    Nc(i,1) = sum(y==i-1);
end
[Ncmax,Ncmaxidx] = max(Nc);
for i=1:N_classes
    weights(i,1) = Ncmax/Nc(i);
end
%parameters_method{end+1} = weights;


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

%% Features discretization
% Nbins = 200;
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
% alternative
% for i=1:D
%     
%     qgrid = [0:0.01:1]';
%     grid = quantile(X(:,i),qgrid);
%     Nbins = length(grid)-1;
%     
%     for j=1:N
%         idx = find(grid>=X(j,i),1,'first');
%         X(j,i) = grid(idx)-(grid(end)-grid(1))/(Nbins*2);
%     end
%     
% end


%%

  
% Initialize global confusion matrix
global_confmat_Te = zeros(N_classes);
global_confmat_Tr = zeros(N_classes);

% Initialize beta ranking vectors
%beta_rankvec_sum = [];
%beta_rankvec_mean = [];
beta_sum_mat = [];

for i = 1:N_it
    setSeed(i);
    fprintf('\nStarting iteration number %.d%\n\n', i );
    
    idx_perm = randperm(N);
    Nk = floor(N/K);
    idxCV=[];
    for k = 1:K
        idxCV(k,:) = idx_perm(1+(k-1)*Nk:k*Nk);
    end
    
    for k = 1:K
        fprintf('\nK-fold %.d/%.d% \n\n', k,K );
        fprintf('\n')
        
        % get k'th subgroup in test, others in train
        idxTe = idxCV(k,:);
        idxTr = idxCV([1:k-1 k+1:end],:);
        idxTr = idxTr(:);
        yTe = y(idxTe);
        XTe = X(idxTe,:);
        yTr = y(idxTr);
        XTr = X(idxTr,:); 
        
        % Fit model using the wrapping function
        model = trainClassifier(method, type_classif, yTr, XTr, parameters_method);
        
        % Features ranking
        if features_ranking
            betas = abs(model);
            beta_sum = sum(betas,2);
            beta_sum = beta_sum(2:end,:);
            beta_sum_mat(:,end+1) = beta_sum;
            %beta_mean = mean(betas,2);
            %[~,beta_rank_sum] = sort(beta_sum,'descend');
            %[~,beta_rank_mean] = sort(beta_mean,'descend');
            %beta_rankvec_sum = [beta_rankvec_sum beta_rank_sum];
            %beta_rankvec_mean = [beta_rankvec_mean beta_rank_mean];
        end
        
        
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
        
        % Update global confusion matrix
        cmTe = confusionmat(yTe,predTe);
        
        % if there is no melting snow in the test set...
        %if size(cmTe,1) < 6
        %    cmTe = [cmTe(:,1:2) [0;0;0;0;0;] cmTe(:,3:5)];
        %    cmTe = [cmTe(1:2,:); [0 0 0 0 0 0]; cmTe(3:5,:)];
        %end
        
        global_confmat_Te = global_confmat_Te + cmTe;
        global_confmat_Tr = global_confmat_Tr + confusionmat(yTr,predTr);
        %         
        %         idx_fail{i,k}=idxTe((predTe~=yTe));
        %         real_labels{i,k}=yTe((predTe~=yTe));
        %         predicted_labels{i,k}=predTe((predTe~=yTe));
        %          
        fprintf('\nTesting BER: %.2f%%\n', BER_Te(i,k) * 100 );
        fprintf('Training BER: %.2f%%\n\n', BER_Tr(i,k) * 100 );

    end
        
end


%% diagnostic
fprintf('*****************************\n');
fprintf('Mean Testing BER    : %.2f%% pm %.2f%%\n', mean(BER_Te(:))*100, std(BER_Te(:))*100);
fprintf('Mean Training BER   : %.2f%% pm %.2f%%\n\n', mean(BER_Tr(:))*100, std(BER_Tr(:))*100);
fprintf('Mean Testing OA     : %.2f%% pm %.2f%%\n', mean(OA_Te(:))*100, std(OA_Te(:))*100);
fprintf('Mean Training OA    : %.2f%% pm %.2f%%\n\n', mean(OA_Tr(:))*100, std(OA_Tr(:))*100);
fprintf('Mean Testing Kappa  : %.2f%% pm %.2f%%\n', mean(kappa_Te(:))*100, std(kappa_Te(:))*100);
fprintf('Mean Training Kappa : %.2f%% pm %.2f%%\n', mean(kappa_Tr(:))*100, std(kappa_Tr(:))*100);
fprintf('*****************************\n');
%%
% confusion matrices
global_confmat_Te = round(global_confmat_Te./sum(global_confmat_Te(:)) * 100,2);
global_confmat_Tr = round(global_confmat_Tr./sum(global_confmat_Tr(:)) * 100,2);

if display_confmat
    figure;
    heatmap(global_confmat_Tr,labels,labels,1,'Colormap','red','ShowAllTicks',1,'UseLogColorMap',true,'Colorbar',true,'Fontsize',14);
    set(gca,'Fontsize',14);
    title('Train confusion matrix');
    figure;
    heatmap(global_confmat_Te,labels,labels,1,'Colormap','red','ShowAllTicks',1,'UseLogColorMap',true,'Colorbar',true,'Fontsize',14);
    set(gca,'Fontsize',14);
    title('Test confusion matrix');
end


%% features ranking illustration
if features_ranking
    beta_sum_vec = mean(beta_sum_mat,2);
    [beta_weight,beta_id] = sort(beta_sum_vec,'descend');
    beta_weight = beta_weight./sum(beta_weight)*100;
    
    figure; hold on; box on;
    bar(beta_weight);
    set(gca,'xtick',1:D,'xticklabel',int2str(feat_vec(beta_id)'));
    xlabel('Feature ID');
    ylabel('Feature weight [%]');
    set(gca,'Fontsize',14);
    set(gca,'Linewidth',1.5);
    title('Feature Ranking');
  
end

