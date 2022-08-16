% first attempt to automatically classify snowflakes data
% ../events/snowflake_folders
clearvars; %close all;

                % alpha = 0.0001, lambda = 0.01, mini-batch size = 100
                % momentum = 0.9, Max no of iterations = 5000
                %parameters={0.00001,0.01,100,0.9,5000}; 

% USER DEFINED PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K=4; % Number of subsamples in the K-fold cross-validation
N_it=10; % Number of iterations of K-fold cross-validation
normalization_type = 'standardization';
method='baseline'; % can be either 'logistic', 'naive_bayes,'nn','svm','random_forest' or 'baseline'
parameters_method =  {1,1,'linear'};
% best for svm before Hu : parameters_method =  {10000,0.0001,'rbf'};
% default for logistic : {0.001,0,779,0.9,20000}; 
% default for svm-binary (6 classes) : {100,0.1,'rbf'} or {1,1,'linear'}
take_pca=true; % can be either true of false, if true, the PCA of the features is taken
type_classif='multiclass'; % can be binary or multiclass
save_log=false; % Can be true of false


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

dir_col = '/home/praz/Documents/MASC/masclab/events/columns';
col_filenames = dir(fullfile(dir_col,'*.mat'));
col_filenames = {col_filenames.name}';
col_picnames = dir(fullfile(dir_col,'*.png'));
col_picnames = {col_picnames.name}';

dir_dend = '/home/praz/Documents/MASC/masclab/events/dendrites';
dend_filenames = dir(fullfile(dir_dend,'*.mat'));
dend_filenames = {dend_filenames.name}';
dend_picnames = dir(fullfile(dir_dend,'*.png'));
dend_picnames = {dend_picnames.name}';

dir_bulros = '/home/praz/Documents/MASC/masclab/events/rosettes';
bulros_filenames = dir(fullfile(dir_bulros,'*.mat'));
bulros_filenames = {bulros_filenames.name}';
bulros_picnames = dir(fullfile(dir_bulros,'*.png'));
bulros_picnames = {bulros_picnames.name}';

%all_filenames = [fullfile(dir_graupels,graupels_filenames); fullfile(dir_melt,melt_filenames)];
%all_picnames  = [fullfile(dir_graupels,graupels_picnames); fullfile(dir_melt,melt_picnames)];

t_str_start = '20150101000000';
t_str_stop  = '20170101000000';

% small validation test
% for i=1:length(bulros_filenames)
%     
%     picname = bulros_picnames{i};
%     filename = strcat(picname(1:end-3),'mat');
%     tmp = sum(strcmp(bulros_filenames,filename));
%     if tmp ~= 1
%         disp(picname);
%     end
%     
% end

%% load flakes descriptors

% load matrix data
% graupels "g"
[Xg,Xlab_g,Xname_g,Xt_g] = load_processed_data(dir_graupels,t_str_start,t_str_stop);
yg = zeros(size(Xg,1),1); yg(yg==0) = 1;
% aggregates "a"
[Xa,Xlab_a,Xname_a,Xt_a] = load_processed_data(dir_agg,t_str_start,t_str_stop);
ya = zeros(size(Xa,1),1); ya(ya==0) = 2;
% melting snow "m"
[Xm,Xlab_m,Xname_m,Xt_m] = load_processed_data(dir_melt,t_str_start,t_str_stop);
ym = zeros(size(Xm,1),1); ym(ym==0) = 3;
% small crystals "s"
[Xs,Xlab_s,Xname_s,Xt_s] = load_processed_data(dir_small,t_str_start,t_str_stop);
ys = zeros(size(Xs,1),1); ys(ys==0) = 4;
% columns and needles "c"
[Xc,Xlab_c,Xname_c,Xt_c] = load_processed_data(dir_col,t_str_start,t_str_stop);
yc = zeros(size(Xc,1),1); yc(yc==0) = 5;
% dendrites and plates "d"
[Xd,Xlab_d,Xname_d,Xt_d] = load_processed_data(dir_dend,t_str_start,t_str_stop);
yd = zeros(size(Xd,1),1); yd(yd==0) = 6;    
labels = {'graup','agg','melt','smal','col','dend'};

% merge to one matrix
X = [Xg; Xa; Xm; Xs; Xc; Xd];
y = [yg; ya; ym; ys; yc; yd];
% Xname = [Xname_g; Xname_a; Xname_m; Xname_s; Xname_c; Xname_d];

% hogs features
[hog_g,hog2_g,hogname_g] = load_hog_features(dir_graupels,t_str_start,t_str_stop);
[hog_a,hog2_a,hogname_a] = load_hog_features(dir_agg,t_str_start,t_str_stop);
[hog_m,hog2_m,hogname_m] = load_hog_features(dir_melt,t_str_start,t_str_stop);
[hog_s,hog2_s,hogname_s] = load_hog_features(dir_small,t_str_start,t_str_stop);
[hog_c,hog2_c,hogname_c] = load_hog_features(dir_col,t_str_start,t_str_stop);
[hog_d,hog2_d,hogname_d] = load_hog_features(dir_dend,t_str_start,t_str_stop);

hog = [hog_g; hog_a; hog_m; hog_s; hog_c; hog_d];
hog2 = [hog2_g; hog2_a; hog2_m; hog2_s; hog2_c; hog2_d];
hogname = [hogname_g; hogname_a; hogname_m; hogname_s; hogname_c; hogname_d];

%remove empty columns
idx_null = find(sum(hog)==0);
hog(:,idx_null) = [];
idx_null = find(sum(hog2)==0);
hog2(:,idx_null) = [];

hog = [hog hog2];
[coeff,score,latent] = pca(hog);
%[coeff2,score2,latent2] = pca(hog2);

X = [X score(:,1:20)];% score2(:,1:10)];

N = size(X,1);
D = size(X,2);
%parameters_method{1} = ones(D,1);
N_classes = length(unique(y));

%%
%X_glcm = X(:,78:end);
%X = X(:,1:77);
%D = size(X,2);

%[coeff,score,latent] = pca(X_glcm);
%X = [X score(:,1:3)];
%D = size(X,2);

% remove features with only zero entries in glcm
% tmp_sum = sum(X_glcm);
% idx2remove = find(tmp_sum==0);
% X_glcm(:,idx2remove) = [];
% D_glcm = size(X,2);

% normalization and feature transformation
%idx1 = [1,2,3,4,5,6,7,8,9,10,11,13,14,15,16,28,29,31,39,47,48,49,50,51,52,55,69,70,72,74,75,76,77,78,79,81,83,84,86,87,88,89,90,92];
%idx2 = [27,32,33,34,35,38,50,51,52,55];
%idx3 = [80,91,95,96,97];%[19,20,22,25,50,53,63,64,67];


%Xlab_a{49:52} = {};
%X_old = X;
%X2 = X.^2;
%X3 = sqrt(abs(X));
%X = [X X2];

D = size(X,2);

for i=1:D
      
    skew = skewness(X(:,i));
    
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
        
    elseif strcmp(normalization_type,'rescaling')
        
        tmp_min = min(X(:,i));
        tmp_max = max(X(:,i));
        X(:,i) = (X(:,i) - tmp_min)/(tmp_max - tmp_min);
        
    end
    
end

%%


% PCA
%[coeff,score,latent] = pca(X);
%X = score;

% randomize 
% rand_idx = randperm(N);
% X = X(rand_idx,:);
% y = y(rand_idx);


% Initialize global confusion matrix
global_confmat_Te = zeros(N_classes);
global_confmat_Tr = zeros(N_classes);


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
        
        % Get predictions
        [predTe,~] = predictClassifier(method,type_classif,model,XTe);
        [predTr,~] = predictClassifier(method,type_classif,model,XTr);    
           
        % Compute accuracy
        BER_Te(i,k) = computeBER(yTe,predTe);
        BER_Tr(i,k) = computeBER(yTr,predTr);
        OA_Te(i,k) = computeOA(yTe,predTe);
        OA_Tr(i,k) = computeOA(yTr,predTr);
        kappa_Te(i,k) = computeKAPPA(yTe,predTe);
        kappa_Tr(i,k) = computeKAPPA(yTr,predTr);
        
        % Update global confusion matrix
        global_confmat_Te = global_confmat_Te + confusionmat(yTe,predTe);
        global_confmat_Tr = global_confmat_Tr + confusionmat(yTr,predTr);
        
        idx_fail{i,k}=idxTe((predTe~=yTe));
        real_labels{i,k}=yTe((predTe~=yTe));
        predicted_labels{i,k}=predTe((predTe~=yTe));
         
        fprintf('\nTesting BER: %.2f%%\n', BER_Te(i,k) * 100 );
        fprintf('Training BER: %.2f%%\n\n', BER_Tr(i,k) * 100 );
        
    end
        
%         if take_pca % Perform PCA on training set
%             % First normalize testing set with training centroids
%             for m=1:size(XTr,2)
%                 XTe(:,m)=XTe(:,m)-mean(XTr(:,m));
%             end
%             % Then take (centered) PCA of training set
%             [coeff,score,latent] = pca(XTr);
%             % Assign PCA components to training set
%             XTr=score; 
%             % Multiply testing set by PCA coefficient (change of basis)
%             XTe=XTe*coeff;
%         end
%      
%         % Handle some special cases for some methods
%         if strcmp(method, 'random_forest') && pca
%            XTr=XTr(:,1:100); % Take only 100 first components 
%            XTe=XTe(:,1:100); % Take only 100 first components 
%         elseif strcmp(method,'naive_bayes')
%             idx_ok=var(XTr)>1E-6; % We have to remove features with zero variance
%             XTr=XTr(:,idx_ok);
%             XTe=XTe(:,idx_ok);
%         elseif strcmp(method,'nn') && ~ take_pca
%             % In this case we need to normalize train and test
%              [XTr, mu, sigma] = zscore(XTr);
%              XTe=normalize(XTe,mu,sigma);
%         end
%            
%         % Fit model using the wrapping function
%         model = trainClassifier(method, type_classif, yTr, XTr, parameters_method); 
%         
%         % Get predictions
%         [predTe,~] = predictClassifier(method,type_classif,model,XTe);
%         [predTr,~] = predictClassifier(method,type_classif,model,XTr);    
%     
%         global_confusion_mat=global_confusion_mat+confusionmat(int8(yTe),int8(predTe));
%         
%         predErrTe(i,k) = computeBER(yTe,predTe);
%         predErrTr(i,k) = computeBER(yTr,predTr);
%         
%         if predErrTe(i,k)<best_score
%             best_model=model;
%             if take_pca % Also keep track of the change of basis matrix
%                 best_coeff=coeff;
%             else
%                 best_coeff=[];
%             end
%         end
%         
%         idx_fail{i,k}=idxTe((predTe~=yTe));
%         real_labels{i,k}=yTe((predTe~=yTe));
%         predicted_labels{i,k}=predTe((predTe~=yTe));
%         
%         fprintf('\nTesting error: %.2f%%\n\n', predErrTe(i,k) * 100 );
%         fprintf('\nTraining error: %.2f%%\n\n', predErrTr(i,k) * 100 );
%     end
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

% confusion matrices
% normalization by lines
for i=1:size(global_confmat_Te,1)
    global_confmat_Te(i,:) = round(global_confmat_Te(i,:)./sum(global_confmat_Te(i,:)) * 100,2);
    global_confmat_Tr(i,:) = round(global_confmat_Tr(i,:)./sum(global_confmat_Tr(i,:)) * 100,2);
end

% global normalization
%global_confmat_Te = round(global_confmat_Te./sum(global_confmat_Te(:)) * 100,2);
%global_confmat_Tr = round(global_confmat_Tr./sum(global_confmat_Tr(:)) * 100,2);
figure;
heatmap(global_confmat_Tr,labels,labels,1,'Colormap','red','ShowAllTicks',1,'UseLogColorMap',true,'Colorbar',true);
title('Train confusion matrix');
set(gca,'Fontsize',14);
figure;
heatmap(global_confmat_Te,labels,labels,1,'Colormap','red','ShowAllTicks',1,'UseLogColorMap',true,'Colorbar',true);
title('Test confusion matrix');
set(gca,'Fontsize',14);


%% visualization

if false

    idx_fail_all = [];
    for i=1:size(idx_fail,1)

        idx_fail_all = [idx_fail_all [idx_fail{i,:}]];

    end

    idx_fail_all = sort(idx_fail_all);
    occurence = 1;
    idx_fail_final = [];
    idx_fail_occurence = [];
    for i=2:length(idx_fail_all)

            if idx_fail_all(i) == idx_fail_all(i-1)

                occurence = occurence + 1;

            else

                idx_fail_final(end+1) = idx_fail_all(i-1);
                idx_fail_occurence(end+1) = occurence;
                occurence = 1;

            end       

    end

    if idx_fail_all(end) == idx_fail_all(end-1)

        idx_fail_occurence(end) = idx_fail_occurence(end) + 1;

    else

        idx_fail_final(end+1) = idx_fail_all(end);
        idx_fail_occurence(end+1) = 1;

    end

    % sort by occurence
    [idx_fail_occurence,idx_order] = sort(idx_fail_occurence,'descend');
    idx_fail_final = idx_fail_final(idx_order);

    figure;
    for i=1:length(idx_fail_final)

        clf();
        imgname = all_picnames{idx_fail_final(i)};
        img = imread(imgname);
        imshow(img);
        title(sprintf('Occurrences : %f/%u',idx_fail_occurence,N_it));
        pause;

    end
  
end
    


