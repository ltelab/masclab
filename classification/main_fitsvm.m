% first attempt to automatically classify snowflakes data
% ../events/snowflake_folders
clearvars; %close all;

                % alpha = 0.0001, lambda = 0.01, mini-batch size = 100
                % momentum = 0.9, Max no of iterations = 5000
                %parameters={0.00001,0.01,100,0.9,5000}; 

% USER DEFINED PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K=4; % Number of subsamples in the K-fold cross-validation
N_it=1; % Number of iterations of K-fold cross-validation
normalization_type = 'standardization';
method='svm'; % can be either 'logistic', 'naive_bayes,'nn','svm','random_forest' or 'baseline'
parameters_method =  {1,1,'linear'};
% best for svm before Hu : parameters_method =  {10000,0.0001,'rbf'};
% default for logistic : {0.001,0,779,0.9,20000}; 
% default for svm-binary (6 classes) : {100,0.1,'rbf'} or {1,1,'linear'}
take_pca=true; % can be either true of false, if true, the PCA of the features is taken
type_classif='multiclass'; % can be binary or multiclass
save_log=false; % Can be true of false


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% path to training data
dir_data = '/home/praz/Documents/MASC/masclab/events/mixed_sample_final';
data_filenames = dir(fullfile(dir_data,'20*.mat'));
data_filenames = {data_filenames.name}';
data_picnames = dir(fullfile(dir_data,'*.png'));
data_picnames = {data_picnames.name}';

% time interval to take into consideration
t_str_start = '20150101000000';
t_str_stop  = '20170101000000';

% path to labels
dir_labels = '/home/praz/Documents/MASC/masclab/labelling_flakes';
labels_filename = 'common_all_labelled_6classes.mat';

%% load flakes descriptors

% load the training matrix X
[X,Xlab,Xname,Xt] = load_processed_data(dir_data,t_str_start,t_str_stop);

% load the labels vector y
y = load_labels(dir_data,dir_labels,labels_filename);

% load the hog features
%[hog,hog2,hogname] = load_hog_features(dir_data,t_str_start,t_str_stop);

%remove empty columns
% idx_null = find(sum(hog)==0);
% hog(:,idx_null) = [];
% idx_null = find(sum(hog2)==0);
% hog2(:,idx_null) = [];

% take pca of hog features
%hog_tot = [hog hog2];
%[coeff,score,latent] = pca(hog_tot);
%classif_params.pca_hogcoeff = coeff;
%[coeff2,score2,latent2] = pca(hog2);

% add some pc to the X matrix
%X = [X score(:,1:20)];% score(:,1:10) score2(:,1:10)];

% remove "undetermined" labels
idx_unknown = find(y==0);
y(idx_unknown) = [];
X(idx_unknown,:) = [];
data_picnames(idx_unknown) = [];
data_filenames(idx_unknown) = [];

% 5 = dendrites & plates
%y(y==6) = 5;
% 6 = columsn & needles
%y(y==7 | y==8) = 6;


% some hack
%idx_small = find(X(:,2)<30 & y==3);
%y(idx_small) = 4;

N = size(X,1);
D = size(X,2);
N_classes = length(unique(y));

labels = {'graup','agg','melt','smal','dend','col'};

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

%% Algorithm training & testing
boxconstraint_vec = 1;%[1e-4 1e-3 1e-2 1e-1 1 1e1 1e2 1e3 1e4 1e5];
kernelscale_vec = 1;%[1e-4 1e-3 1e-2 1e-1 1 1e1 1e2 1e3 1e4 1e5];

for m=1:length(boxconstraint_vec)
for n=1:length(kernelscale_vec)    

    
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
        
        % apply 1vsall training
        model = {};

        % loop over each class and apply "one vs all" rule
        y_tmp = yTr;
        for j=1:N_classes

            y_tmp(yTr~=j) = -1;
            y_tmp(yTr==j) = +1;   

            SVMModel = fitcsvm(XTr,y_tmp);%,'BoxConstraint',boxconstraint_vec(m),'KernelFunction','rbf','KernelScale',kernelscale_vec(n)); %'auto'
            model{j} = SVMModel;

        end 
        
        
        % classify training set
        scores = [];
        for j=1:N_classes
            
            [label,score] = predict(model{j},XTr);
            scores(:,j) = score(:,2);     
            
        end
        [~,predTr] = max(scores,[],2);
        
        % classify testing set
        scores = [];
        for j=1:N_classes
            
            [label,score] = predict(model{j},XTe);
            scores(:,j) = score(:,2);     
            
        end
        [~,predTe] = max(scores,[],2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % apply 1vs1 training
%         model = {};
%         
%         % double loop over the classes
%         y_tmp = yTr;
%         for m=1:N_classes-1
%             for n=(m+1):N_classes
%                 
%                 model{m,n}.idx2keep = find(yTr==m | yTr==n);
%                 y_tmp = yTr(model{m,n}.idx2keep);
%                 y_tmp(y_tmp~=m) = -1;
%                 y_tmp(y_tmp==m) = 1;
%                 model{m,n}.XTr = XTr(model{m,n}.idx2keep,:);
%                 model{m,n}.SVMModel = fitcsvm(model{m,n}.XTr,y_tmp,'BoxConstraint',0.01,'KernelFunction','linear','KernelScale',1);
%                 
%             end
%         end
%         
%         % classifiy training set 
%         % initialisation of a matrix NTrx6 of classvotes
%         classvote = zeros(length(yTr),N_classes);
%         for m=1:N_classes-1
%             for n=(m+1):N_classes
%                 
%                 XTr_tmp = XTr;%(model{m,n}.idx2keep,:);
%                 [label,score] = predict(model{m,n}.SVMModel,XTr_tmp);
%                 %label_ext = zeros(length(yTr),1);
%                 %label_ext(model{m,n}.idx2keep) = label;
%                 label_ext = label;
%                 
%                 classvote(label_ext==1,m) = classvote(label_ext==1,m) + 1;
%                 classvote(label_ext==-1,n) = classvote(label_ext==-1,n) + 1;
%         
%             end
%         end
%         [~,predTr] = max(classvote,[],2);
%         
%         % classify testing set
%         % initialisation of a matrix NTex6 of classvotes
%         classvote = zeros(length(yTe),N_classes);
%         for m=1:N_classes-1
%             for n=(m+1):N_classes
%                 
%                 XTe_tmp = XTe;
%                 [label,score] = predict(model{m,n}.SVMModel,XTe_tmp);
%                 label_ext = label;
%                 
%                 classvote(label_ext==1,m) = classvote(label_ext==1,m) + 1;
%                 %classvote(label_ext==1,n) = classvote(label_ext==1,n) - score(label_ext==1);
%                 classvote(label_ext==-1,n) = classvote(label_ext==-1,n) + 1;
%                 %classvote(label_ext==-1,m) = classvote(label_ext==-1,m) - score(label_ext==-1);
%                 
%             end
%         end
%         [~,predTe] = max(classvote,[],2);
%         
%         % compute centroid of train sample
%         classes=unique(yTr);
%         for j=1:length(classes)
%             centroids(j,:)=median(XTr(yTr==classes(j),:));
%         end
%         
%         
%         idx_equality = find(max(classvote,[],2)<N_classes-1);
%         for j=1:length(idx_equality)
%             idx_duel = find(classvote(idx_equality(j),:)==max(classvote(idx_equality(j),:)));
%             if length(idx_duel)>1
% 
%                 dist2centroid = pdist2(XTe(idx_equality(j),:),centroids);
%                 dist2centroid(idx_duel) = -dist2centroid(idx_duel);
%                 [~,newPred] = min(dist2centroid);
%                 %fprintf('old pred : %u',predTe(idx_equality(j)));
%                 predTe(idx_equality(j)) = newPred;
%                 %fprintf('     new pred : %u\n',predTe(idx_equality(j)));
%             
%             end
%         end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

         % Compute accuracy
         BER_Te(i,k) = computeBER(yTe,predTe);
         BER_Tr(i,k) = computeBER(yTr,predTr);
         OA_Te(i,k) = computeOA(yTe,predTe);
         OA_Tr(i,k) = computeOA(yTr,predTr);
         kappa_Te(i,k) = computeKAPPA(yTe,predTe);
         kappa_Tr(i,k) = computeKAPPA(yTr,predTr);
%         
%         % Update global confusion matrix
         global_confmat_Te = global_confmat_Te + confusionmat(yTe,predTe);
         global_confmat_Tr = global_confmat_Tr + confusionmat(yTr,predTr);
%         
%         idx_fail{i,k}=idxTe((predTe~=yTe));
%         real_labels{i,k}=yTe((predTe~=yTe));
%         predicted_labels{i,k}=predTe((predTe~=yTe));
%          
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

% confusion matrices
global_confmat_Te = round(global_confmat_Te./sum(global_confmat_Te(:)) * 100,2);
global_confmat_Tr = round(global_confmat_Tr./sum(global_confmat_Tr(:)) * 100,2);
figure;
heatmap(global_confmat_Tr,labels,labels,1,'Colormap','red','ShowAllTicks',1,'UseLogColorMap',true,'Colorbar',true);
title('Train confusion matrix');
figure;
heatmap(global_confmat_Te,labels,labels,1,'Colormap','red','ShowAllTicks',1,'UseLogColorMap',true,'Colorbar',true);
title('Test confusion matrix');



%% visualization

if 0

idx_fail_all = [];
real_labels_all = [];
predicted_labels_all = [];
for i=1:size(idx_fail,1)

    idx_fail_all = [idx_fail_all [idx_fail{i,:}]];
    real_labels_all = [real_labels_all [real_labels{i,:}]];
    predicted_labels_all = [predicted_labels_all [predicted_labels{i,:}]];

end


[idx_fail_all,idx_sorted] = sort(idx_fail_all);
real_labels_all = real_labels_all(idx_sorted);
predicted_labels_all = predicted_labels_all(idx_sorted);

occurence = 1;
idx_fail_final = [];
idx_fail_occurence = [];
idx_fail_real_label = [];
idx_fail_pred_label = [];

for i=2:length(idx_fail_all)

    if idx_fail_all(i) == idx_fail_all(i-1)

        occurence = occurence + 1;

    else

        idx_fail_final(end+1) = idx_fail_all(i-1);
        idx_fail_real_label(end+1) = real_labels_all(i-1);
        idx_fail_pred_label(end+1) = predicted_labels_all(i-1);
        idx_fail_occurence(end+1) = occurence;
        occurence = 1;

    end       

end

if idx_fail_all(end) == idx_fail_all(end-1)

    idx_fail_occurence(end) = idx_fail_occurence(end) + 1;

else

    idx_fail_final(end+1) = idx_fail_all(end);
    idx_fail_real_label(end+1) = real_labels_all(end);
    idx_fail_pred_label(end+1) = predicted_labels_all(end);
    idx_fail_occurence(end+1) = 1;

end


% sort by occurence
[idx_fail_occurence,idx_order] = sort(idx_fail_occurence,'descend');
idx_fail_final = idx_fail_final(idx_order);
idx_fail_real_label = idx_fail_real_label(idx_order);
idx_fail_pred_label = idx_fail_pred_label(idx_order);

figure;
for i=1:length(idx_fail_final)

    clf();
    imgname = data_picnames{idx_fail_final(i)};
    img = imread(imgname);
    subplot(2,3,1);
    imshow(img);
    h_text = subplot(2,3,2:3);
    xl = xlim(h_text); 
    xPos = xl(1) + diff(xl) / 2; 
    yl = ylim(h_text); 
    yPos = yl(1) + diff(yl) / 2; 
    plot_text = text(xPos, yPos, sprintf('Pred : %u. Real : %u. Occurrences : %2.0f/%u',idx_fail_pred_label(i),idx_fail_real_label(i),idx_fail_occurence(i),N_it), 'Parent', h_text);
    set(plot_text, 'HorizontalAlignment', 'center','FontSize',12,'FontWeight','demi');
    set(h_text,'visible','off');  
    h_text2 = subplot(2,3,4:6);
    xl = xlim(h_text); 
    xPos = xl(1) + diff(xl) / 2; 
    yl = ylim(h_text); 
    yPos = yl(1) + diff(yl) / 2; 
    plot_text2 = text(xPos,yPos,sprintf('%s',imgname),'interpreter','none');
    set(plot_text2, 'HorizontalAlignment', 'center','FontSize',12,'FontWeight','demi');
    set(h_text2,'visible','off'); 
    %title(sprintf('Pred : %u. Real : %u. Occurrences : %u/%u',idx_fail_pred_label(i),idx_fail_real_label(i),idx_fail_occurence(i),N_it));
    pause;

end

end

    


