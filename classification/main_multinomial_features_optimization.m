function [beta_weight,beta_id,mean_BER_Te,mean_kappa_Te] = main_multinomial_features_optimization(feat_vec)


% USER DEFINED PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K=4; % Number of subsamples in the K-fold cross-validation
N_it=1; % Number of iterations of K-fold cross-validation
normalization_type = 'standardization';
method='logistic'; % can be either 'logistic', 'naive_bayes,'nn','svm','random_forest' or 'baseline'
parameters_method =  {0.001,0.01,1000,0,10000};
%parameters_method = {1,1,'linear'};
% best for svm before Hu : parameters_method =  {10000,0.0001,'rbf'};
% default for logistic : {0.001,0,779,0.9,20000}; 
% default for svm-binary (6 classes) : {100,0.1,'rbf'} or {1,1,'linear'}
type_classif='multiclass'; % can be binary or multiclass
save_log=false; % Can be true of false

features_ranking = true;
display_confmat = false;
label_version = '1.1';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% path to training data
dir_data = '/home/praz/Documents/MASC/MASC_labelling/mixed_sample_merged';
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
dir_labels = '/home/praz/Documents/MASC/masclab/labelling_GUI/labels_v.1.1';
labels_filename = 'FINAL1.mat';
%dir_labels = '/media/praz/MyData/Massimo_data';
%labels_filename = 'CP_labels_all_sorted.mat';

%% load flakes descriptors

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
%beta_rankvec_sum = [];
%beta_rankvec_mean = [];
beta_sum_mat = [];

fprintf('\nStarting Classification CV with %u features\n',D);

for i = 1:N_it
    %setSeed(randi(1000));
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
    
%     figure; hold on; box on;
%     bar(beta_weight);
%     set(gca,'xtick',1:D,'xticklabel',int2str(beta_id));
%     xlabel('Feature ID');
%     ylabel('Feature weight [%]');
%     set(gca,'Fontsize',12);
    %drawnow;
    
    
end




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

mean_BER_Te = mean(BER_Te(:));
mean_kappa_Te = mean(kappa_Te(:));

end

    


