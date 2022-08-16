% small script to do features optimization
clear all; close all;

% load data
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

% what to optimize ?
classif_subject = 'melting'; % geometry lambda = 1| riming lambda = 0.01 | melting lambda = 

% list of descriptors wanted
feat_vec = [1:1:96]';

% compute new descriptors (if necessary)
process_new_descriptors(dir_data,t_str_start,t_str_stop)

% load the training matrix X
[X,Xlab,Xname,Xt] = load_processed_data(dir_data,t_str_start,t_str_stop,feat_vec);
%feat_vec(101) = [];
%feat_vec([2 6 7 14 15 17 18 26 33 34 37 38 39 41:48 51 53 70 77 79 80 92]) = [];
feat_vec([32 33 37 38 41:48 67 77 79 80]) = [];
dim_ini = length(feat_vec);

%%
if strcmp(classif_subject,'geometry')
    
    % load the labels vector y
    [y,~,~] = load_labels(dir_data,dir_labels,labels_filename,'1.1');

    % remove "undetermined", "not labelled yet" or "ambiguous" labels
    idx_unknown = find(y<=0);
    y(idx_unknown) = [];
    X(idx_unknown,:) = [];
    data_picnames(idx_unknown) = [];
    data_filenames(idx_unknown) = [];
    fprintf('%u samples discarded because labelled as "undetermined" or "ambiguous" \n',length(idx_unknown));

    % classification scheme 
    
    % 1: small particle
    % 2: column & needle
    y(y==3) = 2;
    % 3: sectored plates & dendrites
    y(y==5 | y==6) = 3;
    % 4: remove plates, add aggregates
    idx_plates = find(y==4);
    y(idx_plates) = [];
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
    X(idx_combcrystals,:) = [];
    data_picnames(idx_combcrystals) = [];
    data_filenames(idx_combcrystals) = [];
    % add graupel-like particles into graupel
    %y(y==4 & (yR==4 | yR==5)) = 5;
    labels = {'smal','col','dend','agg','graup','capcol'};

elseif strcmp(classif_subject,'riming')
    
    % load the labels vector y
    [~,y,~] = load_labels(dir_data,dir_labels,labels_filename,'1.1');

    % remove "undetermined", "not labelled yet" or "ambiguous" labels
    idx_unknown = find(y<=0);
    y(idx_unknown) = [];
    X(idx_unknown,:) = [];
    data_picnames(idx_unknown) = [];
    data_filenames(idx_unknown) = [];
    fprintf('%u samples discarded because labelled as "undetermined" or "ambiguous" \n',length(idx_unknown));    
    
elseif strcmp(classif_subject,'melting')
    
    % load the labels vector y
    [~,~,y] = load_labels(dir_data,dir_labels,labels_filename,'1.1');
    y = double(y);
    
    % remove "undetermined", "not labelled yet" or "ambiguous" labels
    idx_unknown = find(y<0);
    y(idx_unknown) = [];
    X(idx_unknown,:) = [];
    data_picnames(idx_unknown) = [];
    data_filenames(idx_unknown) = [];
    fprintf('%u samples discarded because labelled as "undetermined" or "ambiguous" \n',length(idx_unknown));  
    
end




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


%% relieff
if 0
    k_relief = [1:10:400]';
    ranked_vec = [];
    weight_vec = [];
    for i=1:length(k_relief)

        fprintf('%u/%u\n',i,length(k_relief));
        [ranked,weight] = relieff(X,y,k_relief(i));
        ranked_vec = [ranked_vec ranked'];
        weight_vec = [weight_vec weight'];

    end
end


%% forward selection according to relieff features weights
if 0
    load('feat_list_relieff_k100.mat');
    ranked = ranked';
    mean_BER_Te = [];
    mean_kappa_Te = [];
    xaxis4BER = [];
    feat_mat = [];
    dim_max = 30;

    for k=1:dim_max

        fprintf('Start iteration %u with %u features \n',k,40);
        feat_test = ranked(1:k);
        Xtest = X(:,feat_test);
        out = CrossValidation('logistic','multiclass',{0.001,1,1000,0,5000},4,1,Xtest,y,1,0,0,1);
        mean_BER_Te(end+1) = out.BER_Te;
        mean_kappa_Te(end+1) = out.kappa_Te;
        xaxis4BER(end+1) = k;
        feat_mat = [feat_mat [feat_test; zeros(dim_max-length(feat_test),1)]];

    end
end




%% backward elimination based on beta weights
if 0
    k=0;
    mean_BER_Te = [];
    mean_kappa_Te = [];
    xaxis4BER = [];
    feat_mat = feat_vec;
    beta_vec = [];
    while length(feat_vec) > 2

        k=k+1;
        
        fprintf('Start iteration %u with %u features \n',k,length(feat_vec));
        
        % run Nit instance of CV
        out = CrossValidation('logistic','multiclass',{0.001,1,1000,0,5000},4,1,X,y,0,0,0);
        mean_BER_Te(end+1) = out.BER_Te;
        mean_kappa_Te(end+1) = out.kappa_Te;
        xaxis4BER(end+1) = length(feat_vec);
        
        % retrieve beta rank, and eliminate one beta accordingly
        feat_vec = feat_vec(out.beta_id);
        feat_vec(end) = [];
        feat_vec = sort(feat_vec,'ascend');
        X(:,out.beta_id(end)) = [];
        
        % store the new feat_vec in feat_mat
        feat_mat = [feat_mat [feat_vec; zeros(dim_ini-length(feat_vec),1)]];
        
        
        
%         beta_vec(end+1) = beta_id(end);
%         feat_mat = [feat_mat [feat_vec; zeros(dim_ini-length(feat_vec),1)]];
% 
%         fprintf('\nfeature to remove : %u \n',beta_id(end));
%         feat_vec = feat_vec(beta_id);
%         fprintf('feature removed : %u \n',feat_vec(end));
%         %feat_vec = feat_vec(1:end-1);
%         feat_vec(end-1) = [];
%         feat_vec = sort(feat_vec,'ascend');

    end
end

%% backward elimination based on kappa performance (long as hell)
if 0
    k=0;
    mean_BER_Te = [];
    mean_kappa_Te = [];
    xaxis4BER = [];
    feat_mat = feat_vec;
    beta_vec = [];
    dim_ini = length(feat_vec);
    
    while length(feat_vec)>2
        
        k=k+1;
        fprintf('\n Start iteration %u with %u features \n',k,length(feat_vec));
        
        % try to remove one feature, compute and store kappa
        pool_BER = [];
        pool_kappa = [];
        for i=1:length(feat_vec)
            
            fprintf('Try removing feature %u...\n',feat_vec(i));
            feat_vec_test = feat_vec;
            feat_vec_test(i) = [];
            Xtest = X(:,feat_vec_test);
            out = CrossValidation('logistic','multiclass',{0.001,1,1000,0,5000},4,10,Xtest,y,1,0,0);
            pool_BER(end+1) = out.BER_Te;
            pool_kappa(end+1) = out.kappa_Te;
            
        end
        
        [mean_kappa_Te(end+1),idx_best] = max(pool_kappa);
        mean_BER_Te(end+1) = pool_BER(idx_best);
        feat_vec(idx_best) = [];
        xaxis4BER(end+1) = length(feat_vec);
        feat_mat = [feat_mat [feat_vec; zeros(dim_ini-length(feat_vec),1)]];
        
    end
               
end







%% forward selection based on best kappa/RMSE
if 1
    k=0;
    mean_BER_Te = [];
    mean_kappa_Te = [];
    mean_RMSE_Te = [];
    mean_softOA_Te = [];
    mean_OA_Te = [];
    xaxis4BER = [];
    feat_vec_pool = feat_vec;
    feat_vec_ini = feat_vec;
    feat_vec_current = [];
    dim_stop = 50;
    feat_mat = zeros(dim_stop,1);
    
    while k<dim_stop

        k = k+1;
        fprintf('Iteration %u : selecting the %uth feature...\n\n',k,k);
        kappa_test = [];
        BER_test = [];
        RMSE_test = [];
        softOA_test = [];
        OA_test = [];

        parfor i=1:length(feat_vec_pool)

            fprintf('Trying adding feature %u...',feat_vec_pool(i));
            feat_vec_test = [feat_vec_current; feat_vec_pool(i)];
            feat_vec_test = sort(feat_vec_test,'ascend');
            Xtest = X(:,feat_vec_test);
            % best params for geometry, tol = 10 {0.0001,1,1000,0,10000}
            % out = CrossValidation('logistic','multiclass',{0.0001,1,1000,0,10000},4,1,Xtest,y,0,0,0,1);
            % best params for riming, tol = 10 {0.0001,0.01,1000,0,10000};
            % best params for melting, tol = 0.1 {{0.001,1,1000,0,10000}
            out = CrossValidation('logistic','binary',{0.001,1,1000,0,10000},4,1,Xtest,y,0,0,0,1);
            kappa_test(i) = out.kappa_Te;
            BER_test(i) = out.BER_Te;
            RMSE_test(i) = out.RMSE_Te;
            softOA_test(i) = out.softOA_Te;
            OA_test(i) = out.OA_Te;
            fprintf('   Done!\n');

        end
        
        if strcmp(classif_subject,'geometry');
        
            [val,idx] = max(kappa_test);% max(kappa_test);
            feat_vec_current = [feat_vec_current; feat_vec_pool(idx)];
            feat_vec_pool(idx) = [];
            mean_BER_Te(end+1) = BER_test(idx);
            mean_kappa_Te(end+1) = kappa_test(idx);
            mean_RMSE_Te(end+1) = RMSE_test(idx);
            mean_softOA_Te(end+1) = softOA_test(idx);
            mean_OA_Te(end+1) = OA_test(idx);
            fprintf('Feature %u was added ! New HSS : %2.2f%%\n',feat_vec_current(end),val*100);
        
        elseif strcmp(classif_subject,'riming');
            
            [val,idx] = min(RMSE_test);
            feat_vec_current = [feat_vec_current; feat_vec_pool(idx)];
            feat_vec_pool(idx) = [];
            mean_BER_Te(end+1) = BER_test(idx);
            mean_kappa_Te(end+1) = kappa_test(idx);
            mean_RMSE_Te(end+1) = RMSE_test(idx);
            mean_softOA_Te(end+1) = softOA_test(idx);
            mean_OA_Te(end+1) = OA_test(idx);
            fprintf('Feature %u was added ! New RMSE : %2.2f\n',feat_vec_current(end),val);
            
        elseif strcmp(classif_subject,'melting');
            
            [val,idx] = max(kappa_test);
            feat_vec_current = [feat_vec_current; feat_vec_pool(idx)];
            feat_vec_pool(idx) = [];
            mean_BER_Te(end+1) = BER_test(idx);
            mean_kappa_Te(end+1) = kappa_test(idx);
            %mean_RMSE_Te(end+1) = RMSE_test(idx);
            mean_softOA_Te(end+1) = softOA_test(idx);
            mean_OA_Te(end+1) = OA_test(idx);
            fprintf('Feature %u was added ! New OA : %2.2f%%\n',feat_vec_current(end),val*100);
            
        end
        
        feat_mat = [feat_mat [feat_vec_current; zeros(dim_stop-length(feat_vec_current),1)]];
        
    end
    
    xaxis4BER = [1:1:length(mean_BER_Te)];
    
end
    



%% backward elimination based on "merit"
if 0
    k=0;
    mean_BER_Te = [];
    mean_kappa_Te = [];
    xaxis4BER = [];
    feat_mat = feat_vec;
    beta_vec = [];
    merit_evo = [];
    feat_evo = [];
    method_ff = 'MDL';
    method_fc = 'MDL';    %'centroidBER';

    while length(feat_vec) > 20

        fprintf('%u\n',length(feat_vec));
        k=k+1; 
        merit_vec = [];

        if k==1
            old_merit = compute_correlation_merit(X,y,method_fc,method_ff);
            merit_evo(end+1) = old_merit;
            feat_evo(end+1) = length(feat_vec);
        end

        % loop over the remaining features : compute merit
        parfor i=1:length(feat_vec)
            X_test = X;
            X_test(:,i) = [];
            merit_vec(i) = compute_correlation_merit(X_test,y,method_fc,method_ff); % Pearson Kendall Spearman
        end

        % remove the feature which is leading to the highest merit
        [~,idx2remove] = max(merit_vec);
        X(:,idx2remove) = [];
        feat_vec(idx2remove) = [];
        old_merit = merit_vec(idx2remove);
        merit_vec(idx2remove) = [];
        merit_evo(end+1) = old_merit;
        feat_evo(end+1) = length(feat_vec);

        % every X iterations, we run a classification
        if mod(length(feat_vec),10)==0
            [~,~,mean_BER_Te(end+1),mean_kappa_Te(end+1)] = main_multinomial_features_optimization(feat_vec);
            xaxis4BER(end+1) = length(feat_vec);
        end


    end

    figure; 
    subplot(211); hold on; grid on; box on;
    stem(feat_vec,merit_vec-merit_evo(end-1),':diamondr');
    xlabel('feature ID');
    ylabel('merit diff');
    subplot(212); hold on; grid on; box on;
    plot(feat_evo,merit_evo,'o-');
    xlabel('# features');
    ylabel('merit');


    figure;
    subplot(211); hold on; grid on; box on;
    plot(xaxis4BER,mean_BER_Te,'k.-');
    subplot(212); hold on; grid on; box on;
    plot(xaxis4BER,mean_kappa_Te,'k.-');


end
    

        
        
        
        
        
        
        
    