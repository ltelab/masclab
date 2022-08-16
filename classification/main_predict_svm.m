% main predict
clear all; %close all;

input_type = 'folder'; %  folder,
blurry_threshold = 9;
classifier = 'fitsvm_7.59BER_noHOG.mat';
load(classifier); %model
%load('fitsvm_5.73BER.mat');
%dir_data = '/home/praz/Documents/MASC/masclab/events/dendrites'; % data directory
dir_data = '/media/praz/MyData/MASC/processed_APRES3/GOOD/2016.01.28/01';
t_str_start = '20150101000000';
t_str_stop = '20170101000000';
normalization_type = 'standardization';
save_results = true;
prediction_scheme = {'graupel','aggregate','melting snow','small particle','dendrite','column'};

fprintf('Classifying snowflakes in : %s...\n',dir_data);

%% load flakes descriptors

% process new descriptors (symetry features, ...) if necessary
process_new_descriptors(dir_data,t_str_start,t_str_stop);

% load the training matrix X
[X,Xlab,Xname,Xt] = load_processed_data(dir_data,t_str_start,t_str_stop);

% load the hog feature
%[hog,hog2,hogname] = load_hog_features(dir_data,t_str_start,t_str_stop);
%hog_tot = [hog hog2];
% retrieve princ. comp. from training PCA
%hog_tot_PCs = hog_tot * classif_params.pca_hogcoeff;

%remove empty columns
%idx_null = find(sum(hog)==0);
%hog(:,idx_null) = [];
%idx_null = find(sum(hog2)==0);
%hog2(:,idx_null) = [];

% take pca of hog features
%hog_tot = [hog hog2];
%[coeff,score,latent] = pca(hog_tot);
%[coeff2,score2,latent2] = pca(hog2);

% add some pc to the X matrix
X = [X];% hog_tot_PCs(:,1:20)];% score(:,1:20)];% score(:,1:10) score2(:,1:10)]; 
N = size(X,1);
D = size(X,2);
labels = {'graup','agg','melt','smal','dend','col'};
N_classes = 6;

%% Features transformation
% according to classif_params

fprintf('Assembling processed data in X matrix and transforming features...');

D = size(X,2);

for i=1:D
         
    if classif_params.skew(i) > 1
        
        X(:,i) = log(abs(X(:,i)+1));
        
    elseif classif_params.skew(i) > 0.75
        
        X(:,i) = sqrt(abs(X(:,i)));
        
    elseif classif_params.skew(i) < -1
        
        X(:,i) = exp(X(:,i));
        
    elseif classif_params.skew(i) < -0.75
        
        X(:,i) = X(:,i).^2;
        
    end

    if strcmp(normalization_type,'standardization') 
        
        X(:,i) = (X(:,i) - classif_params.mean(i))/classif_params.std(i);
        
    elseif strcmp(normalization_type,'rescaling')
        
        disp('warning: rescaling not implemented yet!');
        
    end
    
end

fprintf('  Done!\n'); 

%% Prediction

fprintf('Classifying images according to %s classifier...',classifier);

scores = [];
for j=1:N_classes

    [label,score] = predict(model{j},X);
    scores(:,j) = score(:,2);     

end
[~,pred] = max(scores,[],2);

fprintf('   Done!\n');

%% saving
if save_results
    
    fprintf('Saving results in roi structures...');
    N_fine = 0;
    N_blurry = 0;
    
    for i=1:numel(Xname)
        load(fullfile(dir_data,Xname{i}));
        if roi.xhi > blurry_threshold
            roi.label_ID = pred(i);
            roi.label_name = prediction_scheme{pred(i)};
            N_fine = N_fine + 1;
        else
            roi.label_ID = 0;
            roi.label_name = 'blurry';
            N_blurry = N_blurry + 1;
        end
            
        save(fullfile(dir_data,Xname{i}),'roi');
    end

    N_tot = N_fine + N_blurry;
    fprintf('   Done!\n');
    fprintf('***** %u images classified ***** %u (%u%%) good classifications ***** %u (%u%%) blurry images *****\n',N_tot,N_fine,round(N_fine/N_tot*100),N_blurry,round(N_blurry/N_tot*100));
    
end




