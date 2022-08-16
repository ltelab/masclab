clearvars; close all;

                % alpha = 0.0001, lambda = 0.01, mini-batch size = 100
                % momentum = 0.9, Max no of iterations = 5000
                %parameters={0.00001,0.01,100,0.9,5000}; 

% USER DEFINED PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K=4; % Number of subsamples in the K-fold cross-validation
N_it=10; % Number of iterations of K-fold cross-validation
normalization_type = 'standardization';
method='logistic'; % can be either 'logistic', 'naive_bayes,'nn','svm','random_forest' or 'baseline'
parameters_method =  {0.0001,0.01,1000,0,10000};
%parameters_method = {1,1,'linear'};
% best for svm before Hu : parameters_method =  {10000,0.0001,'rbf'};
% default for logistic : {0.001,0,779,0.9,20000}; 
% default for svm-binary (6 classes) : {100,0.1,'rbf'} or {1,1,'linear'}
type_classif='multiclass'; % can be binary or multiclass
save_log=false; % Can be true of false

continuous_riming = true;
features_ranking = false;
display_confmat = true;
display_RMSE_histo = true;
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
%for geometry
%feat_vec = [1 2 3 4 5 9 10 16 19 25 28 29 34 37 49 50 53 54 55 59 63 64 65 70 75 83 84 87]; 
%for riming
%feat_vec = [1 3 27 28 30 34 36 37 53 54 55 56 58 62 63 64 65 74 76 83 84 87];
%all except absurd
%feat_vec = [1:1:96]';
%feat_vec([32 33 37 38 41:48 67 77 79 80]) = [];
% selection of 30 from feat. ranking analysis
%feat_vec = [4 5 7 9 15 22 24 25 28 29 39 53 54 56 58 60 61 62 71 74 75 78 84 85 86 87 88 93 94 95];

Dwanted = 25;
F=load('out_riming_10it_rand_FS_CV4_RMSE_D80.mat');
feat_vec = F.feat_mat(:,Dwanted+1);
feat_vec = feat_vec(feat_vec~=0);
feat_vec = sort(feat_vec,'ascend');

% compute new descriptors (if necessary)
%process_new_descriptors(dir_data,t_str_start,t_str_stop)

% load the training matrix X
[X,Xlab,Xname,Xt] = load_processed_data(dir_data,t_str_start,t_str_stop,feat_vec);

% load the labels vector y
[y2,y,~] = load_labels(dir_data,dir_labels,labels_filename,'1.1');

% remove "undetermined", "not labelled yet" or "ambiguous" labels
idx_unknown = find(y<=0);
y(idx_unknown) = [];
X(idx_unknown,:) = [];
data_picnames(idx_unknown) = [];
data_filenames(idx_unknown) = [];
fprintf('%u samples discarded because riming labelled as "undetermined" or "ambiguous" \n',length(idx_unknown));

N = size(X,1);
D = size(X,2);
N_classes = length(unique(y));
%legends = {'none','moderately rimed','heavily rimed','graupel-like','graupel'};
legends = {'1','2','3','4','5'};
labels = legends;

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

%%
    
% Initialize global confusion matrix
global_confmat_Tr = zeros(N_classes);

% Initialize beta ranking vectors
beta_rankvec_sum = [];
beta_rankvec_mean = [];

% Initialize error vectors 
err_Tr_vec = [];


yTr = y;
XTr = X; 
    
% Fit model using the wrapping function
model = trainClassifier(method, type_classif, yTr, XTr, parameters_method);
          
% Features ranking
if features_ranking
    beta_sum = sum(abs(model),2);
    beta_mean = mean(abs(model),2);
    [~,beta_rank_sum] = sort(beta_sum,'descend');
    [~,beta_rank_mean] = sort(beta_mean,'descend');
    beta_rankvec_sum = [beta_rankvec_sum beta_rank_sum];
    beta_rankvec_mean = [beta_rankvec_mean beta_rank_mean];
end

% Predict
[predTr,scoresTr] = predictClassifier(method,type_classif,model,XTr); 
        
if continuous_riming
    scoreTr = scoresTr(:,1)*1 + scoresTr(:,2)*2 + scoresTr(:,3)*3 + scoresTr(:,4)*4 + scoresTr(:,5)*5;
end

% Compute accuracy
BER_Tr = computeBER(yTr,predTr);
OA_Tr = computeOA(yTr,predTr);
kappa_Tr = computeKAPPA(yTr,predTr);
% for riming
softOA_Tr = compute_softOA_riming(yTr,predTr);
if strcmp(method,'logistic')
    RMSE_Tr = computeRMSE(yTr,scoreTr);
else
    RMSE_Tr = -1;            
end

% Update global confusion matrix
global_confmat_Tr = confusionmat(yTr,predTr);
        
% update vector of errors
if strcmp(method,'logistic')
    err_Tr_vec = scoreTr - yTr;
end
        
     
fprintf('Training BER: %.2f%%\n\n', BER_Tr * 100 );
fprintf('Training kappa: %.2f%%\n\n', kappa_Tr * 100 );
fprintf('Training soft OA: %.2f%%\n\n', softOA_Tr * 100 );
fprintf('Training RMSE: %.2f\n\n', RMSE_Tr );


% save stuff in a structure
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
classifier.softOA = softOA_Tr;
classifier.RMSE_Tr = RMSE_Tr;
classifier.feat_vec = feat_vec;

%% illustration

if display_confmat
    figure;
    heatmap(global_confmat_Tr,legends,legends,1,'Colormap','red','ShowAllTicks',1,'UseLogColorMap',true,'Colorbar',true,'Fontsize',14);
    set(gca,'Fontsize',14);
    title('Train confusion matrix');
end

% error histograms
if display_RMSE_histo
    figure;
    hold on; box on;
    histogram(err_Tr_vec);
    set(gca,'Ytick',[]);
    set(gca,'Xlim',[-2 2]);
    xlabel('True label - Prediction');
    set(gca,'Fontsize',14);
    title('Error distribution');
end

% boxplots
scatt = 0.15;
scatt_vec = -scatt +(scatt+scatt).*rand(length(yTr),1);
yTr_scatt = yTr + scatt_vec;

% training set boxplot
pause(0.5);
grey_c = [89 89 89]/255;
red_c = [240 59 32]/255;
figure; hold on; box on; grid on;
% plot all the datapoints
plot(yTr_scatt,scoreTr,'k.','markersize',5,'MarkerEdgeColor',grey_c);
% generate the boxplot
bh = boxplot(scoreTr,yTr,'symbol','','positions',1:1:5,'widths',0.5,'colors','k');
% fill the boxes with color
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),'b','FaceAlpha',.5);
end
% redraw the boxplot over it
bh = boxplot(scoreTr,yTr,'symbol','','positions',1:1:5,'widths',0.5,'colors','k');
% thicken the outline of the boxplot
set(bh(:,:),'linewidth',2);
% set the color of the median line
set(bh(6,:),'linewidth',2,'Color',red_c);
% set the axis
axis([0.5 5.5 0.5 5.5]);
set(gca,'YTick',[1 2 3 4 5]);
set(gca,'XTick',[1 2 3 4 5]);
xlabel('true riming degrees (discrete)');
ylabel('predicted riming degree (continuous 1~5)');
set(gca,'fontsize',12);
title('Train sample');

