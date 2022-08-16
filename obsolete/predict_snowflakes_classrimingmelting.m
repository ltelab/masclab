% % function to predict the snowflake class, riming degree and melting snow
% % at the same time.
% % WARNING: automatically assumes that the normalization
% function predict_snowflakes_classrimingmelting(dir_data,t_str_start,t_str_stop,classifier,classifer_riming,classifier_melting,save_results,blurry_threshold)
% 
%     if nargin == 7
%         blurry_threshold = 0;
%     end
%     
%     load(classifier);
%     normalization_type = classifier.normalization;
%     prediction_scheme = classifier.N_labels;
% 
%     % exemple
%     % dir_data = '/home/praz/Documents/MASC/MASC_labelling/mixed_sample_part2';
%     % t_str_start = '20150101000000';
%     % t_str_stop = '20170101000000';
%     % classifier = 'logit_test.mat';
% 
% 
%     fprintf('Classifying snowflakes type in : %s...\n',dir_data);
%     %% load flakes descriptors
% 
%     % process new descriptors (symetry features, ...) if necessary
%     process_new_descriptors(dir_data,t_str_start,t_str_stop);
% 
%     % load the training matrix X
%     [X,~,Xname,~] = load_processed_data(dir_data,t_str_start,t_str_stop);
% 
%     %% Features transformation
%     % according to classif_params
% 
%     fprintf('Assembling processed data in X matrix and transforming features...');
% 
%     D = size(X,2);
% 
%     for i=1:D
% 
%         if classifier.normalization_params.skew(i) > 1
% 
%             X(:,i) = log(abs(X(:,i)+1));
% 
%         elseif classifier.normalization_params.skew(i) > 0.75
% 
%             X(:,i) = sqrt(abs(X(:,i)));
% 
%         elseif classifier.normalization_params.skew(i) < -1
% 
%             X(:,i) = exp(X(:,i));
% 
%         elseif classifier.normalization_params.skew(i) < -0.75
% 
%             X(:,i) = X(:,i).^2;
% 
%         end
% 
%         if strcmp(normalization_type,'standardization') 
% 
%             X(:,i) = (X(:,i) - classifier.normalization_params.mean(i))/classifier.normalization_params.std(i);
% 
%         elseif strcmp(normalization_type,'rescaling')
% 
%             disp('warning: rescaling not implemented yet!');
% 
%         end
% 
%     end
% 
%     fprintf('  Done!\n'); 
% 
%     %% Prediction
%     fprintf('Classifying images according to %s classifier...',classifier.type);  
%     [pred,scores] = predictClassifier(classifier.type,classifier.type_classif,classifier.model,X);    
%     fprintf('   Done!\n');
% 
%     %% saving
%     if save_results
% 
%         fprintf('Saving results in roi structures...');
%         N_fine = 0;
%         N_blurry = 0;
% 
%         for i=1:numel(Xname)
%             load(fullfile(dir_data,Xname{i}));
%             if roi.xhi > blurry_threshold
%                 roi.label_ID = pred(i);
%                 roi.label_name = prediction_scheme{pred(i)};
%                 roi.label_probs = scores(i,:);
%                 N_fine = N_fine + 1;
%             else
%                 roi.label_ID = -9;
%                 roi.label_name = 'blurry';
%                 roi.label_probs = scores(i,:);
%                 N_blurry = N_blurry + 1;
%             end
% 
%             save(fullfile(dir_data,Xname{i}),'roi');
%         end
% 
%         N_tot = N_fine + N_blurry;
%         fprintf('   Done!\n');
%         fprintf('***** %u images classified ***** %u (%u%%) good classifications ***** %u (%u%%) blurry images *****\n',N_tot,N_fine,round(N_fine/N_tot*100),N_blurry,round(N_blurry/N_tot*100));
% 
%     end
% 
% end
% 
% 
