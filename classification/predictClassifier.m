function [pred_class, scores] = predictClassifier( method, type_classif, fitted_model, X )
%predict_classifier : classifies data points by using a trained
%                     classification model
%   Inputs:
%     method       : the classification method, can be either
%                   ('svm','baseline','nn','random_forest','naive_bayes' 
%                    'or 'logistic')
%     type_classif : 'binary' or 'multiclass'
%     fitted_model : the trained model, output of train_classifier
%     X            : the set of features (matrix size NxD)

%   Outputs:
%      pred_class  : the predicted classes that correspond to input X
%      scores      : the scores that correspond to input X, the meaning of
%                    these scores depend on the model but are usually
%                    probabilities to be in one class or the other

% Credit
%--------
% Christophe Praz & Daniel Wolfensberger, December 2015


switch method
    case 'nn' 
        if strcmp(type_classif,'multiclass')
            fitted_model.testing = 1;
            fitted_model = nnff(fitted_model, X, zeros(size(X,1), fitted_model.size(end)));
            fitted_model.testing = 0;
            % predict on the test set
            scores = fitted_model.a{end};
            % get the most likely class
            [~,pred_class] = max(scores,[],2);
        elseif strcmp(type_classif,'binary')
            fitted_model.testing = 1;
            fitted_model = nnff(fitted_model, X, zeros(size(X,1), fitted_model.size(end)));
            fitted_model.testing = 0;
            % predict on the test set
            scores=fitted_model.a{end};
            pred_class = round(scores);
        end
    case 'logistic'
        if strcmp(type_classif,'multiclass')
            tX=[ones(size(X,1),1),X];
            scores=exp(tX*fitted_model);
            scores=bsxfun (@rdivide, scores, sum(scores,2)); % Divide by the sum of exp to normalize probabilities to 1
            [~,pred_class] = max(scores,[],2);  
        elseif strcmp(type_classif,'binary')
            tX=[ones(size(X,1),1),X];
            scores=sigmoid(tX*fitted_model);
            pred_class=round(scores);
        end
    case 'random_forest'
        [pred_class, scores]=predict(fitted_model,X);
        pred_class=cellfun(@str2num,pred_class);
    case 'baseline'
        % Compute distances from all points to centroids
        dist=pdist2(X,fitted_model);
        
        [~,pred_class] = min(dist,[],2);
        if strcmp(type_classif,'binary')
            % We defined scores as the ratio of distances between class 1 and 2
            % for the baseline
            scores=dist(:,1)./dist(:,2);
            pred_class=pred_class-1; % Get from matlab index (1,2) to class index (0,1)
        else
            scores=-1; % We haven't defined any scores for the multiclass case
        end
    case 'naive_bayes'
        pred_class=predict(fitted_model, X);
        scores=posterior(fitted_model,X);
    case 'svm'
        if strcmp(type_classif,'binary')
            % compute prediction kernel
            if strcmp(fitted_model.kernel,'linear')
                K_pred = linear_kernel(X,fitted_model.X_SV);
            elseif strcmp(fitted_model.kernel,'rbf')
                K_pred = rbf_kernel(X,fitted_model.X_SV,fitted_model.gamma);
            elseif strcmp(fitted_model.kernel,'lorentz')
                K_pred = lorentz_kernel(X,fitted_model.X_SV,fitted_model.gamma);
            end
            % compute scores
            scores = K_pred*(fitted_model.alphas_SV.*fitted_model.y_SV) + fitted_model.beta0;
            % fit into {-1;1} 
            pred_class=zeros(length(scores),1);
            pred_class(scores >= 0) = 1;
            pred_class(scores <= 0) = 0;  
        elseif strcmp(type_classif,'multiclass') % we use 1vsall here, see commented code for 1vs1 at the bottom
            % for loop over the classes to compute score
            for i=1:length(fitted_model.y_SV)
                % this time, kernel has to be recomputed for each class
                if strcmp(fitted_model.kernel,'linear')
                    K_pred = linear_kernel(X,fitted_model.X_SV{i});
                elseif strcmp(fitted_model.kernel,'rbf')
                    K_pred = rbf_kernel(X,fitted_model.X_SV{i},fitted_model.gamma);
                elseif strcmp(fitted_model.kernel,'lorentz')
                    K_pred = lorentz_kernel(X,fitted_model.X_SV(i),fitted_model.gamma);
                end
                scores(:,i) = K_pred*(fitted_model.alphas_SV{i}.*fitted_model.y_SV{i}) + fitted_model.beta0{i};       
            end
            % get the most likely class (= higher score value)
            [~,pred_class] = max(scores,[],2);
        end
end


end

function sigma = sigmoid(x)
    for i=1:length(x)
        if sign(x(i)) == -1
            sigma(i,1) = exp(x(i))/(1+exp(x(i)));
        else
            sigma(i,1) = 1/(1+exp(-x(i)));
        end
    end
end


% code for SVM 1vs1 classification (obsolete)
%
%         it = 1;
%         %X_ref = X;
%         % initialisation of a matrix Nx4 of classvotes
%         classvote = zeros(size(X,1),4);
%         % for loop over each instance of 1vs1 (i vs j)
%         for i=1:3
%             for j=(i+1):4
%                 % retrieve corresponding X
%                 %X = X_ref(fitted_model.idx{it},:);
%                 % compute kernel
%                 if strcmp(fitted_model.kernel,'linear')
%                     K_pred = linear_kernel(X,fitted_model.X_SV{it});
%                 elseif strcmp(fitted_model.kernel,'rbf')
%                     K_pred = rbf_kernel(X,fitted_model.X_SV{it},fitted_model.gamma);
%                 elseif strcmp(fitted_model.kernel,'lorentz')
%                     K_pred = lorentz_kernel(X,fitted_model.X_SV{it},fitted_model.gamma);
%                 end
%                 scores = K_pred*(fitted_model.alphas_SV{it}.*fitted_model.y_SV{it}) + fitted_model.beta0{it};
%                 %scores_extended = zeros(size(X,1),1);
%                 %scores_extended(scores_extended(fitted_model.idx{it})) = scores;
%                 %idx_iplus = find(scores_extended >= 0);
%                 classvote(scores>=0,i) = classvote(scores>=0,i) + 1;
%                 %idx_jplus = find(scores_extended < 0);
%                 classvote(scores<0,j) = classvote(scores<0,j) + 1;
%                 
%                 it = it + 1;
%                 
%             end
  