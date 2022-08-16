function [ model ] = trainClassifier( method, type_classif, y, X, parameters )
%train_classifier Trains a classification algorithm using vector of labels
%y and set of features X. 

%   Inputs:
%     method       : the classification method, can be either
%                   ('svm','baseline','nn','random_forest','naive_bayes' 
%                    'or 'logistic')
%     type_classif : 'binary' or 'multiclass'
%     y            : the vector of labels  (vector size Nx1)
%     X            : the set of features (matrix size NxD)
%     parameters   : (Optional) parameters for the classifier, given as a 
%                    cell array, if not provided, the best ones we found 
%                    will be used

%   Outputs:
%           model  : the fitted model that can be used as input to 
%                    the predict_classifier function

% Credit
%--------
% Christophe Praz & Daniel Wolfensberger, December 2015

    take_defaults = isempty(parameters); % Check if parameters have been defined
    
    switch method
        case 'nn'
            % Train neural network
            if take_defaults
                parameters={[200],50,size(X,2)}; % hid.lay + hid.neurons, epochs, samples mini-batch
            end
            model=trainNN(y,X,parameters{:});
            
        case 'random_forest'
            % Train 100 decision trees
            if take_defaults
                % 100 trees, take all variables
                parameters={100,'NVarToSample', (size(X,2))}; 
            end
            model=TreeBagger(parameters{1},X,y,parameters{2:end});
        case 'logistic'
            % Train either penalized binomial or multinomial logistic
            % regression,
            
            if take_defaults
                % alpha = 0.0001, lambda = 0.01, mini-batch size = 100
                % momentum = 0.9, Max no of iterations = 5000
                % parameters={0.00001,0.01,100,0.9,5000}; 
            end
            if strcmp(type_classif, 'binary')
                model=trainLogReg_binomial(y,X,parameters{:});
            elseif strcmp(type_classif, 'multiclass')
 
                model=trainLogReg_multinomial(y,X,parameters{:});
            end
        case 'baseline'
            % No parameters
            % Model consists simply of the centroids of every class
            classes=unique(y);
            for i=1:length(classes)
                centroids(i,:)=median(X(y==classes(i),:));
            end
            model=centroids;
        case 'naive_bayes'
            % No parameters (we could not use any other distribution)
            model=fitNaiveBayes(X,y);
        case 'svm'
            if take_defaults
                parameters={100,1e-6,'rbf'};
            end
            if strcmp(type_classif, 'binary')
                % Define y as [-1,1] and as double
                y=double(y);
                y(y==0) = -1;
                model=SVMbinary(y,X,parameters{:}); 
            else
                model=SVM1vsall(y,X,parameters{:});
            end         
    end

end

