function model = SVM1vs1(y,X,C,gamma,kernel)
%SVM1vsall    : train a SVM model for multiclass classification by using
%                 the 1-vs-all methodology
%   Inputs:
%     y       : the vector of labels
%     X       : the matrix of features
%     C       : the margin of the SVM
%     gamma   : the parameter of the kernel (not used for linear kernels)
%     kernel  : the type of kernel (can be linear, rbf or lorentz)

%   Outputs:
%      model  :  set of trained SVM models consisting of their support 
%                vectors and their parameters, to be used with 
%                predict_classifier

% Credit
%--------
% Christophe Praz & Daniel Wolfensberger, December 2015

    % convert y into doubles
    y = double(y);

    % check number of classes
    n_classes = length(unique(y));
    
    % double loop over each class, keep only concerned labels
    % For 1vs1, 1 kernel has to computed for each 1vs1 "confrontation"
    y_ref = y;
    X_ref = X;
    % increment determining the iteration
    it = 1;
    for i=1:(n_classes-1)
        for j=(i+1):n_classes

            fprintf('compute 1vs1 for class %u vs %u .... \n',i,j);
            
            % retrieve indexes belonging to class i or j and store it in
            % order to use it for prediction (computing appropriate X_test)
            model.idx{it} = find(y_ref==i | y_ref==j);
            y = y_ref(model.idx{it});
            y(y~=i) = -1;
            y(y==i) = +1;
            X = X_ref(model.idx{it},:);

            % compute kernel
            if strcmp(kernel,'linear')        
                K = linear_kernel(X,X);
            elseif strcmp(kernel,'rbf')
                K = rbf_kernel(X,X,gamma);
            elseif strcmp(kernel,'lorentz');
                K = lorentz:kernel(X,X,gamma);
            else
                fprintf('Error : kernel unkown \n');
            end       

            % compute model parameters
            [alphas, model.beta0{it}] = SMO(K, y, C);

            % keep only parameters correspondig to non-zero alpha entries (support
            % vectors)
            SV_inds = find(alphas>0);
            model.X_SV{it} = X(SV_inds,:);
            model.y_SV{it} = y(SV_inds,:);
            model.alphas_SV{it} = alphas(SV_inds);
            
            it = it + 1;
     
        end
        
    end
    
    model.kernel = kernel;
    model.gamma = gamma;

end