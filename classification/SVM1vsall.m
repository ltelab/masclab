function model = SVM1vsall(y,X,C,gamma,kernel)
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

    % create kernel (kernel remains the same for each instance of 1vsall)
    if strcmp(kernel,'linear')        
        K = linear_kernel(X,X);
    elseif strcmp(kernel,'rbf')
        K = rbf_kernel(X,X,gamma);
    elseif strcmp(kernel,'lorentz');
        K = lorentz:kernel(X,X,gamma);
    else
        fprintf('Error : kernel unkown \n');
    end

    % check number of classes
    n_classes = length(unique(y));
    
    % loop over each class and apply "one vs all" rule
    y_ref = y;
    for i=1:n_classes
        
        y(y_ref~=i) = -1;
        y(y_ref==i) = +1;
        
        % compute model parameters
        [alphas, model.beta0{i}] = SMO(K, y, C);

        % keep only parameters correspondig to non-zero alpha entries (support
        % vectors)
        SV_inds = find(alphas>0);
        model.X_SV{i} = X(SV_inds,:);
        model.y_SV{i} = y(SV_inds,:);
        model.alphas_SV{i} = alphas(SV_inds);
        
        %fprintf('compute 1vsall for class %u .... \n',i);
        
    end
    
    model.kernel = kernel;
    model.gamma = gamma;

end