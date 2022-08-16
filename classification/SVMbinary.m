function model = SVMbinary(y,X,C,gamma,kernel)
%SVM1binary   : train a SVM model for binary classification
%
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

    % create kernel
    if strcmp(kernel,'linear')        
        K = linear_kernel(X,X);
    elseif strcmp(kernel,'rbf')
        K = rbf_kernel(X,X,gamma);
    elseif strcmp(kernel,'lorentz')
        K = lorentz_kernel(X,X,gamma);
    else
        fprintf('Error : kernel unkown \n');
    end

    % compute model parameters
    [alphas, model.beta0] = SMO(K, y, C);

    % keep only parameters correspondig to non-zero alpha entries (support
    % vectors)
    SV_inds = find(alphas>0);
    model.X_SV = X(SV_inds,:);
    model.y_SV = y(SV_inds,:);
    model.alphas_SV = alphas(SV_inds);
    model.kernel = kernel;
    model.gamma = gamma;

end