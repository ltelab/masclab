% Build a gaussian kernel
function [ K ] = rbf_kernel(X1, X2, gamma )

    if nargin == 2 
        
        gamma = 1;
        
    end

    K = exp(-gamma.*(pdist2(X1,X2)).^2);

end