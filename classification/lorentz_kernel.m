function [ K ] = lorentz_kernel(X1, X2, gamma )
% Build a kernel based on Lorentz function.

    if nargin == 2 
        
        gamma = 1;
        
    end

    K = exp(-gamma.*pdist2(X1,X2));

end