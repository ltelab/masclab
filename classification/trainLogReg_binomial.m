function [beta,converged] = trainLogReg_binomial(y,X,alpha,lambda,batchSize,momentum, maxIters,weights,beta,tol)
%trainLogReg_binomial Trains (penalized) binomial logistic regression
%  Inputs:
%  y          : the vector of labels
%  X          : the set of features (matrix size NxD)
%  alpha      : the learning rate in the gradient descent method
%  lambda     : the penalization parameter
%  momentum   : the momentum (inertia) used in the gradient descent
%  maxIters   : the  maximum number of iterations (default = 1000)
%  beta       : the initial guess of the beta parameters (default = zeros)
%  tol        : the tolerance in the convergence criterion (default = 1e-6)


%   Outputs:
%      beta      : Optimized vector of parameters beta
%      converged : A flag that tells if the gradient descent has converged


% Credit
%--------
% Christophe Praz & Daniel Wolfensberger, December 2015

    % Get tX matrix
    tX=[ones(length(y),1),X];
    
    % Specify default tolerance if needed
    if nargin < 10 
        tol = 0.01;       
    end
    
    % Specify default beta0 if needed (vector of zeros)
    if nargin < 9
        beta = zeros(size(tX,2),1);
    end
    
    % Specify default weights in the logistic cost
    if nargin < 8
        weights = ones(2,1);
    end

    % Specify default maxIters if needed
    if nargin < 7
        maxIters=10000;
    end
    
    L_old = 0; % Initalize cost    
    step=0; % Initialize step in gradient descent
    converged = 0; % Initialize convergence flag
    
    for k = 1:maxIters % Until max number of iterations is reached
        
        % Compute gradient
        g=calculateBatchGradientPenLogistic(y,tX,beta,lambda,batchSize,weights); 

        step=momentum*step+alpha*g; % Step of gradient descent
        
        beta = beta + step; % Update beta with gradient descent
        if mod(k,100)==0 % To sav time we compute the cost only every 100 iterations
            L = computeCostPenLogistic(y,tX,beta,lambda,weights); % Compute cost
            %disp(L);
        
            if abs(L-L_old) < tol % Check for convergence
                converged = 1;
                fprintf('convergence threshold has been reached after %u iterations ! \n',k);
                break;
            end

             L_old = L; % Update cost
        end
    end
   
    
end

% Method to compute the gradient of logistic regression
function g = calculateBatchGradientPenLogistic(y,tX,beta, lambda, batchSize,weights)
    if nargin<6
        weights = ones(2,1);
    end
    if batchSize > length(y)
        batchSize = length(y);
    end
    perm=randperm(length(y),batchSize);  % Take mini-batch
    %g = -(tX(perm,:)'*(sigmoid(tX(perm,:)*beta) - y(perm,:))+lambda*[0;beta(2:end)]);
    g = -(tX(perm,:)'*(weights(1)*sigmoid(tX(perm,:)*beta) -weights(2)*y(perm,:) -(weights(1)-weights(2))*y(perm,:).*sigmoid(tX(perm,:)*beta))+lambda*[0;beta(2:end)]);
end

function sigma = sigmoid(x)
    for i=1:length(x)  % To guarantee numerical stability choose one of two
        % equivalent functions depending on sign of x.
        if sign(x(i)) == -1 
            sigma(i,1) = exp(x(i))/(1+exp(x(i)));
        else
            sigma(i,1) = 1/(1+exp(-x(i)));
        end
    end   
end

function L = computeCostPenLogistic(y,tX,beta,lambda,weights) % Compute cost of logistic regression
    if nargin<5
        weights = ones(2,1);
    end
    L = 0;
    for i=1:length(y)
        %L = L + y(i)*tX(i,:)*beta - log(1 + exp(tX(i,:)*beta));
        L = L + weights(2)*y(i)*tX(i,:)*beta - weights(1)*log(1 + exp(tX(i,:)*beta)) + (weights(1)-weights(2))*y(i)*log(1 + exp(tX(i,:)*beta));
    end 
    L = -L +lambda*sum(beta.^2) ; % Take negative of log-likelihood and add regularization
end



