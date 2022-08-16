function [beta,converged] = trainLogReg_multinomial(y,X,alpha,lambda, batchSize,momentum,maxIters,weights,beta,tol)
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
    
    % get classes
    classes=unique(y);

    % Specify default tolerance if needed
    if nargin < 10  
        tol = 10;
        % 0.1 : last value for snow habit
        % 10 : former value for snow habit
    end
    
    % Specify default beta0 if needed (vector of zeros)
    if nargin < 9
        beta = zeros(size(tX,2),length(classes));
    end

    % Specify default weights in the logistic cost
    if nargin < 8
        weights = ones(length(classes),1);
    end
        
    % Specify default maxIters if needed
    if nargin < 7
        maxIters=10000;
    end

    L_old = 0; % Initalize cost    
    step=0; % Initialize step in gradient descent
    converged = 0; % Initialize convergence flag
    itCost = 100;
    
    for k = 1:maxIters % Until max number of iterations is reached
        
        % Compute gradient
        g=calculateBatchGradientPenLogistic(y,tX,beta,lambda, batchSize,classes,weights); 

        step=momentum*step - alpha*g;  % Step of gradient descent
        
        beta = beta + step; % Update beta with gradient descent
        %if k < 10000
        
%         elseif k < 1000000
%             itCost = 10;
%             batchSize = 500;
%         else 
%             itCost = 1;
%             batchSize = 100;
%         end
        if mod(k,itCost)==0 % To sav time we compute the cost only every itCost iterations
            L = computeCostPenLogistic(y,tX,beta,lambda,classes,weights); % Compute cost
            %disp(L);
            if norm(L-L_old) < tol % Check for convergence
                converged = 1;
                %fprintf('convergence threshold has been reached after %u iterations ! \n',k);
                break;
            end

             L_old = L; % Update cost
        end
    end
    
end
% 

function g = calculateBatchGradientPenLogistic(y,tX,beta, lambda, batchSize, classes, weights)
    if nargin<7
        weights=ones(length(classes),1);
    end

    if batchSize > length(y)
        batchSize = length(y);
    end
    perm=randperm(length(y),batchSize);
    tX=tX(perm,:);
    y=y(perm,:);

    sum_probs=sum(exp(tX*beta),2);
    
    g = zeros(size(tX,2),length(classes));

    for j=1:length(classes)
        % somehow the performance are just better by commenting the lambda part of the gradient -.-
        g(:,j) = (weights(j)*tX'*((exp(tX*beta(:,j))./sum_probs)-(y==classes(j)))  + 2*lambda*[0;beta(2:end,j)]);%    +lambda*sign([0;beta(2:end,j)]));%  % ) % 
%         for i=1:length(y)
%             g(:,j) = g(:,j) + tX(i,:)'*exp(tX(i,:)*beta(:,j))/sum(exp(tX(i,:)*beta));
%             if y(i) == classes(j)
%                 g(:,j) = g(:,j) - tX(i,:)';
%             end
%         end
%         g(:,j) = g(:,j) + 2*lambda*[0;beta(2:end,j)];
     end
end


function L = computeCostPenLogistic(y,tX,beta,lambda, classes, weights) % Compute cost of logistic regression
    L = 0;
    
    if nargin<6
        weights = ones(length(classes),1);
    end
    
    % first term of the cost
    for i=1:length(y)
        for j=1:length(classes)
            if(y(i)==classes(j))
                %L = L - (tX(i,:)*beta(:,j)+log(sum(exp(tX(i,:)*beta))));
                L = L - weights(j)*tX(i,:)*beta(:,j); %+ log(sum(exp(tX(i,:)*beta))); %weights*log
            end
        end
    end 
    %second term of the cost
    for i=1:length(y)
       for j=1:length(classes)
           L = L + weights(j)*log(exp(tX(i,:)*beta(:,j)));
       end
    end
    % penalization
    L = L+lambda*sum(beta(:).^2); % Take negative of log-likelihood and add regularization term
end





