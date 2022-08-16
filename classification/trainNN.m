function [ nnOutput ] = trainNN( y,X, nNeurons,nEpochs,nBatchSize)
%trainNN Trains a neural network using the DeepLearning Toolbox
%   Inputs:
%  y          : the vector of labels
%  X          : the set of features (matrix size NxD)
%  nNeurons   : number of neurons in every hidden layer as a vector
%  nBatchSize : size of the mini-batch (in samples) used for the backpropa
%               gation algorithm

%   Outputs:
%      nnOutput    : A trained neural network


% Credit
%--------
% Christophe Praz & Daniel Wolfensberger, December 2015

    nClassesOutputs=length(unique(y));
    if nClassesOutputs==2 % Binary case
        nClassesOutputs=1;  
    end

    nn = nnsetup([size(X,2) nNeurons nClassesOutputs]);
    opts.numepochs =  nEpochs;   %  Number of full sweeps through data
    opts.batchsize = nBatchSize;  %  Take a mean gradient step over this many samples

    % if == 1 => plots trainin error as the NN is trained
    opts.plot               = 0;

    nn.learningRate = 2;

    % this neural network implementation requires number of samples to be a
    % multiple of batchsize, so we remove some for this to be true.
    numSampToUse = opts.batchsize * floor( size(X) / opts.batchsize);
    X = X(1:numSampToUse,:);
    y = y(1:numSampToUse);

    % prepare labels for NN
    if nClassesOutputs>1
        y=dummyvar(double(y));
    end
    [nnOutput, L] = nntrain(nn, X, y, opts);
    
end

