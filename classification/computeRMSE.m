function rmse = computeRMSE( real, pred )
% compute Root Mean Square Error of prediction 
% (valid metric for continuous data)

    rmse = sqrt(1/length(real) * sum((real-pred).^2));

end

