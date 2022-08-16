function OA = computeOA( real, pred )
% compute Overall Accuracy of Classification

    OA = sum(real==pred)/length(real);

end