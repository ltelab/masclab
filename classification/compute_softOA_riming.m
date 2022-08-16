function softOA = compute_softOA_riming( real, pred )
% compute a soft version of OA for rimng degree

    softOA = (sum(pred==real) + sum(pred-1==real) + sum(pred+1==real))/length(real);

end