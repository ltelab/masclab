function kappa = computeKAPPA( real, pred )
% compute Cohen's Kappa

    p_obs = sum(real==pred)/length(real);
    cm = confusionmat(real,pred);
    user_acc = sum(cm,1);
    prod_acc = sum(cm,2);
    p_est = user_acc * prod_acc/(sum(cm(:))^2);
    kappa = (p_obs - p_est)/(1 - p_est);
    
end