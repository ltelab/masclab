% calculate "scale invariant" central moments
function eta_pq = normalized_central_moments(im,xnorm,ynorm,p,q)
    
    mu_pq = sum(sum((xnorm.^p).*(ynorm.^q).*im));
    mu_00 = sum(sum(im)); %mu_00
    % normalise moments for scale invariance
    eta_pq = mu_pq/(mu_00^(1+(p+q)/2));
    
end