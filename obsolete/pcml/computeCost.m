function L = computeCost(y,tX,beta)
    
    e = y - tX*beta;
    L = 1/(2*length(y)) * e' * e;


end