% LeastSquares : compute the regression coefficients vector beta using the
% normal equation.
function beta = leastSquares(y,tX)

    beta = (tX'*tX)\(tX'*y);
    % we used Matlab mldivide solver (or "\"), more efficient and robust to
    % ill-conditionning than inv(tX'*tX)*(tX'*y). 

end