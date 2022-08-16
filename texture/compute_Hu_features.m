function Hu = compute_Hu_features(im)
    
% define a co-ordinate system for image
im = double(im);
[height,width] = size(im);
xgrid = repmat([1:height]',1,width);
ygrid = repmat([1:width],height,1);

% compute "textural centroid"
M_00 = sum(im(:));
M_10 = sum(sum((xgrid.^1).*(ygrid.^0).*im));
M_01 = sum(sum((xgrid.^0).*(ygrid.^1).*im));
xc = M_10/M_00;
yc = M_01/M_00;

% normalize coordinate system by substracting centroid
xnorm = xgrid - xc;
ynorm = ygrid - yc;

% normalized central moments computation
%eta_01 = normalized_central_moments(im,xnorm,ynorm,0,1);
%eta_10 = normalized_central_moments(im,xnorm,ynorm,1,0);
eta_11 = normalized_central_moments(im,xnorm,ynorm,1,1);
eta_02 = normalized_central_moments(im,xnorm,ynorm,0,2);
eta_20 = normalized_central_moments(im,xnorm,ynorm,2,0);
eta_21 = normalized_central_moments(im,xnorm,ynorm,2,1);
eta_12 = normalized_central_moments(im,xnorm,ynorm,1,2);
eta_22 = normalized_central_moments(im,xnorm,ynorm,2,2);
eta_03 = normalized_central_moments(im,xnorm,ynorm,0,3);
eta_30 = normalized_central_moments(im,xnorm,ynorm,3,0);
%eta_13 = normalized_central_moments(im,xnorm,ynorm,1,3);
%eta_31 = normalized_central_moments(im,xnorm,ynorm,3,1);
%eta_23 = normalized_central_moments(im,xnorm,ynorm,2,3);
%eta_32 = normalized_central_moments(im,xnorm,ynorm,3,2);
%eta_33 = normalized_central_moments(im,xnorm,ynorm,3,3);


% calculate first 8 hu moments of order 3
Hu.h1 = eta_20 + eta_02;
Hu.h2 = (eta_20 - eta_02)^2 + 4*eta_11;
Hu.h3 = (eta_30 - 3*eta_12)^2 + (eta_03 - 3*eta_21)^2;
Hu.h4 = (eta_30 + eta_12)^2 + (eta_03 + eta_21)^2;
Hu.h5 = (eta_30 - 3*eta_12)*(eta_30 + eta_12)*((eta_30 + eta_12)^2 - 3*(eta_21 + eta_03)^2) + (3*eta_21 - eta_03)*(eta_21 + eta_03)*(3*(eta_30 + eta_12)^2 - (eta_03 + eta_21)^2);
Hu.h6 = (eta_20 - eta_02)*((eta_30 + eta_12)^2 - (eta_21 + eta_03)^2) + 4*eta_11*(eta_30 + eta_12)*(eta_21 + eta_03);
Hu.h7 = (3*eta_21 - eta_03)*(eta_30 + eta_12)*((eta_30 + eta_12)^2 - 3*(eta_21 + eta_03)^2) - (eta_30 - 3*eta_12)*(eta_21 + eta_03)*(3*(eta_30 + eta_12)^2 - (eta_03 + eta_21)^2);
Hu.h8 = eta_11*((eta_30 + eta_12)^2 - (eta_03 + eta_21)^2) - (eta_20 - eta_02)*(eta_30 + eta_12)*(eta_21 + eta_03);
 
end