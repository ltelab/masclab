%load example
clear all; close all;
imdir = '/home/praz/Documents/MASC/masclab/events/aggregates';
datafiles = dir(fullfile(imdir,'*.mat'));
datafiles = {datafiles.name}'; 
load(fullfile(imdir,datafiles{1}));   

% define a co-ordinate system for image
im = double(roi.data);
[height,width] = size(roi.data);
xgrid = repmat([1:height]',1,width);
ygrid = repmat([1:width],height,1);

% compute centroid 
s = regionprops(roi.bw_mask_filled,'centroid');
if length(s)~=1
    disp('Warning ! There migth be 2 ROIs in that bw_mask_filled');
end

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


%calculate first 8 hu moments of order 3
I_one   = eta_20 + eta_02;
I_two   = (eta_20 - eta_02)^2 + 4*eta_11;
I_three = (eta_30 - 3*eta_12)^2 + (eta_03 - 3*eta_21)^2;
I_four  = (eta_30 + eta_12)^2 + (eta_03 + eta_21)^2;
I_five  = (eta_30 - 3*eta_12)*(eta_30 + eta_12)*((eta_30 + eta_12)^2 - 3*(eta_21 + eta_03)^2) + (3*eta_21 - eta_03)*(eta_21 + eta_03)*(3*(eta_30 + eta_12)^2 - (eta_03 + eta_21)^2);
I_six   = (eta_20 - eta_02)*((eta_30 + eta_12)^2 - (eta_21 + eta_03)^2) + 4*eta_11*(eta_30 + eta_12)*(eta_21 + eta_03);
I_seven = (3*eta_21 - eta_03)*(eta_30 + eta_12)*((eta_30 + eta_12)^2 - 3*(eta_21 + eta_03)^2) - (eta_30 - 3*eta_12)*(eta_21 + eta_03)*(3*(eta_30 + eta_12)^2 - (eta_03 + eta_21)^2);
I_eight = eta_11*((eta_30 + eta_12)^2 - (eta_03 + eta_21)^2) - (eta_20 - eta_02)*(eta_30 + eta_12)*(eta_21 + eta_03);



% 
% [x_bar, y_bar] = centerOfMass(image,xgrid,ygrid);
% 
% % normalize coordinate system by subtracting mean
% xnorm = x_bar - xgrid;
% ynorm = y_bar - ygrid;
% 
% % central moments
% mu_11 = central_moments( image ,xnorm,ynorm,1,1);
% mu_20 = central_moments( image ,xnorm,ynorm,2,0);
% mu_02 = central_moments( image ,xnorm,ynorm,0,2);
% mu_21 = central_moments( image ,xnorm,ynorm,2,1);
% mu_12 = central_moments( image ,xnorm,ynorm,1,2);
% mu_03 = central_moments( image ,xnorm,ynorm,0,3);
% mu_30 = central_moments( image ,xnorm,ynorm,3,0);


% % calculate "scale invariant" central moments
% function eta_pq = normalized_central_moments(im,xnorm,ynorm,p,q)
%     
%     mu_pq = sum(sum((xnorm.^p).*(ynorm.^q).*im));
%     mu_00 = sum(sum(im)); %mu_00
%     % normalise moments for scale invariance
%     eta_pq = mu_pq/(mu_00^(1+(p+q)/2));
%     
% end