% whole GLCM as a vector investigation

%load example
clear all; close all;
imdir = '/home/praz/Documents/MASC/masclab/events/aggregates';
datafiles = dir(fullfile(imdir,'*.mat'));
datafiles = {datafiles.name}'; 
load(fullfile(imdir,datafiles{1}));  

im = roi.data;

% compute "1-pixel distance" glcm
NumLevels = 256;
% angle = 0 deg
glcm_0  = graycomatrix(im,'NumLevels',NumLevels,'Offset',[0 1],'Symmetric',true);
% angle = 45 deg
glcm_45  = graycomatrix(im,'NumLevels',NumLevels,'Offset',[-1 1],'Symmetric',true);
% angle = 90 deg
glcm_90  = graycomatrix(im,'NumLevels',NumLevels,'Offset',[-1 0],'Symmetric',true);
% angle = 135 deg
glcm_135 = graycomatrix(im,'NumLevels',NumLevels,'Offset',[-1 -1],'Symmetric',true);

% keep only the important part of the co-occurence matrix, ie. the part
% which describes the snowflake and not the background 
% /!\ holes in the sf are therefore neglected
glcm_0   = glcm_0(17:end,17:end);
glcm_45  = glcm_45(17:end,17:end);
glcm_90  = glcm_90(17:end,17:end);
glcm_135 = glcm_135(17:end,17:end);

glcm_all = zeros([size(glcm_0),4]);

glcm_all(:,:,1) = glcm_0;
glcm_all(:,:,2) = glcm_45;
glcm_all(:,:,3) = glcm_90;
glcm_all(:,:,4) = glcm_135;

results = GLCM_Features1(glcm_all);










% % put them together
% glcm_tot = glcm_0 + glcm_45 + glcm_90 + glcm_135;
% 
% % normalization
% glcm_tot = glcm_tot ./ sum(glcm_tot(:));
% 
% % keep only the upper triangular part (matrix is symetric)
% triangle_up_mask = true(size(glcm_tot));
% triangle_up_mask = triu(triangle_up_mask);
% 
% glcm_vec = glcm_tot(triangle_up_mask);



    % keep only the important part of the co-occurence matrix, ie. the part
    % which describes the snowflake and not the background 
    % /!\ holes in the sf are therefore neglected
    %glcm_0   = glcm_0(17:end,17:end);
    %glcm_45  = glcm_45(17:end,17:end);
    %glcm_90  = glcm_90(17:end,17:end);
    %glcm_135 = glcm_135(17:end,17:end);

    % compute 4 Haralick feature
    % stats_0 = graycoprops(glcm_0,{'Contrast','Correlation','Energy','Homogeneity'});
    % stats_45 = graycoprops(glcm_45,{'Contrast','Correlation','Energy','Homogeneity'});
    % stats_90 = graycoprops(glcm_90,{'Contrast','Correlation','Energy','Homogeneity'});
    % stats_135 = graycoprops(glcm_135,{'Contrast','Correlation','Energy','Homogeneity'});
