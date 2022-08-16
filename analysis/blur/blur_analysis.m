% investigate the blind deconvolution
clear all; close all;

% img_names = {'2015.06.20_09.29.52_flake_47260_cam_0.png','2015.06.20_09.33.11_flake_47426_cam_0.png','2015.06.20_09.51.15_flake_48415_cam_2.png',...
%     '2015.06.20_09.51.26_flake_48425_cam_1.png','2015.06.20_09.52.19_flake_48468_cam_2.png','2015.06.20_09.53.42_flake_48509_cam_2.png',...
%     '2015.06.20_09.53.44_flake_48511_cam_1.png','2015.06.20_09.54.58_flake_48553_cam_2.png','2015.06.20_09.56.59_flake_48612_cam_2.png',...
tmp = dir('*.png');
img_names = {tmp.name};
    
img = {}; img_blur = {}; 
sigma_blur = 0.2;
for i=1:numel(img_names)
    img{i} = imread(img_names{i});
    tmp_mean = mean(size(img{i}));
    sigma_blur = tmp_mean/25;
    img_blur{i} = imgaussfilt(img{i},sigma_blur);
    
    diff = (img{i}-img_blur{i});
    stdv(i) = std2(diff);
    meanv(i) = mean(diff(:));
    
    figure;
    subplot(121);
    imshow(img{i});
    title(sprintf('mean : %2.2f',meanv(i)));
    subplot(122);
    imshow(img_blur{i});
    title(sprintf('std : %2.2f',stdv(i)));
end







 
% figure;
% nlines = 5;
% k=1;
% for i=1:nlines
%     subplot(nlines,2,k);
%     imshow(img{i});
%     subplot(nlines,2,k+1);
%     imshow(img_blur{i});
%     k=k+2;
% end

%PSF_ini = fspecial('gaussian',3,10);
%img_deblur = deconvblind(img,PSF_ini);

% img_blur = imgaussfilt(img,2);
% 
% figure;
% subplot(121);
% imshow(img);
% subplot(122);
% imshow(img_blur);