% hog investigation 

clear all; close all;
imdir = '/home/praz/Documents/MASC/masclab/events/dendrites';
datafiles = dir(fullfile(imdir,'*.mat'));
datafiles = {datafiles.name}'; 
load(fullfile(imdir,datafiles{60})); 

im = roi.data;

% width and height of resized image (in pixels)
width = 300;
height = 300;
nBins = 30;
nOrientations = 8;

% ratio = width/max(size(im));
% im_alt1 = imresize(im,ratio,'bicubic');
% im_alt1 = imrotate(im_alt1,-roi.orientation,'bicubic');

im_alt2 = imrotate(im,-roi.orientation,'bicubic');
ratio = width/max(size(im_alt2));
im_alt2 = imresize(im_alt2,ratio,'bicubic');
im = im_alt2;

if size(im,1) > height || size(im,2) > width
    
    max_dim = max(size(im));
    ratio = width/(max_dim+1);
    im = imresize(im,ratio);
end

n = 1;
while size(im,2) < width
   new_column = zeros(size(im,1),1);
   if mod(n,2) == 1
       im = [im new_column];
   else
       im = [new_column im];
   end
   n = n + 1;
end

n = 1;
while size(im,1) < height
   new_line = zeros(1,size(im,2));
   if mod(n,2) == 1
       im = [im; new_line];
   else
       im = [new_line; im];
   end
   n = n + 1;
end


% figure; imshow(roi.data);
% figure; imshow(im);
% figure; imshow(imrotate(im,-roi.orientation,'bicubic'));

hog_features = hog(single(im)/255,nBins,nOrientations,1);

figure;
subplot(1,3,1);
imshow(roi.data);
subplot(1,3,2);
imshow(im);
subplot(1,3,3);
imshow( hogDraw(hog_features) ); colormap gray;
axis off; colorbar off;






% function refresh_image(handles,roi)
% pos =  getpixelposition(handles.axes1);%,'Position');
% width = floor(pos(3));
% height = floor(pos(4));
% if roi.width+2 >= width || roi.height+2 >= height
%     imshow(roi.data);
% else
%     im = uint8(zeros(width,height));
%     im(height/2-floor(roi.height/2)+1:height/2+ceil(roi.height/2),width/2-floor(roi.width/2)+1:width/2+ceil(roi.width/2)) = roi.data;
%     imshow(im);
% end
% if isempty(handles.labelName_list{handles.current_idx})
%     title('not labelled yet.');
% else
%     title(sprintf('current label : %s.',handles.labelName_list{handles.current_idx}));
% end   