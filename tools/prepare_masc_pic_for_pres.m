% small script loading all images in the "dir" folder and change their size
% to an unique one (adding a black background, not scaling the snowflake)
% if the input image is already larger than the desired output, no
% modification is done
clear all; close all;


% input parameters
in_dir = '/home/praz/Desktop/lte_on_enac1files/commun1/Christophe/MASC/snowflakes_selection/melting_snow_raw';
out_dir = '/home/praz/Desktop/lte_on_enac1files/commun1/Christophe/MASC/snowflakes_selection/melting_snow';
if isempty(dir(out_dir))
   mkdir(out_dir);
end
% final dimensions of the image [pixels]
width = 300; 
height = 300;

% load pictures
file_list = dir(fullfile(in_dir,'*.png'));
file_only_list = {file_list.name};
file_list = fullfile(in_dir,file_only_list);

for i=1:length(file_list)
   im = imread(file_list{i});
   if length(size(im)) > 2
       im = rgb2gray(im);
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
   
   if size(im,1) == 300 && size(im,2) == 300
       imwrite(im,fullfile(out_dir,strcat('squared_',file_only_list{i})),'png','BitDepth', 8);
   end
       
        
end

