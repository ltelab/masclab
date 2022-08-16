% 3d reconstruction tests
clearvars; close all;

im_mid = imread('triplets/2015.02.23_20.13.27_flake_18829_cam_1.png');
w_mid = size(im_mid,2);
h_mid = size(im_mid,1);

im_left = imread('triplets/2015.02.23_20.13.27_flake_18829_cam_0.png');
w_left = size(im_left,2);
h_left = size(im_left,1);

im_right = imread('triplets/2015.02.23_20.13.27_flake_18829_cam_2.png');
w_right = size(im_right,2);
h_right = size(im_right,1);

% square centering images
max_dim = max([h_left,h_mid,h_right,w_left,w_mid,w_right]);
im_mid_sq = square_center(im_mid,max_dim);
im_left_sq = square_center(im_left,max_dim);
im_right_sq = square_center(im_right,max_dim);

% silhouettes
bw_mid = logical(im_mid_sq > 0);
bw_left = logical(im_left_sq > 0);
bw_right = logical(im_right_sq > 0);

% the cube
hull_mid = ones(max_dim,max_dim,100);
%hull_mid(10:100,10:100,10:100) = 2;

% extrusion of the middle image
for k=1:size(hull_mid,3)
    slice = hull_mid(:,:,k);
    slice(bw_mid==0) = 0;
    hull_mid(:,:,k) = slice;
end
hull_mid(:,:,1) = 0;
hull_mid(:,:,end) = 0;

% rotate middle image -36 degree
theta = 36;
M_rot = [cos(theta), -sin(theta), 0; ...
    sin(theta), cos(theta), 0; ...
    0, 0, 1];
hull_mid_swapped_dim = permute(hull_mid,[3,2,1]);
hull_left_swapped_dim = imrotate(hull_mid_swapped_dim,36);
hull_left = ipermute(hull_left_swapped_dim,[3,2,1]);
%resize slice by slice
for k=1:size(hull_left,3)
    hull_left_rescaled(:,:,k) = imresize(hull_left(:,:,k),[max_dim,max_dim]);
end

% extrusion of the left image
for k=1:size(hull_left_rescaled,3)
    slice = hull_left_rescaled(:,:,k);
    slice(bw_left==0) = 0;
    hull_left_rescaled(:,:,k) = slice;
end


figure; grid on; box on;
%[X,Y,Z] = meshgrid(1:max_dim,1:max_dim,1:max_dim);
p = patch(isosurface(hull_mid,0.9999));
%isonormals(X,Y,Z,hull_mid,p)
p.FaceColor = 'red';
p.EdgeColor = 'none';
daspect([1,1,1]);
view(3); axis([0 max_dim 0 max_dim 0 max_dim]);
camlight ;
lighting gouraud;
xlabel('X - width');
ylabel('Y - height');
zlabel('Z - depth');
set(gca,'Ydir','reverse');
%set(gca,'Ydir','reverse');

figure; grid on; box on;
%[X,Y,Z] = meshgrid(1:max_dim,1:max_dim,1:max_dim);
p = patch(isosurface(hull_left_rescaled,0.9999));
%isonormals(X,Y,Z,hull_mid,p)
p.FaceColor = 'red';
p.EdgeColor = 'none';
daspect([1,1,1]);
view(3); axis([0 max_dim 0 max_dim 0 max_dim]);
camlight ;
lighting gouraud;
xlabel('X - width');
ylabel('Y - height');
zlabel('Z - depth');
set(gca,'Ydir','reverse');
%set(gca,'Ydir','reverse');



figure;
subplot(2,3,1);
imshow(im_left_sq);
subplot(2,3,2);
imshow(im_mid_sq);
subplot(2,3,3);
imshow(im_right_sq);
subplot(2,3,4);
imshow(bw_left);
subplot(2,3,5);
imshow(bw_mid);
subplot(2,3,6);
imshow(bw_right);