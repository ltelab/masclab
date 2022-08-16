% skeleton analysis exploration

clear all; close all;

path = '/home/praz/Documents/MASC/sample_haralick/';
pic_names = dir(fullfile(path,'*.png'));
pic_names = {pic_names.name}';

for i=1:length(pic_names)

    data = imread(fullfile(path,pic_names{i}));
    bw = zeros(size(data)); 
    bw(data>15) = 1;

    % smooth the silhouette
    bw_close = bwmorph(bw,'close');
    bw_close_open = bwmorph(bw_close,'open');
% 
     bw_skel = bwmorph(bw,'skel',Inf);
    bw_skel_close = bwmorph(bw_close,'skel',Inf);
    %bw_skel_open_close = bwmorph(bw_open_close,'skel',Inf);
    bw_skel_close_open = bwmorph(bw_close_open,'skel',Inf);
    
    % compute a "kind of" simplified skeleton
    bw_thin = bwmorph(bw_close_open,'thin',Inf);
     
    % remove potential leftover pixels
    all_roi = regionprops(bw_thin,'Area','PixelIdxList');
    all_areas = [all_roi.Area];
    [~,idx] = max(all_areas);
    bw_thin(:) = 0;
    bw_thin(all_roi(idx).PixelIdxList) = 1;
    
    % compute associated features
    perim = sum(sum(bwperim(bw_close_open)));
    area = sum(bw_close_open(:)); 
    perim_skel = sum(bw_thin(:));
   
    skel.p_ratio = perim_skel/perim;
    skel.A_ratio = perim_skel/area;
    skel.N_ends = sum(sum(bwmorph(bw_thin,'endpoints')));
    skel.N_branches = sum(sum(bwmorph(bw_thin,'branchpoints')));

    % se = strel('disk',1);        
    % bw_eroded = imerode(bw_close_open,se);



    % to do : dilate and erode before computing the skeleton

    figure;
    subplot(411);
    data_plot = data;
    data_plot(bw_skel==1) = 255;
    imshow(data_plot);
    title('raw mask skeleton');
    subplot(412);
    data_plot = data;
    data_plot(bw_skel_close_open==1) = 255;
    imshow(data_plot);
    title('opening + closing');
    subplot(413);
    data_plot = data;
    data_plot(bw_thin==1) = 255;
    imshow(data_plot);
    title('opening + closing + thining');
%     h_text = subplot(414);
%     xl = xlim(h_text); 
%     xPos = xl(1) + diff(xl) / 2; 
%     yl = ylim(h_text); 
%     yPos = yl(1) + diff(yl) / 2; 
%     plot_text = text(xPos, yPos, sprintf('R_{p} raw = %2.2f \n R_{p} smooth = %2.2f \n R_{A} raw = %2.2f \n R_{A} smooth = %2.2f \n',R1,R2,R3,R4), 'Parent', h_text);
%     set(plot_text, 'HorizontalAlignment', 'center','FontSize',12,'FontWeight','demi');
%     set(h_text,'visible','off');       


end

% figure;
% subplot(321);
% imshow(bw);
% title('mask');
% subplot(322);
% imshow(bw_open);
% title('mask opened');
% subplot(323);
% imshow(bw_close);
% title('mask closed');
% subplot(324);
% imshow(bw_open_close);
% title('mask opened + closed');
% subplot(325);
% imshow(bw_close_open);
% title('mask closed + opened');




