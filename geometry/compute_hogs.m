function hog_feats = compute_hogs(im,orientation,width,height,nBins,nOrientations,crop,illustration)

    % rotate the image
    im = imrotate(im,-orientation,'bicubic');
    
    % default values 
    if nargin == 2
        width = 300;
        height = 300;
        nBins = 30;
        nOrientations = 8;
        crop = 1;
        illustration = 0;
    end
    
    % image rescaling 
    ratio = width/max(size(im));
    im = imresize(im,ratio,'bicubic');
    
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

    hog_feats = hog(single(im)/255,nBins,nOrientations,crop);   
    
    if illustration
        figure;
        subplot(1,2,1);
        imshow(im);
        subplot(1,2,2);
        imshow( hogDraw(hog_feats) ); colormap gray;
        axis off; colorbar off;
    end
    
    hog_feats = hog_feats(:)';

end
