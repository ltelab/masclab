function im_out = square_center(im_in,max_dim)

    if size(im_in,1) > max_dim || size(im_in,1) > max_dim
        fprintf('error : input image is larger than max_dim, original image returned \n');
        im_out = im_in;
        return;
    end

    n = 1;
    im_out = im_in;
       while size(im_out,2) < max_dim
           new_column = zeros(size(im_out,1),1);
           if mod(n,2) == 1
               im_out = [im_out new_column];
           else
               im_out = [new_column im_out];
           end
           n = n + 1;
       end

    n = 1;
    while size(im_out,1) < max_dim
       new_line = zeros(1,size(im_out,2));
       if mod(n,2) == 1
           im_out = [im_out; new_line];
       else
           im_out = [new_line; im_out];
       end
       n = n + 1;
    end
    
    
end