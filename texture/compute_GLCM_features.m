function glcm_vec = compute_GLCM_features(im,NumLevels)

    % compute "1-pixel distance" glcm
    % angle = 0 deg
    glcm_0  = graycomatrix(im,'NumLevels',NumLevels,'Offset',[0 1],'Symmetric',true);
    % angle = 45 deg
    glcm_45  = graycomatrix(im,'NumLevels',NumLevels,'Offset',[-1 1],'Symmetric',true);
    % angle = 90 deg
    glcm_90  = graycomatrix(im,'NumLevels',NumLevels,'Offset',[-1 0],'Symmetric',true);
    % angle = 135 deg
    glcm_135 = graycomatrix(im,'NumLevels',NumLevels,'Offset',[-1 -1],'Symmetric',true);

    % put them together
    glcm_tot = glcm_0 + glcm_45 + glcm_90 + glcm_135;


    % normalization
    % need to check if this is necessary !
    glcm_tot = glcm_tot ./ sum(glcm_tot(:));

    % keep only the upper triangular part (matrix is symetric)
    triangle_up_mask = true(size(glcm_tot));
    triangle_up_mask = triu(triangle_up_mask);

    % put everything into one vector 
    glcm_vec = glcm_tot(triangle_up_mask);
    
 
end