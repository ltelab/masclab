% === Compute a set of Haralick features based on snowflake picture ======
%
% Compute basic Haralick features from the snowflake image (Contrast,
% Correlation, Energy, Homogeneity).
%
% h is the output structure containing the 4 features.
%
% Author : Christophe Praz (christophe.praz@epfl.ch)
% Last update : November 2015
% ========================================================================
% path = '/home/praz/Documents/MASC/sample_haralick/';
% pic_names = dir(fullfile(path,'*.png'));
% pic_names = {pic_names.name}';
function h = haralick_props(im,illustration)

    if nargin == 1
        
        illustration = false;
        
    end 
    
    pix_value_thresh = 2;

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
    glcm_0   = glcm_0(pix_value_thresh:end,pix_value_thresh:end);
    glcm_45  = glcm_45(pix_value_thresh:end,pix_value_thresh:end);
    glcm_90  = glcm_90(pix_value_thresh:end,pix_value_thresh:end);
    glcm_135 = glcm_135(pix_value_thresh:end,pix_value_thresh:end);
    
    
    % piece of code to compute the whole haralick features set
    % put everything in 1 matrix
%     glcm_all = zeros([size(glcm_0),4]);
%     glcm_all(:,:,1) = glcm_0;
%     glcm_all(:,:,2) = glcm_45;
%     glcm_all(:,:,3) = glcm_90;
%     glcm_all(:,:,4) = glcm_135;
%     h = GLCM_Features1(glcm_all);
%     return;
    
    
    

    % compute 4 Haralick feature
    stats_0 = graycoprops(glcm_0,{'Contrast','Correlation','Energy','Homogeneity'});
    stats_45 = graycoprops(glcm_45,{'Contrast','Correlation','Energy','Homogeneity'});
    stats_90 = graycoprops(glcm_90,{'Contrast','Correlation','Energy','Homogeneity'});
    stats_135 = graycoprops(glcm_135,{'Contrast','Correlation','Energy','Homogeneity'});

    % average them
    h.Contrast = mean([stats_0.Contrast,stats_45.Contrast,stats_90.Contrast,stats_135.Contrast]);
    h.Correlation = mean([stats_0.Correlation,stats_45.Correlation,stats_90.Correlation,stats_135.Correlation]);
    h.Energy = mean([stats_0.Energy,stats_45.Energy,stats_90.Energy,stats_135.Energy]);
    h.Homogeneity = mean([stats_0.Homogeneity,stats_45.Homogeneity,stats_90.Homogeneity,stats_135.Homogeneity]);
    
    % standard deviation displayed on plot
    h.Contrast_std = std([stats_0.Contrast,stats_45.Contrast,stats_90.Contrast,stats_135.Contrast]);
    h.Correlation_std = std([stats_0.Correlation,stats_45.Correlation,stats_90.Correlation,stats_135.Correlation]);
    h.Energy_std = std([stats_0.Energy,stats_45.Energy,stats_90.Energy,stats_135.Energy]);
    h.Homogeneity_std = std([stats_0.Homogeneity,stats_45.Homogeneity,stats_90.Homogeneity,stats_135.Homogeneity]);
    
    if illustration
  
        % plot I : co-occurence matrix
        figure;
        subplot(221);
        contourf(glcm_0,linspace(0,5,25),'Edgecolor','none');
        axis equal; 
        colormap(hot); colorbar; %caxis([0 7]);
        title('GLCM \theta = 0');
        subplot(222);
        contourf(glcm_45,linspace(0,5,25),'Edgecolor','none');
        axis equal;
        colormap(hot); colorbar;
        title('GLCM \theta = 45');
        subplot(223);
        contourf(glcm_90,linspace(0,5,25),'Edgecolor','none');
        axis equal;
        colormap(hot); colorbar;
        title('GLCM \theta = 90');
        subplot(224);
        contourf(glcm_135,linspace(0,5,25),'Edgecolor','none');
        axis equal;
        colormap(hot); colorbar;
        title('GLCM \theta = 135');

        % plot II : snowflake image + haralick features values
        figure;
        subplot(131);
        imshow(im);
        h_text = subplot(1,3,2:3);
        xl = xlim(h_text); 
        xPos = xl(1) + diff(xl) / 2; 
        yl = ylim(h_text); 
        yPos = yl(1) + diff(yl) / 2; 
        plot_text = text(xPos, yPos, sprintf('Contrast : %2.4f +- %2.4f \n Correlation : %2.4f +- %2.4f \n Energy : %f +- %f \n Homogeneity : %2.4f +- %2.4f \n', ...
            h.Contrast,h.Contrast_std,h.Correlation,h.Correlation_std,h.Energy,h.Energy_std,h.Homogeneity,h.Homogeneity_std), 'Parent', h_text);
        set(plot_text, 'HorizontalAlignment', 'center','FontSize',12,'FontWeight','demi');
        set(h_text,'visible','off'); 
        
    end

    
 
    
end