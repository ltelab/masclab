% Haralick features exploration

clear all; close all;

path = '/home/praz/Documents/MASC/sample_haralick/';
pic_names = dir(fullfile(path,'*.png'));
pic_names = {pic_names.name}';

for i=1:length(pic_names)

    data = imread(fullfile(path,pic_names{i}));
    bw = zeros(size(data)); 
    bw(data>0) = 1;

    % compute "1-pixel distance" glcm
    % angle = 0 deg
    glcm_0  = graycomatrix(data,'NumLevels',256,'Offset',[0 1],'Symmetric',true);
    % angle = 45 deg
    glcm_45  = graycomatrix(data,'NumLevels',256,'Offset',[-1 1],'Symmetric',true);
    % angle = 90 deg
    glcm_90  = graycomatrix(data,'NumLevels',256,'Offset',[-1 0],'Symmetric',true);
    % angle = 135 deg
    glcm_135 = graycomatrix(data,'NumLevels',256,'Offset',[-1 -1],'Symmetric',true);

    % keep only the important part of the co-occurence matrix, ie. the part
    % which describes the snowflake and not the background /!\ assuming that
    % there is no holes in the snowflake
    glcm_0   = glcm_0(17:end,17:end);
    glcm_45  = glcm_45(17:end,17:end);
    glcm_90  = glcm_90(17:end,17:end);
    glcm_135 = glcm_135(17:end,17:end);

    % compute 4 Haralick feature
    stats_0 = graycoprops(glcm_0,{'Contrast','Correlation','Energy','Homogeneity'});
    stats_45 = graycoprops(glcm_45,{'Contrast','Correlation','Energy','Homogeneity'});
    stats_90 = graycoprops(glcm_90,{'Contrast','Correlation','Energy','Homogeneity'});
    stats_135 = graycoprops(glcm_135,{'Contrast','Correlation','Energy','Homogeneity'});

    % average them
    Contrast = mean([stats_0.Contrast,stats_45.Contrast,stats_90.Contrast,stats_135.Contrast]);
    Contrast_std = std([stats_0.Contrast,stats_45.Contrast,stats_90.Contrast,stats_135.Contrast]);

    Correlation = mean([stats_0.Correlation,stats_45.Correlation,stats_90.Correlation,stats_135.Correlation]);
    Correlation_std = std([stats_0.Correlation,stats_45.Correlation,stats_90.Correlation,stats_135.Correlation]);

    Energy = mean([stats_0.Energy,stats_45.Energy,stats_90.Energy,stats_135.Energy]);
    Energy_std = std([stats_0.Energy,stats_45.Energy,stats_90.Energy,stats_135.Energy]);

    Homogeneity = mean([stats_0.Homogeneity,stats_45.Homogeneity,stats_90.Homogeneity,stats_135.Homogeneity]);
    Homogeneity_std = std([stats_0.Homogeneity,stats_45.Homogeneity,stats_90.Homogeneity,stats_135.Homogeneity]);

    % figure;
    % imshow(data);
    % figure;
    % subplot(221);
    % contourf(glcm_0,'Edgecolor','none');
    % axis equal;
    % colormap(hot); colorbar;
    % subplot(222);
    % contourf(glcm_45,'Edgecolor','none');
    % axis equal;
    % colormap(hot); colorbar;
    % subplot(223);
    % contourf(glcm_90,'Edgecolor','none');
    % axis equal;
    % colormap(hot); colorbar;
    % subplot(224);
    % contourf(glcm_135,'Edgecolor','none');
    % axis equal;
    % colormap(hot); colorbar;

    figure;
    subplot(131);
    imshow(data);
    h_text = subplot(1,3,2:3);
    xl = xlim(h_text); 
    xPos = xl(1) + diff(xl) / 2; 
    yl = ylim(h_text); 
    yPos = yl(1) + diff(yl) / 2; 
    plot_text = text(xPos, yPos, sprintf('Contrast : %2.4f +- %2.4f \n Correlation : %2.4f +- %2.4f \n Energy : %f +- %f \n Homogeneity : %2.4f +- %2.4f \n', ...
        Contrast,Contrast_std,Correlation,Correlation_std,Energy,Energy_std,Homogeneity,Homogeneity_std), 'Parent', h_text);
    set(plot_text, 'HorizontalAlignment', 'center','FontSize',12,'FontWeight','demi');
    set(h_text,'visible','off');  

end



