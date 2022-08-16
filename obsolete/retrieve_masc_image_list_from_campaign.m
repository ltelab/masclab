% script to retrieve a list of all masc images, id and cams for a campaign
clearvars; close all;

datadir = '/data2/MASC/Davos_winter_2015-16';
img_list = dir(fullfile(datadir,'**','*.png'));
img_list = {img_list.name}';

%%
ntot = numel(img_list);
data = nan(ntot,8);
%data.datestr = cell(numel(img_list),1);
%data.cam_id = nan(numel(img_list),1);
%data.flake_id = nan(numel(img_list),1);

for i=1:ntot
    
    fprintf('%u/%u \n',i,ntot);
    
    img_current = img_list{i};
    data(i,1:6) = datevec(img_current,'yyyy.mm.dd_HH.MM.SS');
    
    tmp_cam_id = img_current(end-4);
    data(i,7)  = str2num(tmp_cam_id);
    
    tmp_flake_id = img_current(27:end-10);
    data(i,8)    = str2num(tmp_flake_id);
    
end

legend = 'col. 1-6 = datevec; col. 7 = cam ID; col.8 = flake ID';
    

