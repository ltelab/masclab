% script to create datasets of blowing snow/real snowfall particles from
% APRES3 2015-2016 campaign
clear all; close all;

load('classif_datastruct_APRES3_triplet_withgraupfix.mat');
t_all_flakes = datetime(data.Yt,'ConvertFrom','datenum');


% pure blowing snow
tmp = load('timesteps_blowing_snow_APRES32015-16.txt');
bs_min_str = cell(numel(tmp),1);
for i=1:length(tmp)
    bs_min_str{i} = num2str(tmp(i));
end
bs_min_num = datenum(bs_min_str,'yyyymmddHHMMSS');
bs_min = datetime(bs_min_num,'ConvertFrom','datenum');

idx_all = [];
for i=1:numel(bs_min)
    
    idx_tmp = find(t_all_flakes >= bs_min(i) & t_all_flakes <= bs_min(i) + minutes);
    if ~isempty(idx_tmp)
        idx_all(end+1,1) = idx_tmp;
    end
    
end





% quiet events : 22.11.2015

% windy events : 2015.11.11 8h-24h (wind >= 60 km/h)