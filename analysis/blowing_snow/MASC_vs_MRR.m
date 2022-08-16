% retrieve MASC #images and MRR precip rate near ground (300-600m)
% in order to identify pure blowing snow / mixed / pure precipitation
% events during APRES3

%clear all; close all;

MASC_datadir = '/ltedata/APRES3_2015/MASC/Raw_data';
MRR_datadir = '/ltedata/APRES3_2015/MRR/Data_Netcdf_format/ProcessedData';

t_vec_start = [2015 11 11 00 00 00];
t_vec_stop = [2016 01 30 00 00 00]; 
min_inter = 10;

dt_start = datetime(t_vec_start);
dt_stop = datetime(t_vec_stop);
dt_vec = [dt_start:minutes(min_inter):dt_stop]';


MASC_data = dir(fullfile(MASC_datadir,'**','*.png'));

%%
MASC = [];
MASC.name = {MASC_data.name}';
%MASC.dt = cell(numel(MASC.name),1);

for i=1:numel(MASC.name)
    MASC.dt(i,1) = datetime(MASC.name{i}(1:16),'InputFormat','yyyy.MM.dd_HH.mm');
end

N_Masc = zeros(numel(dt_vec)-1,1);
for i=1:(numel(dt_vec)-1)
    idx = find(MASC.dt >= dt_vec(i) & MASC.dt < dt_vec(i+1));
    N_Masc(i) = numel(idx);
end
        
%%
MRR_file = dir(fullfile(MRR_datadir,'**','*.nc'));
MRR_name = {MRR_file.name}';
MRR_folder = {MRR_file.folder}';

%%

MRR_dt = [];
MRR_rr = [];

for i=1:numel(MRR_folder)

    MRR_current = fullfile(MRR_folder{i},MRR_name{i});

    t = ncread(MRR_current,'Time');
    h = ncread(MRR_current,'H'); % h(1) = 100, h(2) = 200, ...
    rr = ncread(MRR_current,'RR');

    if exist('tmp_dt','var')
        clear tmp_dt
    end
    tmp_rr_vec = [];
    tmp_rr = [];
    for j=1:size(t,2)
        tmp_dt(j,1) = datetime(t(:,j)','InputFormat','yyMMddHHmmSS');
        tmp_rr_vec = rr(3:5,j);
        tmp_rr_vec(tmp_rr_vec == -99900) = NaN;
        tmp_rr(j,1) = nanmean(tmp_rr_vec)*10/3600; % convert mm/h into mm (10sec timestep)
    end
    
    MRR_dt = [MRR_dt; tmp_dt];
    MRR_rr = [MRR_rr; tmp_rr];
   
end

RR_MRR = zeros(numel(dt_vec)-1,1);
for i=1:(numel(dt_vec)-1)
    idx = find(MRR_dt >= dt_vec(i) & MRR_dt < dt_vec(i+1));
    RR_MRR(i) = nansum(MRR_rr(idx));
end




