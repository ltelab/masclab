% some folders :
%'/home/praz/Documents/MASC/sample_snow_20150620/sample_cool/processed_05-Feb-2016/DATA/GOOD';
%'/media/praz/Masc-Data/APRES3_2015/PROCESSED_20160119/DATA/GOOD';


function [Xt,xhi,Xname,Xname_only] = load_processed_time(dirname,t_str_start,t_str_stop)

    fprintf('Load processed data time footprint...');

    tmin = datenum(t_str_start,'yyyymmddHHMMSS');
    tmax = datenum(t_str_stop,'yyyymmddHHMMSS');

    file_list = dir(fullfile(dirname,'20*.mat'));
    % for Massimo data
    if isempty(file_list)
        file_list = dir(fullfile(dirname,'ICE*.mat'));
    end
    file_only_list = {file_list.name};
    file_list = fullfile(dirname,file_only_list);

    % loop over the snowflakes and save each feature in a vector of the main
    % structure data
    n_flakes = 0;
    Xt = [];
    xhi = [];
    Xname = {};
    Xname_only = {};
    for i=1:length(file_list)

        load(file_list{i});
        % load only pictures within the time interval
        if roi.tnum >= tmin && roi.tnum <= tmax
            
            n_flakes = n_flakes + 1;
            Xt(i) = roi.tnum;
            xhi(i) = roi.xhi;
            % the full name!!!
            Xname{i} = file_list{i};
            Xname_only{i} = file_only_list{i};

        end

    end
    fprintf('   Done! %u images found. \n',n_flakes);
    Xt = Xt';
    xhi = xhi';
    Xname = Xname';
    Xname_only = Xname_only';
    
end
