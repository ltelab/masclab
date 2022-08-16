function [X_hog, X_hog_2, Xname] = load_hog_features(dirname,t_str_start,t_str_stop)

    %dirname = '/home/praz/Documents/MASC/sample_snow_20150620/sample_cool/processed_05-Feb-2016/DATA/GOOD';
    %'/media/praz/Masc-Data/APRES3_2015/PROCESSED_20160119/DATA/GOOD';

    tmin = datenum(t_str_start,'yyyymmddHHMMSS');
    tmax = datenum(t_str_stop,'yyyymmddHHMMSS');

    file_list = dir(fullfile(dirname,'20*.mat'));
    file_only_list = {file_list.name};
    file_list = fullfile(dirname,file_only_list);

    % loop over the snowflakes and retrieve hog features
    n_flakes = 0;
    X_hog = [];
    X_hog_2 = [];
    for i=1:length(file_list)

        load(file_list{i});

        % load only pictures within the time interval
        if roi.tnum >= tmin && roi.tnum <= tmax
            
            n_flakes = n_flakes + 1;     
            Xname{i} = file_only_list{i};
            X_hog(end+1,:) = roi.hog1;
            X_hog_2(end+1,:) = roi.hog2;

        end

    end
    fprintf('flakes found : %u \n',n_flakes);
    Xname = Xname';
    
end
