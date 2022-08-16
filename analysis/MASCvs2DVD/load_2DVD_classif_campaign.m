% load 2DVD classification campaign for a period of time
function dvd = load_2DVD_classif_campaign(tstr_start,tstr_stop,datadir)

    if nargin < 3
        datadir = '/data2/2DVD/Classified';
    end
    
    tstart = datenum(tstr_start,'yyyymmddHHMMSS');
    tstop = datenum(tstr_stop,'yyyymmddHHMMSS');
    tday = datenum([0 0 1 0 0 0]);
    
    dvd.t = [];
    dvd.N = [];
    dvd.classif = [];
    dvd.label = [];
    ini_label = true;
    
    while tstart <= tstop+tday % +tday to include a 1 day margin
        
        dvd_tmp = load_2DVD_classif(datestr(tstart,'yyyymmddHHMMSS'));
        
        if ~isempty(dvd_tmp.t)
            dvd.t = [dvd.t; dvd_tmp.t];
            dvd.N = [dvd.N; dvd_tmp.N];
            dvd.classif = [dvd.classif; dvd_tmp.classif];
            if ini_label
                dvd.label = dvd_tmp.label;
                ini_label = false;
            end
        end
        
        tstart = tstart+tday;
        
    end
    
    % final check to filter data out of the given time interval
    idx = find(dvd.t<datenum(tstr_start,'yyyymmddHHMMSS') | dvd.t>tstop);
    dvd.t(idx) = [];
    dvd.N(idx) = [];
    dvd.classif(idx) = [];
    
    
              
end


% load 2DVD classification data for 1 day
function dvd = load_2DVD_classif(tstr,datadir)

    if nargin < 2
        datadir = '/data2/2DVD/Classified';
    end

    tnum = datenum(tstr,'yyyymmddHHMMSS');
    tvec = datevec(tnum);
    tvec(4:6) = 0;
    tnum = datenum(tvec);

    cyear = num2str(tvec(1));
    cDOY = DOY(tnum);
    cfile = sprintf('V%s%03u_distribution_descriptors_60.txt',cyear(3:4),cDOY);

    try
        data = importdata(fullfile(datadir,cfile));
        %dvd.t = data.data(:,1);
        n = length(data.data(:,1));
        dvd.t = tnum + datenum([zeros(n,5) data.data(:,2)]); 
        dvd.N = data.data(:,3);
        dvd.classif = data.data(:,236);
        dvd.label = {'AG','DE','GR','SP','CO','MS','R','RIM'};

    catch ME
        dvd.t = [];
        dvd.N = [];
        dvd.classif = [];
        dvd.label = [];
        fprintf('Error in load_2DVD_classif : %s (%s) \n',ME.message,cfile);

    end
    
end



