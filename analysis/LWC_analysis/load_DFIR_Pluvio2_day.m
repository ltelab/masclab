% load DFIR_Pluvio2_day
function pluvio = load_DFIR_Pluvio2_day(daystr)

    datadir = '/home/praz/Documents/MASC/masclab/LWC_analysis/Pluvio2DFIR_data';
    filename = strcat('WFJ',daystr,'.R2Pluvio2DFIR');
    
    try
        data = importdata(fullfile(datadir,filename));
    catch
        fprintf('Error: unable to open file %s\n',filename);
        pluvio = [];
        return;
    end

    try

        % remove header in textdata
        data.textdata = data.textdata(17:end,1:2);
        date = data.textdata(:,1);
        hour = data.textdata(:,2);
        pluvio.t_str = strcat(date,{' '},hour);
        pluvio.t = datetime(pluvio.t_str,'InputFormat','yyyyMMdd HH:mm:SS');
        pluvio.accuNRT = data.data(:,3); % in mm
        pluvio.accuNRT_total = data.data(:,4); % in mm
        pluvio.bucketRT = data.data(:,5); % in mm, raw instrument output
        
        if length(pluvio.t_str) ~= length(pluvio.accuNRT)

            fprintf('Warning : time vector and accumulation do not have the same size, something went wrong. File: %s\n',filename);

        end

    catch 

        fprintf('Error: unable to correctly read data from %s\n',filename);
        pluvio = [];
        return;

    end

end


