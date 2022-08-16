%% function thygan = load_wfj_thygan(filepath,tsrt_start,tstr_stop,illustration)
%
% Input  :  tsrt_start   : from (yyyymmddHHMMSS)
%           tsrt_stop    : to (yyyymmddHHMMSS)
%           illustration : plot of T and RH over the period
%
% Output :  thygan       : structure containing tnum, T and RH
%
% -- created on 24.06.2015, Praz Christophe
function thygan = load_wfj_thygan(tstr_start,tstr_stop,illustration)

filepath = '/home/praz/Documents/MASC/Davos_ground_data/slf-sma/15WFJ_METEO_thygan.dat';
% tstr_start = '20150620000000';
% tstr_stop = '20150620235959';
% illustration = true;

try
    
    tmp_data = importdata(filepath);
    
    
catch
    
    fprintf('Error (in load_wfj_thygan) : impossible to load data from %s \n',filepath);
    fprintf('Function returns an empty structure thygan = [] \n');
    
end

tmp_tnum = datenum(tmp_data.textdata(4:end),'dd.mm.yyyy HH:MM'); % 3 lines of headers

if nargin == 1
    
    % if only 1 argument, we do a linear interpol between last meas. before
    % and first meas. after
    tnum_val = datenum(tstr_start,'yyyymmddHHMMSS');
    idx_inf = find(tmp_tnum <= tnum_val,1,'last');
    idx_sup = find(tmp_tnum >= tnum_val,1,'first');
    
    if ~isempty(idx_inf) && ~isempty(idx_sup) && idx_sup-idx_inf == 1
    
        thygan.T = tmp_data.data(idx_inf,1) + (tnum_val - tmp_tnum(idx_inf))/(tmp_tnum(idx_sup) - tmp_tnum(idx_inf)) * (tmp_data.data(idx_sup,1) - tmp_data.data(idx_inf,1));
        thygan.RH = tmp_data.data(idx_inf,2) + (tnum_val - tmp_tnum(idx_inf))/(tmp_tnum(idx_sup) - tmp_tnum(idx_inf)) * (tmp_data.data(idx_sup,2) - tmp_data.data(idx_inf,2));
        thygan.tnum = tnum_val;
        
    elseif ~isempty(idx_inf) && ~isempty(idx_sup) && idx_sup-idx_inf == 0
        
        thygan.T = tmp_data.data(idx_inf,1);
        thygan.RH = tmp_data.data(idx_inf,2);
        thygan.tnum = tnum_val;
        
    else
             
        thygan.T = NaN;
        thygan.RH = NaN;
        thygan.tnum = tnum_val;
        fprintf('Warning : no reliable thygan measurements found for %s \n',datestr(tnum_val));
        
    end
    
    
else

    % select period of interest
    tnum_start = datenum(tstr_start,'yyyymmddHHMMSS');
    tnum_stop = datenum(tstr_stop,'yyyymmddHHMMSS');
    idx_start = find(tmp_tnum >= tnum_start,1,'first');
    idx_stop = find(tmp_tnum <= tnum_stop,1,'last');

    % select data of interest
    thygan.tnum = tmp_tnum(idx_start:idx_stop);
    thygan.T = tmp_data.data(idx_start:idx_stop,1);
    thygan.RH = tmp_data.data(idx_start:idx_stop,2);

    % plot if desired
    if illustration == true

        figure;
        subplot(2,1,1);
        plot(thygan.tnum,thygan.T,'r-');
        datetick('x');
        xlabel('Time');
        ylabel('T [C]');
        subplot(2,1,2);
        plot(thygan.tnum,thygan.RH,'b-');
        datetick('x');
        xlabel('Time');
        ylabel('Rel. Hum. [%]');

    end

end

end


