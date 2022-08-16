% function to load atmospheric variables during a measurement campaign
% atm variables are retrieved from IDAWEB for Davos
%
% written by  : Christophe Praz
% last update : October 2016
function atm = load_atm_vars(filepath,tstr_start,tstr_stop,illustration)

% function arguments
% filepath = '/home/praz/Documents/MASC/masclab/atm_vars/order_45574_data.txt';
% tstr_start = '20160101000000';
% tstr_stop = '20160131235959';
% illustration = true;

try   
    tmp_data = importdata(filepath);  
catch   
    fprintf('Error in load_atm_vars : impossible to load data from %s \n',filepath);
    fprintf('Function returns an empty structure atm = [] \n');
    atm = [];
    return
end

tmp_tnum = datenum(num2str(tmp_data.data(:,1)),'yyyymmddHHMM');

% case 1 : t_start ~= t_stop and we load the whole interval
if ~strcmp(tstr_start,tstr_stop)
    idx = find(tmp_tnum>=datenum(tstr_start,'yyyymmddHHMMSS') & tmp_tnum<=datenum(tstr_stop,'yyyymmddHHMMSS'));  
    
    if ~isempty(idx)
        atm.tnum = tmp_tnum(idx);
        atm.gust3s = tmp_data.data(idx,2);
        atm.gust1s = tmp_data.data(idx,3);
        atm.p = tmp_data.data(idx,4);
        atm.pred = tmp_data.data(idx,5);
        atm.T = tmp_data.data(idx,6);
        atm.precip = tmp_data.data(idx,7);
        atm.precip = atm.precip.*6; % from mm/10min -> mm/h
        atm.RH = tmp_data.data(idx,8);
        atm.winv = tmp_data.data(idx,9);
        atm.windir = tmp_data.data(idx,10);
    
    else 
        fprintf('Error in load_atm_vars : no data found in the entered time interval \n');
        fprintf('Function returns an empty structure atm = []\n');
        atm = [];
        
    end

% case 2 : t_start == t_stop and we load/interpolate the vars at that time
else
    tnum_val = datenum(tstr_start,'yyyymmddHHMMSS');
    idx_inf = find(tmp_tnum <= tnum_val,1,'last');
    idx_sup = find(tmp_tnum >= tnum_val,1,'first');
    tnum_1h = datenum([0 0 0 1 0 0]);
    
    if isempty(idx_inf) || isempty(idx_sup)
        fprintf('Error in load_atm_vars : entered time is outside of the time interval loaded \n');
        fprintf('Function returns an empty structure atm = [] \n');
        atm = [];
        return
        
    elseif idx_inf==idx_sup
        atm.tnum = tmp_tnum(idx_inf);
        atm.gust3s = tmp_data.data(idx_inf,2);
        atm.gust1s = tmp_data.data(idx_inf,3);
        atm.p = tmp_data.data(idx_inf,4);
        atm.pred = tmp_data.data(idx_inf,5);
        atm.T = tmp_data.data(idx_inf,6);
        atm.precip = tmp_data.data(idx_inf,7);
        atm.precip = atm.precip*6;
        atm.RH = tmp_data.data(idx_inf,8);
        atm.winv = tmp_data.data(idx_inf,9);
        atm.windir = tmp_data.data(idx_inf,10);
        
    elseif abs(tmp_tnum(idx_sup)-tmp_tnum(idx_inf)>=tnum_1h)
        fprintf('Error in load_atm_vars : entered time had to be interpolated between two measurements separated by more than 1 hour \n');
        fprintf('Function returns an empty structure atm = [] \n');
        
    else
        fprintf('Atm variables interpolated between %s and %s \n',datestr(tmp_tnum(idx_inf)),datestr(tmp_tnum(idx_sup)));
        atm.tnum = tnum_val;
        atm.gust3s = tmp_data.data(idx_inf,2) + (tnum_val - tmp_tnum(idx_inf))/(tmp_tnum(idx_sup) - tmp_tnum(idx_inf)) * (tmp_data.data(idx_sup,2) - tmp_data.data(idx_inf,2));
        atm.gust1s = tmp_data.data(idx_inf,3) + (tnum_val - tmp_tnum(idx_inf))/(tmp_tnum(idx_sup) - tmp_tnum(idx_inf)) * (tmp_data.data(idx_sup,3) - tmp_data.data(idx_inf,3));
        atm.p = tmp_data.data(idx_inf,4) + (tnum_val - tmp_tnum(idx_inf))/(tmp_tnum(idx_sup) - tmp_tnum(idx_inf)) * (tmp_data.data(idx_sup,4) - tmp_data.data(idx_inf,4));
        atm.pred = tmp_data.data(idx_inf,5) + (tnum_val - tmp_tnum(idx_inf))/(tmp_tnum(idx_sup) - tmp_tnum(idx_inf)) * (tmp_data.data(idx_sup,5) - tmp_data.data(idx_inf,5));
        atm.T = tmp_data.data(idx_inf,6) + (tnum_val - tmp_tnum(idx_inf))/(tmp_tnum(idx_sup) - tmp_tnum(idx_inf)) * (tmp_data.data(idx_sup,6) - tmp_data.data(idx_inf,6));
        atm.precip = tmp_data.data(idx_inf,7) + (tnum_val - tmp_tnum(idx_inf))/(tmp_tnum(idx_sup) - tmp_tnum(idx_inf)) * (tmp_data.data(idx_sup,7) - tmp_data.data(idx_inf,7));
        atm.precip = atm.precip*6;
        atm.RH = tmp_data.data(idx_inf,8) + (tnum_val - tmp_tnum(idx_inf))/(tmp_tnum(idx_sup) - tmp_tnum(idx_inf)) * (tmp_data.data(idx_sup,8) - tmp_data.data(idx_inf,8));
        atm.winv = tmp_data.data(idx_inf,9) + (tnum_val - tmp_tnum(idx_inf))/(tmp_tnum(idx_sup) - tmp_tnum(idx_inf)) * (tmp_data.data(idx_sup,9) - tmp_data.data(idx_inf,9));
        atm.windir = tmp_data.data(idx_inf,10) + (tnum_val - tmp_tnum(idx_inf))/(tmp_tnum(idx_sup) - tmp_tnum(idx_inf)) * (tmp_data.data(idx_sup,10) - tmp_data.data(idx_inf,10));
        
    end
    
end

% plot if desired
if illustration == true

    figure;
    subplot(3,1,1);
    plot(atm.tnum,atm.T,'r-');
    datetick('x');
    xlabel('Time [UT]');
    ylabel('T [C]');
    subplot(3,1,2);
    plot(atm.tnum,atm.RH,'b-');
    datetick('x');
    xlabel('Time [UT]');
    ylabel('Rel. Hum. [%]');
    subplot(3,1,3);
    plot(atm.tnum,atm.winv,'k-');
    datetick('x');
    xlabel('Time [UT]');
    

end

end
    


% % % function thygan = load_wfj_thygan(filepath,tsrt_start,tstr_stop,illustration)
% % 
% % Input  :  tsrt_start   : from (yyyymmddHHMMSS)
% %           tsrt_stop    : to (yyyymmddHHMMSS)
% %           illustration : plot of T and RH over the period
% % 
% % Output :  thygan       : structure containing tnum, T and RH
% % 
% % -- created on 24.06.2015, Praz Christophe
% % function thygan = load_wfj_thygan(tstr_start,tstr_stop,illustration)
% % 
% % filepath = '/home/praz/Documents/MASC/Davos_ground_data/slf-sma/15WFJ_METEO_thygan.dat';
% % tstr_start = '20150620000000';
% % tstr_stop = '20150620235959';
% % illustration = true;
% % 
% % try
% %     
% %     tmp_data = importdata(filepath);
% %     
% %     
% % catch
% %     
% %     fprintf('Error (in load_wfj_thygan) : impossible to load data from %s \n',filepath);
% %     fprintf('Function returns an empty structure thygan = [] \n');
% %     
% % end
% % 
% % tmp_tnum = datenum(tmp_data.textdata(4:end),'dd.mm.yyyy HH:MM'); % 3 lines of headers
% % 
% % if nargin == 1
% %     
% %     if only 1 argument, we do a linear interpol between last meas. before
% %     and first meas. after
% %     tnum_val = datenum(tstr_start,'yyyymmddHHMMSS');
% %     idx_inf = find(tmp_tnum <= tnum_val,1,'last');
% %     idx_sup = find(tmp_tnum >= tnum_val,1,'first');
% %     
% %     if ~isempty(idx_inf) && ~isempty(idx_sup) && idx_sup-idx_inf == 1
% %     
% %         thygan.T = tmp_data.data(idx_inf,1) + (tnum_val - tmp_tnum(idx_inf))/(tmp_tnum(idx_sup) - tmp_tnum(idx_inf)) * (tmp_data.data(idx_sup,1) - tmp_data.data(idx_inf,1));
% %         thygan.RH = tmp_data.data(idx_inf,2) + (tnum_val - tmp_tnum(idx_inf))/(tmp_tnum(idx_sup) - tmp_tnum(idx_inf)) * (tmp_data.data(idx_sup,2) - tmp_data.data(idx_inf,2));
% %         thygan.tnum = tnum_val;
% %         
% %     elseif ~isempty(idx_inf) && ~isempty(idx_sup) && idx_sup-idx_inf == 0
% %         
% %         thygan.T = tmp_data.data(idx_inf,1);
% %         thygan.RH = tmp_data.data(idx_inf,2);
% %         thygan.tnum = tnum_val;
% %         
% %     else
% %              
% %         thygan.T = NaN;
% %         thygan.RH = NaN;
% %         thygan.tnum = tnum_val;
% %         fprintf('Warning : no reliable thygan measurements found for %s \n',datestr(tnum_val));
% %         
% %     end
% %     
% %     
% % else
% % 
% %     select period of interest
% %     tnum_start = datenum(tstr_start,'yyyymmddHHMMSS');
% %     tnum_stop = datenum(tstr_stop,'yyyymmddHHMMSS');
% %     idx_start = find(tmp_tnum >= tnum_start,1,'first');
% %     idx_stop = find(tmp_tnum <= tnum_stop,1,'last');
% % 
% %     select data of interest
% %     thygan.tnum = tmp_tnum(idx_start:idx_stop);
% %     thygan.T = tmp_data.data(idx_start:idx_stop,1);
% %     thygan.RH = tmp_data.data(idx_start:idx_stop,2);
% % 
% %     plot if desired
% %     if illustration == true
% % 
% %         figure;
% %         subplot(2,1,1);
% %         plot(thygan.tnum,thygan.T,'r-');
% %         datetick('x');
% %         xlabel('Time');
% %         ylabel('T [C]');
% %         subplot(2,1,2);
% %         plot(thygan.tnum,thygan.RH,'b-');
% %         datetick('x');
% %         xlabel('Time');
% %         ylabel('Rel. Hum. [%]');
% % 
% %     end
% % 
% % end
% % 
% % end


