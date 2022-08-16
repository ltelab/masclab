% assuming a string format like XXX_cam_id_YYY.ZZZ
function id = get_cam_id(string)

    try
        idx_start = strfind(string,'cam_') + 4; 
        string = string(idx_start:end);
        idx_stop_1 = strfind(string,'_') - 1;
        idx_stop_2 = strfind(string,'.') - 1;
        if isempty(idx_stop_1)
            idx_stop = idx_stop_2;
        elseif isempty(idx_stop_2)
            idx_stop = idx_stop_1;
        else
            idx_stop = min(idx_stop_1(1).idx_stop_2(i));
        end     
        string = string(1:idx_stop);
        id = str2num(string);

    catch
        fprintf('Error in get_snowflake_id : could not retrieve flake ID from string : %s \n',string);
        id = NaN;      
    end

end
        