% assuming a string format like XXX_flake_id_YYY.ZZZ
function id = get_snowflake_id(string)

    try
        idx_start = strfind(string,'flake_') + 6; 
        string = string(idx_start:end);
        idx_stop = strfind(string,'_') - 1;
        string = string(1:idx_stop);
        id = str2num(string);

    catch
        fprintf('Error in get_snowflake_id : could not retrieve flake ID from string : %s \n',string);
        id = NaN;      
    end

end
        