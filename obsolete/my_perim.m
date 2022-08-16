function P = my_perim(x_vec,y_vec)
    P = 0;
    for i=1:(length(x_vec)-1)
        if x_vec(i+1) ~= x_vec(i) && y_vec(i+1) ~= y_vec(i)
            P = P + sqrt(2);
        else
            P = P+1;
        end
    end
    
    if x_vec(end) ~= x_vec(1) && y_vec(end) ~= y_vec(1)
        P = P + sqrt(2);
    else
        P = P+1;
    end

end