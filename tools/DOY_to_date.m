
function [year, month, day] = DOY_to_date(year, doy)

monthlength = [31 28 31 30 31 30 31 31 30 31 30 31];

% sum(monthlength)
% length(monthlength)

if mod(year, 100) == 0 && ~mod(year,400) == 0
    monthlength(2) = 28;
elseif mod(year,4) == 0
    monthlength(2) = 29;
end

% monthlength(2)

month_counter = 1;

while doy > monthlength(month_counter)
    doy = doy - monthlength(month_counter);
    month_counter = month_counter+1;
end

month = month_counter;
day = doy;