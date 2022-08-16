% function doy = DOY(year, month, day)
function doy = DOY(varargin)

if nargin == 3
    year = varargin{1};
    month = varargin{2};
    day = varargin{3};
elseif nargin == 1;
    [year, month, day, ~, ~, ~] = datevec(varargin{1});
end

%from wikipedia:
monthlength = [31 28 31 30 31 30 31 31 30 31 30 31];

if mod(year, 100) == 0 && ~mod(year,400) == 0
    monthlength(2) = 28;
elseif mod(year,4) == 0
    monthlength(2) = 29;
end

% monthlength(2)

if month ~= 1
monthdays = sum(monthlength(1:(month-1)));
else
    monthdays = 0;
end

doy = monthdays+day;
