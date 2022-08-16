% === mascpaths() =======================================================
%
% add local folder and sub-folders to MATLAB paths
%
% this function has to be run first before using masclab
%
% Author : Tim Garrett, University of Utah
%
% =======================================================================
function mascpaths()

    try
    %addpath(genpath('pwd'));
        addpath(genpath(pwd));
        
    catch
        
    end

end
