% === Parameters describing the Masc setup ==========
%
% LTE cameras specs :
% 5 MP (2448 X 2048) each
% lens : 12 mm FD
% 
% res [mm/pix] = fovmat / width
%
% Adapted from Tim Garrett, University of Utah. 
% Author : Christophe Praz (christophe.praz@epfl.ch)
% Last update : August 2015
% ===================================================
function cam_params = cam_params()
    
    % min time between 2 shots from the same camera
    cam_params.interarrivaltime = 0.5;
    
    % horizontal field of view of the 3 cameras [mm ?]
    cam_params.fovmat = [75 75 75];
    

end

%%%% FIELD OF VIEW FOR EACH LENS IN THE HORIZONTAL DIRECTION  %%%%
%%%% TRUST, BUT VERIFY!!!! %%%%%%
%%% 12 mm lens = 47 mm FOV on old MASC, 44 mm on new MASC with 1.2 MP
%%% 12 mm lens = 75 mm FOV with 5 MP !!!!
%%% 16 mm lens = 33 mm FOV with 1.2 MP
%%% 35 mm lens = 22 mm FOV with 5 MP
%%% 25 mm lens = 33 mm FOV with 5 MP
%%% 16 mm lens = 63 mm FOV with 5 MP