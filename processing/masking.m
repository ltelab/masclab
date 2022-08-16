% === Create a dark background for MASC images ======
%
% data = masking(flake_data,flake_cam,process)
%
% inputs
%
%   flake_data : MASC picture (2d matrix)
%   flake_cam  : id of the camera
%   process    : struct of process_params
%
% output
%
%   data       : MASC picture masked
%
%
% Adapted from Tim Garrett, University of Utah. 
% Author : Christophe Praz (christophe.praz@epfl.ch)
% Last update : August 2015
% ===================================================
function data = masking(flake_data,flake_cam,process)

        data = flake_data;
        
        % remove clutter according to the discard matrix
        if flake_cam == process.camera_order(1) % = led à gauche
            
            data(:,1:process.discardmat(3)) = 0;
            
        elseif flake_cam == process.camera_order(3) % = led à droite
            
            data(:,end-process.discardmat(4):end) = 0;
            
        end
          
        if ismember(flake_cam,process.camera_order(1:3))
            data(1:process.discardmat(1),:) = 0;
            data(end-process.discardmat(2):end,:) = 0;
        end
       
        % Remove clutter below a certain lum. treshold
        data(data<=process.backthresh) = 0;
          
end



%         if MASCtype == 0;     
%             
%             background = find(flakebwin<40);
%             backthresh = quantile(flakebwin(background),0.995) + 2;
%             backmin = quantile(flakebwin(background),0.01);
%             flakebwin(1:topdiscard,:) = backmin; %Remove top clutter
%             flakebwin(find(flakebwin<=backthresh)) = backmin; %Smooth background so background edges aren't detected
% 
%         elseif MASCtype == 1
%             
%             flakebwin(find(flakebwin<=backthresh)) = mean(flakebwin(flakebwin<=backthresh)); %Smooth background so background edges aren't detected
%         
%         
%         elseif MASCtype == 2
%             
%             if cam == 1 %Only central camera has top and bottom clutter.
%                 flakebwin(1:topdiscard,:) = backthresh; %Remove top clutter
%                 flakebwin(vert-botdiscard:vert,:) = backthresh; %Remove bottom clutter
%                 flakebwin(:,1:leftdiscard) = backthresh; %Remove left clutter
%                 flakebwin(:,horz-rightdiscard:horz) = backthresh; %Remove right clutter
%             end
%             
%             flakebwin(find(flakebwin<=backthresh)) = mean(flakebwin(flakebwin<=backthresh)); %Smooth background so background edges aren't detected
%         
%         elseif MASCtype > 2; %Commercial cameras
%             
%                 flakebwin(1:topdiscard,:) = backthresh; %Remove top clutter
%                 flakebwin(vert-botdiscard:vert,:) = backthresh; %Remove bottom clutter
%                 flakebwin(:,1:leftdiscard) = backthresh; %Remove left clutter
%                 flakebwin(:,horz-rightdiscard:horz) = backthresh; %Remove right clutter
%             
%             flakebwin(find(flakebwin<=backthresh)) = mean(flakebwin(flakebwin<=backthresh)); %Smooth background so background edges aren't detected
%             
%         end
% 
%         flakebwout = flakebwin;