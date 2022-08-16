% === brighten a masc snowflake picture ===================================
%
% This function modify and brigthen an input picture using CLAHE histogram
% equalization algorithm
%
% Author : Christophe Praz (christophe.praz@epfl.ch)
% Last update : August 2015
% ========================================================================
function imageout = brightening(imagein)

    tile_dim = 8;
    [nlines,ncols] = size(imagein);
    
    if nlines < tile_dim && ncols < tile_dim
        
        imagemod = uint8(zeros(tile_dim));
        imagemod(1:nlines,1:ncols) = imagein;
        imageout = adapthisteq(imagemod);
        imageout = imageout(1:nlines,1:ncols);
        
    elseif nlines < tile_dim
        
        imagemod = uint8(zeros(tile_dim,ncols));
        imagemod(1:nlines,1:ncols) = imagein;
        imageout = adapthisteq(imagemod);
        imageout = imageout(1:nlines,1:ncols);
        
    elseif ncols < tile_dim

        imagemod = uint8(zeros(nlines,tile_dim));
        imagemod(1:nlines,1:ncols) = imagein;
        imageout = adapthisteq(imagemod);
        imageout = imageout(1:nlines,1:ncols);
        
    else
        
        imageout = adapthisteq(imagein);
        
    end
    
end
    
%%%%%%%%%%%% old code
%                 try
%                     
%                     imageout = adapthisteq(imagein);
%                     
%                 catch
%                     
%                     disp('WARNING : could not brighten the snowflake picture (area probably too small).');
%                     
%                     imageout = imagein;
%                     
%                 end
%%%%%%%%%%%%%%%%%%%%
            

