% === Detect edges and filled areas in a MASC picture ======
%
% data_out = edge_detection(data_in,illustration)
%
% inputs
%
%   data_in      : MASC picture (2d matrix)
%   illustration : boolean whether you want to display the process
%   (optional parameter)
%
% output
%
%   data_out     : MASC picture edged and eroded (BW)
%
%
% Author : Christophe Praz (christophe.praz@epfl.ch)
% Last update : August 2015
% =========================================================
function data_out = edge_detection(data_in,strel_size,illustration)

    if nargin < 3
        
        illustration = false;
        
    end
    
    if nargin < 2
        
        strel_size = 5; %default value

    end

    % Temporary extend mask in order to avoid boarder effects when applying morphological dilation/erosion
    extended = extend_BW_mask(data_in,10*strel_size);

    % Detecting edges
    se0 = strel('line', strel_size, 0); %horz floor(1.5*process.linefill)/flake.imfinfo.XRes
    se90 = strel('line', strel_size, 90); %vert % 9.79 floor(1.5*process.linefill)/flake.imfinfo.XRes
    edged = edge(extended,'Sobel',0.008); %'log' is better but significantly slower       

    
    % Dilate and erode in order to blurry internal complexity
    dilated = imdilate(edged, [se0 se90]);
    filled = imfill(dilated,'holes');
    eroded = imerode(filled,[se0 se90]); 
    
    % Reduce the size of the mask to the initial one again
    data_out = shrink_BW_mask(eroded,10*strel_size); 
    
        
    % Illustration if desired
    if illustration 
        
            figure(6);
            subplot(3,2,1);
            imshow(data_in);
            title('raw');

            subplot(3,2,2);
            imshow(extended);
            title('extended');

            subplot(3,2,3);
            title('edges');
            imshow(edged);

            subplot(3,2,4);
            imshow(dilated);
            title('dilated');

            subplot(3,2,5);
            imshow(filled);
            title('filled');

            subplot(3,2,6);
            imshow(data_out);
            title('final');
            
    end
 
end