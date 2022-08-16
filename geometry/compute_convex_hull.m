% compute convex hull (x,y) coordinates, solidity and convexity of a shape
function hull = compute_convex_hull(x,y,perim,illustration)
    
    % load example
%     clear all; close all;
%     illustration = true;
%     imdir = '/home/praz/Documents/MASC/masclab/events/columns';
%     datafiles = dir(fullfile(imdir,'*.mat'));
%     datafiles = {datafiles.name}'; 
%     load(fullfile(imdir,datafiles{j}));   
%     x = roi.x;
%     y = roi.y;
    
    try
        k = convhull(x,y);
    catch
        hull = NaN;
        disp('Warning : convhull failed, hull = NaN');
        return;
    end
    
    hull.xh = x(k);
    hull.yh = y(k);
    hull.solidity = length(x)/polyarea(hull.xh,hull.yh);
    % hull perimeter computation
    hull.perim = 0;
    for i=1:length(hull.xh)
        p1 = [hull.xh(i),hull.yh(i)];
        if i==length(hull.xh)
            p2 = [hull.xh(1),hull.yh(1)];
        else
            p2 = [hull.xh(i+1),hull.yh(i+1)];
        end
        hull.perim = hull.perim+pdist2(p1,p2);
    end
    hull.convexity = hull.perim/perim;
  
    if illustration
       figure; hold on; box on;
       plot(x,y,'kx');
       plot(hull.xh,hull.yh,'b--');
       title(sprintf('solidity=%2.2f, convexity=%2.2f',hull.solidity,hull.convexity));
    end 

 end