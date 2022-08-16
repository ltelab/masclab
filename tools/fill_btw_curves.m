% x, y1 and y2 are single column vectors
% y1 : lower bound
% y2 : upper bound
% y3 : middle line
% the function can handles non-continuous time-series (NaN) and will fill
% only the continuous parts
function fill_btw_curves(x,y1,y2,y3,C,alpha)
 
    while ~isempty(x) && sum(isnan(y1)) ~= numel(y1)
    
        idx_start = find(~isnan(y1),1,'first');
        idx_stop = find(isnan(y1(idx_start:end)),1,'first')-1;
        
        if idx_stop ~= idx_start

            xd = x(idx_start:idx_stop);
            y1d = y1(idx_start:idx_stop);
            y2d = y2(idx_start:idx_stop);
            y3d = y3(idx_start:idx_stop);
            
            fill([xd; flipud(xd)], [y1d; flipud(y2d)],C,'edgecolor','none','facealpha',alpha);
            plot(xd, y3d, 'Color', C, 'linewidth', 2);
            
        else
           
            xd = x(idx_start);
            y1d = y1(idx_start);
            y2d = y2(idx_start);
            y3d = y3(idx_start);
            
            plot([xd xd],[y1d y2d],'Color',C);
            plot(xd, y3d,'ko','MarkerFaceColor',C);
            

        end

        x(1:idx_stop) = [];
        y1(1:idx_stop) = [];
        y2(1:idx_stop) = [];
        y3(1:idx_stop) = [];
  
    end
       

end