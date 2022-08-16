function Rect = drawRectangles(feat_id,height,facecolor)
    for i=1:length(feat_id)
        rectangle('Position',[feat_id(i)-0.5 height 1 1],'Curvature',[0.5,1],'Facecolor',facecolor,'Edgecolor','k');
    end
end