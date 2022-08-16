function rect = compute_rectangularity(x,y,perim,illustration)
    
    % load example
%     imdir = '/home/praz/Documents/MASC/masclab/events/dendrites';
%     datafiles = dir(fullfile(imdir,'*.mat'));
%     datafiles = {datafiles.name}'; 
%     load(fullfile(imdir,datafiles{i}));   
%     x = roi.x;
%     y = roi.y;
    
    % use John d'Errico's bounding objects function
    [rectx,recty,rectA,rectP] = minboundrect(x,y);
    P1 = [rectx(1) recty(1)];
    P2 = [rectx(2) recty(2)];
    P3 = [rectx(3) recty(3)];
    dist1 = pdist2(P1,P2);
    dist2 = pdist2(P2,P3);
    
    if dist1 > dist2
        vec = P2-P1;
        vec_0 = [1;0];
        theta = acos(dot(vec,vec_0)/(norm(vec)*norm(vec_0)));
        if theta > pi/2
            theta = pi - theta;
        end 
    else
        vec = P2-P3;
        vec_0 = [1;0];
        theta = acos(dot(vec,vec_0)/(norm(vec)*norm(vec_0)));
        if theta > pi/2
            theta = pi - theta;
        end  
    end
    rect.theta = theta;
    
    rect.width = min(dist1,dist2);
    rect.length = max(dist1,dist2);
    rect.A_ratio = length(x)/rectA;
    rect.p_ratio = rectP/perim;
    rect.aspect_ratio = rect.width/rect.length;
    rect.eccentricity = sqrt(1-rect.width/rect.length);
    rect.rectx = rectx;
    rect.recty = recty;
    
     if illustration
        figure; hold on; box on;
        plot(x,y,'kx');
        line(rectx,recty);
        set(gca,'Ydir','reverse','Fontsize',14); 
        axis equal;
        title(sprintf('w=%2.1f, l=%2.1f, Ar=%2.2f, Pr=%2.2f, theta=%2.1f',rect.width,rect.length,rect.A_ratio,rect.p_ratio,rect.theta*180/pi));
     end 

end