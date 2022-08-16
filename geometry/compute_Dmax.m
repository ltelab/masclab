function [Dmax,theta,A,B] = compute_Dmax(xh,yh,illustration)
    
    % compute the smallest convex hull
    %k = convhull(xP,yP);
    
    % compute the maximal distance between two vertices of the convex hull
    % this represent for sure the maximal distance between 2 points
    % belonging to the particle
    % old version : Dmax = max(pdist([xP(k),yP(k)]));
    Dmax = 0;
    for i=1:length(xh)-1
        for j=i+1:length(xh)
            tmp_dist = sqrt((xh(i)-xh(j))^2 + (yh(i)-yh(j))^2);
            if tmp_dist > Dmax
                Dmax = tmp_dist;
                idxA = i;
                idxB = j;
            end
        end
    end
        
%     for i=1:length(k)-1
%         for j=i+1:length(k)
%             tmp_dist = sqrt((xP(k(i))-xP(k(j)))^2 + (yP(k(i))-yP(k(j)))^2);
%             if tmp_dist > Dmax
%                 Dmax = tmp_dist;
%                 idxA = k(i);
%                 idxB = k(j);
%             end
%         end
%     end
    
    % retrieve the 2 extremities of Dmax and compute angle formed with
    % horizontal 
    A = [xh(idxA); yh(idxA)];
    B = [xh(idxB); yh(idxB)];
    vec = B-A;
    vec_0 = [1;0];
    theta = acos(dot(vec,vec_0)/(norm(vec)*norm(vec_0)));
    if theta > pi/2
        theta = pi - theta;
    end
    
    % add - sign to this angle if we have to rotate counterclockwise to align the
    % main direction to the horizontal
    if (A(1) < B(1) && A(2) < B(2)) || (A(1) > B(1) && A(2) > B(2))
        theta = -theta;
    end
    
    % transform in degrees
    theta = theta*180/pi;
     
    if illustration
        
        figure; hold on; axis equal; box on;
        plot(xh,yh,'k--');
        %plot(xP(k),yP(k),'b--');
        plot([A(1) B(1)],[A(2) B(2)],'r--');
        xlabel('x [pix]');
        ylabel('y [pix]');
        set(gca,'Ydir','reverse','Fontsize',14); 
        axis equal;
        title(sprintf('D_{max} = %2.2f pixels. theta = %2.1f deg.',Dmax, theta));
        
    end

end
