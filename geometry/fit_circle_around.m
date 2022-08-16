function circle=fit_circle_around(xP,yP,xh,yh,illustration)

    stepsize = 0.1;
    tol = 1;

    x0_new = mean(xh);
    y0_new = mean(yh);
    
    it = 0;
    r_old = inf;

    % heuristic algorithm start
    while true     
        
        % try to update the circle
        xd = xh - x0_new;
        yd = yh - y0_new;
        r_vec = sqrt(xd.^2+yd.^2);
        r_new = max(r_vec);
        idx_support = find(r_vec >= r_new - tol);
        n_support = length(idx_support);

        % if the new circle is better, we update radius and center
        if r_new < r_old

            normal = stepsize*[mean(xd(idx_support)), mean(yd(idx_support))]/norm([mean(xd(idx_support)),mean(yd(idx_support))]);
            x0_old = x0_new;
            y0_old = y0_new;
            x0_new = x0_new + normal(1);
            y0_new = y0_new + normal(2);
            it = it + 1;
            if illustration
                fprintf('iteration %u : %u support vertices, r = %2.2f, r_diff = %2.2f \n',it,n_support,r_new,r_old-r_new);
            end
            r_old = r_new;
            
        else
            
            r = r_old;
            x0 = x0_old;
            y0 = y0_old;
            break;
            
        end
   
    end
    
    
    circle = struct( ...
     'r',r,...
     'A',pi*r^2,...
     'X0',x0,...
     'Y0',y0,...
     'status','ok' );
      
    % illustration
    if illustration
        
       t = linspace(0,2*pi);
       xt = x0 + r.*cos(t);
       yt = y0 + r.*sin(t);
       
       figure; hold on; box on; axis equal;
       plot(xP,yP,'k.');
       plot(xt,yt,'r-');
       set(gca,'Ydir','reverse'); 
       title(sprintf('Area = %2.2f',pi*r^2));
       
    end
    
end