% Least Squares applied the pixels of the perimeter to fit.
% x,y = perimeter coordinates
function ellipse = fit_ellipse_ls(x,y,im_width,im_height,illustration)
  
    mean_x = mean(x);
    mean_y = mean(y);
    x = x-mean_x;
    y = y-mean_y;

    f = -1;
    X = [x.^2, x.*y, y.^2, x, y];
    A = sum(X)/(X'*X);
    [a,b,c,d,e] = deal(A(1),A(2),A(3),A(4),A(5));

    % ellipse equation :
    % A(1)*x^2 + A(2)*x*y + A(3)*y^2 + A(4)*x + A(5)*y + F = 0
    % let us transform this into a paramteric (plotable) form
    M0 = [f, d/2, e/2; ...
        d/2, a b/2; ...
        e/2, b/2, c];

    M = [a, b/2; ...
        b/2, c];

    %eigenvalues of M
    eigM = eig(M);
    if abs(eigM(1) - a) <= abs(eigM(1) -c)
        lambda1 = eigM(1);
        lambda2 = eigM(2); 
    else
        lambda1 = eigM(2);
        lambda2 = eigM(1);
    end

    % parameters of the parametric form
    u = sqrt(-det(M0)/(det(M)*lambda1)); % semi grand-axe
    v = sqrt(-det(M0)/(det(M)*lambda2)); % semi petit-axe
    
%     if u~=v
%         tmp_max = max(u,v);
%         tmp_min = min(u,v);
%         u = tmp_max;
%         v = tmp_min;
%     end
    
    
    h = (b*e - 2*c*d)/(4*a*c - b^2); % x-position of the center
    k = (b*d - 2*a*e)/(4*a*c - b^2); % y-position of the center
    theta = acot((a-c)/b)/2; % the orientation (tilt from x-axis)
    
    
    % reshift the whole ellipse to the original place
    h = h + mean_x;
    k = k + mean_y;

    % discretization of the ellipse outline on a 100 pts raster (100 should be enough but can be increased)
    t = linspace(0,2*pi);
    xt = h + cos(theta).*u.*cos(t) - sin(theta).*v.*sin(t);
    yt = k + sin(theta).*u.*cos(t) + cos(theta).*v.*sin(t);
    
    % compute the mask of the ellipse fit on the same raster as the original mask of the particle 
    try
        BW_fit_area = poly2mask(xt,yt,im_height,im_width);
    catch
        disp('Warning : impossible to rasterize ellipse fitted using LS method, something went wrong...');
        BW_fit_area = NaN;
    end
    
    ellipse = struct( ...
    'a',max(u,v),...
    'b',min(u,v),...
    'theta',theta,...
    'X0',h + mean_x,...
    'Y0',h + mean_y,...
    'X0_centered',h,...
    'Y0_centered',k,...
    'BW_fit_area',BW_fit_area,...
    'status','ok' );

%     if ellipse.a < ellipse.b
%         tmp =ellipse.a;
%         ellipse.a = ellipse.b;
%         ellipse.b = tmp;
%     end

    if illustration == 1



        % 2 extremities of the semi petit-axe
        P1 = [h - sin(theta)*v, k + cos(theta)*v]; 
        P2 = [h + sin(theta)*v, k - cos(theta)*v];

        % 2 extremitites of the semi grand-axe
        P3 = [h + cos(theta)*u, k + sin(theta)*u];
        P4 = [h - cos(theta)*u, k - sin(theta)*u];

        figure; hold on; axis equal;
        plot(xt,yt,'r-');
        plot(x+mean_x,y+mean_y,'kx');
        plot(h,k,'bx');
        plot(P1(1),P1(2),'b.');
        plot(P2(1),P2(2),'b.');
        plot(P3(1),P3(2),'b.'); 
        plot(P4(1),P4(2),'b.');
        plot([P4(1) P3(1)],[P4(2) P3(2)],'b--');
        plot([P2(1) P1(1)],[P2(2) P1(2)],'b--');
        set(gca,'Ydir','reverse');
        
        % illustration
        figure;  
        imshow(BW_fit_area); hold on;
        plot(xt,yt,'r-');
        plot(x+mean_x,y+mean_y,'kx');
        plot(h,k,'bx');
        plot(P1(1),P1(2),'b.');
        plot(P2(1),P2(2),'b.');
        plot(P3(1),P3(2),'b.'); 
        plot(P4(1),P4(2),'b.');
        plot([P4(1) P3(1)],[P4(2) P3(2)],'b--');
        plot([P2(1) P1(1)],[P2(2) P1(2)],'b--');        
        
        

    end

end