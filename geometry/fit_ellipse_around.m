function ellipse=fit_ellipse_around(x,y,xh,yh,theta,illustration)

    % center the points around the origin
    centroid = [mean(x), mean(y)];
    xc = x - centroid(1);
    yc = y - centroid(2);
    
    % retrieve the tilt from the PCA
    %[coeff,~,~] = pca([xc,yc]); 
    %theta = acos(coeff(1,1));
    %if coeff(2,1) < 0 
    %   theta = -theta;
    %end
    
    theta = -theta;
    
    % untilt the data
    xc0 = cos(-theta).*xc - sin(-theta).*yc;
    yc0 = sin(-theta).*xc + cos(-theta).*yc;

    % instead of working on the whole data sample, it is easier and faster to
    % keep only the data points forming a polygon encompassing the data sample
    %k = convhull(x,y);
    %xh = x(k);
    %yh = y(k); 
    xhc = xh - centroid(1);
    yhc = yh - centroid(2);
    
    % untilt the hull
    xhc0 = cos(-theta).*xhc - sin(-theta).*yhc;
    yhc0 = sin(-theta).*xhc + cos(-theta).*yhc;

    % lets find the ellipse with the minimum area encompassing the hull in
    % function of d, the distance between the 2 focal points 
    maxd = max(xc0) - min(xc0); % guess for maximal focal distance
    d = 0:1:maxd;
    s = zeros(length(d),1);
    A = zeros(length(d),1);
    for i=1:length(d)
        s(i) = max(sqrt((xhc0-0.5*d(i)).^2 + yhc0.^2) + sqrt((xhc0+0.5*d(i)).^2 + yhc0.^2));
        A(i) = pi*s(i)/4 * sqrt(s(i)^2-d(i)^2);
    end
    [~,idx] = min(A);
    s = s(idx);
    d = d(idx);

    % from the optimal s and d, we retrieve a and b
    a = s/2;
    b = sqrt((s/2)^2 - (d/2)^2);
    
    ellipse = struct( ...
    'a',a,...
    'b',b,...
    'theta',-theta,...
    'X0',centroid(1),...
    'Y0',centroid(2),...
    'status','ok' );

    % illustration plot if desired
    if illustration == 1
        t = linspace(0,2*pi);
        xt = centroid(1) + cos(theta).*a.*cos(t) - sin(theta).*b.*sin(t);
        yt = centroid(2) + sin(theta).*a.*cos(t) + cos(theta).*b.*sin(t);

        figure; hold on; axis equal;
        plot(x,y,'k.');
        plot(xh,yh,'b--');
        plot(xt,yt,'r-');
        set(gca,'Ydir','reverse');
    end

end
