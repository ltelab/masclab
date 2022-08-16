function ellipse=fit_ellipse_inside(x,y,xP,yP,theta,illustration)

    % center the points around the origin
    centroid = [mean(x), mean(y)];
    xPc = xP - centroid(1);
    yPc = yP - centroid(2);
    xc = x - centroid(1);
    yc = y - centroid(2);
    
    theta = -theta;
    
    % retrieve the tilt from the PCA
%     [coeff,~,~] = pca([xc,yc]); 
%     theta = acos(coeff(1,1));
%     if coeff(2,1) < 0 
%         theta = -theta;
%     end 

    
    % untilt the data
    xPc0 = cos(-theta).*xPc - sin(-theta).*yPc;
    yPc0 = sin(-theta).*xPc + cos(-theta).*yPc;
    xc0 = cos(-theta).*xc - sin(-theta).*yc;
    
    % lets find the ellipse with the maximal area inscribed within the
    % given perimeter
    maxd = max(xc0) - min(xc0); % guess for maximal focal distance
    d = 0:1:maxd;
    s = zeros(length(d),1);
    A = zeros(length(d),1);
    for i=1:length(d)
        s(i) = min(sqrt((xPc0-0.5*d(i)).^2 + yPc0.^2) + sqrt((xPc0+0.5*d(i)).^2 + yPc0.^2));
        A(i) = pi*s(i)/4 * sqrt(s(i)^2-d(i)^2);
    end
    [~,idx] = max(A);
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
        plot(xt,yt,'r-');
        set(gca,'Ydir','reverse');
    end

end