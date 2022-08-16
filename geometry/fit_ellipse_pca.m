% fit_ellipse_pca
% method based on pca
function ellipse = fit_ellipse_pca(x,y,Q,illustration)

% Q is the only free parameter : the quantile used to compute ratio a/b of the
% ellipse from the data points
%Q = 0.90;

centroid = [mean(x), mean(y)];
x_center = x - centroid(1);
y_center = y - centroid(2);

[coeff,~,~] = pca([x_center,y_center]); 

% retrieve the tilt
theta = acos(coeff(1,1));
%if coeff(2,1) < 0 
if theta > pi/2
    theta = pi-theta;
end


% I take all my data and I tilt them to match pc1-pc2 with x-y axis
xr = cos(-theta).*x_center - sin(-theta).*y_center;
yr = sin(-theta).*x_center + cos(-theta).*y_center;

% take a/b ratio as (x-axis maximal distance)/(x-axis maximal distance ratio)
a_min = quantile(xr,1-Q); a_max = quantile(xr,Q);
a = 0.5 * (a_max-a_min);
b_min = quantile(yr,1-Q); b_max = quantile(yr,Q);
b = 0.5 * (b_max-b_min);
axis_ratio = a/b;

% compute final a,b fixing area(ellipse) = area(snowflake)
a = sqrt(axis_ratio * length(x) / pi);
b = a / axis_ratio;

ellipse = struct( ...
    'a',a,...
    'b',b,...
    'theta',theta,...
    'X0',centroid(1),...
    'Y0',centroid(2),...
    'status','ok' );

theta = -0.5161;

% illustration if desired
if illustration == 1
    
    % figure of the fit
    t = linspace(0,2*pi);
    xt = centroid(1) + cos(theta).*a.*cos(t) - sin(theta).*b.*sin(t);
    yt = centroid(2) + sin(theta).*a.*cos(t) + cos(theta).*b.*sin(t);
    figure; hold on; axis equal;
    plot(x,y,'b.');
    plot(xt,yt,'r-');
    legend('data','fit');
    %set(gca,'Ydir','reverse');

    %figure of the overlap flake/ellipse area
    mask = ones(max(y),max(x));
    [y_im,x_im] = find(mask > -1);
    y_im = y_im - centroid(2);
    x_im = x_im - centroid(1);
    x_imr = cos(-theta).*x_im - sin(-theta).*y_im;
    y_imr = sin(-theta).*x_im + cos(-theta).*y_im;

    [y_flake,x_flake] = deal(y,x);
    y_flake = y_flake - centroid(2);
    x_flake = x_flake - centroid(1);
    x_flaker = cos(-theta).*x_flake - sin(-theta).*y_flake;
    y_flaker = sin(-theta).*x_flake + cos(-theta).*y_flake;

    overlap = zeros(size(x_imr));
    ellipse_only = zeros(size(x_imr));
    flake_only = zeros(size(x_imr));
    nothing = zeros(size(x_imr));

    for i=1:length(x_imr)

        is_in_flake =  ~isempty(find(x_flaker==x_imr(i),1)) && ~isempty(find(y_flaker==y_imr(i),1));
        is_in_ellipse = x_imr(i)^2/a^2 + y_imr(i)^2/b^2 <= 1;

        if is_in_flake && is_in_ellipse

            overlap(i) = 1;

        elseif is_in_flake

            flake_only(i) = 1;

        elseif is_in_ellipse 

            ellipse_only(i) = 1;

        else

            nothing(i) = 1;

        end

    end   
    
%     % compute accuracy and F1 score based on confusion matrix of the overlap
%     accuracy = (sum(overlap(:)) + sum(nothing(:))) / length(x_imr);
%     F1 = 2*sum(overlap(:)) / (2*sum(overlap(:)) + sum(ellipse_only(:)) + sum(flake_only(:)));
%     
%     figure; hold on; axis equal;
%     plot(x_imr(overlap==1),y_imr(overlap==1),'g.');
%     plot(x_imr(flake_only==1),y_imr(flake_only==1),'m.');
%     plot(x_imr(ellipse_only==1),y_imr(ellipse_only==1),'r.');
%     plot(x_imr(nothing==1),y_imr(nothing==1),'k.');
%     title(sprintf('accuracy : %1.2f \t F1 : %1.2f \n',accuracy,F1));
    
    
end


