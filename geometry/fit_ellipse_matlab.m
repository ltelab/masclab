% fit_ellipse_matlab
% Matlab image processing toolbox function 
function ellipse = fit_ellipse_matlab(mask,illustration)

    E = regionprops(mask,'Area','MajorAxisLength', ...
    'Orientation','MinorAxisLength','Centroid'); 

    if length(E) > 1
        idx = find(max([E.Area]));
        E = E(idx);
    end
    
    [y,x] = find(mask==1);

    a = E.MajorAxisLength/2;
    b = E.MinorAxisLength/2;
    theta = E.Orientation;%*pi/180;
    centroid(1) = E.Centroid(1);
    centroid(2) = E.Centroid(2);

    ellipse = struct( ...
    'a',a,...
    'b',b,...
    'theta',theta,...
    'X0',centroid(1),...
    'Y0',centroid(2),...
    'status','ok' );

    if illustration
        
        
    % figure of the fit
    t = linspace(0,2*pi);
    xt = centroid(1) + cos(theta).*a.*cos(t) - sin(theta).*b.*sin(t);
    yt = centroid(2) + sin(theta).*a.*cos(t) + cos(theta).*b.*sin(t);
    figure; hold on; axis equal;
    plot(x,y,'b.');
    plot(xt,yt,'r-');
    legend('data','fit');
    set(gca,'Ydir','reverse');
    title(sprintf('theta = %2.1f',abs(ellipse.theta)));

    %figure of the overlap flake/ellipse area
%     [y_im,x_im] = find(mask > -1);
%     y_im = y_im - centroid(2);
%     x_im = x_im - centroid(1);
%     x_imr = cos(theta).*x_im - sin(theta).*y_im;
%     y_imr = sin(theta).*x_im + cos(theta).*y_im;
%     
%     [y_flake,x_flake] = deal(y,x);
%     y_flake = y_flake - centroid(2);
%     x_flake = x_flake - centroid(1);
%     x_flaker = cos(theta).*x_flake - sin(theta).*y_flake;
%     y_flaker = sin(theta).*x_flake + cos(theta).*y_flake;
% 
%     overlap = zeros(size(x_imr));
%     ellipse_only = zeros(size(x_imr));
%     flake_only = zeros(size(x_imr));
%     nothing = zeros(size(x_imr));
% 
%     for i=1:length(x_imr)
% 
%         is_in_flake =  ~isempty(find(x_flaker==x_imr(i),1)) && ~isempty(find(y_flaker==y_imr(i),1));
%         is_in_ellipse = x_imr(i)^2/a^2 + y_imr(i)^2/b^2 <= 1;
% 
%         if is_in_flake && is_in_ellipse
% 
%             overlap(i) = 1;
% 
%         elseif is_in_flake
% 
%             flake_only(i) = 1;
% 
%         elseif is_in_ellipse 
% 
%             ellipse_only(i) = 1;
% 
%         else
% 
%             nothing(i) = 1;
% 
%         end
% 
%     end   
    
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

end
