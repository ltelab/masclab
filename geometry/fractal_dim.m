% === Compute snowflake fractal dimension ================================
%
% Compute image fractal dimension based on the box-counting algorithm.
% Fractal dimension D is computed as the slope of the log-log relationship
% between size of the box and number of boxes required to cover the area of
% the particle.
%
% edit : added nbox_max parameter allowing to consider only box size
% below 2^(nbox-1).
%
% edit: remove 2 last points in (nbox,rbox) (keep 4 points min)
%
% Author : Christophe Praz (christophe.praz@epfl.ch)
% Last update : November 2015
% ========================================================================
function F = fractal_dim(im,nbox_max,illustration)

    if nargin == 2
        
        illustration = false;
        
    elseif nargin == 1
        
        illustration = false;
        nbox_max = inf;
        
    end
    
    
    % boxcount
    [nbox,rbox] = boxcount(im);
    if length(nbox) > nbox_max
        nbox = nbox(1:nbox_max);
        rbox = rbox(1:nbox_max);
    % mod: we want to avoid the 2 last points, but still keep 4 points min. for the fit
    elseif length(nbox) > 5
        nbox = nbox(1:end-2);
        rbox = rbox(1:end-2);
    elseif length(nbox)==5
        nbox = nbox(1:end-1);
        rbox = rbox(1:end-1);
    end
    nbox = nbox';
    rbox = rbox';
    log_nbox = log(nbox);
    log_rbox = log(rbox);
    
    % least Squares fit
    beta = leastSquares(log_nbox,[ones(length(nbox),1) log_rbox]);
    F = -beta(2);
    
    % illustration
    if illustration
        
        xx = linspace(min(rbox),max(rbox),2);
        F1 = figure; 
        subplot(121);
        imshow(im);
        title('flake img');
        axis equal;
        subplot(122); hold on; box on; grid on;
        plot(rbox,nbox,'ks-');
        plot(xx,exp(beta(1)).*xx.^beta(2),'r--');
        v = axis;
        axis([min(rbox)/2 max(rbox)*2 v(3) v(4)]);
        xlabel('box size');
        ylabel('N boxes');
        title(sprintf('D = %2.2f',F));
        legend('data','LS fit');
        set(gca,'XScale','log','YScale','log');
        %Extract axes handles of all subplots from the figure
        axesHandles = findobj(get(F1,'Children'), 'flat','Type','axes');
        %Set the axis property to square
        axis(axesHandles,'square');
    
    end

end 
