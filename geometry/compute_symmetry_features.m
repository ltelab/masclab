function S = compute_symmetry_features(mask,Dmax,Req,illustration)

    % computation of the centroid
    s = regionprops(mask,'centroid');
    
    % computation of the perimeter (sorted by link)
    B = bwboundaries(mask);
    if length(B) > 1
        disp('Warning : symmetry routine detected more than 1 particle in the filled mask image !');
        disp('----> results are not reliable !');
    end
    B = B{1};
    xp = B(:,2);
    yp = B(:,1);
    
    % center the image (makes the calculations easier)
    xp = xp - s.Centroid(1);
    yp = yp - s.Centroid(2);
   
    % computation of the maximum centroid2boundary distance for all
    % directions
    angles = 0:1:359;
    for i=1:length(angles)
        xl = Dmax * cosd(angles(i));
        yl = Dmax * sind(angles(i));
        [allx,ally] = polyxpoly([0 xl],[0 yl],xp,yp);
        allnorm = sqrt(allx.^2+ally.^2);
        if ~isempty(allnorm)
            distances(i) = max(allnorm);
        else
            distances(i) = Req;
        end
    end
    
    % standardize distance values so that fft are on the same scale
    % it also removes the huge peak centered at 0 (corresponding to the
    % mean radius)
    tmp_mean = mean(distances);
    tmp_std = std(distances);
    distances = (distances-tmp_mean)/tmp_std;
    
    S.mean = tmp_mean;
    S.std = tmp_std;
    
    % fft stuff
    Fs = 1; % sampling frequency
    T = 1/Fs; % sampling period
    L = length(angles); % length of signal
    t = (0:L-1)*T; % time vector
    f_ = fft(distances);
    % two-sided spectrum P2
    P2 = abs(f_/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    fvec = 10*Fs*(0:(L/2))/L;
    
    S.P0 = P1(1);
    S.P1 = P1(2);
    S.P2 = P1(3);
    S.P3 = P1(4);
    S.P4 = P1(5);
    S.P5 = P1(6);
    S.P6 = P1(7);
    S.P7 = P1(8);
    S.P8 = P1(9);
    S.P9 = P1(10);
    S.P10 = P1(11);
    
    
    
    if illustration
        
        % shit for illustration (candidacy)
        % create bw mask
        bw = zeros(size(mask)); 
        bw(mask>0) = 1;  

        % smooth the silhouette
        bw_close = bwmorph(mask,'close');
        bw_close_open = bwmorph(bw_close,'open');

        % compute simple pseudo-skeleton
        bw_thin = bwmorph(bw_close_open,'thin',Inf);

        % remove potential leftover pixels
        all_roi = regionprops(bw_thin,'Area','PixelIdxList');   

        figure;
        %subplot(3,2,1);
        imshow(mask);
        %subplot(3,2,2); hold on; box on;
        figure; hold on;
        plot(xp,yp,'k-');
        plot(0,0,'kx');
        set(gca,'Ydir','reverse');
        axis equal;
        %subplot(3,2,3:4); hold on; box on;
        figure; hold on;
        plot(angles,distances,'r-');
        xlabel('Angle [degree]');
        ylabel('Standardized distance');
        set(gca,'Fontsize',14);
        %subplot(3,2,5:6); hold on; box on;
        figure; hold on;
        plot(fvec,P1);
        xlabel('Frequency');
        ylabel('Spectrum Power');
        set(gca,'Fontsize',14);
        
        figure;
        im_plot = bw_close_open;
        im_plot(bw_thin==1) = 0;
        imshow(im_plot);  
       
    end

end