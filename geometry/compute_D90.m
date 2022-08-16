function Dout = compute_D90(bw_mask_filled,Dmax_ini,Dmax_angle,verbose,illustration)
    
    % from the maximum diameter of a particle Dmax, compute the maximum
    % diameter in the perpendicular direction Dmax_90. Also recompute
    % Dmax_0 after image rotation for verification purposes.
    %
    % output : Dout.Dmax_0, Dout.Dmax_90, Dout.AR, Dout.angle
    
    % computation of Dmax_0 for angle rotation - candidate 1
    tmp_im = imrotate(bw_mask_filled,-Dmax_angle); % -angle -> clockwise
    Dmax_0 = sum(tmp_im,1);
    Dmax_0 = Dmax_0(Dmax_0>0);
    Dmax_0 = length(Dmax_0);
    
    % computation of Dmax_0 for angle rotation - candidate 2
%     tmp_im_C2 = imrotate(bw_mask_filled,Dmax_angle); % +angle -> counterclockwise
%     Dmax_0_C2 = sum(tmp_im_C2,1);
%     Dmax_0_C2 = Dmax_0_C2(Dmax_0_C2>0);
%     Dmax_0_C2 = length(Dmax_0_C2);
    
    % select best candidate and Dmax_90 verification
    err = (Dmax_0-Dmax_ini)/Dmax_ini * 100;
    %err_C2 = (Dmax_0_C2-Dmax_ini)/Dmax_ini * 100;
    
    Dmax_90 = sum(tmp_im,2);
    Dmax_90 = Dmax_90(Dmax_90 > 0);
    Dmax_90 = length(Dmax_90);
%     Dmax_90_C2 = sum(tmp_im_C2,2);
%     Dmax_90_C2 = Dmax_90_C2(Dmax_90_C2 > 0);
%     Dmax_90_C2 = length(Dmax_90_C2);
    
%     if Dmax_90_C1 > Dmax_0_C1 && Dmax_90_C2 > Dmax_0_C2
%         
%         fprintf('ERROR : D90 is larger than D0 for both angle options !!!');
%         alternative = 1;
%         
%     elseif Dmax_90_C1 > Dmax_0_C1
%         
%         alternative = 2;
%         
%     elseif Dmax_90_C2 > Dmax_0_C2
%         
%         alternative = 1;
%         
%     elseif abs(err_C1) > abs(err_C2)
%         
%         alternative = 2;
%         
%     else
%         
%         alternative = 1;
%         
%     end
%     
%     if alternative == 1
%         
%         tmp_im = tmp_im_C1;
%         err = err_C1;
%         Dmax_0 = Dmax_0_C1;
%         Dmax_90 = Dmax_90_C1;        
%     
%     elseif alternative == 2
%         
%         tmp_im = tmp_im_C2;
%         err = err_C2;
%         Dmax_0 = Dmax_0_C2;
%         Dmax_90 = Dmax_90_C2;
%         
%     end
    
    % output structure
    Dout = struct('Dmax_0',Dmax_0,'Dmax_90',Dmax_90,'AR',Dmax_90/Dmax_0,'angle',Dmax_angle);
    
    if verbose
        if abs(err) < 5
            fprintf('Error in Dmax_90 is %2.1f%% \n',err);
        else
            fprintf('WARNING : Error in Dmax_90 is %2.1f%% \n',err);
        end  
    end
    
    
      
    if illustration
        
        figure; hold on; axis equal; box on;
        subplot(121);
        imshow(bw_mask_filled);
        xlabel('original');
        subplot(122);
        imshow(tmp_im);
        xlabel(sprintf('after rotation by %2.1f degrees',Dmax_angle));
        
    end

end