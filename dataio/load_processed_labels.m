% some folders :
%'/home/praz/Documents/MASC/sample_snow_20150620/sample_cool/processed_05-Feb-2016/DATA/GOOD';
%'/media/praz/Masc-Data/APRES3_2015/PROCESSED_20160119/DATA/GOOD';


function [X,Xlab,Xname,Xt,Xfullprob_label,Xfullprob_riming] = load_processed_labels(dirname,t_str_start,t_str_stop,load_fullprob)

    if nargin < 4
        load_fullprob = 0;
    end

    fprintf('Load processed labels...');
    tmin = datenum(t_str_start,'yyyymmddHHMMSS');
    tmax = datenum(t_str_stop,'yyyymmddHHMMSS');
    file_list = dir(fullfile(dirname,'20*.mat'));
    % for Massimo data
    if isempty(file_list)
        file_list = dir(fullfile(dirname,'ICE*.mat'));
    end
    % for Vanderbilt triplets
    if isempty(file_list)
        file_list = dir(fullfile(dirname,'FLAKE*.mat'));
    end
    file_only_list = {file_list.name};
    file_list = fullfile(dirname,file_only_list);

    % loop over the snowflakes and save each feature in a vector of the main
    % structure data
    n_flakes = 0;
    X = [];
    Xlab = {};
    Xname = {};
    Xt = [];
    Xfullprob_label = [];
    Xfullprob_riming = [];

    
    for i=1:length(file_list)

        load(file_list{i});
        % load only pictures within the time interval
        if roi.tnum >= tmin && roi.tnum <= tmax

            j = 1;
            n_flakes = n_flakes + 1;
            
            Xname{i} = file_only_list{i};
            %disp(Xname{i});
            Xt(i) = roi.tnum;
           
            X(i,j) = roi.label_ID; Xlab{j} = 'label ID'; j=j+1;
            X(i,j) = roi.area; Xlab{j} = 'area'; j=j+1;
            X(i,j) = roi.Dmax; Xlab{j} = 'Dmax'; j=j+1;
            X(i,j) = roi.Dmean; Xlab{j} = 'Dmean'; j=j+1;
            X(i,j) = roi.complex; Xlab{j} = 'complex'; j=j+1;
            X(i,j) = roi.xhi; Xlab{j} = 'xhi'; j=j+1;
            if ~isempty(roi.fallspeed) 
                X(i,j) = roi.fallspeed; Xlab{j} = 'fallspeed'; j=j+1;
            else
                X(i,j) = NaN; Xlab{j} = 'fallspeed'; j=j+1;
            end
            X(i,j) = roi.eq_radius; Xlab{j} = 'eq_radius'; j=j+1;
            X(i,j) = roi.E.b/roi.E.a; Xlab{j} = 'aspect_ratio'; j=j+1;
            X(i,j) = roi.E.theta; Xlab{j} = 'orientation'; j=j+1;
            X(i,j) = max(roi.label_probs); Xlab{j} = 'label prob'; j=j+1;
            
            % new info
            X(i,j) = roi.riming_ID; Xlab{j} = 'riming ID'; j=j+1;
            X(i,j) = max(roi.riming_probs); Xlab{j} = 'riming prob'; j=j+1;
            X(i,j) = 1*roi.riming_probs(1) + 2*roi.riming_probs(2) + 3*roi.riming_probs(3) + 4*roi.riming_probs(4) + 5*roi.riming_probs(5); Xlab{j} = 'cont. riming ID'; j=j+1;
            X(i,j) = compute_riming_idx(1*roi.riming_probs(1) + 2*roi.riming_probs(2) + 3*roi.riming_probs(3) + 4*roi.riming_probs(4) + 5*roi.riming_probs(5)); Xlab{j} = 'riming index'; j=j+1;
            X(i,j) = roi.melting_ID; Xlab{j} = 'melting ID'; j=j+1;
            X(i,j) = roi.melting_probs; Xlab{j} = 'melting prob'; j=j+1;
            
            % more stuff to see repartition of particles on screen
            X(i,j) = roi.cam; Xlab{j} = 'cam'; j=j+1;
            X(i,j) = roi.centroid(1); Xlab{j} = 'centroid Xpos'; j=j+1;
            X(i,j) = roi.centroid(2); Xlab{j} = 'centroid Ypos'; j=j+1;
            X(i,j) = 0*roi.riming_probs(1) + 0.15*roi.riming_probs(2) + 0.5*roi.riming_probs(3) + 0.85*roi.riming_probs(4) + 1*roi.riming_probs(5); Xlab{j} = 'old bs riming index'; j=j+1;
            
            
            % stuff related to particle sharpness (beta)
            X(i,j) = roi.lap; Xlab{j} = 'lap'; j=j+1;
            X(i,j) = roi.wavs; Xlab{j} = 'wavs'; j=j+1;
            X(i,j) = roi.hist_entropy; Xlab{j} = 'histE'; j=j+1;
            X(i,j) = roi.local_std7; Xlab{j} = 'std7'; j=j+1;
            X(i,j) = roi.mean_intens; Xlab{j} = 'mean_intens'; j=j+1;
            X(i,j) = roi.max_intens; Xlab{j} = 'max_intens'; j=j+1;
            
            % number of particles on screen to potentially identify blowing snow
            X(i,j) = roi.n_roi; Xlab{j} = 'number of ROIs detected'; j=j+1;
            
            % particle roundness / area ratio
            X(i,j) = roi.roundness; Xlab{j} = 'area ratio/ roundness'; j=j+1;
            
            % new ARs
            X(i,j) = roi.D90.AR; Xlab{j} = 'AR D90'; j=j+1;
            X(i,j) = roi.Rect.aspect_ratio; Xlab{j} = 'AR Rect'; j=j+1;
            X(i,j) = roi.E_out.b/roi.E_out.a; Xlab{j} = 'AR E_out'; j=j+1;
            
            % new orientations
            X(i,j) = roi.Dmax_theta; Xlab{j} = 'orientation Dmax'; j=j+1;
            X(i,j) = roi.Rect.theta; Xlab{j} = 'orientation Rect'; j=j+1;
           
            
            
            % tmp fractal dimension (for paper)
            % X(i,j) = roi.F
            
            % new xhi test
            %X(i,j) = roi.blur_idx.xhi2; Xlab{j} = 'xhi2'; j=j+1;
            %X(i,j) = roi.blur_idx.xhi4; Xlab{j} = 'xhi4'; j=j+1;
            
            % ALWAYS KEEP flake ID as the last entry in X!
            X(i,j) = roi.id; Xlab{j} = 'flake ID'; j=j+1;
            
            
            % load full label probabilities (for triplet merging)
            if load_fullprob
                Xfullprob_label(i,:) = roi.label_probs;
                Xfullprob_riming(i,:) = roi.riming_probs;
            end
            
        end

    end
    fprintf('   Done! %u labelled images found.\n',n_flakes);
    Xlab = Xlab';
    Xname = Xname';
    Xt = Xt';
    
end
