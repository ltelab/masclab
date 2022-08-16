% ========== useful tool to load all processed and classified .mat files into a single data structure ==========
%
% inputs :
% ========
% data_dir    : path to the folder where the processed and classified MASC .mat files are located
% tstart_str  : start time string 'yyyymmddHHMMSS' (MASC images recorded before this will be ignored)
% tstop_str   : stop time string 'yyyymmddHHMMSS' (MASC images recorded before this will be ignored)
%
% outputs :
% =========
% data        : a matlab structure containing different fields
%   data.X contains classification data for each MASC image separatly (each line of the matrix = 1 MASC image)
%   data.Y contains classification data after merging the probabilities on the 3 views of the same triplet (each line = 1 triplet)
%   data.Xlab / data.Ylab lists what the different columns of data.X/Y correspond to
%   data.Xt/Yt contains a datenum vector containing the timestamps of data.X/data.Y
%   data.Xname/Yname contains a string vector containing the name of the images
%
% Author : Christophe Praz (christophe.praz@epfl.ch)
% Last update : March 2018
% ============================================================================================================
function data = merge_predictions_for_campaign(data_dir,tstart_str,tstop_str)
%data_dir = '/media/praz/MyData/MASC/processed_APRES3/processed_single/GOOD';
%tstart_str = '20151112000000';
%tstop_str = '20151113000000';

    tstart_vec = datevec(tstart_str,'yyyymmddHHMMSS');
    tstart = datenum(tstart_vec);
    tstop_vec  = datevec(tstop_str,'yyyymmddHHMMSS');
    tstop  = datenum(tstop_vec);
    res = 33.5; % 1 pixel = res [mu]

    % load all data and store in large matrices (X)
    X = [];
    Xname = {};
    Xt = [];
    Xfullprob_label = [];
    Xfullprob_riming = [];
    dir_list = uploaddirs(data_dir,tstart_vec,tstop_vec);
    perc_old = 0;
    for i=1:numel(dir_list)

        fprintf('%u / %u\n',i,numel(dir_list));
        
        % if necessary, we process new descriptors (new blur index)
        % process_new_descriptors(dir_list{i},tstart_str,tstop_str);
           
        [X_tmp,Xlab_tmp,Xname_tmp,Xt_tmp,Xfullprob_label_tmp,Xfullprob_riming_tmp] = load_processed_labels(dir_list{i},tstart_str,tstop_str,1);
        if ~isempty(X_tmp)
            X = [X; X_tmp];
            Xname = [Xname; Xname_tmp];
            Xt = [Xt; Xt_tmp];
            Xfullprob_label = [Xfullprob_label; Xfullprob_label_tmp];
            Xfullprob_riming = [Xfullprob_riming; Xfullprob_riming_tmp];
            Xlab = Xlab_tmp;
        end
        perc_current = floor(100*i/numel(dir_list));
        if perc_current > perc_old
           fprintf('%u %% done \n',perc_current);
           perc_old = perc_current;
        end

    end
    
    % if there is no particle found within the time interval, return an empty structure
    if isempty(X)
        data.X = [];
        data.Xt = [];
        data.Xlab = {};
        data.Xname = {};
        data.Xfullprob_label = [];
        data.Y = [];
        data.Yt = [];
        data.Ylab = {};
        data.Yname = {};
        return;
    end

    % creation of new  snowflake IDs starting at 1 and ensuring that each ID is unique
    diffID = [0;diff(X(:,end))];
    newID = zeros(length(diffID),1);
    newID(1) = 1;
    for i=2:length(diffID)
        if diffID(i) == 0
            newID(i) = newID(i-1);
        else
            newID(i) = newID(i-1) + 1;
        end
    end
    X(:,end+1) = newID;
    Xlab{end+1} = 'new uniqueID';

    % merging the triplets found in X to create a new merged matrix Y
    % for each flake ID, we find all views and average stuff over it
    fprintf('\nAveraging the classification over the different views of a same snowflake...');
    uniqueID = unique(X(:,end));
    Y = zeros(length(uniqueID),size(X,2)+1);
    Y(Y==0) = NaN;
    Yfullprob_label = zeros(length(uniqueID),size(Xfullprob_label,2));
    Yt = zeros(length(uniqueID),1);
    Yname = cell(length(uniqueID),1);

    for i=1:length(uniqueID)
        idx_views = find(X(:,end)==uniqueID(i));
        % label and label probability matrix
        if numel(idx_views) > 1
            tmp_probs = nansum(Xfullprob_label(idx_views,:));
            Yfullprob_label(i,:) = tmp_probs./sum(tmp_probs);
            [~,Y(i,1)] = nanmax(Yfullprob_label(i,:));
            % exclusive classification - not used!
%             if length(unique(X(idx_views,1))) == 1
%                 Y(i,1) = unique(X(idx_views,1));
%             else
%                 Y(i,1) = NaN;
%             end
            
        else
            Yfullprob_label(i,:) = Xfullprob_label(idx_views,:);
            Y(i,1) = X(idx_views,1);
        end
        % riming index
        Y(i,15) = nanmean(X(idx_views,15));
        % melting ID and probs
        Y(i,17) = nanmean(X(idx_views,17));
        if Y(i,17) > 0.5
            Y(i,16) = 1;
        else
            Y(i,16) = 0;
        end
        % other (averageable) descriptors
        if numel(idx_views) > 1
            Y(i,[2,3,4,5,6,7,8,9,10,14,19,20,21:27,29]) = nanmean(X(idx_views,[2,3,4,5,6,7,8,9,10,14,19,20,21:27,29]));
        else
            Y(i,[2,3,4,5,6,7,8,9,10,14,19,20,21:27,29]) = X(idx_views,[2,3,4,5,6,7,8,9,10,14,19,20,21:27,29]);
        end
        % Dmax (3)
        Y(i,3) = nanmax(X(idx_views,3));
        % complexity (5)
        Y(i,5) = nanmax(X(idx_views,5));
        % AR (9)
        Y(i,9) = nanmin(X(idx_views,9));
        % angle (10)
        Y(i,10) = nanmin(abs(X(idx_views,10)));
        
        
        Y(i,end-2) = X(idx_views(1),end-1);
        % number of cams found for each triplet
        Y(i,end) = numel(idx_views);
        % copy Xname and Xt
        Yt(i) = Xt(idx_views(1));
        tmp_name = Xname{idx_views(1)};
        Yname{i} = tmp_name(1:end-10);

    end
    % add the newlabelID
    Y(:,end-1) = uniqueID;
    % copy Xlab in Ylab
    Ylab = Xlab;
    Ylab{end+1} = 'number of views';
    fprintf(' Done! \n');

    data.X = X;
    data.Xt = Xt;
    data.Xlab = Xlab;
    data.Xname = Xname;
    data.Xfullprob_label = Xfullprob_label;
    data.Y = Y;
    data.Yt = Yt;
    data.Ylab = Ylab;
    data.Yname = Yname;
    
end

   


     
    


