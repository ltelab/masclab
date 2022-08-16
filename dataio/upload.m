% === upload MASc fallspeeds and image filenames found in dirname ========
%
% This function reads in the text files found in dirname and outputs a list
% of all images found in pic.files, each with a timestamp pic.time_vec, and
% a picture id pic.id, and a camera id pic.cam, and a fallspeed and
% fallspeed id pic.fallspeed and pic.fallid
%
% If the fallspeed information is not available, the function returns a
% vector of nan.
%
% Author : Christophe Praz (christophe.praz@epfl.ch)
% Last update : November 2015 (add granularity within an hour, as defined
% in label_params
% 
% ========================================================================
function pic = upload(dirname,label_params)

    % Import  all imagelist and time from imgInfo.txt
    try
        [imginfo delim nheaderlines] = importdata(fullfile(dirname,'imgInfo.txt')); 
        [imgrowN imgcolN] = size(imginfo.textdata);  
   
        % Read picture data from imgInfo.txt
        pic.time_vec = zeros(imgrowN-nheaderlines,6);
        pic.id = zeros(imgrowN-nheaderlines,1);
        pic.cam = zeros(imgrowN-nheaderlines,1);
        pic.files = cell(imgrowN-nheaderlines,1);

        for i = 1:length(pic.id)

            pic.id(i) = str2double(imginfo.textdata{i+nheaderlines,1});
            pic.cam(i) = str2double(imginfo.textdata{i+nheaderlines,2});
            day = imginfo.textdata{i+nheaderlines,3};
            timestamp = imginfo.textdata{i+nheaderlines,4};
            pic.files{i} = imginfo.textdata{i+nheaderlines,5};
            temp = cell2mat({day});
            temp2 = cell2mat({timestamp});

            % some bricolage to extract time_vec
            d = regexp( temp, '\.','split' );
            t = regexp( temp2, '\:','split' );                 
            picyr = str2double(d{3});
            picmo = str2double(d{1});
            picdd = str2double(d{2});
            pichh = str2double(t{1});
            picmm = str2double(t{2});
            picss = str2double(t{3});
            pic.time_vec(i,:) = [picyr picmo picdd pichh picmm picss]; %Matrix of timestamps 

        end
        
    catch
        pic.time_vec = [];
        pic.id = [];
        pic.cam = [];
        fprintf('ERROR : unable to open imgInfo.txt in %s.\n',dirname);
        return;
    end
    
    % Keep only the images within the time interval defined in label_params
    
    pic.time_num = datenum(pic.time_vec);
    t_min = datenum(label_params.starthr_vec);
    t_max = datenum(label_params.endhr_vec);
    idx2keep = find(pic.time_num >= t_min & pic.time_num <= t_max);
    if length(idx2keep) ~= length(pic.time_vec)
        pic.time_vec = pic.time_vec(idx2keep,:);
        pic.id = pic.id(idx2keep);
        pic.cam = pic.cam(idx2keep);
        pic.files = pic.files(idx2keep);
        pic.time_num = pic.time_num(idx2keep);
    end
        
    pic.time_str = datestr(pic.time_vec);

    % Import fallspeed data
    try
        [datainfo, delim, nheaderlines]= importdata(fullfile(dirname,'dataInfo.txt'));
        [datarowN, datacolN] = size(datainfo.data);

        % Read picture data from dataInfo.txt
        pic.fallspeed = datainfo.data; % fallspeed
        pic.fallid = zeros(datarowN-nheaderlines,1); %fallspeed id
        for i = 1:length(pic.fallspeed);
            pic.fallid(i) = str2double(datainfo.textdata{i+nheaderlines,1});
        end
        
        % Keep only the fallspeed within the time interval defined in
        % label_params
        [~,idx2keep,~] = intersect(pic.fallid,pic.id);
        pic.fallid = pic.fallid(idx2keep);
        pic.fallspeed = pic.fallspeed(idx2keep);
        
    catch 
        % /!\ TO CHECK ON MONDAY
        %pic.fallspeed = NaN*ones(length(pic.id)/3,1);
        pic.fallid = unique(pic.id);
        pic.fallspeed = NaN*ones(length(pic.fallid),1);
        fprintf('ERROR : unable to open dataInfo.txt in %s.\n',dirname);
        return;
    end
    
end
    


