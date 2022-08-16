function [classId,rimingId,isMelting] = load_labels(dir_data,dir_labels,labels_filename,version)

    data_filenames = dir(fullfile(dir_data,'20*.mat'));
    % for Massimo data
    if isempty(data_filenames)
        data_filenames = dir(fullfile(dir_data,'ICE*.mat'));
    end
    data_filenames = {data_filenames.name}';
    %data_picnames = dir(fullfile(dir_data,'*.png'));
    %data_picnames = {data_picnames.name}';
    
    if strcmp(version,'1.0')

        load(fullfile(dir_labels,labels_filename));

        [~,idx_sorted] = sort(label.idx_rnd);
        flakenames_sorted = label.flakes(idx_sorted);
        classId = label.id(idx_sorted);
        
        % there is no rimingId and isMelting in version 1.0
        rimingId = nan(size(classId));
        isMelting = nan(size(classId));

        % small verifications
        if isequal(flakenames_sorted,data_filenames)
            disp('Labels correctly loaded and associated to the training sample');
        else
            disp(length(flakenames_sorted));
            disp(length(data_filenames));
            disp('Error : unable to associate labels with the given training sample');
            classId = NaN;
            rimingId = NaN;
            isMelting = NaN;
        end
        
    elseif strcmp(version,'1.1')
        
        load(fullfile(dir_labels,labels_filename));
        
        [~,idx_sorted] = sort(label.idx_rnd);
        flakenames_sorted = label.flakes(idx_sorted);
        classId = label.classId(idx_sorted);
        rimingId = label.rimingId(idx_sorted);
        isMelting = label.isMelting(idx_sorted);
        
         % small verifications
        if isequal(flakenames_sorted,data_filenames)
            disp('Labels correctly loaded and associated to the training sample');
        else
            disp(length(flakenames_sorted));
            disp(length(data_filenames));
            disp('Error : unable to associate labels with the given training sample');
            classId = NaN;
            rimingId = NaN;
            isMelting = NaN;
        end
        
    end
    
end
