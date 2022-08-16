function create_proc_params_file(pathname,label,process)

    filename = 'proc_params.txt';
    fileID = fopen(fullfile(pathname,filename),'w');
    
    fprintf(fileID,'Processing parameters associated with data in this folder \n\n\n');
    
    
    fprintf(fileID,'creation date          : %s\n',date);
    fprintf(fileID,'data processed         : %s\n',label.campaigndir);
    fprintf(fileID,'starting time          : %s\n',datestr(label.starthr_vec));
    fprintf(fileID,'ending time            : %s\n',datestr(label.endhr_vec));
    fprintf(fileID,'processed in parallel  : %u\n\n',process.parallel);
    fprintf(fileID,'backtresh limit        : %u\n',process.backthresh);
    fprintf(fileID,'size min               : %u\n',process.sizemin);
    fprintf(fileID,'min area               : %u\n',process.min_area);
    fprintf(fileID,'min brightness         : %2.2f\n',process.minbright);
    fprintf(fileID,'max intens thresh      : %2.2f\n',process.max_intensthresh);
    fprintf(fileID,'min hole area          : %u\n',process.min_hole_area);
    fprintf(fileID,'discardmat [t,b,l,r]   : [%u %u %u %u]\n\n',process.discardmat(1),process.discardmat(2),process.discardmat(3),process.discardmat(4));
    fprintf(fileID,'use triplet algo       : %u\n',process.use_triplet_algo);
    if process.use_triplet_algo
        fprintf(fileID,'matching tol. pix.     : %u\n',process.matching_tol_pix);
        fprintf(fileID,'matching tol. percent. : %u\n',process.matching_tol_percent);   
    end
    
    fclose(fileID);


end