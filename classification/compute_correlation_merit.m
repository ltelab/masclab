% compute the heuristics "merit" of a feature subset. the higher the merit,
% the better this subset is supposed to explain the output variable y
%
% method based on M. Hall thesis (1998)
% 
function merit = compute_correlation_merit(X,y,method_fc,method_fs)

    N = size(X,1);
    D = size(X,2);
    classes = unique(y);
    N_classes = length(classes);
    
    if strcmp(method_fc,'MDL')
        
        MDL = [];
        
        for k=1:D

            feat = X(:,k);
            class = X(:,1);


            N = size(X,1);
            D = size(X,2);
            C = length(unique(class));
            classes = unique(class);
            nc = [];

            for i=1:C
                nc(i) = sum(class==classes(i));
            end

            % Prior MDL
            p1 = log_facto(N);
            p2 = 0;
            for i=1:C
                p2 = p2 + log_facto(nc(i));
            end
            p3 = log_facto(N+C-1);
            p4 = log_facto(N);
            p5 = log_facto(C-1);
            prior_MDL = p1-p2+p3-p4-p5;

            % Post MDL
            p1 = 0;
            p2 = 0;
            p3 = 0;
            p4 = 0;
            p5 = 0;

            ufeat = unique(feat);

            for j=1:length(ufeat)

                nj = sum(feat==ufeat(j));
                p1 = p1 + log_facto(nj);
                for i=1:C
                    ncj(i) = sum(feat==ufeat(j) & class==classes(i));
                    p2 = p2 + log_facto(ncj(i));
                end
                p3 = p3 + log_facto(nj+C-1);
                p4 = p4 + log_facto(nj);
                p5 = p5 + log_facto(C-1);

            end
            post_MDL = p1-p2+p3-p4-p5;

            MDL(k) = (prior_MDL - post_MDL)/N;
            %MDL(k) = MDL(k)/(prior_MDL/N);

        end
        
        fc_score = MDL;
        
    end
        
    
    if strcmp(method_ff,'MDL')
        
        MDL = [];
        
        for k=1:D
            for l=1:D
            
                if k==l
                    MDL(k,l) = NaN;
                
                else
                    
                    feat = X(:,k);
                    class = X(:,l);

                    N = size(X,1);
                    D = size(X,2);
                    C = length(unique(class));
                    classes = unique(class);
                    nc = [];

                    for i=1:C
                        nc(i) = sum(class==classes(i));
                    end

                    % Prior MDL
                    p1 = log_facto(N);
                    p2 = 0;
                    for i=1:C
                        p2 = p2 + log_facto(nc(i));
                    end
                    p3 = log_facto(N+C-1);
                    p4 = log_facto(N);
                    p5 = log_facto(C-1);
                    prior_MDL = p1-p2+p3-p4-p5;

                    % Post MDL
                    p1 = 0;
                    p2 = 0;
                    p3 = 0;
                    p4 = 0;
                    p5 = 0;

                    ufeat = unique(feat);

                    for j=1:length(ufeat)

                        nj = sum(feat==ufeat(j));
                        p1 = p1 + log_facto(nj);
                        for i=1:C
                            ncj(i) = sum(feat==ufeat(j) & class==classes(i));
                            p2 = p2 + log_facto(ncj(i));
                        end
                        p3 = p3 + log_facto(nj+C-1);
                        p4 = p4 + log_facto(nj);
                        p5 = p5 + log_facto(C-1);

                    end
                    post_MDL = p1-p2+p3-p4-p5;

                    MDL(k,l) = (prior_MDL - post_MDL)/N;
                    %MDL(k) = MDL(k)/(prior_MDL/N);
                
                end
                
            end
        end
        
        ff_score = MDL(:);
        
        
    end
    
    
    
    
    % compute feature-class score
    if strcmp(method_fc,'centroidOA') || strcmp(method_fc,'centroidBER')
        
        for i=1:D
            feat = X(:,i);
            for j=1:N_classes
                centroids(j,1) = median(feat(y==classes(j)));
            end
            all_dist = pdist2(feat,centroids);
            [~,pred_class] = min(all_dist,[],2);
            pred_class = classes(pred_class);
            if strcmp(method_fc,'centroidOA')
                fc_score(i) = computeOA(y,pred_class);
            elseif strcmp(method_fc,'centroidBER')
                fc_score(i) = 1-computeBER(y,pred_class);
            end
        end
        
    elseif strcmp(method_fc,'intra-inter_var')
        
        for i=1:D
            feat = X(:,i);
            for j=1:N_classes
                idx_class = find(y==classes(j));
                class_std(j) = std(feat(idx_class));
                class_mean(j) = mean(feat(idx_class));
            end
            % intraclass mean std
            avg_intra_std = mean(class_std);
            % std of interclass mean
            std_inter_mean = std(class_mean);
            % score
            fc_score(i) = std_inter_mean/avg_intra_std;
              
        end
        
        fc_score = fc_score./max(fc_score);
        
        
    end
    
    
    % compute feature-feature score
    if strcmp(method_ff,'Paerson') || strcmp(method_ff,'Spearman')
        C = corr(X,'type',method_fs);
        C = triu(C,1);
        C = C(:);
        ff_score = C(C~=0);
    end

    % compute merit based on fc_score and ff_score
    %up = D*mean(fc_score)
    %down = sqrt(D + D*(D-1)*mean(ff_score))
    merit = D*mean(fc_score)/ sqrt(D + D*(D-1)*mean(ff_score));
    
end