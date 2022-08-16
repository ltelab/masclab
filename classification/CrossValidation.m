function out = CrossValidation(method,type_classif,parameters_method,K,N_it,X,y,random,verbose,illustration,use_weights,apply_feat_transfo,feat_vec)

    [N,D] = size(X);
    N_classes = length(unique(y));
    y = double(y);
    
    % computing the weights, if we want to penalize the cost
    if use_weights
        if strcmp(type_classif,'binary')
            for i=1:N_classes
                Nc(i,1) = sum(y==(i-1));
            end
        else
            for i=1:N_classes
                Nc(i,1) = sum(y==i);
            end
        end
        [Ncmax,~] = max(Nc);
        for i=1:N_classes
            weights(i,1) = Ncmax/Nc(i);
        end
        parameters_method{end+1} = weights;
       
    end
    

    % Initialize global confusion matrix
    global_confmat_Te = zeros(N_classes);
    global_confmat_Tr = zeros(N_classes);

    % Initialize beta ranking vectors
    beta_sum_mat = [];

    for i = 1:N_it
        
        if ~random
            setSeed(i);
        else
            setSeed(randi(1000));
        end
        
        if verbose
            fprintf('\nStarting iteration number %.d%\n\n', i );
        end
        idx_perm = randperm(N);
        Nk = floor(N/K);
        idxCV=[];
        for k = 1:K
            idxCV(k,:) = idx_perm(1+(k-1)*Nk:k*Nk);
        end

        scoreTr_all = [];
        scoreTe_all = [];
        RiTr = [];
        RiTe = [];
        yTr_all = [];
        yTe_all = [];
        
        for k = 1:K
            if verbose
                fprintf('\nK-fold %.d/%.d% \n\n', k,K );
                fprintf('\n')
            end
            
            % get k'th subgroup in test, others in train
            idxTe = idxCV(k,:);
            idxTr = idxCV([1:k-1 k+1:end],:);
            idxTr = idxTr(:);
            yTe = y(idxTe);
            XTe = X(idxTe,:);
            yTr = y(idxTr);
            XTr = X(idxTr,:);
            
            % create synthetic training data
            
%             fprintf('size of XTr before = %u \n',numel(yTr));
%             
%             for j=1:numel(yTr)
%                 
%                 Xsample = XTr(j,:);
%                 ysample = yTr(j);
%                 for s=1:20
%                     Xnewsample = zeros(size(Xsample));
%                     sigma_sampling = 0.1;
%                     for l=1:numel(Xsample)
%                         Xnewsample(l) = normrnd(Xsample(l),sigma_sampling);
%                     end
%                     XTr(end+1,:) = Xnewsample;
%                     yTr(end+1) = ysample;
%                 end
%                 
%             end
%             
%             fprintf('size of XTr after = %u \n',numel(yTr));
            
            
            % if apply_feat_transfo, we do it inside the CV, on the train data and apply it to test data
            if apply_feat_transfo
                for j=1:D
                    skew = skewness(XTr(:,j));

                    if skew > 1

                        XTr(:,j) = log(abs(XTr(:,j)+1));
                        XTe(:,j) = log(abs(XTe(:,j)+1));

                    elseif skew > 0.75

                        XTr(:,j) = sqrt(abs(XTr(:,j)));
                        XTe(:,j) = sqrt(abs(XTe(:,j)));

                    elseif skew < -1

                        XTr(:,j) = exp(XTr(:,j));
                        XTe(:,j) = exp(XTe(:,j));

                    elseif skew < -0.75

                        XTr(:,j) = XTr(:,j).^2;
                        XTe(:,j) = XTe(:,j).^2;

                    end

                    tmp_mean = mean(XTr(:,j));
                    tmp_std  = std(XTr(:,j));
                    XTr(:,j) = (XTr(:,j) - tmp_mean)/tmp_std;
                    XTe(:,j) = (XTe(:,j) - tmp_mean)/tmp_std;

                end
            end
            

            % Fit model using the wrapping function
            model = trainClassifier(method, type_classif, yTr, XTr, parameters_method);

            % Features ranking
            if strcmp(method,'logistic')
                betas = abs(model);
                beta_sum = sum(betas,2);
                beta_sum = beta_sum(2:end,:);
                beta_sum_mat(:,end+1) = beta_sum;
            end

            % Predict
            [predTe,scoresTe] = predictClassifier(method,type_classif,model,XTe);
            [predTr,scoresTr] = predictClassifier(method,type_classif,model,XTr);  
            
            % Compute continuous score (for riming)
            if strcmp(method,'logistic') && strcmp(type_classif,'multiclass')
                scoreTr = scoresTr(:,1)*1 + scoresTr(:,2)*2 + scoresTr(:,3)*3 + scoresTr(:,4)*4 + scoresTr(:,5)*5;
                scoreTe = scoresTe(:,1)*1 + scoresTe(:,2)*2 + scoresTe(:,3)*3 + scoresTe(:,4)*4 + scoresTe(:,5)*5;
                scoreTr_all = [scoreTr_all; scoreTr];
                scoreTe_all = [scoreTe_all; scoreTe];
            end
            
            if strcmp(type_classif,'multiclass') && strcmp(method,'logistic')
                % Compute continuous riming index
                RiTr = [RiTr; 0.15*scoresTr(:,2) + 0.5*scoresTr(:,3) + 0.85*scoresTr(:,4) + 1*scoresTr(:,5)];
                RiTe = [RiTe; 0.15*scoresTe(:,2) + 0.5*scoresTe(:,3) + 0.85*scoresTe(:,4) + 1*scoresTe(:,5)];
                % Keep yTr/yTe for plots
                yTr_all = [yTr_all; yTr];
                yTe_all = [yTe_all; yTe];
            end

            % Compute accuracy
            BER_Te(i,k) = computeBER(yTe,predTe);
            BER_Tr(i,k) = computeBER(yTr,predTr);
            OA_Te(i,k) = computeOA(yTe,predTe);
            OA_Tr(i,k) = computeOA(yTr,predTr);
            kappa_Te(i,k) = computeKAPPA(yTe,predTe);
            kappa_Tr(i,k) = computeKAPPA(yTr,predTr);
            softOA_Te(i,k) = compute_softOA_riming(yTe,predTe);
            softOA_Tr(i,k) = compute_softOA_riming(yTr,predTr);
            if strcmp(method,'logistic') && strcmp(type_classif,'multiclass')
                RMSE_Te(i,k) = computeRMSE(yTe,scoreTe);
                RMSE_Tr(i,k) = computeRMSE(yTr,scoreTr);
            end
            
            % Compute Error Rate per class
            for c=1:N_classes
                if strcmp(type_classif,'binary')
                    idx_c = find(yTe==c-1);
                elseif strcmp(type_classif,'multiclass')   
                    idx_c = find(yTe==c);
                end
                ER_Te(i,k,c) = sum(yTe(idx_c)~=predTe(idx_c))/length(idx_c);
            end
                
            % Update global confusion matrix
            global_confmat_Te = global_confmat_Te + confusionmat(yTe,predTe);
            global_confmat_Tr = global_confmat_Tr + confusionmat(yTr,predTr);
            
            if verbose
                fprintf('\nTesting BER: %.2f%%\n', BER_Te(i,k) * 100 );
                fprintf('Training BER: %.2f%%\n\n', BER_Tr(i,k) * 100 );
            end

        end

    end


    %% diagnostic and illustration
    if verbose
        fprintf('*****************************\n');
        fprintf('Mean Testing BER    : %.2f%% pm %.2f%%\n', mean(BER_Te(:))*100, std(BER_Te(:))*100);
        fprintf('Mean Training BER   : %.2f%% pm %.2f%%\n\n', mean(BER_Tr(:))*100, std(BER_Tr(:))*100);
        fprintf('Mean Testing OA     : %.2f%% pm %.2f%%\n', mean(OA_Te(:))*100, std(OA_Te(:))*100);
        fprintf('Mean Training OA    : %.2f%% pm %.2f%%\n\n', mean(OA_Tr(:))*100, std(OA_Tr(:))*100);
        fprintf('Mean Testing Kappa  : %.2f%% pm %.2f%%\n', mean(kappa_Te(:))*100, std(kappa_Te(:))*100);
        fprintf('Mean Training Kappa : %.2f%% pm %.2f%%\n', mean(kappa_Tr(:))*100, std(kappa_Tr(:))*100);
        fprintf('*****************************\n');
    end
    
    
    % confusion matrices
    global_confmat_Te = round(global_confmat_Te./sum(global_confmat_Te(:)) * 100,2);
    global_confmat_Tr = round(global_confmat_Tr./sum(global_confmat_Tr(:)) * 100,2);
    
    % beta ranking
    beta_sum_vec = mean(beta_sum_mat,2);
    [beta_weight,beta_id] = sort(beta_sum_vec,'descend');
    beta_weight = beta_weight./sum(beta_weight)*100;

    if illustration
        
        if strcmp(type_classif,'multiclass')
            %labels = {'SP','CC','PC','AG','GR','CPC'};
            labels = {'Agg','Col','Gra','Ros','Sph','Oth'};
        else
            labels = {'dry','melting'};
        end
        % confmat
        figure('units','pixels','Position',[100 100 762 638]);
        heatmap(global_confmat_Tr',labels,labels,1,'Colormap','red','ShowAllTicks',1,'UseLogColorMap',true,'Colorbar',true,'Fontsize',14);
        set(gca,'Fontsize',14);
        title('Train confusion matrix');
        figure('units','pixels','Position',[100 100 762 638]);
        heatmap(global_confmat_Te',labels,labels,1,'Colormap','red','ShowAllTicks',1,'UseLogColorMap',true,'Colorbar',true,'Fontsize',14);
        set(gca,'Fontsize',14);
        title('Test confusion matrix');

        % beta weight
        figure('units','pixels','Position',[100 100 1015 544]); hold on; box on;
        bar(beta_weight);
        set(gca,'xtick',1:D,'xticklabel',int2str(feat_vec(beta_id)));
        xlabel('Feature ID');
        ylabel('Feature importance [%]');
        set(gca,'Fontsize',14);
        set(gca,'Linewidth',1.5);
        set(gca,'XLim',[0 D+1]);
        %title('Feature Ranking');
        
        
        %% boxplots
        if 0% strcmp(type_classif,'multiclass') && strcmp(method,'logistic')

            scatt = 0.15;
            scatt_vec = -scatt +(scatt+scatt).*rand(length(yTr_all),1);
            yTr_scatt = yTr_all + scatt_vec;
            scatt_vec = -scatt +(scatt+scatt).*rand(length(yTe_all),1);
            yTe_scatt = yTe_all + scatt_vec;
            

            % training set boxplot
            pause(0.5);
            grey_c = [89 89 89]/255;
            red_c = [240 59 32]/255;
            figure; hold on; box on; grid on;
            % plot all the datapoints
            plot(yTr_scatt,RiTr,'k.','markersize',5,'MarkerEdgeColor',grey_c);
            % generate the boxplot
            bh = boxplot(RiTr,yTr_all,'symbol','','positions',1:1:5,'widths',0.5,'colors','k');
            % fill the boxes with color
            h = findobj(gca,'Tag','Box');
            for j=1:length(h)
                patch(get(h(j),'XData'),get(h(j),'YData'),'b','FaceAlpha',.5);
            end
            % redraw the boxplot over it
            bh = boxplot(RiTr,yTr_all,'symbol','','positions',1:1:5,'widths',0.5,'colors','k');
            % thicken the outline of the boxplot
            set(bh(:,:),'linewidth',2);
            % set the color of the median line
            set(bh(6,:),'linewidth',2,'Color',red_c);
            % set the axis
            %axis([0.5 5.5 0.5 5.5]);
            %set(gca,'YTick',[1 2 3 4 5]);
            %set(gca,'XTick',[1 2 3 4 5]);
            xlabel('true riming degrees (discrete)');
            ylabel('predicted riming degree (continuous 1~5)');
            set(gca,'fontsize',12);
            title('Train sample');

            % test set boxplot
            figure; hold on; box on; grid on;
            % plot all the datapoints
            plot(yTe_scatt,RiTe,'k.','markersize',5,'MarkerEdgeColor',grey_c);
            % generate the boxplot
            bh = boxplot(RiTe,yTe_all,'symbol','','positions',1:1:5,'widths',0.5,'colors','k');
            % fill the boxes with color
            h = findobj(gca,'Tag','Box');
            for j=1:length(h)
                patch(get(h(j),'XData'),get(h(j),'YData'),'r','FaceAlpha',.5);
            end
            % redraw the boxplot over it
            bh = boxplot(RiTe,yTe_all,'symbol','','positions',1:1:5,'widths',0.5,'colors','k');
            % thicken the outline of the boxplot
            set(bh(:,:),'linewidth',2);
            % set the color of the median line
            set(bh(6,:),'linewidth',2,'Color',red_c);
            % set the axis
            axis([0.5 5.5 -0.1 1.1]);
            %set(gca,'YTick',[1 2 3 4 5]);
            set(gca,'XTick',[1 2 3 4 5]);%'XTickLabel',{'0','0.15','0.5','0.85','1'});
            set(gca,'YTick',[0 0.15 0.5 0.85 1]);
            set(gca,'Ydir','reverse');
            xlabel('Labelled R_d');
            ylabel('Predicted R_i');
            set(gca,'fontsize',14);
            
            figure;
            hold on; box on;
            histogram(yTe_all - scoreTe_all);
            set(gca,'Ytick',[]);
            set(gca,'Xlim',[-2 2]);
            xlabel('True label - Prediction');
            set(gca,'Fontsize',14);
            title('Error distribution');
            
            trueRi_Tr = 0.5*(sin(pi/4*(yTr_all-3))+1);
            trueRi_Te = 0.5*(sin(pi/4*(yTe_all-3))+1);
            
            figure;
            hold on; box on;
            histogram(trueRi_Te - RiTe);   
            %set(gca,'Ytick',[]);
            %set(gca,'Xlim',[-2 2]);
            xlabel('Labelled R_i - predicted R_i');
            ylabel('# of test samples');
            set(gca,'Fontsize',14);  
            fprintf('Riming index error mean : %2.2f \n',mean(trueRi_Te - RiTe));
            fprintf('Riming index error standard deviation : %2.2f \n',std(trueRi_Te - RiTe));
            fprintf('Interquartile range : %2.4f \n', quantile(trueRi_Te - RiTe,0.75)-quantile(trueRi_Te - RiTe,0.25));
        
        end

    end
    
    % output structure
    out.BER_Te = mean(BER_Te(:));
    out.BER_Tr = mean(BER_Tr(:));
    out.OA_Te = mean(OA_Te(:));
    out.OA_Tr = mean(OA_Tr(:));
    out.kappa_Te = mean(kappa_Te(:));
    out.kappa_Tr = mean(kappa_Tr(:));
    out.softOA_Te = mean(softOA_Te(:));
    out.softOA_Tr = mean(softOA_Tr(:));
    if strcmp(method,'logistic') && strcmp(type_classif,'multiclass')
        out.RMSE_Te = mean(RMSE_Te(:));
        out.RMSE_Tr = mean(RMSE_Tr(:));
    else
        out.RMSE_Te = NaN;
        out.RMSE_Tr = NaN;
    end
    out.beta_id = beta_id;
    out.beta_weight = beta_weight;
    out.ER_Te = ER_Te;
    
end