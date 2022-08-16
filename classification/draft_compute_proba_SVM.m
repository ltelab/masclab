clear all; close all;

% load data
load('XyTrain_labelCP_1.mat');
N = size(X,1);
D = size(X,2);
N_classes = length(unique(y));

% define train/test samples
trainRatio = 0.9;
setSeed(1);
idx_perm = randperm(N);
idx_Tr = idx_perm(1:floor(N*trainRatio));
idx_Te = idx_perm(floor(N*trainRatio)+1:end);
XTr = X(idx_Tr,:);
yTr = y(idx_Tr);
XTe = X(idx_Te,:);
yTe = y(idx_Te);


% train SVM 1vs all
model = {};
y_tmp = yTr;

for j=1:N_classes

    y_tmp(yTr~=j) = -1;
    y_tmp(yTr==j) = +1;   

    SVMModel = fitcsvm(XTr,y_tmp,'BoxConstraint',1,'KernelFunction','linear','IterationLimit',1e4,'Verbose',0); %'auto'
    % update the model to fit posterior probabilities
    model{j} = fitPosterior(SVMModel);

end 


% classify train sample
scoresTr = [];
for j=1:N_classes

    [label,score] = predict(model{j},XTr);
    scoresTr(:,j) = score(:,2); 
end
[confindexTr,predTr] = max(scoresTr,[],2);

% compute "probabilities"
tmp_sum = sum(scoresTr,2);
for j=1:size(scoresTr,1)
    
    probTr(j,:) = scoresTr(j,:)./tmp_sum(j);
    
end

       
% classify test sample
scoresTe = [];
for j=1:N_classes

    [label,score] = predict(model{j},XTe);
    scoresTe(:,j) = score(:,2);   
    
end
[confindexTe,predTe] = max(scoresTe,[],2);

% compute "probabilities
tmp_sum = sum(scoresTe,2);
for j=1:size(scoresTe,1)
    
    probTe(j,:) = scoresTe(j,:)./tmp_sum(j);
    
end

% accuracy
BER_Te = computeBER(yTe,predTe);
BER_Tr = computeBER(yTr,predTr);
OA_Te = computeOA(yTe,predTe);
OA_Tr = computeOA(yTr,predTr);
kappa_Te = computeKAPPA(yTe,predTe);
kappa_Tr = computeKAPPA(yTr,predTr);
fprintf('*****************************\n');
fprintf('Testing BER    : %.2f%%\n', BER_Te*100);
fprintf('Training BER   : %.2f%%\n\n', BER_Tr*100);
fprintf('Testing Kappa  : %.2f%%\n', kappa_Te*100);
fprintf('Training Kappa : %.2f%%\n', kappa_Tr*100);
fprintf('*****************************\n');

%%
% path to train sample to retrieve images
% dir_data = '/home/praz/Documents/MASC/masclab/events/mixed_sample_final';
% data_picnames = dir(fullfile(dir_data,'*.png'));
% data_picnames = {data_picnames.name}';

% display some diagnostics /!\ labels have to be sorted !!!
idx_fail_Te = idx_Te(predTe~=yTe);
idx_fail_Te_pred_label = predTe(predTe~=yTe);
idx_fail_Te_real_label = yTe(predTe~=yTe);
idx_fail_Te_confindex = confindexTe(predTe~=yTe);
idx_fail_Te_prob = probTe(predTe~=yTe,:);

figure;
for i=1:length(idx_fail_Te)
    
    clf();
    imgname = data_picnames{idx_fail_Te(i)};
    img = imread(imgname);
    subplot(3,1,1);
    imshow(img);
    h_text = subplot(3,1,2:3);
    xl = xlim(h_text);
    xPos = xl(1) + diff(xl) / 2; 
    yl = ylim(h_text); 
    yPos = yl(1) + diff(yl) / 2; 
    plot_text = text(xPos, yPos, sprintf('Real : %u\n Pred : %u\n conf. index : %2.2f\n graupel : %2.2f\n aggregates : %2.2f\n melting snow : %2.2f\n small particle : %2.2f\n dendrite-like : %2.2f\n column-like : %2.2f\n', ...
        idx_fail_Te_real_label(i),idx_fail_Te_pred_label(i),idx_fail_Te_confindex(i),idx_fail_Te_prob(i,1),idx_fail_Te_prob(i,2),idx_fail_Te_prob(i,3),idx_fail_Te_prob(i,4),idx_fail_Te_prob(i,5),idx_fail_Te_prob(i,6)), 'Parent', h_text);
    set(plot_text, 'HorizontalAlignment', 'center','FontSize',12,'FontWeight','demi');
    set(h_text,'visible','off');
    pause;
    
    
end






% figure;
% for i=1:length(idx_fail_final)
% 
%     clf();
%     imgname = data_picnames{idx_fail_final(i)};
%     img = imread(imgname);
%     subplot(2,3,1);
%     imshow(img);
%     h_text = subplot(2,3,2:3);
%     xl = xlim(h_text); 
%     xPos = xl(1) + diff(xl) / 2; 
%     yl = ylim(h_text); 
%     yPos = yl(1) + diff(yl) / 2; 
%     plot_text = text(xPos, yPos, sprintf('Pred : %u. Real : %u. Occurrences : %2.0f/%u',idx_fail_pred_label(i),idx_fail_real_label(i),idx_fail_occurence(i),N_it), 'Parent', h_text);
%     set(plot_text, 'HorizontalAlignment', 'center','FontSize',12,'FontWeight','demi');
%     set(h_text,'visible','off');  
%     h_text2 = subplot(2,3,4:6);
%     xl = xlim(h_text); 
%     xPos = xl(1) + diff(xl) / 2; 
%     yl = ylim(h_text); 
%     yPos = yl(1) + diff(yl) / 2; 
%     plot_text2 = text(xPos,yPos,sprintf('%s',imgname),'interpreter','none');
%     set(plot_text2, 'HorizontalAlignment', 'center','FontSize',12,'FontWeight','demi');
%     set(h_text2,'visible','off'); 
%     %title(sprintf('Pred : %u. Real : %u. Occurrences : %u/%u',idx_fail_pred_label(i),idx_fail_real_label(i),idx_fail_occurence(i),N_it));
%     pause;
% 
% end

