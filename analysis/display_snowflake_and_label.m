% small script to display some flakes and their label
clear all; close all;
%dir_data = '/media/praz/MyData/MASC/processed_Davos_test_LTE/GOOD/2015.06.20/11';
%dir_data = '/media/praz/MyData/Massimo_data/ICE2_04_01_2015_08_08/DATA';

% validation data set
dir_data = '/home/praz/Documents/MASC/MASC_labelling/mixed_sample_part2';
is_validation = true;
load('CP_all_labelled_sorted_part2.mat');
%%%%%%%%%%%%%%%%%%%%%

file_list = dir(fullfile(dir_data,'20*.mat'));
% for Massimo data
if isempty(file_list)
    file_list = dir(fullfile(dir_data,'ICE*.mat'));
end
file_only_list = {file_list.name};
file_list = fullfile(dir_data,file_only_list);

if is_validation
    tmp = ismember(label.flakes,file_only_list);
    label.flakes(~tmp) = [];
    label.names(~tmp) = [];
    label.id(~tmp) = [];
    label.idx_rnd(~tmp) = [];
end


t_str_start = '20150101000000';
t_str_stop = '20170101000000';
[X,Xlab,Xname,Xt] = load_processed_labels(dir_data,t_str_start,t_str_stop);

slices(1,1) = sum(X==1);
slices(2,1) = sum(X==2);
slices(3,1) = sum(X==3);
slices(4,1) = sum(X==4);
slices(5,1) = sum(X==5);
slices(6,1) = sum(X==6);
%slices(7,1) = sum(X==0);
labels = {'grau','agg','melt','small','dend','col'};

figure;
pie(slices,labels);


%figure;
while true

    i = randi(numel(file_list));
    load(file_list{i});
    
    clf();
    subplot(3,3,[1 2 4 5 7 8]);
    imshow(roi.data);
    title(sprintf('Label : %s,  xhi : %2.2f',roi.label_name,roi.xhi)); 
    
    h_text = subplot(3,3,[3 6 9]);
    xl = xlim(h_text); 
    xPos = xl(1) + diff(xl) / 2; 
    yl = ylim(h_text); 
    yPos = yl(1) + diff(yl) / 2; 
    % text to display
    [probs_sorted,idx] = sort(roi.label_probs,'descend');
    labels_sorted = labels(idx);
    plot_text = text(xPos, yPos, sprintf('%s : %2.1f%% \n %s : %2.1f%% \n %s : %2.1f%% \n %s : %2.1f%% \n %s : %2.1f%% \n %s : %2.1f%%',...
        labels_sorted{1},probs_sorted(1)*100,labels_sorted{2},probs_sorted(2)*100,labels_sorted{3},probs_sorted(3)*100,...
        labels_sorted{4},probs_sorted(4)*100,labels_sorted{5},probs_sorted(5)*100,labels_sorted{6},probs_sorted(6)*100), 'Parent', h_text);
    set(plot_text, 'HorizontalAlignment', 'center','FontSize',12,'FontWeight','demi');
    set(h_text,'visible','off');       
    pause;

end

% check the labels of a validation data set
if is_validation
    
    y = label.id;
    y(y==6) = 5;
    y(y==7) = 6;    
    
    yPred = X;
    yPred(y==0) = [];
    y(y==0) = [];
 
    % Compute accuracy
    BER = computeBER(y,yPred);
    OA = computeOA(y,yPred);
    kappa = computeKAPPA(y,yPred);
    
end



