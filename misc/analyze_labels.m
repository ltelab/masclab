% small script to analyze the similarities and discrepancies between labels
% made by different observers


%% small code to sort labels that are in a random order
if true
    clear all; close all;
    dir_labels = '/media/praz/MyData/Massimo_data';
    load(fullfile(dir_labels,'CP_labels_all.mat'));
    label_rnd = label;
    label_sorted = label_rnd;
    label_sorted.flakes(label_rnd.idx_rnd) = label_rnd.flakes;
    label_sorted.names(label_rnd.idx_rnd) = label_rnd.names;
    label_sorted.id(label_rnd.idx_rnd) = label_rnd.id;
    label_sorted.idx_rnd(label_rnd.idx_rnd) = label_rnd.idx_rnd;
    label = label_sorted;
end

%% main part
clear all; close all;
dir_labels = '/home/praz/Documents/MASC/masclab/labelling_flakes';
load(fullfile(dir_labels,'CP_all_labelled_sorted3.mat'));
l_CP = label;
load(fullfile(dir_labels,'JG_all_labelled_sorted.mat'));
l_JG = label;
load(fullfile(dir_labels,'DW_all_labelled_sorted.mat'));
l_DW = label;


% create a set of unanimous labels
if true
l_CP.id(l_CP.id==6) = 5;
l_JG.id(l_JG.id==6) = 5;
l_DW.id(l_DW.id==6) = 5;
l_CP.id(l_CP.id==7) = 6;
l_JG.id(l_JG.id==7) = 6;
l_DW.id(l_DW.id==7 | l_DW.id==8) = 6;

label.id = l_DW.id;
%label.id(l_CP.id ~= l_JG.id) = 0;
%label.id(l_CP.id ~= l_DW.id) = 0;
label.id(l_JG.id ~= l_DW.id) = 0;
end
% 67.36% all in common
% 86.23% CP-DW
% 73.43% CP-JG
% 70.43% DW-JG


%%

% pie plot
slices_labels = {'grau','agg','melt','small','dend','plate','col','bulros','other'};
spaces = {' ',' ',' ',' ',' ',' ',' ',' ',' '};
perc_symbols = {'%','%','%','%','%','%','%','%','%'};
%CP
slices_CP(1) = sum(l_CP.id==1);
slices_CP(2) = sum(l_CP.id==2);
slices_CP(3) = sum(l_CP.id==3);
slices_CP(4) = sum(l_CP.id==4);
slices_CP(5) = sum(l_CP.id==5);
slices_CP(6) = sum(l_CP.id==6);
slices_CP(7) = sum(l_CP.id==7);
slices_CP(8) = sum(l_CP.id==8);
slices_CP(9) = sum(l_CP.id==0);
for i=1:9
    perc_CP(i) = round(100*sum(slices_CP(i))/length(l_CP.id));
end
slices_labels_CP = strcat(slices_labels,spaces,num2strs(perc_CP),perc_symbols);
%JG
slices_JG(1) = sum(l_JG.id==1);
slices_JG(2) = sum(l_JG.id==2);
slices_JG(3) = sum(l_JG.id==3);
slices_JG(4) = sum(l_JG.id==4);
slices_JG(5) = sum(l_JG.id==5);
slices_JG(6) = sum(l_JG.id==6);
slices_JG(7) = sum(l_JG.id==7);
slices_JG(8) = sum(l_JG.id==8);
slices_JG(9) = sum(l_JG.id==0);
for i=1:9
    perc_JG(i) = round(100*sum(slices_JG(i))/length(l_JG.id));
end
slices_labels_JG = strcat(slices_labels,spaces,num2strs(perc_JG),perc_symbols);
%DW
slices_DW(1) = sum(l_DW.id==1);
slices_DW(2) = sum(l_DW.id==2);
slices_DW(3) = sum(l_DW.id==3);
slices_DW(4) = sum(l_DW.id==4);
slices_DW(5) = sum(l_DW.id==5);
slices_DW(6) = sum(l_DW.id==6);
slices_DW(7) = sum(l_DW.id==7 | l_DW.id==8);
slices_DW(8) = 0;
slices_DW(9) = sum(l_DW.id==0);
for i=1:9
    perc_DW(i) = round(100*sum(slices_DW(i))/length(l_DW.id));
end
slices_labels_DW = strcat(slices_labels,spaces,num2strs(perc_DW),perc_symbols);
% illustration
figure;
subplot(1,3,1);
pie(slices_CP,slices_labels_CP);
title('CP');
subplot(1,3,2);
pie(slices_JG,slices_labels_JG);
title('JG');
subplot(1,3,3);
pie(slices_DW,slices_labels_DW);
title('DW');

% merging plates and dendrites
l_CP.id(l_CP.id==6) = 5;
l_JG.id(l_JG.id==6) = 5;
l_DW.id(l_DW.id==6) = 5;
% Daniel: merge bullets and columns
l_DW.id(l_DW.id==8) = 7;

% facultative: remove undetermined
%idx = find(l_CP.id==0 | l_JG.id==0 | l_DW.id==0);
%l_CP.id(idx) = [];
%l_JG.id(idx) = [];
%l_DW.id(idx) = [];
% new labels
confmat_labels = {'other','grau','agg','melt','small','dend','col'};

% compute OA, kappa, BER
cpjg_ber = computeBER(l_CP.id,l_JG.id);
cpjg_OA  = computeOA(l_CP.id,l_JG.id);
cpjg_kap = computeKAPPA(l_CP.id,l_JG.id);

cpdw_ber = computeBER(l_CP.id,l_DW.id);
cpdw_OA  = computeOA(l_CP.id,l_DW.id);
cpdw_kap = computeKAPPA(l_CP.id,l_DW.id);

jgdw_ber = computeBER(l_JG.id,l_DW.id);
jgdw_OA  = computeOA(l_JG.id,l_DW.id);
jgdw_kap = computeKAPPA(l_JG.id,l_DW.id);



% CP - JG
cpjg = confusionmat(l_CP.id,l_JG.id);
figure;
heatmap(cpjg,confmat_labels,confmat_labels,1,'Colormap','red','ShowAllTicks',1,'UseLogColorMap',true,'Colorbar',true);
xlabel('Jacopo');
ylabel('Christophe');
title(sprintf('OA : %2.1f%%   BER : %2.1f%%   HSS : %2.1f%%',cpjg_OA*100,cpjg_ber*100,cpjg_kap*100));

% CP - DW
cpdw = confusionmat(l_CP.id,l_DW.id);
figure;
heatmap(cpdw,confmat_labels,confmat_labels,1,'Colormap','red','ShowAllTicks',1,'UseLogColorMap',true,'Colorbar',true);
xlabel('Daniel');
ylabel('Christophe');
title(sprintf('OA : %2.1f%%   BER : %2.1f%%   HSS : %2.1f%%',cpdw_OA*100,cpdw_ber*100,cpdw_kap*100));

% JG - DW
jgdw = confusionmat(l_JG.id,l_DW.id);
figure;
heatmap(jgdw,confmat_labels,confmat_labels,1,'Colormap','red','ShowAllTicks',1,'UseLogColorMap',true,'Colorbar',true);
xlabel('Daniel');
ylabel('Jacopo');
title(sprintf('OA : %2.1f%%   BER : %2.1f%%   HSS : %2.1f%%',jgdw_OA*100,jgdw_ber*100,jgdw_kap*100));








