% small script to analyze the similarities and discrepancies between labels
% made by different observers

%% tmp script to remove 2 wrong entries in Jacopo labels
if false
    clear tmp_vec
    idx2rmv = [];
    for i=1:length(l_JG.flakes)
        tmp_vec = strfind(l_CP.flakes,l_JG.flakes{i});
        tmp_vec = tmp_vec(:);
        tmp_vec = tmp_vec(~cellfun(@isempty, tmp_vec));
        if isempty(tmp_vec)
            idx2rmv(end+1) = i;
        end
    end
end


%% small code to sort labels that are in a random order
if false
    clear all; close all;
    dir_labels = '/home/praz/Documents/MASC/masclab/labelling_GUI/labels_v.1.1';
    load(fullfile(dir_labels,'CP_v11_part5.mat'));
    label_rnd = label;
    label_sorted = label_rnd;
    label_sorted.flakes(label_rnd.idx_rnd) = label_rnd.flakes;
    label_sorted.className(label_rnd.idx_rnd) = label_rnd.className;
    label_sorted.classId(label_rnd.idx_rnd) = label_rnd.classId;
    label_sorted.rimingName(label_rnd.idx_rnd) = label_rnd.rimingName;
    label_sorted.rimingId(label_rnd.idx_rnd) = label_rnd.rimingId;
    label_sorted.isMelting(label_rnd.idx_rnd) = label_rnd.isMelting;
    label_sorted.idx_rnd(label_rnd.idx_rnd) = label_rnd.idx_rnd;
    label = label_sorted;
end

%% small code to merge sorted labels part1 and part2
if true
    clear all; close all;
    dir_labels = '/home/praz/Documents/MASC/masclab/labelling_GUI/labels_v.1.1';
    load(fullfile(dir_labels,'FINAL2_corrected.mat'));
    label_p1 = label;
    load(fullfile(dir_labels,'CP_v11_part5_s.mat'));
    label_p2 = label;
    % merging
    clear label;
    label.flakes = [label_p1.flakes; label_p2.flakes];
    label.className = [label_p1.className; label_p2.className];
    label.classId = [label_p1.classId; label_p2.classId];
    label.rimingName = [label_p1.rimingName; label_p2.rimingName];
    label.rimingId = [label_p1.rimingId; label_p2.rimingId];
    label.isMelting = [label_p1.isMelting; label_p2.isMelting];
    label.idx_rnd = [label_p1.idx_rnd label_p2.idx_rnd];
end

%% small code to merge FINAL1 with part 3 and part 4
if false
    clear all; close all;
    dir_labels = '/home/praz/Documents/MASC/masclab/labelling_GUI/labels_v.1.1';
    load(fullfile(dir_labels,'FINAL1.mat'));
    label_p1 = label;
    load(fullfile(dir_labels,'CP_v11_part3_s.mat'));
    label_p2 = label;
    load(fullfile(dir_labels,'CP_v11_part4_s.mat'));
    label_p3 = label;
    % merging
    clear label;
    label.flakes = [label_p1.flakes; label_p2.flakes; label_p3.flakes];
    label.className = [label_p1.className; label_p2.className; label_p3.className];
    label.classId = [label_p1.classId; label_p2.classId; label_p3.classId];
    label.rimingName = [label_p1.rimingName; label_p2.rimingName; label_p3.rimingName];
    label.rimingId = [label_p1.rimingId; label_p2.rimingId; label_p3.rimingId];
    label.isMelting = [label_p1.isMelting; label_p2.isMelting; label_p3.isMelting];
    label.idx_rnd = [label_p1.idx_rnd label_p2.idx_rnd label_p3.idx_rnd];
end

%% main part
clear all; close all;
dir_labels = '/home/praz/Documents/MASC/masclab/labelling_GUI/labels_v.1.1';
load(fullfile(dir_labels,'CP_v11_full_s.mat'));
l_CP = label;
load(fullfile(dir_labels,'AB_v11_full_s.mat'));
l_AB = label;
load(fullfile(dir_labels,'JG_v11_full_s.mat'));
l_JG = label;
load(fullfile(dir_labels,'DW_v11_full_s.mat'));
l_DW = label;
Ntot = length(l_CP.idx_rnd);

%% raw histograms 
figure(1);
subplot(221);
histogram(l_CP.classId);
axis([-1 12 0 800]);
title('Christophe');
subplot(222);
histogram(l_AB.classId);
axis([-1 12 0 800]);
title('Alexis');
subplot(223);
histogram(l_DW.classId);
axis([-1 12 0 800]);
title('Daniel');
subplot(224);
histogram(l_JG.classId);
axis([-1 12 0 800]);
title('Jacopo');

figure(2);
subplot(221);
histogram(l_CP.rimingId);
axis([-1 6 0 920]);
title('Christophe');
subplot(222);
histogram(l_AB.rimingId);
axis([-1 6 0 920]);
title('Alexis');
subplot(223);
histogram(l_DW.rimingId);
axis([-1 6 0 920]);
title('Daniel');
subplot(224);
histogram(l_JG.rimingId);
axis([-1 6 0 920]);
title('Jacopo');

%% take labels following majority rule

% some pre filtering
% merge columns and needles
l_CP.classId(l_CP.classId==3) = 2;
l_AB.classId(l_AB.classId==3) = 2;
l_DW.classId(l_DW.classId==3) = 2;
l_JG.classId(l_JG.classId==3) = 2;
% merge sectored plate & dendrite
l_CP.classId(l_CP.classId==6) = 5;
l_AB.classId(l_AB.classId==6) = 5;
l_DW.classId(l_DW.classId==6) = 5;
l_JG.classId(l_JG.classId==6) = 5;
% merge capped columns with columns
% l_CP.classId(l_CP.classId==10) = 2;
% l_AB.classId(l_AB.classId==10) = 2;
% l_DW.classId(l_DW.classId==10) = 2;
% l_JG.classId(l_JG.classId==10) = 2;
% merge graupel-like with graupel
% l_CP.classId(l_CP.rimingId==4) = 11;
% l_AB.classId(l_AB.rimingId==4) = 11;
% l_DW.classId(l_DW.rimingId==4) = 11;
% l_JG.classId(l_JG.rimingId==4) = 11;


% snowflake class
all_classId = [l_CP.classId l_AB.classId l_DW.classId l_JG.classId];
[M,F,C] = mode(all_classId,2);
for i=1:length(C)
    n_cand(i,1) = length(C{i});
end
idx_bad = find(n_cand>1);
Nok = Ntot - length(idx_bad);
idx_good = find(F>=3);
Ngood = length(idx_good);

% create new labels according to majority rule
l_final.flakes = l_CP.flakes;
l_final.idx_rnd = 1:1:Ntot;
l_final.classId = M;
l_final.classId(idx_bad) = -2;
l_final.className = cell(Ntot,1);
l_final.className(l_final.classId == -2) = {'ambiguous'};
l_final.className(l_final.classId == -1) = {'not labelled yet.'};
l_final.className(l_final.classId == 0) = {'undetermined'};
l_final.className(l_final.classId == 1) = {'small'};
l_final.className(l_final.classId == 2) = {'column'};
l_final.className(l_final.classId == 3) = {'needle'};
l_final.className(l_final.classId == 4) = {'plate'};
l_final.className(l_final.classId == 5) = {'sec_plate'};
l_final.className(l_final.classId == 6) = {'dendrite'};
l_final.className(l_final.classId == 7) = {'aggregate'};
l_final.className(l_final.classId == 8) = {'combin_col'};
l_final.className(l_final.classId == 9) = {'combin_planar'};
l_final.className(l_final.classId == 10) = {'combin_mixte'};
l_final.className(l_final.classId == 11) = {'graupel'};

%riming degree
all_rimingId = [l_CP.rimingId l_AB.rimingId l_DW.rimingId l_JG.rimingId];
[M,F,C] = mode(all_rimingId,2);
for i=1:length(C)
    n_cand_riming(i,1) = length(C{i});
end
idx_bad_riming = find(n_cand_riming>1);
Nok_riming = Ntot - length(idx_bad_riming);
idx_good_riming = find(F>=3);
Ngood_riming = length(idx_good_riming);

% create new labels according to average rule
all_rimingId(all_rimingId == 0) = NaN;
l_final.rimingId = nanmean(all_rimingId,2);
l_final.rimingId(isnan(l_final.rimingId)) = 0;
l_final.rimingId = round(l_final.rimingId);
l_final.rimingName = cell(Ntot,1);
l_final.rimingName(l_final.rimingId == -1) = {'not labelled yet.'};
l_final.rimingName(l_final.rimingId == 0) = {'undetermined'};
l_final.rimingName(l_final.rimingId == 1) = {'none'};
l_final.rimingName(l_final.rimingId == 2) = {'m_rimed'};
l_final.rimingName(l_final.rimingId == 3) = {'h_rimed'};
l_final.rimingName(l_final.rimingId == 4) = {'graupel_like'};
l_final.rimingName(l_final.rimingId == 5) = {'graupel'};

%melting snow
all_isMelting = [l_CP.isMelting l_AB.isMelting l_DW.isMelting l_JG.isMelting];
[M,F,C] = mode(all_isMelting,2);
for i=1:length(C)
    n_cand_melting(i,1) = length(C{i});
end
idx_bad_melting = find(n_cand_melting>1);
Nok_melting = Ntot - length(idx_bad_melting);
idx_good_melting = find(F>=3);
Ngood_melting = length(idx_good_melting);

% create new labels according to majority rule (2-2 = not melting)
l_final.isMelting = M;

fprintf('Analyzing the coherence between labels... \n');
fprintf('SNOWFLAKE CLASS : \n');
fprintf('%2.1f%% of the dataset has a predominant label (%u/%u) \n',Nok/Ntot*100,Nok,Ntot);
fprintf('%2.1f%% of the dataset has at least 3 identical labels (%u/%u) \n\n',Ngood/Ntot*100,Ngood,Ntot);
fprintf('RIMING DEGREE : \n');
fprintf('%2.1f%% of the dataset has a predominant label (%u/%u) \n',Nok_riming/Ntot*100,Nok_riming,Ntot);
fprintf('%2.1f%% of the dataset has at least 3 identical labels (%u/%u) \n\n',Ngood_riming/Ntot*100,Ngood_riming,Ntot);
fprintf('MELTING SNOW : \n');
fprintf('%2.1f%% of the dataset has a predominant label (%u/%u) \n',Nok_melting/Ntot*100,Nok_melting,Ntot);
fprintf('%2.1f%% of the dataset has at least 3 identical labels (%u/%u) \n\n',Ngood_melting/Ntot*100,Ngood_melting,Ntot);


%% Removing undetermined
idx = find(l_CP.classId == 0 | l_AB.classId ==0 | l_JG.classId == 0 | l_DW.classId == 0);
l_CP.classId(idx) = [];
l_AB.classId(idx) = [];
l_DW.classId(idx) = [];
l_JG.classId(idx) = [];
% merge columns and needles
l_CP.classId(l_CP.classId==3) = 2;
l_AB.classId(l_AB.classId==3) = 2;
l_DW.classId(l_DW.classId==3) = 2;
l_JG.classId(l_JG.classId==3) = 2;
% merge sectored plates and stellar crystals (5-6)
l_CP.classId(l_CP.classId==6) = 5;
l_AB.classId(l_AB.classId==6) = 5;
l_DW.classId(l_DW.classId==6) = 5;
l_JG.classId(l_JG.classId==6) = 5;
% merge comb of planar crystals & aggregates
% l_CP.classId(l_CP.classId==9) = 7;
% l_AB.classId(l_AB.classId==9) = 7;
% l_DW.classId(l_DW.classId==9) = 7;
% l_JG.classId(l_JG.classId==9) = 7;


fprintf('\n SNOWFLAKE CLASS : Removing particles classified at least once as undetermined... %u particles removed (%2.1f%%) \n\n',length(idx),length(idx)/Ntot*100);

idx = find(l_CP.rimingId == 0 | l_AB.rimingId ==0 | l_JG.rimingId == 0 | l_DW.rimingId == 0);
l_CP.rimingId(idx) = [];
l_AB.rimingId(idx) = [];
l_DW.rimingId(idx) = [];
l_JG.rimingId(idx) = [];
fprintf('\n RIMING DEGREE : Removing particles classified at least once as undetermined... %u particles removed (%2.1f%%) \n\n',length(idx),length(idx)/Ntot*100);

% BER confmat
label_ppl = {'Christophe','Alexis','Daniel','Jacopo'};
BER = zeros(4);
BER(1,2) = computeBER(l_CP.classId,l_AB.classId); BER(2,1) = BER(1,2);
BER(1,3) = computeBER(l_CP.classId,l_DW.classId); BER(3,1) = BER(1,3);
BER(1,4) = computeBER(l_CP.classId,l_JG.classId); BER(4,1) = BER(1,4);
BER(2,3) = computeBER(l_AB.classId,l_DW.classId); BER(3,2) = BER(2,3);
BER(2,4) = computeBER(l_AB.classId,l_JG.classId); BER(4,2) = BER(2,4);
BER(3,4) = computeBER(l_DW.classId,l_JG.classId); BER(4,3) = BER(3,4);
figure(3);
heatmap(BER,label_ppl,label_ppl,1,'Colormap','red','ShowAllTicks',1,'UseLogColorMap',false,'Colorbar',true,'MinColorValue',0,'MaxColorValue',0.5);
title('Balanced Error Rate between different sets of labels');

% HSS confmat
HSS = ones(4);
HSS(1,2) = computeKAPPA(l_CP.classId,l_AB.classId); HSS(2,1) = HSS(1,2);
HSS(1,3) = computeKAPPA(l_CP.classId,l_DW.classId); HSS(3,1) = HSS(1,3);
HSS(1,4) = computeKAPPA(l_CP.classId,l_JG.classId); HSS(4,1) = HSS(1,4);
HSS(2,3) = computeKAPPA(l_AB.classId,l_DW.classId); HSS(3,2) = HSS(2,3);
HSS(2,4) = computeKAPPA(l_AB.classId,l_JG.classId); HSS(4,2) = HSS(2,4);
HSS(3,4) = computeKAPPA(l_DW.classId,l_JG.classId); HSS(4,3) = HSS(3,4);
figure(4);
heatmap(round(HSS,2),label_ppl,label_ppl,1,'Colormap','red','ShowAllTicks',1,'UseLogColorMap',false,'Colorbar',true,'MinColorValue',0.7,'MaxColorValue',1);
title('Heide Skill Score between different sets of labels');
set(gca,'Fontsize',12);

% OA confmat
OA = ones(4);
OA(1,2) = computeOA(l_CP.classId,l_AB.classId); OA(2,1) = OA(1,2);
OA(1,3) = computeOA(l_CP.classId,l_DW.classId); OA(3,1) = OA(1,3);
OA(1,4) = computeOA(l_CP.classId,l_JG.classId); OA(4,1) = OA(1,4);
OA(2,3) = computeOA(l_AB.classId,l_DW.classId); OA(3,2) = OA(2,3);
OA(2,4) = computeOA(l_AB.classId,l_JG.classId); OA(4,2) = OA(2,4);
OA(3,4) = computeOA(l_DW.classId,l_JG.classId); OA(4,3) = OA(3,4);
figure(5);
heatmap(OA,label_ppl,label_ppl,1,'Colormap','red','ShowAllTicks',1,'UseLogColorMap',false,'Colorbar',true,'MinColorValue',0.5,'MaxColorValue',1);
title('Overall Accuracy between different sets of labels');


if false
    % Snwoflake Class 1-1 comparison
    fprintf('\n\n SNOWFLAKE CLASS \n\n');
    fprintf('***** Christophe - Alexis comparison ***** \n');
    ber = computeBER(l_CP.classId,l_AB.classId);
    OA = computeOA(l_CP.classId,l_AB.classId);
    HSS = computeKAPPA(l_CP.classId,l_AB.classId);
    fprintf('   BER : %2.2f%% \n',ber*100);
    fprintf('   HSS : %2.2f%% \n',HSS*100);
    fprintf('   OA  : %2.2f%% \n',OA*100);
    fprintf('****************************************** \n\n');

    fprintf('***** Christophe - Daniel comparison ***** \n');
    ber = computeBER(l_CP.classId,l_DW.classId);
    OA = computeOA(l_CP.classId,l_DW.classId);
    HSS = computeKAPPA(l_CP.classId,l_DW.classId);
    fprintf('   BER : %2.2f%% \n',ber*100);
    fprintf('   HSS : %2.2f%% \n',HSS*100);
    fprintf('   OA  : %2.2f%% \n',OA*100);
    fprintf('****************************************** \n\n');

    fprintf('***** Christophe - Jacopo comparison ***** \n');
    ber = computeBER(l_CP.classId,l_JG.classId);
    OA = computeOA(l_CP.classId,l_JG.classId);
    HSS = computeKAPPA(l_CP.classId,l_JG.classId);
    fprintf('   BER : %2.2f%% \n',ber*100);
    fprintf('   HSS : %2.2f%% \n',HSS*100);
    fprintf('   OA  : %2.2f%% \n',OA*100);
    fprintf('****************************************** \n\n');

    fprintf('***** Alexis - Daniel comparison ***** \n');
    ber = computeBER(l_AB.classId,l_DW.classId);
    OA = computeOA(l_AB.classId,l_DW.classId);
    HSS = computeKAPPA(l_AB.classId,l_DW.classId);
    fprintf('   BER : %2.2f%% \n',ber*100);
    fprintf('   HSS : %2.2f%% \n',HSS*100);
    fprintf('   OA  : %2.2f%% \n',OA*100);
    fprintf('****************************************** \n\n');

    fprintf('***** Alexis - Jacopo comparison ***** \n');
    ber = computeBER(l_AB.classId,l_JG.classId);
    OA = computeOA(l_AB.classId,l_JG.classId);
    HSS = computeKAPPA(l_AB.classId,l_JG.classId);
    fprintf('   BER : %2.2f%% \n',ber*100);
    fprintf('   HSS : %2.2f%% \n',HSS*100);
    fprintf('   OA  : %2.2f%% \n',OA*100);
    fprintf('****************************************** \n\n');

    fprintf('***** Daniel - Jacopo comparison ***** \n');
    ber = computeBER(l_DW.classId,l_JG.classId);
    OA = computeOA(l_DW.classId,l_JG.classId);
    HSS = computeKAPPA(l_DW.classId,l_JG.classId);
    fprintf('   BER : %2.2f%% \n',ber*100);
    fprintf('   HSS : %2.2f%% \n',HSS*100);
    fprintf('   OA  : %2.2f%% \n',OA*100);
    fprintf('****************************************** \n\n');

    % Riming Degree 1-1 comparison
    fprintf('\n\n RIMING DEGREE \n\n');
    fprintf('***** Christophe - Alexis comparison ***** \n');
    ber = computeBER(l_CP.rimingId,l_AB.rimingId);
    OA = computeOA(l_CP.rimingId,l_AB.rimingId);
    HSS = computeKAPPA(l_CP.rimingId,l_AB.rimingId);
    fprintf('   BER : %2.2f%% \n',ber*100);
    fprintf('   HSS : %2.2f%% \n',HSS*100);
    fprintf('   OA  : %2.2f%% \n',OA*100);
    fprintf('****************************************** \n\n');

    fprintf('***** Christophe - Daniel comparison ***** \n');
    ber = computeBER(l_CP.rimingId,l_DW.rimingId);
    OA = computeOA(l_CP.rimingId,l_DW.rimingId);
    HSS = computeKAPPA(l_CP.rimingId,l_DW.rimingId);
    fprintf('   BER : %2.2f%% \n',ber*100);
    fprintf('   HSS : %2.2f%% \n',HSS*100);
    fprintf('   OA  : %2.2f%% \n',OA*100);
    fprintf('****************************************** \n\n');

    fprintf('***** Christophe - Jacopo comparison ***** \n');
    ber = computeBER(l_CP.rimingId,l_JG.rimingId);
    OA = computeOA(l_CP.rimingId,l_JG.rimingId);
    HSS = computeKAPPA(l_CP.rimingId,l_JG.rimingId);
    fprintf('   BER : %2.2f%% \n',ber*100);
    fprintf('   HSS : %2.2f%% \n',HSS*100);
    fprintf('   OA  : %2.2f%% \n',OA*100);
    fprintf('****************************************** \n\n');

    fprintf('***** Alexis - Daniel comparison ***** \n');
    ber = computeBER(l_AB.rimingId,l_DW.rimingId);
    OA = computeOA(l_AB.rimingId,l_DW.rimingId);
    HSS = computeKAPPA(l_AB.rimingId,l_DW.rimingId);
    fprintf('   BER : %2.2f%% \n',ber*100);
    fprintf('   HSS : %2.2f%% \n',HSS*100);
    fprintf('   OA  : %2.2f%% \n',OA*100);
    fprintf('****************************************** \n\n');

    fprintf('***** Alexis - Jacopo comparison ***** \n');
    ber = computeBER(l_AB.rimingId,l_JG.rimingId);
    OA = computeOA(l_AB.rimingId,l_JG.rimingId);
    HSS = computeKAPPA(l_AB.rimingId,l_JG.rimingId);
    fprintf('   BER : %2.2f%% \n',ber*100);
    fprintf('   HSS : %2.2f%% \n',HSS*100);
    fprintf('   OA  : %2.2f%% \n',OA*100);
    fprintf('****************************************** \n\n');

    fprintf('***** Daniel - Jacopo comparison ***** \n');
    ber = computeBER(l_DW.rimingId,l_JG.rimingId);
    OA = computeOA(l_DW.rimingId,l_JG.rimingId);
    HSS = computeKAPPA(l_DW.rimingId,l_JG.rimingId);
    fprintf('   BER : %2.2f%% \n',ber*100);
    fprintf('   HSS : %2.2f%% \n',HSS*100);
    fprintf('   OA  : %2.2f%% \n',OA*100);
    fprintf('****************************************** \n\n');

    % Melting Snow 1-1 comparison
    fprintf('\n\n MELTING SNOW \n\n');
    fprintf('***** Christophe - Alexis comparison ***** \n');
    ber = computeBER(l_CP.isMelting,l_AB.isMelting);
    OA = computeOA(l_CP.isMelting,l_AB.isMelting);
    HSS = computeKAPPA(l_CP.isMelting,l_AB.isMelting);
    fprintf('   BER : %2.2f%% \n',ber*100);
    fprintf('   HSS : %2.2f%% \n',HSS*100);
    fprintf('   OA  : %2.2f%% \n',OA*100);
    fprintf('****************************************** \n\n');

    fprintf('***** Christophe - Daniel comparison ***** \n');
    ber = computeBER(l_CP.isMelting,l_DW.isMelting);
    OA = computeOA(l_CP.isMelting,l_DW.isMelting);
    HSS = computeKAPPA(l_CP.isMelting,l_DW.isMelting);
    fprintf('   BER : %2.2f%% \n',ber*100);
    fprintf('   HSS : %2.2f%% \n',HSS*100);
    fprintf('   OA  : %2.2f%% \n',OA*100);
    fprintf('****************************************** \n\n');

    fprintf('***** Christophe - Jacopo comparison ***** \n');
    ber = computeBER(l_CP.isMelting,l_JG.isMelting);
    OA = computeOA(l_CP.isMelting,l_JG.isMelting);
    HSS = computeKAPPA(l_CP.isMelting,l_JG.isMelting);
    fprintf('   BER : %2.2f%% \n',ber*100);
    fprintf('   HSS : %2.2f%% \n',HSS*100);
    fprintf('   OA  : %2.2f%% \n',OA*100);
    fprintf('****************************************** \n\n');

    fprintf('***** Alexis - Daniel comparison ***** \n');
    ber = computeBER(l_AB.isMelting,l_DW.isMelting);
    OA = computeOA(l_AB.isMelting,l_DW.isMelting);
    HSS = computeKAPPA(l_AB.isMelting,l_DW.isMelting);
    fprintf('   BER : %2.2f%% \n',ber*100);
    fprintf('   HSS : %2.2f%% \n',HSS*100);
    fprintf('   OA  : %2.2f%% \n',OA*100);
    fprintf('****************************************** \n\n');

    fprintf('***** Alexis - Jacopo comparison ***** \n');
    ber = computeBER(l_AB.isMelting,l_JG.isMelting);
    OA = computeOA(l_AB.isMelting,l_JG.isMelting);
    HSS = computeKAPPA(l_AB.isMelting,l_JG.isMelting);
    fprintf('   BER : %2.2f%% \n',ber*100);
    fprintf('   HSS : %2.2f%% \n',HSS*100);
    fprintf('   OA  : %2.2f%% \n',OA*100);
    fprintf('****************************************** \n\n');

    fprintf('***** Daniel - Jacopo comparison ***** \n');
    ber = computeBER(l_DW.isMelting,l_JG.isMelting);
    OA = computeOA(l_DW.isMelting,l_JG.isMelting);
    HSS = computeKAPPA(l_DW.isMelting,l_JG.isMelting);
    fprintf('   BER : %2.2f%% \n',ber*100);
    fprintf('   HSS : %2.2f%% \n',HSS*100);
    fprintf('   OA  : %2.2f%% \n',OA*100);
    fprintf('****************************************** \n\n');
end



% % facultative: remove undetermined
% %idx = find(l_CP.id==0 | l_JG.id==0 | l_DW.id==0);
% %l_CP.id(idx) = [];
% %l_JG.id(idx) = [];
% %l_DW.id(idx) = [];
% % new labels
% confmat_labels = {'other','grau','agg','melt','small','dend','col'};

% % compute OA, kappa, BER
% cpjg_ber = computeBER(l_CP.id,l_JG.id);
% cpjg_OA  = computeOA(l_CP.id,l_JG.id);
% cpjg_kap = computeKAPPA(l_CP.id,l_JG.id);
% 
% cpdw_ber = computeBER(l_CP.id,l_DW.id);
% cpdw_OA  = computeOA(l_CP.id,l_DW.id);
% cpdw_kap = computeKAPPA(l_CP.id,l_DW.id);
% 
% jgdw_ber = computeBER(l_JG.id,l_DW.id);
% jgdw_OA  = computeOA(l_JG.id,l_DW.id);
% jgdw_kap = computeKAPPA(l_JG.id,l_DW.id);
% 
% 
% 
% % CP - JG
% cpjg = confusionmat(l_CP.id,l_JG.id);
% figure;
% heatmap(cpjg,confmat_labels,confmat_labels,1,'Colormap','red','ShowAllTicks',1,'UseLogColorMap',true,'Colorbar',true);
% xlabel('Jacopo');
% ylabel('Christophe');
% title(sprintf('OA : %2.1f%%   BER : %2.1f%%   HSS : %2.1f%%',cpjg_OA*100,cpjg_ber*100,cpjg_kap*100));
% 
% % CP - DW
% cpdw = confusionmat(l_CP.id,l_DW.id);
% figure;
% heatmap(cpdw,confmat_labels,confmat_labels,1,'Colormap','red','ShowAllTicks',1,'UseLogColorMap',true,'Colorbar',true);
% xlabel('Daniel');
% ylabel('Christophe');
% title(sprintf('OA : %2.1f%%   BER : %2.1f%%   HSS : %2.1f%%',cpdw_OA*100,cpdw_ber*100,cpdw_kap*100));
% 
% % JG - DW
% jgdw = confusionmat(l_JG.id,l_DW.id);
% figure;
% heatmap(jgdw,confmat_labels,confmat_labels,1,'Colormap','red','ShowAllTicks',1,'UseLogColorMap',true,'Colorbar',true);
% xlabel('Daniel');
% ylabel('Jacopo');
% title(sprintf('OA : %2.1f%%   BER : %2.1f%%   HSS : %2.1f%%',jgdw_OA*100,jgdw_ber*100,jgdw_kap*100));








