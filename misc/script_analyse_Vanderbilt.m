% small script to analyse Vanderbilt data
clear all; close all;

load('prediction/data/newFormat/classif_datastruct_Vanderbilt_triplets_provided_withgraupfix_cleaned.mat');

%% time series

% create a structure masc in the required shape for load_MASC_classif_2
masc.Xt = data.Yt;
masc.xhi = data.Y(:,6);
masc.label_ID = data.Y(:,1);
masc.riming_idx = data.Y(:,15);
masc.melting = data.Y(:,16);
% time serie
TMP1 = load_MASC_classif_2(masc,[2015 07 27 00 00 00],[2015 07 27 13 30 00], 8.5, 10, 10, 0, true); % window size, window shift

TMP2 = load_MASC_classif_2(masc,[2015 08 04 00 00 00],[2015 08 04 09 00 00], 8.5, 10, 10, 0, true);