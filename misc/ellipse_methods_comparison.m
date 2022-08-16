clear all; close all;
% col 1 : least square fit on flake centered using perimeter
% col 2 : least square fit on flake centered using all data
% col 3 : 

load('ellipse_F1_accuracy_comparison.mat');

figure; hold on;
plot(accuracy_block(:,1),'b-');
plot(accuracy_block(:,2),'r-');
plot(accuracy_block(:,3),'k');
plot(accuracy_block(:,4),'y');

mean(accuracy_block(:,1))
mean(accuracy_block(:,2))
mean(accuracy_block(:,3))
mean(accuracy_block(:,4))