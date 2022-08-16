% analyze_optimization_output
clear all; close all;


% forward selection greedy
a = load('out_melting_1it_norand_FS_CV4_lambda0.0001.mat'); 
b = load('out_melting_1it_norand_FS_CV4_lambda0.01.mat');
c = load('out_melting_1it_norand_FS_CV4_lambda1.mat');
d = load('out_melting_1it_norand_FS_CV4_lambda1_noovtest.mat');
% e = load('opt_out_1_0ini_1it_norand_1lambda_CV3.mat');
blues = [221,218,226; ...
        189,201,225; ...
        116,169,207; ...
        43,140,190; ...
        4,90,141];
blues = blues./255;
% % backward elimination beta
% f = load('opt_out_1_88ini_1it_norand_1lambda.mat');
% g = load('opt_out_1_88ini_1it_norand2_1lambda.mat');
% h = load('opt_out_1_88ini_1it_rand_1lambda.mat');
% reds = [252,141,89; ...
%         227,74,51; ...
%         179,0,0];
% reds = reds./255;
% different lambda
%i = load('out_riming_1it_norand_FS_CV4_lambda0.0001.mat');
%j = load('out_riming_1it_norand_FS_CV4_lambda0.01.mat');
% %j = load('opt_out_1_0ini_10it_norand_1lambda_CV3.mat');
% greens = [116,196,118; ...
%             49,163,84];
% greens = greens./255;
% % backward elimination greedy
% k = load('opt_out_1_88ini_1it_norand_1lambda_greedy.mat');
% purples = [122,1,119]./255;
% 

%overview
figure;
subplot(211); grid on; hold on; box on;
plot(a.xaxis4BER,a.mean_BER_Te*100,'k.-','Color',blues(1,:),'linewidth',1);
plot(b.xaxis4BER,b.mean_BER_Te*100,'k.-','Color',blues(2,:),'linewidth',1);
plot(c.xaxis4BER,c.mean_BER_Te*100,'k.-','Color',blues(3,:),'linewidth',1);
plot(d.xaxis4BER,d.mean_BER_Te*100,'k.-','Color',blues(4,:),'linewidth',1);
%plot(i.xaxis4BER,i.mean_BER_Te*100,'k.-','Color',greens(1,:),'linewidth',1);
%plot(j.xaxis4BER,j.mean_BER_Te*100,'k.-','Color',greens(2,:),'linewidth',1);
ylabel('BER [%]');
subplot(212); grid on; hold on; box on;
plot(a.xaxis4BER,a.mean_kappa_Te*100,'k.-','Color',blues(1,:),'linewidth',1);
plot(b.xaxis4BER,b.mean_kappa_Te*100,'k.-','Color',blues(2,:),'linewidth',1);
plot(c.xaxis4BER,c.mean_kappa_Te*100,'k.-','Color',blues(3,:),'linewidth',1);
plot(d.xaxis4BER,d.mean_kappa_Te*100,'k.-','Color',blues(4,:),'linewidth',1);
%plot(i.xaxis4BER,i.mean_kappa_Te*100,'k.-','Color',greens(1,:),'linewidth',1);
%plot(j.xaxis4BER,j.mean_kappa_Te*100,'k.-','Color',greens(2,:),'linewidth',1);
ylabel('HSS [%]');
% subplot(223); grid on; hold on; box on;
% plot(a.xaxis4BER,a.mean_softOA_Te*100,'k.-','Color',blues(1,:),'linewidth',1);
% plot(b.xaxis4BER,b.mean_softOA_Te*100,'k.-','Color',blues(2,:),'linewidth',1);
% plot(c.xaxis4BER,c.mean_softOA_Te*100,'k.-','Color',blues(3,:),'linewidth',1);
%plot(d.xaxis4BER,d.mean_softOA_Te*100,'k.-','Color',blues(4,:),'linewidth',1);
%plot(i.xaxis4BER,i.mean_softOA_Te*100,'k.-','Color',greens(1,:),'linewidth',1);
%plot(j.xaxis4BER,j.mean_softOA_Te*100,'k.-','Color',greens(2,:),'linewidth',1);
% ylabel('soft OA [%]');
% subplot(224); grid on; hold on; box on;
% plot(a.xaxis4BER,a.mean_RMSE_Te,'k.-','Color',blues(1,:),'linewidth',1);
% plot(b.xaxis4BER,b.mean_RMSE_Te,'k.-','Color',blues(2,:),'linewidth',1);
% plot(c.xaxis4BER,c.mean_RMSE_Te,'k.-','Color',blues(3,:),'linewidth',1);
%plot(d.xaxis4BER,d.mean_RMSE_Te,'k.-','Color',blues(4,:),'linewidth',1);
%plot(i.xaxis4BER,i.mean_RMSE_Te,'k.-','Color',greens(1,:),'linewidth',1);
%plot(j.xaxis4BER,j.mean_RMSE_Te,'k.-','Color',greens(2,:),'linewidth',1);
%ylabel('RMSE');




% BER/kappa evolution
%figure;
%hold on; box on; grid on;
% % forward selection greedy
%plot(a.xaxis4BER,a.mean_BER_Te*100,'k.-','Color',blues(1,:),'linewidth',1);
% plot(b.xaxis4BER,b.mean_BER_Te*100,'k.-','Color',blues(2,:),'linewidth',1);
% plot(c.xaxis4BER,c.mean_BER_Te*100,'k.-','Color',blues(3,:),'linewidth',1);
% plot(d.xaxis4BER,d.mean_BER_Te*100,'k.-','Color',blues(4,:),'linewidth',1);
% plot(e.xaxis4BER,e.mean_BER_Te*100,'k.-','Color',blues(5,:),'linewidth',1);
% % backward elimination beta
% plot(f.xaxis4BER,f.mean_BER_Te*100,'k.-','Color',reds(1,:),'linewidth',1);
% plot(g.xaxis4BER,g.mean_BER_Te*100,'k.-','Color',reds(2,:),'linewidth',1);
% plot(h.xaxis4BER,h.mean_BER_Te*100,'k.-','Color',reds(3,:),'linewidth',1);
% % relieff
% plot(i.xaxis4BER,i.mean_BER_Te*100,'k.-','Color',greens(1,:),'linewidth',1);
% plot(j.xaxis4BER,j.mean_BER_Te*100,'k.-','Color',greens(2,:),'linewidth',1);
% % backward elimination 
% plot(k.xaxis4BER,k.mean_BER_Te*100,'k.-','Color',purples(1,:),'linewidth',1);
%xlabel('# features used');
%ylabel('BER [%]');
% axis([0 80 5 15]);
% legend('greedy FS CV4 - constant split #1','greedy FS CV4 - constant split #2','greedy FS CV4 - random split','greedy FS CV3 - random split','greedy FS CV3 - constant split #3',...
%     'BE with weights - constant split #1','BE with weights - constant split #2','BE with weights - random split',...
%     'FS Filter method 1 : Merit + MDL', 'FS Filter method 2 : ReliefF',...
%     'greedy BE CV 4 - constant split #1');
% set(gca,'Fontsize',14);
% figure; hold on; box on; grid on;
% % forward selection greedy
% plot(a.xaxis4BER,a.mean_kappa_Te*100,'k.-','Color',blues(1,:),'linewidth',1);
% plot(b.xaxis4BER,b.mean_kappa_Te*100,'k.-','Color',blues(2,:),'linewidth',1);
% plot(c.xaxis4BER,c.mean_kappa_Te*100,'k.-','Color',blues(3,:),'linewidth',1);
% plot(d.xaxis4BER,d.mean_kappa_Te*100,'k.-','Color',blues(4,:),'linewidth',1);
% plot(e.xaxis4BER,e.mean_kappa_Te*100,'k.-','Color',blues(5,:),'linewidth',1);
% % backward elimination beta
% plot(f.xaxis4BER,f.mean_kappa_Te*100,'k.-','Color',reds(1,:),'linewidth',1);
% plot(g.xaxis4BER,g.mean_kappa_Te*100,'k.-','Color',reds(2,:),'linewidth',1);
% plot(h.xaxis4BER,h.mean_kappa_Te*100,'k.-','Color',reds(3,:),'linewidth',1);
% % relieff
% plot(i.xaxis4BER,i.mean_kappa_Te*100,'k.-','Color',greens(1,:),'linewidth',1);
% plot(j.xaxis4BER,j.mean_kappa_Te*100,'k.-','Color',greens(2,:),'linewidth',1);
% % backward elimination 
% plot(k.xaxis4BER,k.mean_kappa_Te*100,'k.-','Color',purples(1,:),'linewidth',1);
% xlabel('# features used');
% ylabel('HSS [%]');
% axis([0 80 90 95]);
% set(gca,'Fontsize',14);
% 
nf = 30;
idx_a = find(a.xaxis4BER==nf);
a_feat_id = a.feat_mat(:,idx_a+1);
a_feat_id(a_feat_id==0) = [];
idx_b = find(b.xaxis4BER==nf);
b_feat_id = b.feat_mat(:,idx_b+1);
b_feat_id(b_feat_id==0) = [];
idx_c = find(c.xaxis4BER==nf);
c_feat_id = c.feat_mat(:,idx_c+1);
c_feat_id(c_feat_id==0) = [];
idx_d = find(d.xaxis4BER==nf);
d_feat_id = d.feat_mat(:,idx_d+1);
d_feat_id(d_feat_id==0) = [];
% idx_e = find(e.xaxis4BER==nf);
% e_feat_id = e.feat_mat(:,idx_e+1);
% e_feat_id(e_feat_id==0) = [];
% idx_f = find(f.xaxis4BER==nf);
% f_feat_id = f.feat_mat(:,idx_f);
% f_feat_id(f_feat_id==0) = [];
% idx_g = find(g.xaxis4BER==nf);
% g_feat_id = g.feat_mat(:,idx_g);
% g_feat_id(g_feat_id==0) = [];
% idx_h = find(h.xaxis4BER==nf);
% h_feat_id = h.feat_mat(:,idx_h);
% h_feat_id(h_feat_id==0) = [];
% idx_i = find(i.xaxis4BER==nf);
% i_feat_id = i.feat_mat(:,idx_i);
% i_feat_id(i_feat_id==0) = [];
% idx_j = find(j.xaxis4BER==nf);
% j_feat_id = j.feat_mat(:,idx_j);
% j_feat_id(j_feat_id==0) = [];
% idx_k = find(k.xaxis4BER==nf);
% k_feat_id = k.feat_mat(:,idx_k+1);
% k_feat_id(k_feat_id==0) = [];
% 
% 
% % presence of features
figure; hold on; box on; %grid on;
drawRectangles(a_feat_id,0,blues(1,:));
drawRectangles(b_feat_id,1,blues(2,:));
drawRectangles(c_feat_id,2,blues(3,:));
drawRectangles(d_feat_id,3,blues(4,:));
% drawRectangles(e_feat_id,4,blues(5,:));
% drawRectangles(f_feat_id,5,reds(1,:));
% drawRectangles(g_feat_id,6,reds(2,:));
% drawRectangles(h_feat_id,7,reds(3,:));
% drawRectangles(i_feat_id,4,greens(1,:));
% drawRectangles(j_feat_id,5,greens(2,:));
% drawRectangles(k_feat_id,10,purples(1,:));
axis([0 97 -0.2 11.2]);
set(gca,'Xtick',[0:1:96]);
set(gca,'Xticklabel',{'',1,'','','',5,'','','','',10,'','','','',15,'','','','',20,'','','','',25,'','','','',30,'','','','',35,'','','','',40,'','','','',45,'','','','',50,'','','','',55,'','','','',60,...
    '','','','',65,'','','','',70,'','','','',75,'','','','',80,'','','','',85,'','','','',90,'','','','',95,''});
set(gca,'Ytick',[]);
xlabel('Feature ID (1-96)');
title('Remaining features, D = 20');
set(gca,'Fontsize',14);
grid on;








