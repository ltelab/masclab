% analyze_optimization_output
clear all; close all;

% generating the figures for the manuscript
figw = 750;
figh = 500;

% hydrometeor type
% load('out_geo_10it_rand_FS_CV4_HSS_D80.mat');
% n=50;
% xx = 1:1:n;
% myfig = figure('Visible','on','units','pixels','position',[800 600 figw figh]);
% hold on; grid on; box on;
% plot(xx,mean_kappa_Te,'linewidth',1.25);
% axis([0 50 0.5 1]);
% addaxislabel(1,'HSS');
% addaxis(xx,mean_BER_Te,[0 0.5],'linewidth',1.25);
% addaxislabel(2,'BER');
% legend('HSS','BER');
% plot3([25 25],[0.5 1],[-0.1 -0.1],'k--');
% xlabel('# of features considered');
% box on;
% set(gca,'Fontsize',14);
% 
% % riming degree
% load('out_riming_10it_rand_FS_CV4_RMSE_D80.mat');
% n=50;
% xx = 1:1:n;
% myfig = figure('Visible','on','units','pixels','position',[800 400 figw figh]);
% hold on; grid on; box on;
% plot(xx,mean_kappa_Te,'linewidth',1.25);
% axis([0 50 0.5 1]);
% addaxislabel(1,'HSS');
% addaxis(xx,mean_BER_Te,[0 0.5],'linewidth',1.25);
% addaxislabel(2,'BER');
% legend('HSS','BER');
% plot3([25 25],[0.5 1],[-0.1 -0.1],'k--');
% xlabel('# of features considered');
% box on;
% set(gca,'Fontsize',14);
% 
% 
% % melting snow
% load('out_melting_10it_rand_FS_CV4_HSS_80D.mat');
% n=50;
% xx = 1:1:n;
% myfig = figure('Visible','on','units','pixels','position',[800 200 figw figh]);
% hold on; grid on; box on;
% plot(xx,mean_kappa_Te,'linewidth',1.25);
% axis([0 50 0.5 1]);
% addaxislabel(1,'HSS');
% addaxis(xx,mean_BER_Te,[0 0.5],'linewidth',1.25);
% addaxislabel(2,'BER');
% legend('HSS','BER');
% plot3([25 25],[0.5 1],[-0.1 -0.1],'k--');
% xlabel('# of features considered');
% box on;
% set(gca,'Fontsize',14);

% all on the same graph
color = hsv(6);
type_color = [49,130,189]/255;
rim_color = [116,196,118]/255;
melt_color = [252,141,89]./255;
a = load('out_geo_10it_rand_FS_CV4_HSS_D80.mat');
b = load('out_riming_10it_rand_FS_CV4_RMSE_D80.mat');
c = load('out_melting_10it_rand_FS_CV4_HSS_80D.mat');
n=50;
xx = 1:1:n;
myfig = figure('Visible','on','units','pixels','position',[800 600 figw figh]);
hold on; grid on; box on;
plot([25 25],[0.5 1],'k--','color',[150,150,150]/255);
p1 = plot(xx,a.mean_kappa_Te,'color',type_color,'linewidth',1.5);
addaxislabel(1,'HSS');
addaxis(xx,a.mean_BER_Te,[0 0.5],'--','color',type_color,'linewidth',1.5);
addaxislabel(2,'BER');
addaxis(xx,b.mean_BER_Te,[0 0.5],'--','color',rim_color,'linewidth',1.5);
addaxis(xx,c.mean_BER_Te,[0 0.5],'--','color',melt_color,'linewidth',1.5);
p2 = plot(xx,b.mean_kappa_Te,'color',rim_color,'linewidth',1.5);
p3 = plot(xx,c.mean_kappa_Te,'color',melt_color,'linewidth',1.5);
axis([0 50 0.5 1]);
legend([p1 p2 p3],{'type','riming','melting'});
box on;
set(gca,'Fontsize',16);



if 0

% first script to see BER/kappa evolution
if 0
    load('optimization_output_5_1it_notmerged_fmeasplus_noNEW.mat');
    nfeat_vec = dim_ini:-1:dim_ini-length(mean_BER_Te)+1;
    figure;
    subplot(2,1,1); hold on; box on; grid on;
    plot(nfeat_vec,mean_BER_Te*100,'k.-');
    ylabel('BER [%]');
    subplot(2,1,2); hold on; box on; grid on;
    plot(nfeat_vec,mean_kappa_Te*100,'k.-');
    xlabel('# features');
    ylabel('HSS [%]');
end

% second script
if 1
    % forward selection greedy
    a = load('w_costWeights/out_geo_1it_rand_FS_CV10_OA_weird.mat');%load('opt_out_1_0ini_1it_norand_1lambda.mat'); 
    b = load('w_costWeights/out_geo_3it_rand_FS_CV10_OA.mat');%load('opt_out_1_0ini_1it_norand2_1lambda.mat');
    c = load('w_costWeights/out_geo_2it_rand_FS_CV4_OA.mat');%load('opt_out_1_0ini_1it_rand_1lambda.mat'); 
    d = load('w_costWeights/out_geo_10it_rand_FS_CV4_OA.mat');%load('opt_out_1_0ini_1it_rand_1lambda_CV3.mat');
    e = load('w_costWeights/out_geo_20it_rand_FS_CV4_OA.mat');%load('opt_out_1_0ini_1it_norand_1lambda_CV3.mat');
    blues = [221,218,226; ...
            189,201,225; ...
            116,169,207; ...
            43,140,190; ...
            4,90,141];
    blues = blues./255;
    % backward elimination beta
    f = load('BEST_CONF_weights/out_geo_1it_rand1_FS_CV4_HSS_80D.mat');%load('opt_out_1_88ini_1it_norand_1lambda.mat');
    g = load('BEST_CONF_weights/out_geo_1it_seed1_FS_CV4_HSS_80D.mat');%load('opt_out_1_88ini_1it_norand2_1lambda.mat');
    h = load('BEST_CONF_weights/out_geo_10it_rand_FS_CV4_HSS_D80.mat');%load('opt_out_1_88ini_1it_rand_1lambda.mat');
    reds = [252,141,89; ...
            227,74,51; ...
            179,0,0];
    reds = reds./255;
    % filter methods
    i = load('NEW_Weights/out_geo_1it_seed1_FS_CV4_OA.mat');
    j = load('NEW_Weights/out_geo_3it_seed1_FS_CV4_OA.mat');
    %i = load('opt_out_relieff_k100.mat');
    %j = load('opt_out_relieff_k100_anotherseed.mat');
    %i = load('out_geo_1it_norand_FS_CV4_BER.mat');
    %j = load('opt_out_1_0ini_10it_rand_1lambda_CV4.mat');
    greens = [82,82,82; ...
            37,37,37];
    greens = greens./255;
    % backward elimination greedy
    k = load('opt_out_1_88ini_1it_norand_1lambda_greedy.mat');
    purples = [122,1,119]./255;
    % based on OA
    l = load('out_geo_1it_norand_FS_CV4_OA.mat');
    m = load('NEW_Weights/out_geo_1it_seed1_FC_CV4_OA_68D.mat');%load('out_geo_1it_norand_FS_CV4_OA2.mat');
    n = load('NEW_Weights/out_geo_1it_seed2_FC_CV4_OA_68D.mat');%load('out_geo_1it_norand_FS_CV4_OA3.mat');
    o = load('NEW_Weights/out_geo_1it_rand_FC_CV4_OA_68D.mat');%load('out_geo_1it_rand_FS_CV4_OA.mat');
    p = load('NEW_Weights/out_geo_1it_rand2_FC_CV4_OA_68D.mat');%load('out_geo_10it_rand_FS_CV4_OA.mat');
    q = load('NEW_Weights/out_geo_10it_seed1_FS_CV4_OA.mat');
    fgreens = [161,217,155; ...
        116,196,118;
        65,171,93;
        35,139,69;
        0,109,44
        0,68,27];
    fgreens = fgreens./255;
    

    % BER/kappa evolution
    figure;
    hold on; box on; grid on;
    % forward selection greedy
%     plot(a.xaxis4BER,a.mean_BER_Te*100,'k.-','Color',blues(1,:),'linewidth',1);
%     plot(b.xaxis4BER,b.mean_BER_Te*100,'k.-','Color',blues(2,:),'linewidth',1);
%     plot(c.xaxis4BER,c.mean_BER_Te*100,'k.-','Color',blues(3,:),'linewidth',1);
%     plot(d.xaxis4BER,d.mean_BER_Te*100,'k.-','Color',blues(4,:),'linewidth',1);
%     plot(e.xaxis4BER,e.mean_BER_Te*100,'k.-','Color',blues(5,:),'linewidth',1);
    % backward elimination beta
    plot(f.xaxis4BER,f.mean_BER_Te*100,'k.-','Color',reds(1,:),'linewidth',1);
    plot(g.xaxis4BER,g.mean_BER_Te*100,'k.-','Color',reds(2,:),'linewidth',1);
    plot(h.xaxis4BER,h.mean_BER_Te*100,'k.-','Color',reds(3,:),'linewidth',1);
    % relieff
%     plot(i.xaxis4BER,i.mean_BER_Te*100,'k.-','Color',greens(1,:),'linewidth',1);
%     plot(j.xaxis4BER,j.mean_BER_Te*100,'k.-','Color',greens(2,:),'linewidth',1);
    % backward elimination 
%     plot(k.xaxis4BER,k.mean_BER_Te*100,'k.-','Color',purples(1,:),'linewidth',1);
    % FS based on OA
%     plot(l.xaxis4BER,l.mean_BER_Te*100,'k.-','Color',fgreens(1,:),'linewidth',1);
    plot(m.xaxis4BER,m.mean_BER_Te*100,'k.-','Color',fgreens(2,:),'linewidth',1);
    plot(n.xaxis4BER,n.mean_BER_Te*100,'k.-','Color',fgreens(3,:),'linewidth',1);
    plot(o.xaxis4BER,o.mean_BER_Te*100,'k.-','Color',fgreens(4,:),'linewidth',1);
    plot(p.xaxis4BER,p.mean_BER_Te*100,'k.-','Color',fgreens(5,:),'linewidth',1);
%     plot(q.xaxis4BER,q.mean_BER_Te*100,'k.-','Color',fgreens(6,:),'linewidth',1);
    xlabel('# features used');
    ylabel('BER [%]');
    axis([0 50 5 15]);
    legend('greedy FS CV4 - constant split #1','greedy FS CV4 - constant split #2','greedy FS CV4 - random split #1','greedy FS CV4 - random split #2');
%     legend('greedy FS CV4 - constant split #1','greedy FS CV4 - constant split #2','greedy FS CV4 - random split','greedy FS CV3 - random split','greedy FS CV3 - constant split #3',...
%         'BE with weights - constant split #1','BE with weights - constant split #2','BE with weights - random split',...
%         'FS Filter method 1 : Merit + MDL', 'FS Filter method 2 : ReliefF',...
%         'greedy BE CV 4 - constant split #1','FS CV4 OA norand #1','FS CV4 OA norand#2','FS CV4 OA norand #3','FS CV4 OA rand','FS CV4 OA rand 10it');
    set(gca,'Fontsize',14);
    figure; hold on; box on; grid on;
    % forward selection greedy
%     plot(a.xaxis4BER,a.mean_kappa_Te*100,'k.-','Color',blues(1,:),'linewidth',1);
%     plot(b.xaxis4BER,b.mean_kappa_Te*100,'k.-','Color',blues(2,:),'linewidth',1);
%     plot(c.xaxis4BER,c.mean_kappa_Te*100,'k.-','Color',blues(3,:),'linewidth',1);
%     plot(d.xaxis4BER,d.mean_kappa_Te*100,'k.-','Color',blues(4,:),'linewidth',1);
%     plot(e.xaxis4BER,e.mean_kappa_Te*100,'k.-','Color',blues(5,:),'linewidth',1);
    % backward elimination beta
    plot(f.xaxis4BER,f.mean_kappa_Te*100,'k.-','Color',reds(1,:),'linewidth',1);
    plot(g.xaxis4BER,g.mean_kappa_Te*100,'k.-','Color',reds(2,:),'linewidth',1);
    plot(h.xaxis4BER,h.mean_kappa_Te*100,'k.-','Color',reds(3,:),'linewidth',1);
    % relieff
%     plot(i.xaxis4BER,i.mean_kappa_Te*100,'k.-','Color',greens(1,:),'linewidth',1);
%     plot(j.xaxis4BER,j.mean_kappa_Te*100,'k.-','Color',greens(2,:),'linewidth',1);
    % backward elimination 
%     plot(k.xaxis4BER,k.mean_kappa_Te*100,'k.-','Color',purples(1,:),'linewidth',1);
    % FS based on OA
%     plot(l.xaxis4BER,l.mean_kappa_Te*100,'k.-','Color',fgreens(1,:),'linewidth',1);
    plot(m.xaxis4BER,m.mean_kappa_Te*100,'k.-','Color',fgreens(2,:),'linewidth',1);
    plot(n.xaxis4BER,n.mean_kappa_Te*100,'k.-','Color',fgreens(3,:),'linewidth',1);
    plot(o.xaxis4BER,o.mean_kappa_Te*100,'k.-','Color',fgreens(4,:),'linewidth',1);
    plot(p.xaxis4BER,p.mean_kappa_Te*100,'k.-','Color',fgreens(5,:),'linewidth',1);
%     plot(q.xaxis4BER,q.mean_kappa_Te*100,'k.-','Color',fgreens(6,:),'linewidth',1);
    
    xlabel('# features used');
    ylabel('HSS [%]');
    axis([0 50 90 95]);
    set(gca,'Fontsize',14);
    
    nf = 30;
%     idx_a = find(a.xaxis4BER==nf);
%     a_feat_id = a.feat_mat(:,idx_a+1);
%     a_feat_id(a_feat_id==0) = [];
%     idx_b = find(b.xaxis4BER==nf);
%     b_feat_id = b.feat_mat(:,idx_b+1);
%     b_feat_id(b_feat_id==0) = [];
%     idx_c = find(c.xaxis4BER==nf);
%     c_feat_id = c.feat_mat(:,idx_c+1);
%     c_feat_id(c_feat_id==0) = [];
%     idx_d = find(d.xaxis4BER==nf);
%     d_feat_id = d.feat_mat(:,idx_d+1);
%     d_feat_id(d_feat_id==0) = [];
%     idx_e = find(e.xaxis4BER==nf);
%     e_feat_id = e.feat_mat(:,idx_e+1);
%     e_feat_id(e_feat_id==0) = [];
    idx_f = find(f.xaxis4BER==nf); 
    f_feat_id = f.feat_mat(:,idx_f);
    f_feat_id(f_feat_id==0) = [];
    idx_g = find(g.xaxis4BER==nf);
    g_feat_id = g.feat_mat(:,idx_g);
    g_feat_id(g_feat_id==0) = [];
    idx_h = find(h.xaxis4BER==nf);
    h_feat_id = h.feat_mat(:,idx_h);
    h_feat_id(h_feat_id==0) = [];
%     idx_i = find(i.xaxis4BER==nf);
%     i_feat_id = i.feat_mat(:,idx_i);
%     i_feat_id(i_feat_id==0) = [];
%     idx_j = find(j.xaxis4BER==nf);
%     j_feat_id = j.feat_mat(:,idx_j);
%     j_feat_id(j_feat_id==0) = [];
%     idx_k = find(k.xaxis4BER==nf);
%     k_feat_id = k.feat_mat(:,idx_k+1);
%     k_feat_id(k_feat_id==0) = [];
%     idx_l = find(l.xaxis4BER==nf);
%     l_feat_id = l.feat_mat(:,idx_l+1);
%     l_feat_id(l_feat_id==0) = [];
    idx_m = find(m.xaxis4BER==nf);
    m_feat_id = m.feat_mat(:,idx_m+1);
    m_feat_id(m_feat_id==0) = [];
    idx_n = find(n.xaxis4BER==nf);
    n_feat_id = n.feat_mat(:,idx_n+1);
    n_feat_id(n_feat_id==0) = [];
    idx_o = find(o.xaxis4BER==nf);
    o_feat_id = o.feat_mat(:,idx_o+1);
    o_feat_id(o_feat_id==0) = [];
    idx_p = find(p.xaxis4BER==nf);
    p_feat_id = p.feat_mat(:,idx_p+1);
    p_feat_id(p_feat_id==0) = [];
%     idx_q = find(q.xaxis4BER==nf);
%     q_feat_id = q.feat_mat(:,idx_q+1);
%     q_feat_id(q_feat_id==0) = [];
    
    
    % presence of features
    figure; hold on; box on; %grid on;
%     drawRectangles(a_feat_id,0,blues(1,:));
%     drawRectangles(b_feat_id,1,blues(2,:));
%     drawRectangles(c_feat_id,2,blues(3,:));
%     drawRectangles(d_feat_id,3,blues(4,:));
%     drawRectangles(e_feat_id,4,blues(5,:));
    drawRectangles(f_feat_id,9,reds(1,:));
    drawRectangles(g_feat_id,10,reds(2,:));
    drawRectangles(h_feat_id,11,reds(3,:));
%     drawRectangles(i_feat_id,8,greens(1,:));
%     drawRectangles(j_feat_id,9,greens(2,:));
%     drawRectangles(k_feat_id,10,purples(1,:));
%     drawRectangles(l_feat_id,11,fgreens(1,:));
    drawRectangles(m_feat_id,12,fgreens(2,:));
    drawRectangles(n_feat_id,13,fgreens(3,:));
    drawRectangles(o_feat_id,14,fgreens(4,:));
    drawRectangles(p_feat_id,15,fgreens(5,:));
%     drawRectangles(q_feat_id,16,fgreens(6,:));
    axis([0 97 8.8 16.2]);
    set(gca,'Xtick',[0:1:96]);
    set(gca,'Xticklabel',{'',1,'','','',5,'','','','',10,'','','','',15,'','','','',20,'','','','',25,'','','','',30,'','','','',35,'','','','',40,'','','','',45,'','','','',50,'','','','',55,'','','','',60,...
        '','','','',65,'','','','',70,'','','','',75,'','','','',80,'','','','',85,'','','','',90,'','','','',95,''});
    set(gca,'Ytick',[]);
    xlabel('Feature ID (1-96)');
    title('Remaining features, D = 30');
    set(gca,'Fontsize',14);
    grid on;
    
%     xgrid = 0:1:120;
%     c_feat_id_3x = [c_feat_id;c_feat_id;c_feat_id];
%     histogram(c_feat_id_3x,xgrid,'FaceColor','r','FaceAlpha',1);    
%     b_feat_id_2x = [b_feat_id;b_feat_id];
%     histogram(b_feat_id_2x,xgrid,'FaceColor','b','FaceAlpha',1);
%     histogram(a_feat_id,xgrid,'FaceColor','k','FaceAlpha',1);
    
    
    
end

end





