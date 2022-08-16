% script to generate data txt file for Vanderbilt people
clear all; close all;

load('classif_datastruct_Vanderbilt_triplets_provided_withgraupfix_cleaned.mat');

%% blur analysis
if 0

c = hsv(6);
labels = {'small','col','plan','agg','grau','col+plan'};

for i=1:6

    idx = find(data.X(:,1)==i); 
    
    % xhi
    [N_xhi,X_xhi] = hist(data.X(idx,6),[6:0.2:12]);   
    figure(10); 
    subplot(121); hold on; box on; grid on;
    plot(X_xhi,N_xhi,'Color',c(i,:),'linewidth',2);
    subplot(122); hold on; box on; grid on;
    plot(X_xhi,N_xhi./sum(N_xhi),'Color',c(i,:),'linewidth',2);
    %xhi2
    [N_xhi2,X_xhi2] = hist(data.X(idx,28),[0:0.2:50]);
    figure(11);
    subplot(121); hold on; box on; grid on;
    plot(X_xhi2,N_xhi2,'Color',c(i,:),'linewidth',2);
    subplot(122); hold on; box on; grid on;
    plot(X_xhi2,N_xhi2./sum(N_xhi2),'Color',c(i,:),'linewidth',2);
    % scatter xhi - xhi2
    figure(12); hold on; box on; grid on;
    plot(data.X(idx,6),data.X(idx,28),'ko','MarkerfaceColor',c(i,:));
    
end

figure(10);
subplot(121);
xlabel('quality index');
ylabel('count');
set(gca,'Fontsize',14);
subplot(122);
xlabel('quality index');
ylabel('normalized count density');
legend(labels);
set(gca,'Fontsize',14);

figure(11);
subplot(121);
xlabel('alternative quality index');
ylabel('count');
set(gca,'Xlim',[0 10]);
set(gca,'Fontsize',14);
subplot(122);
xlabel('alternative quality index');
ylabel('normalized count density');
legend(labels);
set(gca,'Xlim',[0 10]);
set(gca,'Fontsize',14);

end


%% write data 
if 1
    filename = 'Vanderbilt_MASC_classification_clipped_image_dataset.txt';
    fileID = fopen(filename,'w');
    fprintf(fileID,'# timestamp \t hydrometeor type [1-6] \t riming degree [0-1] \t Dmax [pix] \t mean brightness [0-1] \t quality index \n');
    fprintf(fileID,'# hydrometeor type IDs : 1=small particle, 2=columnar crystal, 3=planar crystal, 4=aggregate, 5=graupel, 6=combination of columnar and planar crystals \n');
    fprintf(fileID,'# riming degree : continuous number between [0,1]. 0=no presence of riming, 1=graupel-type particle \n');
    for i=1:numel(data.Yt)
        fprintf(fileID,'%s \t %u \t %2.2f \t %6.2f \t %2.2f \t %5.2f \n',datestr(data.Yt(i),'yyyy.mm.dd HH:MM:SS'),data.Y(i,1),data.Y(i,15),data.Y(i,3),data.Y(i,26),data.Y(i,6));
    end

    fclose(fileID);
end

