%% write MASC data in a txt file
clear; close all;

% MASC matlab structure to load
struct_name = '/ltedata/ICEGENESIS_2021/MASC/Proc_data/data_classification_20210115_20210131.mat';
% saving results in ASCII file
filename = '/ltedata/ICEGENESIS_2021/MASC/Proc_data/ICEGENESIS_2021_MASC_classification_20210115_20210131.dat';
% resolution of MASC cameras [mm/pix]
res = 0.033;

fprintf('Load MASC data structure.... ');
load(struct_name);
fprintf(' Done!\n');

fprintf('Write classification results and particle properties in a text file.... ');
Dmax_table = [];
fid = fopen(filename,'w');

% header
fprintf(fid,'# particle-by-particle MASC classification output \n');
fprintf(fid,'# columns header: 1) particle_timestamp_[UTC], 2) particle_ID, 3) camera_ID, 4) particle_type, 5) particle_degree_of_riming, 6) is_melting?, 7) quality_parameter, 8) Dmax_[mm], 9) Area_[mm^2], 10) Aspect_Ratio, 11) Shape_complexity, 12) Orientation, 13) Fallspeed, 14) Mean_brightness_[0-1], 15) Maximum_brightness_[0-1], 16-21) Probabilities_for_hydrometeor_classification \n');
fprintf(fid,'# columns delimiter = tabulation \n');
fprintf(fid,'# particle types : 1 = small particle (SP), 2 = columnar crystal (CC), 3 = planar crystal (PC), 4 = aggregate (AG), 5 = graupel (GR), 6 = combination of columnar and planar crystals (CPC) \n');


for i=1:length(data.Xt)
    fprintf(fid,'%s \t %u \t %u \t %u \t %2.2f \t %u \t',datestr(data.Xt(i),'yyyy.mm.dd HH:MM:SS'),data.X(i,size(data.X,2)),data.X(i,18),data.X(i,1),data.X(i,15),data.X(i,16));
    fprintf(fid,'%4.2f \t %4.2f \t %4.2f \t %4.2f \t %4.2f \t %6.2f \t %4.2f \t %4.2f \t %4.2f \t',data.X(i,6),data.X(i,3)*res,data.X(i,2)*res^2,data.X(i,9),data.X(i,5),data.X(i,10),data.X(i,7),data.X(i,26),data.X(i,27));
    try
        fprintf(fid,'%4.2f \t %4.2f \t %4.2f \t %4.2f \t %4.2f \t %4.2f \t',data.Xfullprob_label(i,1),data.Xfullprob_label(i,2),data.Xfullprob_label(i,3),data.Xfullprob_label(i,4),data.Xfullprob_label(i,5),data.Xfullprob_label(i,6));
    catch
        if i==1
            fprintf('Warning : failed to write fullprob matrix in the file \n');
        end
    end
    
    % loop on the 5 views to save proj Area
%     for j=0:4
%         idx = find(data.X(:,30) == data.Y(i,30) & data.X(:,18) == j);
%         if isempty(idx)
%             fprintf(fid,'%u \t',-1);
%         elseif j<2 % r1 = res of normal cams 234 // r2 = res of add. cams 01
%             fprintf(fid,'%6.2f \t',data.X(idx,2)*(r2^2));
%         else
%             fprintf(fid,'%6.2f \t',data.X(idx,2)*(r1^2));
%         end
%     end
%     
%     % loop on the 5 views to save Dmax
%     for j=0:4
%         idx = find(data.X(:,30) == data.Y(i,30) & data.X(:,18) == j);
%         if isempty(idx)
%             fprintf(fid,'%u \t',-1);
%             Dmax_table(i,j+1) = NaN;
%         elseif j<2 % r1 = res of normal cams 234 // r2 = res of add. cams 01
%             fprintf(fid,'%6.2f \t',data.X(idx,3)*r2);
%             Dmax_table(i,j+1) = data.X(idx,3)*r2;
%         else
%             fprintf(fid,'%6.2f \t',data.X(idx,3)*r1);
%             Dmax_table(i,j+1) = data.X(idx,3)*r1;
%         end
%     end
%     
%     % loop on the 5 views to save AR
%     for j=0:4
%         idx = find(data.X(:,30) == data.Y(i,30) & data.X(:,18) == j);
%         if isempty(idx)
%             fprintf(fid,'%u \t',-1);
%         else
%             fprintf(fid,'%6.2f \t',data.X(idx,9));
%         end
%     end  
%     
%     % loop on the 5 views to save cplx
%     for j=0:4
%         idx = find(data.X(:,30) == data.Y(i,30) & data.X(:,18) == j);
%         if isempty(idx)
%             fprintf(fid,'%u \t',-1);
%         else
%             fprintf(fid,'%6.2f \t',data.X(idx,5));
%         end
%     end  
%     
%     % loop on the 5 views to save orientation
%     for j=0:4
%         idx = find(data.X(:,30) == data.Y(i,30) & data.X(:,18) == j);
%         if isempty(idx)
%             fprintf(fid,'%u \t',-999);
%         else
%             fprintf(fid,'%6.2f \t',data.X(idx,10));
%         end
%     end
%      
%     % loop on the 5 views to save area ratio
%     for j=0:4
%         idx = find(data.X(:,30) == data.Y(i,30) & data.X(:,18) == j);
%         if isempty(idx)
%             fprintf(fid,'%u \t',-1);
%         else
%             fprintf(fid,'%6.2f \t',data.X(idx,29));
%         end
%     end
    
    fprintf(fid,'\n');
    
end
fclose(fid);

fprintf(' Done!\n');

