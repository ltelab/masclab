% some ellipse stats

%% loading data
clear all; close all;

dir_data = '/media/praz/Masc-Data/DATA_meteoswiss/PROCESSED_20151109/DATA/xi_above_11';
dir_img  = '/media/praz/Masc-Data/DATA_meteoswiss/PROCESSED_20151109/IMAGES/xi_above_11';

file_list = dir(fullfile(dir_data,'*.mat'));
file_only_list = {file_list.name};
file_list = fullfile(dir_data,file_only_list);

% create features matrix
for i=1:2%length(file_list)
    
    load(file_list{i});
    
    % retrieve structure of info
    label.name{i} = file_only_list{i};
    label.id = roi.id;
    label.cam = roi.cam;
   
    % dimension related features 
    X(i,1) = roi.area;
    X(i,2) = 0.5 * double(roi.width + roi.height);
    X(i,3) = roi.perim; 
    
    % ellipse related features
    X(i,4) = roi.E_a;
    X(i,5) = roi.E_b;
    
    
    
    load(file_list{i});
    data.name{i} = file_only_list{i};
    data.area(i) = roi.area;
    data.dim(i) = 0.5*double(roi.width + roi.height);
    data.centroid(i,:) = roi.centroid;
    data.mean_intens(i) = roi.mean_intens;
    data.max_intens(i) = roi.max_intens;
    data.lap(i) = roi.lap;
    data.complex(i) = roi.complex;
    data.std(i) = roi.std;
    data.local_std(i) = roi.local_std;
    data.local_std5(i) = roi.local_std5;
    data.local_std7(i) = roi.local_std7;
    
    data.new.lap(i) = roi.new.lap;
    data.new.std(i) = roi.new.std;
    data.new.local_std(i) = roi.new.local_std;
    data.new.local_std5(i) = roi.new.local_std5;
    data.new.local_std7(i) = roi.new.local_std7;
    
    data.T(i) = roi.thygan.T;
    data.RH(i) = roi.thygan.RH;
    data.tnum(i) = roi.tnum;
    
    data.cam(i) = roi.cam;
    data.id(i) = roi.id;
    
    data.E_theta(i) = roi.E.theta;
    data.E_a(i) = roi.E.a;
    if roi.E.b > roi.E.a
        roi.E.b = roi.E.a;
    end
    data.E_b(i) = roi.E.b;
    data.E_in_a(i) = roi.E_in.a;
    if roi.E_in.b > roi.E_in.a
        roi.E_in.b = roi.E_in.a;
    end
    data.E_in_b(i) = roi.E_in.b;
    data.E_out_a(i) = roi.E_out.a;
    if roi.E_out.b > roi.E_out.a
        roi.E_out.b = roi.E_out.a;
    end
    data.E_out_b(i) = roi.E_out.b;
    
%    data.E_ecc(i) = sqrt(1-(data.E_b(i)/data.E_a(i))^2);
%     if data.E_ecc(i) == 0
%         t = linspace(0,2*pi);
%         xt = roi.E.X0 + cos(roi.E.theta).*roi.E.a.*cos(t) - sin(roi.E.theta).*roi.E.b.*sin(t);
%         yt = roi.E.Y0 + sin(roi.E.theta).*roi.E.a.*cos(t) + cos(roi.E.theta).*roi.E.b.*sin(t);
%         %plot(x,y,'b.');
%         figure;
%         imshow(roi.data);
%         hold on;
%         plot(xt,yt,'r-');
%         %legend('data','fit');
%         set(gca,'Ydir','reverse');
%     end
    
    % retrieve all idx (x,y) of pixels within the snowflake
    [data.y{i},data.x{i}] = find(roi.bw_mask_filled ==1);
    data.x_perim{i} = roi.x_perim;
    data.y_perim{i} = roi.y_perim;
    data.perim(i) = length(roi.x_perim);  
    
    % computation of new perimeter, area and complexity
    % good.P2(i) = my_perim(roi.x_perim,roi.y_perim);
    data.A2(i) = data.area(i) - 0.5*length(roi.x_perim);% + 0.5*good.P2(i);
    data.complex2(i) = length(roi.x_perim)/(2*sqrt(pi*data.A2(i)));
    
    % the magic parameter
    data.xhi(i) = log((data.lap(i)+data.new.lap(i))/2 * data.complex2(i) * (data.local_std(i) + data.new.local_std(i))/2 * data.dim(i)); 
    
end

% further computation
data.dim_mm = data.dim * 30 / 1000; %[mm]
data.E_area = pi.*data.E_a.*data.E_b;
data.E_area_in = pi.*data.E_in_a.*data.E_in_b;
data.E_area_out = pi.*data.E_out_a.*data.E_out_b;

data.E_perim = ellipse_perim(data.E_a,data.E_b);
data.E_perim_in = ellipse_perim(data.E_in_a,data.E_in_b);
data.E_perim_out = ellipse_perim(data.E_out_a,data.E_out_b);

data.E_ecc = sqrt(1-(data.E_b./data.E_a).^2);
data.E_ecc_in = sqrt(1-(data.E_in_b./data.E_in_a).^2);
data.E_ecc_out = sqrt(1-(data.E_out_b./data.E_out_a).^2);

data = orderfields(data);


%% dimensions analysis
clearvars -except data; close all;

% AREA
min = 0; max = 12000; nbins = 100; % 100000
bins = linspace(min,max,nbins);
%bins = bins + (max-min)/nbins;
N1 = hist(data.area,bins);
N2 = hist(data.E_area,bins);
N3 = hist(data.E_area_in,bins);
N4 = hist(data.E_area_out,bins);
N5 = hist(data.A2,bins);
figure(1); hold on; box on;
bar(bins,N1,1);
plot(bins,N2,'r-','linewidth',2);
plot(bins,N3,'c-','linewidth',2);
plot(bins,N4,'g-','linewidth',2);
%plot(bins,N5,'c--');
xlabel('area [pixels^2]');
ylabel('occurence');
legend('data Area','Ellipse fit','Inscribed','Circumbscribed');
title('Particle Area Distribution');
v = axis;
axis([min max/2 0 860]);


% PERIMETER
min = 0; max = 1000; nbins = 50; % 4500
bins = linspace(min,max,nbins);
N1 = hist(data.perim,bins);
N2 = hist(data.E_perim,bins);
N3 = hist(data.E_perim_in,bins);
N4 = hist(data.E_perim_out,bins);
figure(2);
 hold on; box on;
bar(bins,N1,1);
plot(bins,N2,'r-','linewidth',2);
plot(bins,N3,'c-','linewidth',2);
plot(bins,N4,'g-','linewidth',2);
v = axis;
axis([min max v(3) v(4)]);
xlabel('perimeter [pixels]');
ylabel('occurence');
legend('data Perim','Ellipse fit','Inscribed','Circumbscribed');
title('Particle Perimeter Distribution');

% ORIENTATION
figure(3);
min = 0; max = 1; nbins = 50;
bins = linspace(min,max,nbins);
N2 = hist(data.E_ecc,bins);
N3 = hist(data.E_ecc_in,bins);
N4 = hist(data.E_ecc_out,bins);
figure(3);
hold on; box on;
plot(bins,N2,'r-','linewidth',2);
plot(bins,N3,'c-','linewidth',2);
plot(bins,N4,'g-','linewidth',2);
v = axis;
axis([min max v(3) v(4)]);
xlabel('eccentricity');
ylabel('occurence');
legend('Ellipse fit','Inscribed','Circumbscribed');
title('Ellipses Eccentricity Distribution');
%histogram(data.E_theta*180/pi);



