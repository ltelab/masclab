function plot_area_curve(x,curve,x_axis_label,y_axis_label,c_dark,c_med,c_light)
%plot nice filled curved to analyze variance reduction for several
%iterations

% input:

% X: M x 1 vector
% error_train: M x Nit matrix of train errors
% error_test: M x Nit matrix of test errors

x = x';

% function to plot filled curves
fill_between_lines = @(X,Y1,Y2,C) fill( [X fliplr(X)],  [Y1 fliplr(Y2)], C, 'edgecolor','none','facealpha',.5);
hold on; box on; grid on;

% color choice
if nargin < 5
    c_dark = [153,0,0]/255;
    c_med = [239,59,44]/255;
    c_light = [251,106,74]/255;
end


% plot first curve in shades of red
q_first = quantile(curve,[0.1,0.25,0.5,0.75,0.9],2);
if size(curve,2) > 1
    fill_between_lines(x,q_first(:,1)',q_first(:,5)',c_light);
    fill_between_lines(x,q_first(:,2)',q_first(:,4)',c_med);
end
p1=plot(x,q_first(:,3)','Color',c_dark,'Linewidth',2);


xlabel(x_axis_label);
ylabel(y_axis_label);
%set(gca,'Fontsize',14);


end

%set(h,'facealpha',.5)

% shades of blue
% 239,243,255
% 198,219,239
% 158,202,225
% 107,174,214
% 66,146,198
% 33,113,181
% 8,69,148

% shades of red
% 254,229,217
% 252,187,161
% 252,146,114
% 251,106,74
% 239,59,44
% 203,24,29
% 153,0,13