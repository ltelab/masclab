function plot_learning_curves(x,first_curve,first_label,second_curve,second_label,x_axis_label,y_axis_label)
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

% if given, plot second curve in shades of blue
if nargin > 3
    q_second = quantile(second_curve,[0.1,0.25,0.5,0.75,0.9],2);
    if size(second_curve,2) > 1
        fill_between_lines(x,q_second(:,1)',q_second(:,5)',[107,174,214]/255);
        fill_between_lines(x,q_second(:,2)',q_second(:,4)',[66,146,198]/255);
    end
    p2=plot(x,q_second(:,3)','Color',[8,69,148]/255,'Linewidth',2);
end

% plot first curve in shades of red
q_first = quantile(first_curve,[0.1,0.25,0.5,0.75,0.9],2);
if size(first_curve,2) > 1
    fill_between_lines(x,q_first(:,1)',q_first(:,5)',[251,106,74]/255);
    fill_between_lines(x,q_first(:,2)',q_first(:,4)',[239,59,44]/255);
end
p1=plot(x,q_first(:,3)','Color',[153,0,0]/255,'Linewidth',2);


xlabel(x_axis_label);
ylabel(y_axis_label);
set(gca,'Fontsize',14);
legend([p1 p2],{first_label, second_label});


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