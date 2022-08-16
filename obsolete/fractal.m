% fractal dims exploration

clear all; close all;

path = '/home/praz/Documents/MASC/sample_fractal/';
pic_names = dir(fullfile(path,'*.png'));
pic_names = {pic_names.name}';

for i=1:length(pic_names)

    im = imread(fullfile(path,pic_names{i}));
    bw = false(size(im));
    bw(im>0) = true;
    P = sum(sum(bwperim(bw)));
    A = sum(im(:)>0);
    Fdim_jac = 2*log(P/4)/log(A); % [1,2]
    [nbox,rbox] = boxcount(im);
    nbox = nbox';
    rbox = rbox';
    log_nbox = log(nbox);
    log_rbox = log(rbox);
    
    beta = leastSquares(log_nbox,[ones(length(nbox),1) log_rbox]);
    Fdim = -beta(2);
    
    xx = linspace(min(rbox),max(rbox),2);
    fig1 = figure; 
    subplot(121);
    imshow(im);
    axis equal;
    subplot(122); hold on; box on; grid on;
    plot(rbox,nbox,'ks-');
    plot(xx,exp(beta(1)).*xx.^beta(2),'r--');
    v = axis;
    axis([min(rbox)/2 max(rbox)*2 v(3) v(4)]);
    xlabel('box size');
    ylabel('N boxes');
    title(sprintf('D = %2.2f',Fdim));
    legend('data','LS fit');
    set(gca,'XScale','log','YScale','log');
    %Extract axes handles of all subplots from the figure
    axesHandles = findobj(get(fig1,'Children'), 'flat','Type','axes');
    %Set the axis property to square
    axis(axesHandles,'square');

end



