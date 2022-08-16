% script to enhance all pictures (CLAHE) found in INDIR and save results
% in OUTDIR

clear all; close all;

wannasave = true;

indir = '/home/praz/Desktop/snowflakes_selection';
outdir = '/home/praz/Desktop/snowflakes_selection';

list_pic_name = dir(fullfile(indir,'*.png'));
list_pic_name = {list_pic_name.name}';

for i=1:length(list_pic_name)
    datain = imread(fullfile(indir,list_pic_name{i}));
    dataout = brightening(datain);
    dataname = list_pic_name{i}; 
    dataname = [dataname(1:end-4),'_CLAHE_enhanced.png'];
    if wannasave
        imwrite(dataout,fullfile(outdir,dataname),'png','BitDepth', 8);
    end
end

