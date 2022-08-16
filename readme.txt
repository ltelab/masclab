Instructions for use of masclab :
=================================
masclab is a selection of Matlab tools for analyzing data collected by the MASC. Inspired from the initial code written and shared by Fallgatter Technologies, the current version is significantly different from the original one. Among others, masclab allows the user to automatically process snowflake images and classify them according the the following scheme.
1) Hydrometeor type (6 classes) : Small Particle (1), Columnar Crystal (2), Planar Crystal (3), Aggregate (4), Graupel (5), Combination of Columnar and Planar Crystals (6).
2) Riming index (0-1) : Continuous index taking value between 0 (no presence of riming) and 1 (fully rimed graupel particle).
3) Melting snow detection (0,1) : boolean index indicating if the snowflake seems to be melting (1) or not (0).
More information on the classification methodology and accuracy can be found in Praz et al., Atmospheric Measurement Techniques, 2017.

Overview :
==========
masclab is composed by subfolders containing functions and scripts written in Matlab. Users do not need to read and understand all the codes, only a few of them are necessary to process MASC raw images and get the outputs of interest: the classification outcome and some microstructural descriptors (Dmax, area, perim, aspect ratio, orientation, fallspeed, ...). Note that the code is still under construction(!) As a result, some pieces of code might be buggy, obsolete and poorly documented. Please consider this last point when using masclab :)

masclab folder structure :
==========================
3d 	        : attempt to do snowflake 3d reconstruction from MASC triplets (under construction, almost nothing implemented yet).
analysis        : "hotchpotch" of scripts used to analyze processed snowflake images and the outcome of the classification.
classification  : scripts and functions used to implement the machine learning algorithm (Multinomial Logistic Regression) used to train and test the classification. Several of these functions are required to predict the hydrometeor type/riming index/melting snow, like for instance the .mat classification models stored in "masclab/classification/logit_trained_models".
dataio          : functions used to read/write MASC data.
geometry        : functions used to calculate diverse geometrical descriptors based on a processed (cropped) MASC image.
labelling_GUI   : the labelling GUI was used during the classification training phase where several experts labelled independantly more than 3500 MASC images. It can be used to label more images, if desired.
misc            : another hotchpotch containing mostly temporary scripts that might be reused at some point.
obsolete        : obsolete functions and pieces of code. They will be deleted in a future version.
prediction      : important functions used to classify MASC processed images.
processing      : the most important routines and functions are in there (see next section for details).
texture         : functions used to calculate diverse textural descriptors based on a processed (cropped) MASC image.
tools           : utilitary functions
userspecs       : important folder containing the user specifications to run masclab (see next section for details).

How to use it :
===============
This section gives a simple procedure to follow step-by-step in order to initalize masclab, process a set of MASC data, classify the processed images and generate some quicklook figures.

1) Open Matlab and move to the masclab folder.
2) Run mascpaths to load the path to the subfolders into Matlab.
3) Open masclab/processing/MASC_process_classify_quicklooks.m : this routine will allow the user to run sequentially MASC_process.m (to process MASC images), make_predictions_for_campaign.m, merge_predictions_for_campaign.m and make_time_series.m. 
4) Before running the script, make sure that you adjusted all the user parameters (including the ones contained in masclab/userspecs/params_files/process_params_DEFAULT.txt)
5) When displaying images is disabled, the code takes approx. 0.2 seconds to process and save 1 MASC picture so enjoy a coffee break :)

The output of MASC_process is structured as followed : for each MASC image (.png) contained in label_params.campaigndir, a cropped image (.png) and a datafile (.mat) are generated in label_params.output_dir. make_predictions_for_campaign.m will edit each datafile (.mat) by adding fields corresponding to the hydrometeor type, degree of riming and melting snow index as well as the associated probabilities. 

Once you have run the code over a whole data folder (for example corresponding to a measurement campaign), it is convenient to store all the important information (hydrometeor type, riming, Dmax, aspect ratio, fallspeed, etc.) in a matrix to make data treatment and analysis easier. This can be done by running masclab/prediction/merge_predictions_for_campaign.m on a whole processed MASC data folder. Be warned though, this can take several hours if the folder contains several 100'000s of MASC images and the resulting data structure (.mat) can be quite heavy (100s of Mo).

A last note :
=============
masclab was coded using MATLAB R2016b including many different toolboxes (the most important ones are the Image Processing Toolbox, Statistics and Machine Learning Toolbox and Parallel Computing Toolbox). If you are using a different version of Matlab and/or missing some of the mentionned toolboxes, there are chances that the codes will not work correctly on your machine. From my experience (having switched between differend versions of Matlab during this project), consequent errors are usually minor and easy to fix. 


Author : Christophe Praz
Last update : March 2018
