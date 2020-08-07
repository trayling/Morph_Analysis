%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code initiates the 2d analyses for all images in a folder, assuming all 3
% channels of the same image are stored in the same folder
% make sure the the current directory is the folder containing matlab codes
% This code  saves the output images and output data 
% Refer to Ex1_ImAnalysis.m for detailed steps
% Reference: Ling et al. 'Pressure-Induced Changes in Astrocyte GFAP, Actin
% and Nuclear Morphology in Mouse Optic Nerve' IOVS 2020 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%% define input parameters
% define resolution of the image
res=0.1038;         %um/pixel
% display figures?
showfig=0;          % 1= yes 0= no

% state file name and image size for saving
dataname='Ex1_results';
savesize=2048;      % select a suitable size of downsampling for saving in files

% Default ImageProcess Parameters
Lwin=500;           % window size for overall contrast enhancement corrects 
                    % larger scale intensity difference
Swin=30;            % base on space between thick processes, increase 
                    % contrast in local window, smaller number will reveal thinner details but may also increase noise 

% structure size for morphological smoothing, 0 for no smoothing, larger
% for larger degree of smoothing
sthick=3;           % sturcture size for thick channel
sthin=0;            % sturcture size for thin channel
snuc=4;             % sturcture size for nuclear channel
spore=0;            % struture size for pores

%other thresholds
skelthr=10;         % threshold for removing spurious branches from skeleton network (for bwmorph)
bini=8;             % initial estimate of average fiber width in pixels
threshpix=300;      % size for boundary smoothing 

% define colour for each channel 'r' 'g' or 'b'
colorthick='g';     % green channel
colorthin='r';      % red channel
colornuc='b';       % blue channel

%% assemble options
options_thick=struct('Lwin',Lwin , 'Swin', Swin, 'threshpix', threshpix, 'sthick', sthick, 'skelthr', skelthr, 'bini',bini);
options_thin = struct('Lwin',Lwin , 'Swin', Swin, 'sthin',sthin);
options_nuc = struct('Lwin',Lwin , 'snuc',0,'threshpix',threshpix);
options_pore=struct('outlim', 0, 'spore',spore);

%% run Image Analysis
combinedata=Ex1_ImAnalysis(options_thick,options_thin,options_nuc,options_pore,colorthick,colorthin,colornuc,res,savesize,showfig);

%% Save output data
save([dataname '.mat'],'combinedata');
