%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code initiates the 2d analyses for all images in a folder, assuming
% this folder has single channel images stored in the same folder
% make sure the the current directory is the folder containing matlab codes
% This code  saves the output images and output data in combinedata in
% Ex2_results.mat

% Function is written by YikTungTracy Ling, Johns Hopkins University (July 2019)
% Reference: Ling, Y. T. T., Pease, M. E., Jefferys, J. L., Kimball, E. C., Quigley, H. A., 
% & Nguyen, T. D. (2020). Pressure-Induced Changes in Astrocyte GFAP, Actin, and Nuclear 
% Morphology in Mouse Optic Nerve. Investigative Ophthalmology & Visual Science, 61(11), 14-14.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%% define input parameters
% define resolution of the image
res=0.1038;         %um/pixel
% display figures?
showfig=1;          % 1= yes 0= no

% state file name and image size for saving
dataname='Ex2_results';
savesize=2048;      % select a suitable size of downsampling for saving in files

% Default ImageProcess Parameters
Lwin=500;           % window size for overall contrast enhancement corrects 
                    % larger scale intensity difference
Swin=30;            % base on space between thick processes, increase 
                    % contrast in local window, smaller number will reveal thinner details but may also increase noise 

% structure size for morphological smoothing, 0 for no smoothing, larger
% for larger degree of smoothing
sthick=3;           % sturcture size for thick channel

%other thresholds
skelthr=10;         % threshold for removing spurious branches from skeleton network (for bwmorph)
bini=8;             % initial estimate of average fiber width in pixels
threshpix=300;      % size for boundary smoothing in automatic boundary detection using thick channel

% define colour for each channel 'r' 'g' or 'b'
colorthick='g';     % green channel

%% assemble options
options=struct('Lwin',Lwin , 'Swin', Swin, 'threshpix', threshpix, 'sthick', sthick, 'skelthr', skelthr, 'bini',bini);

%% get current matlab code folder name
codefolderloc = pwd;

% get name of the current directory (main folder)
strloc=strfind(codefolderloc, 'MatlabCodes');
folderloc = codefolderloc(1:strloc-1);

% define filenames for to create new folders
namd_thickIden=fullfile(folderloc,'Ex2_OutputPics','thickIden');

% make new folder
mkdir(namd_thickIden);

% folder where input images are
splitfolder = fullfile(folderloc,'Ex2_Inputsinglechannel');
d = dir(fullfile(splitfolder, '*.tif'));
%number of images
stacks=size(d,1);

% features from thick network channel
ThickStat(1:stacks)=struct('mask_area',NaN,...  % total area of ROI (pixels^2)
    'mask_AREA_um',NaN,...                      % total area of ROI (um^2)   
    'mask_ar',NaN,...                           % aspect ratio of ROI    
    'thickarea',NaN,...                         % total area of thick signals in ROI (pixels^2)
    'thickIarea_um',NaN,...                     % total area of thick signals in ROI (um^2)
    'thickdensity',NaN,...                      % area fraction (thickarea/mask_area) 
    'thickcenden',NaN,...                       % area fraction in central zone
    'thickperiden',NaN,...                      % area fraction in peripheral+rim zone
    'thickperi90den',NaN,...                    % area fraction in rim zone
    'thickrimden',NaN,...                       % area fraction in peripheral zone
    'bwidth',NaN, ...                           % average beam width 
    'bwistd',NaN, ...                           % standard deviation of beam width
    'bwimax',NaN,...                            % maximum beam width 
    'bwicen',NaN,...                            % average beam width in central zone
    'bwicenstd',NaN,...                         % standard deviation of beam width in central zone
    'bwiperi',NaN,...                           % average beam width in peripheral+rim zone
    'bwiperistd',NaN,...                        % standard deviation of beam width in peripheral +rim zone    
    'bwiperi90',NaN,...                         % average beam width in peripheral zone
    'bwiperi90std',NaN,...                      % standard deviation of beam width in peripheral zone
    'bwirim',NaN,...                            % average beam width in rim zone
    'bwirimstd',NaN,...                         % standard deviation of beam width in rim zone
    'thickkfit',NaN,...                         % network anisotropy (fitted k from von-mises distribution)                            
    'thickori',NaN);                            % circular mean orientation of the network

%% run Image Analysis
for i=1:stacks
    
    % read image
    Im=imread(fullfile(splitfolder, d(i).name)); 
    filename = getfield(d(i),'name');
    idx=strfind(filename,'.tif');
    savefname=filename(1:idx-1);
    [thickImask1,thickIden,histthickI,outline,thickstat]=process_thick(Im,showfig,res,options,colorthick);
    ThickStat(i)=thickstat;
    
    % save image in a smaller size
    thickIden = imresize(thickIden,[savesize savesize]);
    imwrite(thickIden, fullfile(namd_thickIden,[savefname '_out.tif']));
    
end

%% Save output data
combinedata=[struct2table(ThickStat)];
save([dataname '.mat'],'combinedata');
