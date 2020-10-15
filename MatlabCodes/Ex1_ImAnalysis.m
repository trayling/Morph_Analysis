function combinedata=Ex1_ImAnalysis(options_thick,options_thin,options_nuc,options_pore,colorthick,colorthin,colornuc,res,savesize,showfig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This function ImAnalysis processes images with 3 channels consecutively 
%  binarizes grey scale images, smooths them and calculates network properties
% 
% combinedata=Ex1_ImAnalysis(options_thick,options_thin,options_nuc,options_pore,colorthick,colorthin,colornuc,res,savesize,showfig);
% 
% inputs,
%       refer to Ex1_main_ImAnalysis.m
% outputs,
%    combinedata: measured data from all 3 channels for all image
%    slices

% Function is written by YikTungTracy Ling, Johns Hopkins University (July 2019)
% Reference: Ling, Y. T. T., Pease, M. E., Jefferys, J. L., Kimball, E. C., Quigley, H. A., 
% & Nguyen, T. D. (2020). Pressure-Induced Changes in Astrocyte GFAP, Actin, and Nuclear 
% Morphology in Mouse Optic Nerve. Investigative Ophthalmology & Visual Science, 61(11), 14-14.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all 
% close all
warning('off','all')
%%
% get current matlab code folder name
codefolderloc = pwd;

% get name of the current directory (main folder)
strloc=strfind(codefolderloc, 'MatlabCodes');
folderloc = codefolderloc(1:strloc-1);

% define filenames for new folders
namd_thickIbound=fullfile(folderloc,'Ex1_OutputPics','thickbound');
namd_thickIden=fullfile(folderloc,'Ex1_OutputPics','thickden');
namd_thinden=fullfile(folderloc,'Ex1_OutputPics','thinden');
namd_nucden=fullfile(folderloc,'Ex1_OutputPics','Nuclden');

% make new folder
mkdir(namd_thickIbound);
mkdir(namd_thickIden);
mkdir(namd_thinden);
mkdir(namd_nucden);

% folder where each section is split into 3 channels Red = actin, Green =
% GFAP, cell = Nucleus
splitfolder = fullfile(folderloc,'Ex1_InputSplit3channel');
d = dir(fullfile(splitfolder, '*.tif'));
numfile=size(d,1);
stacks=numfile/3;


% preallocate space for measured features
slidenum=NaN(stacks,1);                         % slide number
CombineAREA=NaN(stacks,1);                      % combine area thick thin channel (pixel^2)
Combineden=NaN(stacks,1);                       % combine area/mask_area
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

% features from thin network channel
ThinStat(1:stacks)=struct('thinarea',NaN,...    % total area of thin network in ROI (pixels^2)
    'thinarea_um',NaN,...                       % total area of thin network in ROI (um^2)
    'thindensity',NaN,...                       % area fration (thinarea/mask_area)
    'thinkfit',NaN,...                          % network anisotropy (fitted k from von-mises distribution)
    'thinori',NaN);                             % circular mean orientation of the network

% features from nuclear channel
NucStat(1:stacks)=struct('Nucarea_total',NaN,...% total area of nuclear signals in ROI (pixels^2)  
    'Nucarea_med',NaN,...                       % median area of the clusters
    'Nucarea_mean',NaN,...                      % meam area of the clusters  
    'Nucarea_std',NaN,...                       % standard deviation of the areas of the clusters
    'NucAR_med',NaN,...                         % median aspect ratio of the clusters
    'NucAR_mean',NaN,...                        % mean aspect ratioof the clusters                        
    'NucAR_std',NaN,...                         % standard deviation of the aspect ratios of the clusters
    'nuccount',NaN,...                          % number of clusters over area of ROI   
    'Nucden',NaN,...                            % # of cluster/mask_area    
    'NucAREAden',NaN);                          % area fraction (Nucarea_total/mask_area)        

% features from negative thin network and negative nuclear channel (pores)
PStat(1:stacks)=struct('porecount',NaN,...      % number of individual pore compartments
    'poreareatotal',NaN,...                     % total pore area
    'parea_med',NaN,...                       % % median area of pore compartments
    'parea_mean',NaN,...                      % mean area of pore compartments
    'parea_std',NaN);                         % standard deviations of areas of pore compartments


for i=1:stacks
    f = waitbar(i/stacks,['Analyzing: estimated remaining time:' num2str((stacks-i+1)*3) 'min']);
    % channel c=0 thick c=1 blue c=2 red 
    thickidx=(i-1)*3+1;
    nucidx=(i-1)*3+2;
    thinidx=(i-1)*3+3;
    %read image
    ThickIm=imread(fullfile(splitfolder, d(thickidx).name)); 
    NuclIm=imread(fullfile(splitfolder, d(nucidx).name));
    ThinIm=imread(fullfile(splitfolder, d(thinidx).name));
    
    filename = getfield(d(thickidx),'name');
    idx=strfind(filename,'slide');
    savefname=filename(1:idx+7);
    slidenum(i,1)=str2double(filename(idx+5:idx+6));
    %Allfname{i,1}=savefname;
    
    %display loaded figures
    showrgb(ThinIm,ThickIm,NuclIm,showfig);
    showrgb(ThinIm(3198:3930,2356:3090),ThickIm(3198:3930,2356:3090),zeros(size(NuclIm(3198:3930,2356:3090))),showfig);
    showsingle(ThickIm,'g',showfig);
    showsingle(NuclIm,'b',showfig);
    showsingle(ThinIm,'r',showfig);
    %close all
    if showfig ==1
        msgfig = msgbox('Verify Correct Image inputs','Step 1','modal');
        uiwait(msgfig)
        disp('Start Image processing.');
        close all
    end
    
    %% processing thick network channel first
    % ROI mask will be determined from GFAP channel as well
    [thickImask1,thickIden,histthickI,outline,thickstat]=process_thick(ThickIm,showfig,res,options_thick,colorthick);
    ThickStat(i)=thickstat;
    
    %% process thin network channel 
    [thinoutline,thin1,histthin,thinstat]=process_thin(ThinIm,thickImask1,showfig,res,options_thin,colorthin);
    ThinStat(i)=thinstat;
    %% combine thick and thin network channel
    combinethinthick=thin1 | thickIden;
    CombineAREA(i,1)=sum(sum(combinethinthick));
    Combineden(i,1)=sum(sum(combinethinthick))/thickstat.mask_area;
    %%  process Nucleus channel
    [nucoutline,nuclayer,histnuc,~,nucstat]=process_nuc(NuclIm,thickImask1,showfig,res,options_nuc,colornuc);
    NucStat(i)=nucstat;    
    %% process pores from negative thin network and negative nuclear channel
    [pstat]=process_pore(thickImask1,thin1,nuclayer,showfig,options_pore); 
    PStat(i)=pstat;
    delete(f)
    %% save images
    % image outlining the boundary of ROI on thick network channel
    thickIden = imresize(thickIden,[savesize savesize]);
    imwrite(thickIden, fullfile(namd_thickIden,[savefname num2str(i) '.tif']));
    
    % image from thick network channel
    outline = imresize(outline,[savesize savesize]);
    imwrite(outline, fullfile(namd_thickIbound, [savefname num2str(i) '.tif']));
    
    % images from thin network channel
    thinoutline = imresize(thinoutline,[savesize savesize]);
    imwrite(thinoutline, fullfile(namd_thinden, [savefname num2str(i) '.tif']));
    
    % image from nuclear channel
    nucoutline = imresize(nucoutline,[savesize savesize]);
    imwrite(nucoutline, fullfile(namd_nucden, [savefname num2str(i) '.tif']));
    

    showrgb(histthin,histthickI,histnuc,showfig);
    
end

%% combine data
combinedata=[struct2table(ThickStat),struct2table(ThinStat),struct2table(NucStat),struct2table(PStat)];
end