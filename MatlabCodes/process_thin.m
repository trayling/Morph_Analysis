function [thinIden,thin1,thinstat]=process_thin(thinI,maskI,showfig,res,options,color)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This function process_thin processes images with fine network. It
%  binarizes grey scale images, smooths them and calculates network properties
% 
% [thinIden,thin1,thinstat]=process_thin(thinI,maskI,showfig,res,options,color)
% 
% inputs,
%   thinI :         The original greyscale image with thin network
%   maskI :         The 2D logical/binary image indicating area of interest, 
%                   where pixels in region of interest = 1 and pixels in 
%                   the background =0
%   res:            resolution of the image, um/pixel
%   showfig:        Optional, indicator, 1= show figures when running the 
%                   program, 0= don't show the figures
%   color:          must be 'r','g' or 'b' ,text for indicating whether output
%                   image should be red green or blue in colour
%   options:        Stuct with input options:
%       .Lwin:      Large window size for more global contrast enhancement
%       .Swin:      Small window size for local contrast enhancement
%       .sthin:     radius of disk to be used for morphological opening
%                   larger number will provide larger degree of smoothing
%                   number should be base on the pixels spanning thin beams so that it's not
%                   overlly smoothed
% outputs,
%    thinIden:      2D image with final identified network outlined on the
%                   median filtered image
%    thin1:         Final 2D Binarized image after smoothing  
%                   for calculation of thin network properties
%    thinstat:        struct with measured properties
%     .thinarea       % total area of thin network in ROI (pixels^2)
%     .thinarea_um    % total area of thin network in ROI (um^2)
%     .thindensity    % area fration (thinarea/mask_area)
%     .thinkfit       % network anisotropy (fitted k from von-mises distribution)
%     .thinori        % circular mean orientation of the network
%
% example, 
%   I=imread('samplered.tif');       %  import greyscale image
%   load('samplemask.mat');  % import  ROI mask
%   options_thin = struct('Lwin',500 , 'Swin', 30, 'sthin',0);
%   [thinIden,thin1,thinstat]=process_thin(I,samplemask,1,1,options_thin,'r');
%   figure, imshow(thinIden);
%   % or refer to ImAnalysis.m

% Function is written by YikTungTracy Ling, Johns Hopkins University (July 2019)
% Reference: Ling et al. 'Pressure-Induced Changes in Astrocyte GFAP, Actin
% and Nuclear Morphology in Mouse Optic Nerve' IOVS 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defaultoptions = struct('Lwin',500 , 'Swin', 30, 'sthin',0);

% Process inputs
if(~exist('options','var'))
    options=defaultoptions; 
else
    tags = fieldnames(defaultoptions);
    for i=1:length(tags)
         if (~isfield(options,tags{i}))
             options.(tags{i})=defaultoptions.(tags{i}); 
         end
    end
    if  (length(tags)~=length(fieldnames(options)))
        error('process_thin: unknown option','unknown options found');
    end
end


    %%%%%%%%%%%%%%%%%%% median filter %%%%%%%%%%%%%%%%%%%%%%%%%%
    medfiltthin=medfilt2(thinI);
    %%%%%%%%%%%%%%% contrast enhancement %%%%%%%%%%%%%%%%%%%%%
    hisadapt=adapthisteq(medfiltthin,'NumTiles' ,[options.Lwin options.Lwin]);
    hisadapt=adapthisteq(hisadapt,'NumTiles' ,[options.Swin options.Swin]);
    showsingle(hisadapt,color,showfig);
    %%%%%%%%%%%%%%%%%% otsu binarization %%%%%%%%%%%%%%%%%%%%
    thres=graythresh(hisadapt);
    binthin= imbinarize(hisadapt,thres); 
    binthin=maskI.*binthin;
    %%%%%%%%%%% morphological smoothing %%%%%%%%%%%%%%%%%%%%
    thin1=bwmorph(binthin,'bridge');
    thin1 = bwmorph(thin1,'clean');
    thin1 = bwmorph(thin1,'majority');
    se=strel('disk',options.sthin);
    thin1=imopen(thin1,se);
    
    bound=bwmorph(thin1,'remove');
    clear binthin medfiltthin
    %  figure
    %  imshow(actin1(3198:3930,2356:3090))

    medfiltthin1=hisadapt;
    medfiltthin1(bound==1)=255;
    medfiltthin2=hisadapt;
    medfiltthin2(bound==1)=0;
    % clear actinden
    thinIden(:,:,1)=medfiltthin1;
    thinIden(:,:,2)= medfiltthin2;
    thinIden(:,:,3)=medfiltthin2;
    
    %%%%%%%% calculate orientation anisotropy %%%%%%%%
    [kfit_thin,mean_orient_thin]=network_kap(thin1,5,floor((min(size(thin1))-1)/2),showfig);
    mean_orient_thin=rad2deg(mean_orient_thin);
    
    %%%%%%%%%%% record thinI features %%%%%%%%%%%%%%%%%%%%
    thinarea=sum(sum(thin1));
    mask_area=sum(sum(maskI));
    
    thinstat.thinarea=thinarea;
    thinstat.thinarea_um=thinarea*res^2;
    thinstat.thindensity=thinarea/mask_area;
    thinstat.thinkfit=kfit_thin;
    thinstat.thinori=mean_orient_thin;
    
    if showfig==1
        msgfig = msgbox('Actin Channel Processed','Completed','modal');
        uiwait(msgfig)
        disp('Continue Image processing.');
        close all
    end
end