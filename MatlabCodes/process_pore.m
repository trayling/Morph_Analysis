function [pstat]=process_pore(maskI,beamI,removeI,showfig,options) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This function process_pore measures negative space of beams and removes
%  the nuclear space to calculate pore space that represent axons
% 
% [pstat]=process_pore(maskI,beamI,removeI,showfig,options) ;
% 
% inputs,
%     maskI:        Optional, If there's the image is further seperated
%                   into different regions, a binary mask image could be
%                   supplied to measure average width in the region only.
%                   Mask is applied after skeletonization and beam
%                   measurement to avoid error in area that are cutt of for
%                   width measurement
%     beamI:        The 2D binary image of the beams, where pixels in
%                   region of interest = 1 and pixels in the backgroun =0
%     removeI:      Optional, if the negative space of beams include cell
%                   nulcear that needs to be removed
%     res:          resolution of the image, um/pixel
%     options:      structure: 
%           .outlim: eliminate pores with aspect ratio > x                           
%           .spore:  area smaller than this number of pixels are removed/
%                    not counted as pores
%     showfig:      Optional, indicator, 1= show figures when running the 
%                   program, 0= don't show the figures
% outputs,
%   pstat:                  struct with measured properties
%       .porecount:         %number of seperated pore space
%       .poreareatotal      %total pore area (pixel^2)
%       .parea_mean:        %mean area of the pores (pixel^2)
%       .parea_med:         %median area of the pores (pixel^2)
%       .parea_Std:         %standard deviation of pore areas (pixel^2)
%
% example, 
%   beamI=imbinarize(imread('samplered.tif'));
%   removeI=imbinarize(imread('sampleblue.tif'));
%   load('samplemask.mat');  % import  ROI mask
%   options = struct('outlim', 0, 'spore',0);
%   [pstat]=process_pore(samplemask,beamI,removeI,1,options);
%   or refer to Ex1_ImAnalysis.m for more image preprocessing


% Function is written by YikTungTracy Ling, Johns Hopkins University (July 2019)
% Reference: Ling, Y. T. T., Pease, M. E., Jefferys, J. L., Kimball, E. C., Quigley, H. A., 
% & Nguyen, T. D. (2020). Pressure-Induced Changes in Astrocyte GFAP, Actin, and Nuclear 
% Morphology in Mouse Optic Nerve. Investigative Ophthalmology & Visual Science, 61(11), 14-14.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defaultoptions = struct('outlim', 0, 'spore',0);

% Process inputs
if (~exist('options','var')) 
    options=defaultoptions; 
else
    tags = fieldnames(defaultoptions);
    for i=1:length(tags)
         if (~isfield(options,tags{i}))
             options.(tags{i})=defaultoptions.(tags{i}); 
         end
    end
    if  (length(tags)~=length(fieldnames(options)))
        error('process_blue: unknown option','unknown options found');
    end
end

    % save coordinates where nucleus exist
    % find corresponding labeled region in actin
    % set region to not pore
    pore= imcomplement(beamI).*(maskI);
    %%%%%%%%%%%%%%%%%%%% remove pores that are touching the boundary %%%%%%%%%%%%%%%%%%%%%%%%
    [pore]=removeBoundPore(pore,maskI,showfig);
    %%%%%%%%%%%%%%%%%%%%% removes pores that belongs to nucleus %%%%%%%%%%%%%%%%%%%%%%%%%
    poreL = bwlabel(pore);
    % find nucleus label that overlapped with actin pores
    findNuc=(removeI).*poreL;
    Nuclist=unique([findNuc]);
    Nuclist=Nuclist(2:end); % the first element is 0 no a pore label
    retain = (~ismember(poreL, Nuclist)).*poreL;
    % Convert to binary
    pore=retain>0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %remove pores less than x number of pixels
    pore=bwareaopen(pore,options.spore);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %remove pores with aspect ratio larger than x
    porestats=regionprops(pore,'Area','MajorAxisLength','MinorAxisLength');
    poreAR=[porestats.MajorAxisLength]./[porestats.MinorAxisLength];
    outlim=options.outlim;
    stay=find(poreAR>outlim);
    retain1 = ismember(retain,stay);
    retain1=bwareaopen(retain1,options.spore);
    retain1=bwlabel(retain1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear porestats retain DapiB4open poreAR actinporeL
    if showfig==1
        Lrgb = label2rgb(retain1, 'jet', 'w', 'shuffle');
        figure
        im1=imshow(beamI);
        im1.AlphaData = 0.5;
        hold on
        himage = imshow(Lrgb);
        himage.AlphaData = 0.5;
        hold off
        clear Lrgb
%         msgfig = msgbox('process pores','Step4','modal');
%         uiwait(msgfig)
%         disp('finalizing');
%         close all
    end

    % actinden = imresize(actinden,[savesize savesize]);
    % imwrite(actinden, [namd_actinden savefname num2str(i) '.tif']);
    porestats=regionprops(retain1,'Area');
    %poreAR=[porestats.MajorAxisLength]./[porestats.MinorAxisLength];
    % figure
    % histogram(poreAR)
    porecount=numel(porestats);
    allAreas = [porestats.Area];
    
    pstat.porecount=porecount;
    pstat.poreareatotal=sum(allAreas);
    pstat.parea_med=median(allAreas);
    pstat.parea_mean=mean(allAreas);
    pstat.parea_std=std(allAreas);


end