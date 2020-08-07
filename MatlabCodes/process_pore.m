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
%   or refer to ImAnalysis.m


%  Function is written by YikTungTracy Ling, 
%  Johns Hopkins University (July 2019)
%  Reference: 
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
    %  figure
    %  imshow(actinpore(3198:3930,2356:3090))
    % clear  GFAPmask1 
    % imshow(actinpore)
    poreL = bwlabel(pore);
    % find nucleus label that overlapped with actin pores
    findNuc=(removeI).*poreL;
    % find nucleus label that overlapped with actin pores
    maskoutline1=bwmorph(maskI,'remove');
    maskoutline1 = imdilate(maskoutline1, strel('disk',2));
    findNuc1=maskoutline1.*poreL;
    %end
    Nuclist=unique([findNuc,findNuc1]);
    Nuclist=Nuclist(2:end);
    retain = (~ismember(poreL, Nuclist)).*poreL;
    % Convert to binary
    pore=retain>0;
    pore=bwareaopen(pore,options.spore);
    porestats=regionprops(pore,'Area','MajorAxisLength','MinorAxisLength');
    %actinporecount=numel(porestats);
    %allAreas = [porestats.Area];
    poreAR=[porestats.MajorAxisLength]./[porestats.MinorAxisLength];
    outlim=options.outlim;
    stay=find(poreAR>outlim);
    retain1 = ismember(retain,stay);
    retain1=bwareaopen(retain1,options.spore);
    retain1=bwlabel(retain1);
    
    clear porestats retain DapiB4open poreAR actinporeL
    if showfig==1
        sx=size(maskI,1);
        Lrgb = label2rgb(retain1, 'jet', 'w', 'shuffle');
        figure
        im1=imshow(beamI(sx/3:sx/3+800,sx/3:sx/3+800));
        im1.AlphaData = 0.5;
        hold on
        himage = imshow(Lrgb(sx/3:sx/3+800,sx/3:sx/3+800,:));
        himage.AlphaData = 0.5;
        hold off
        clear Lrgb
        msgfig = msgbox('process pores','Step4','modal');
        uiwait(msgfig)
        disp('finalizing');
        close all
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