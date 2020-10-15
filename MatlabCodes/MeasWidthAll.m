 
function [Avgwid,Stdwid,Maxwid]=MeasWidthAll(beamI,BWI,res,maskI,dilsize,showfig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This function measures an average beam width within a masked area
% 
% [Avgwid,Stdwid,Maxwid]=MeasWidthAll(beamI,BWI,res,maskI,dilsize,showfig);
% 
% inputs,
%     beamI:        The 2D skeleton image of the beams, where pixels in
%                   region of interest = 1 and pixels in the backgroun =0

%     BWI:          The 2D logical/binary image with beams to be measured

%     res:          resolution of the image, um/pixel or any unit /pixel

%     maskI:        If there's the image is further seperated
%                   into different regions, a binary mask image could be
%                   supplied to measure average width in the region only.
%                   Mask is applied after skeletonization and beam
%                   measurement to avoid error in area that are cutt of for
%                   width measurement; otherwise maskI=BWI or maskI=ones(size(BWI))
%
%     dilsize:      Optional, for visualization only, approximate
%                   average beam width(in pixels)
%
%     showfig:      Optional, indicator, 1= show figures when running the 
%                   program, 0= don't show the figures
% outputs,
%   Avgwid:         Average beam width in the masked area
%   Stdwid:         Standard deviation of beam width in the masked area
%   Maxwid:         Maximum beam width
%
% example, 
%   I=imread('samplegreen.tif'); % can be greyscale image
%   BWI=imbinarize(adapthisteq(I),graythresh(I)); % binarize the image with beams
%   SkelI = bwmorph(BWI,'thin',Inf); % Skeletonize the beams
%   load('samplemask.mat');  % import  ROI mask
%   [Avgwid,Stdwid,Maxwid]=MeasWidthAll(SkelI,BWI,1,samplemask,8,1);
%   or
%   
%   can be used after [~,SkelI,~]=MergeNodes(SkelI,BWI,bini) and after more
%   imag preprocessing, refer to process_thick.m


% Function is written by YikTungTracy Ling, Johns Hopkins University (July 2019)
% Reference: Ling, Y. T. T., Pease, M. E., Jefferys, J. L., Kimball, E. C., Quigley, H. A., 
% & Nguyen, T. D. (2020). Pressure-Induced Changes in Astrocyte GFAP, Actin, and Nuclear 
% Morphology in Mouse Optic Nerve. Investigative Ophthalmology & Visual Science, 61(11), 14-14.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
    error('not enough input variables')
elseif nargin == 2
    res = 1;
    maskI=BWI;
    dilsize=8;
    showfig=1;    
elseif nargin == 3
    maskI=BWI;
    dilsize=8;
    showfig=1;
elseif nargin == 4
    dilsize=8;
    showfig=1;
elseif nargin == 5
    showfig=1;
end
%% start function
% distance for from all white pixels to closest black pixel
    distmap=bwdist(imcomplement(BWI))*res*2;
    % only record distance from midpoint of each beam
    distskel=distmap.*single(beamI);
    % visualize the beam width map on binarized image 
    if showfig==1
        distdil=imdilate(distskel,strel('disk',dilsize*2));
        dilatedbin=distdil.*BWI.*maskI;
        dilatedbin(BWI==0)=NaN;
        figure
        imagesc(dilatedbin)
        colormap jet
        colorbar
        % caxis([0 max(distskel,[],'all')])
    end
    
    % imagesc(fliplr(flipud(dilatedbin)))
    % hold on
    % plot(cenx,ceny,'r-','LineWidth',2)
    % hold off

    %calculate average beam width in each zone
    width_mask=distskel.*maskI;
    width_mask(BWI==0)=nan;
    Avgwid=nanmean(width_mask,'all'); 
    Stdwid=nanstd(width_mask,0,'all'); 
    Maxwid=max(width_mask,[],'all');
    
%     %% calculate average beam width by fitting gaussian distribution %%%
%     [muall,sigall]=fitgaus(distskel, Maxwid);

  
end
