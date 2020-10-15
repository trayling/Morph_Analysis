function [maskcen,maskperi,maskperi90,maskrim]=seperate_CPR_zones(mask1,showfig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This function seperate_CPR_zones seperates a masked region in an image
%  into central, peripheral and rim zones (0-50%, 50-90% and 90-100%
%  of the distance from centroid of the region to the outer boundary).
% 
% [maskcen,maskperi,maskperi90,maskrim]=seperate_CPR_zones(mask1);
% 
% inputs,
%   mask1 : The 2D logical/binary image, where pixels in
%           region of interest = 1 and pixels in the backgroun =0
%
% outputs,
%   maskcen,maskperi,maskperi90,maskrim: logical image/masks for 
%   1. central (0-50%)
%   2. peripheral (50-100%)
%   3. Peripheral90 (50-90%); and
%   4. rim 90-100% zones. 
%
% example, 
%   I=imread('samplegreen.tiff'); % can be greyscale image
%   BW=imbinrize(I,graythresh(I)); % binarize the image
%   regionlabel=bwlabel(BW); % label seperate regions
%   Imask=logical(regionlabel==1); % pick 1 region to segment
%   [maskcen,maskperi,maskperi90,maskrim]=seperate_CPR_zones(Imask,1)
%   figure, imshowpair(Imask,maskcen);
%   or refer to AutoBoundary.m
%
% Function is written by YikTungTracy Ling, Johns Hopkins University (July 2019)
% Reference: Ling, Y. T. T., Pease, M. E., Jefferys, J. L., Kimball, E. C., Quigley, H. A., 
% & Nguyen, T. D. (2020). Pressure-Induced Changes in Astrocyte GFAP, Actin, and Nuclear 
% Morphology in Mouse Optic Nerve. Investigative Ophthalmology & Visual Science, 61(11), 14-14.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

if nargin < 2
  showfig = 1;
end
%% start function
    maskstat=regionprops(mask1,'Centroid');
    center=maskstat.Centroid; %[x y]

    B=bwboundaries(mask1);
    B=cell2mat(B);
    [theta,rho] = cart2pol(B(:,2)-center(1),B(:,1)-center(2));
    [cenx,ceny] = pol2cart(theta,rho/2);% 0.5 * radius
    cenx=cenx+center(1);
    ceny=ceny+center(2);
    maskcen=poly2mask(cenx,ceny,size(mask1,1),size(mask1,2));
    maskperi=mask1 & ~maskcen;
    %cenoutline=bwmorph(maskcen,'remove');

    [peri90x,peri90y] = pol2cart(theta,rho*0.9);% 0.9 * radius
    peri90x=peri90x+center(1);
    peri90y=peri90y+center(2);
    maskperi90=poly2mask(peri90x,peri90y,size(mask1,1),size(mask1,2));
    maskperi90=maskperi90& ~maskcen;
    % figure
    % imshow(maskperi90)

    maskrim=maskperi & ~maskperi90;
    if showfig==1
    figure
    imshow(mask1)
    hold on
    plot(cenx,ceny,'r-','LineWidth',2)
    plot(peri90x,peri90y,'r-','LineWidth',2)
    end
end