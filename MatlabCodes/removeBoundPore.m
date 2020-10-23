function [poreI_nobound]=removeBoundPore(poreI,maskI,showfig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This function removeBoundPore removes pores that are connected / touching the mask boundary. 
% 
% [poreI_nobound]=removeBoundPore(poreI,maskI,showfig)
% inputs,
%   poreI:      The 2D binary image with pores = 1, backgroun=0; can use
%               imcomplement() to flip foreground (1) and background(0)
%   maskI:      The 2D binary mask image with area of interest =1, same
%               dimensions as poreI, if mask is the bounding square of the
%               image, can input as maskI=ones(size(poreI));
%   showfig:    indicator, 1= show figures when running the program, 0=
%               don't show the figures 
% onputs,
%   poreI_nobound:     Final 2D Binarized image with pores that are not
%                      touching the mask boundary
%
% example 1, 
%%create an image with grid pores
% data = zeros(100,100);
% data([5:5:end],:) = 1;
% data(:,[5:5:end]) = 1;
% poreI=imcomplement(data);
% [poreI_nobound]=removeBoundPore(poreI,ones(size(poreI)),1);

% Function is written by YikTungTracy Ling, Johns Hopkins University (July 2019)
% Reference: Ling, Y. T. T., Pease, M. E., Jefferys, J. L., Kimball, E. C., Quigley, H. A., 
% & Nguyen, T. D. (2020). Pressure-Induced Changes in Astrocyte GFAP, Actin, and Nuclear 
% Morphology in Mouse Optic Nerve. Investigative Ophthalmology & Visual Science, 61(11), 14-14.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% begin function
    poreL = bwlabel(poreI);
    % find nucleus label that overlapped with actin pores
    %findNuc=(removeI).*poreL;
    % find pores that overlapped with boundary
    maskoutline1=bwmorph(maskI,'remove');
    maskoutline1 = imdilate(maskoutline1, strel('disk',2));
    findpore1=maskoutline1.*poreL;
    porelist=unique([findpore1]); % find all pore labels to be remove
    porelist=porelist(2:end); % the first element is 0 no a pore label
    retain = (~ismember(poreL, porelist)).*poreL;
    % Convert to binary
    poreI_nobound=retain>0;
    % visualize pores
    if showfig == 1
        retainL=bwlabel(retain);
        Lrgb = label2rgb(retainL, 'jet', 'w', 'shuffle');
        figure
        im1=imshow(poreI);
        im1.AlphaData = 0.5;
        hold on
        himage = imshow(Lrgb);
        himage.AlphaData = 0.5;
        hold off
    end
end