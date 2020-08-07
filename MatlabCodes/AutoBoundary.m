function [mask1,foutline,mcen,mperi,mperi90,mrim,marea,mAR,mcen_area,mperi_area,mperi90_area,mrim_area]=AutoBoundary(refI,B1,threshpix,showfig)    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This function AutoBoundary identifies a region of interest from an image
%  by selecting the region with largest area within an image, and extract
%  its boundary as an output. The mask for the region of interest is further
%  seperated into central peripheral and rim zones. 
%  Total area (number of pixels)for each zone will also be calculated as an
%  output.

% [GFAPmask1,outline,GFAPcen,GFAPperi,GFAPperi90,GFAPrim,lam_area,lam_ar,
% lamcen_area,lamperi_area,lamperi90_area,lamrim_area]=AutoBoundary(refI,B1,showfig);
% 
% inputs,
%     refI :      Original 2D greyscale image  
%     B1:         Binarized image
%     showfig:    indicator, 1= show figures when running the program, 0=
%                 don't show the figures
% outputs,
%     mask1:        logical image/masks for region of interest (ROI)
%     moutline:     image with boundary outline overlayed on to original
%                   image, (uint8)
%     mcen:         logical image/masks for central zone (0-50%)
%     mperi:        logical image/masks for peripheral+rim zone (50-100%)
%     mperi90:      logical image/masks for peripheral zone (50-90%)
%     mrim:         logical image/masks for rim zoneI (90-100%)
%     marea:        area of ROI (pixels)
%     mAR:          aspect Ratio of ROI 
%     mcen_area:    area of centrol zone (pixels)
%     mperi_area:   area of peripheral + rim zone (pixels)
%     mperi90_area: area of pheripheral zone (pixels)
%     mrim_area:    area of rim zone (pixels)
%
% example, 
%   refI=imread('samplegreen.tif'); % greyscale image
%   BW=imbinarize(adapthisteq(refI),graythresh(refI)); % binarize the image
%   showfig=1;
%   [mask1,foutline,mcen,mperi,mperi90,mrim,marea,mAR,mcen_area,mperi_area,mperi90_area,mrim_area]=AutoBoundary(refI,BW,300,showfig);
%   figure, imshow(foutline);
%   or refer to Ex1_ImAnalysis for better quality
%
%   Function is written by YikTungTracy Ling, Johns Hopkins University (July 2019)
% Reference: Ling et al. 'Pressure-Induced Changes in Astrocyte GFAP, Actin
% and Nuclear Morphology in Mouse Optic Nerve' IOVS 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Z1 = bwmorph(B1,'majority');
    fill= imfill(Z1,'holes');
    regionlabel=bwlabel(fill);
    if showfig==1
        figure
        imshow(fill)
        Lrgb = label2rgb(regionlabel, 'jet', 'w', 'shuffle');
        figure
        imshow(Lrgb)
    end
   
    %find region with largest area = region of interest = lamina
    findregion=regionprops(regionlabel,'Area');
    area=[findregion.Area];
    maxarea=max(area);
    index=find(area==maxarea);
    regionlabel(regionlabel~=index)=0;
    % figure
    % imshow(regionlabel)

    mask=logical(regionlabel);
    %smooth region boundary
    se = strel('disk',threshpix);
    maskdilate=imclose(mask,se);
    mask1=imfill(maskdilate,'holes');
    if showfig==1
        figure
        imshow(mask1)
    end
    % if mask was manually drawn
%     if isempty(find(indivstack==i))~=1
%         clear mask1
%         load([namd_GFAPbound 'mask' num2str(i) '.mat']);
%     end
    %

    Lamstats = regionprops(mask1,'MinorAxisLength','MajorAxisLength');
    mAR=[Lamstats.MajorAxisLength]/[Lamstats.MinorAxisLength];
    clear Lamstats maxarea findregion area
    % % imshow(mask1)
    maskoutline=bwmorph(mask1,'remove');
    maskoutline = imdilate(maskoutline, strel('disk',5));
    % figure
    % imshow(maskoutline)

    %% seperate central and peripheral GFAP mask
    [mcen,mperi,mperi90,mrim]=seperate_CPR_zones(mask1,showfig);
    
    %% visualize mask 
    foutline(:,:,1)=uint8(maskoutline).*255;
    foutline(:,:,2)= refI;
    foutline(:,:,3)=uint8(zeros(size(maskoutline)));
    if showfig ==1
    figure
    imshow(uint8(foutline))
    end
    % hold on
    % plot(cenx,ceny,'r-','LineWidth',2)
    
    %% record measured features
    marea=sum(sum(mask1));
    mcen_area=sum(sum(mcen));
    mperi_area=sum(sum(mperi));
    mperi90_area=sum(sum(mperi90));
    mrim_area=sum(sum(mrim));

    clear moutline out1 out2 thickout maskoutline GFAPmaskdilate GFAPdilate GFAPerode Z1 Z2  GFAPmask regionlabel fill Lrgb
