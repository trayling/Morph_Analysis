function [beamlen,beamlen_std,bstraightlen,bstraightlen_std,tort,tort_std]=MeasLength(SkelI,bini,showfig)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     This function MeasLength measures contour length, straigh line
%     distance and tortuosity/waviness of the skeletonized beams
%
%     [beamlen,beamlen_std,bstraightlen,bstraightlen_std,tort,tort_std]=MeasLength(skelI,bini,showfig);
% inputs,
%     skelI:        The 2D skeleton image from the network, where pixels in
%                   region of interest = 1 and pixels in the backgroun =0
%     bini:         Optional, indicator, 1= show figures when running the 
%                   program, 0= don't show the figures
%     showfig:      indicator, 1= show figures when running the 
%                   program, 0= don't show the figures  
% outputs,
%     beamlen:      average contour length of all beams (pixels)
%     bstraightlen: average straight line distance of all beams (pixels)
%     tort:         average toruosity/waviness of all beams
% example, 
%     I=imread('samplegreen.tif'); % greyscale image
%     BWI=imbinrize(I,graythresh(I)); % binarize the image with beams
%     load('samplemask.mat');
%     SkelI = bwmorph(lastBW.samplemask,'thin',Inf); % Skeletonize the beams
%     [beamNodeI,beamI]=MergeNodes(SkelI,8);
%     MeasLength(skelI,8,showfig)  
%     or refer to runAnalysis.m


%  Function is written by YikTungTracy Ling, 
%  Johns Hopkins University (July 2019)
% Reference: Ling et al. 'Pressure-Induced Changes in Astrocyte GFAP, Actin
% and Nuclear Morphology in Mouse Optic Nerve' IOVS 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [~,beamI]=MergeNodes(SkelI,bini);
    branch_lb=bwlabel(beamI);
    branch_lb1=branch_lb;
    no_branch=max(max(branch_lb)); % number of identified branches
    %pre-allocate output variables
    branch_distance=NaN(no_branch,1);
    blstats = regionprops(beamI,'Perimeter','Orientation');
       % record pixel distance
       for l=1:no_branch
           [posy,posx]=find(branch_lb==l);
           len=size(posx,1);

           pos=sub2ind(size(beamI),posy,posx);

           % calculate straight line distance
           branch_dist=pdist2([posy(1),posx(1)],[posy(len),posx(len)]);
           branch_distance(l)=branch_dist;
           branch_lb1(branch_lb==l)=branch_dist;
       end
       
       if showfig==1
           figure
           imagesc(branch_lb1)
           colormap jet
           colorbar
           caxis([0 max(distskel,[],'all')])
       end
       
       branch_pixel=[blstats.Perimeter]/2;
       % calculate tortuosity/waviness of each beam
       branch_turtuosity=branch_pixel./branch_dist;
       %record mean and std
       beamlen=mean(branch_pixel);
       beamlen_std=std(branch_pixel);
       bstraightlen=mean(branch_distance);
       bstraightlen_std=std(branch_distance);
       tort=mean(branch_turtuosity);
       tort_std=std(branch_turtuosity);
end
       