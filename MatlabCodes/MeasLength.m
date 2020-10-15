function [totbeamlen,beamlen,beamlen_std,bstraightlen,bstraightlen_std,tort,tort_std]=MeasLength(SkelI,BWI,bini,res,showfig)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     This function MeasLength measures contour length, straigh line
%     distance and tortuosity/waviness of the skeletonized beams
%
%     [beamlen,beamlen_std,bstraightlen,bstraightlen_std,tort,tort_std]=MeasLength(skelI,bini,showfig);
% inputs,
%     skelI:        The 2D skeleton image from the network, where pixels in
%                   region of interest = 1 and pixels in the backgroun =0
%     BWI:          The 2D logical/binary image with beams to be measured
%     bini:         for visualization only, approximate
%                   average beam width(in pixels)
%     res:          resolution (um/pixel or any units/pixel); res=1 if just
%                   want to output pixels
%     showfig:      indicator, 1= show figures when running the 
%                   program, 0= don't show the figures  
% outputs,
%     totbeamlen:   total 
%     beamlen:      average contour length of all beams (pixels)
%     bstraightlen: average straight line distance of all beams (pixels)
%     tort:         average toruosity/waviness of all beams
% example, 
%     I=imread('samplegreen.tif'); % greyscale image
%     BWI=imbinrize(I,graythresh(I)); % binarize the image with beams
%     load('samplemask.mat');
%     SkelI = bwmorph(lastBW.samplemask,'thin',Inf); % Skeletonize the beams
%     [beamNodeI,beamI]=MergeNodes(SkelI,BWI,8,1,1);
%     MeasLength(skelI,8,showfig)  
%     or refer to Ex1_ImAnalysis.m

% Function is written by YikTungTracy Ling, Johns Hopkins University (July 2019)
% Referece: Ling, Y.T.T., Shi, R., Midgett, D.E., Jefferys, J.L., Quigley, H.A. and Nguyen, T.D., 
% 2019. Characterizing the collagen network structure and pressure-induced strains of the human 
% lamina cribrosa. Investigative ophthalmology & visual science, 60(7), pp.2406-2422.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [~,beamI,branch_lb]=MergeNodes(SkelI,bini);
    branch_lb1=branch_lb;
    no_branch=max(max(branch_lb)); % number of identified branches
    %pre-allocate output variables
    branch_distance=NaN(no_branch,1);
    blstats = regionprops(branch_lb,'Perimeter');
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
            branch_lb1=imdilate(branch_lb1.*res,strel('disk',bini*2));
            dilatedbin=branch_lb1.*BWI;
            dilatedbin(BWI==0)=NaN;
            figure
            imagesc(dilatedbin)
            colormap jet
            colorbar
            
            figure
            histogram(branch_distance.*res)
       end
       
       branch_pixel=[blstats.Perimeter]/2;
       % total beam length
       totbeamlen=sum(branch_pixel,'all');
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
       