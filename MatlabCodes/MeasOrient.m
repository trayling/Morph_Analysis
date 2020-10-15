function [BeamOrient,BeamKap]=MeasOrient(SkelI,BWI,bini,showfig)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     This function MeasOrient measures Orientation of each identified
%     branch
%
%     [BeamOrient,BeamKap]=MeasOrient(SkelI,BWI,bini,showfig);
% inputs,
%     skelI:        The 2D skeleton image from the network, where pixels in
%                   region of interest = 1 and pixels in the backgroun =0
%     BWI:          The 2D logical/binary image with beams to be measured
%     bini:         for visualization, approximate
%                   average beam width(in pixels)
%     showfig:      indicator, 1= show figures when running the 
%                   program, 0= don't show the figures  
% outputs,
%     beamOrient:   circular average beam orientation (pixels)
%     BeamKap:      fitted dispersion parameter from semi-circular
%                   von-mises distribution    
% example, 
% % create image with horizontal or vertical lines
%   data = zeros(100,100);
%   data([10:10:end],:) = 1;
%   figure
%   imshow(data)
%   BW=data;
%   SkelI = bwmorph(BW,'thin',Inf); % Skeletonize the beams
%   [beamNodeI,beamI]=MergeNodes(SkelI,1);
%   [beamOrient,BeamKap]=MeasOrient(SkelI,BW,1,1);  

% Function is written by YikTungTracy Ling, Johns Hopkins University (July 2019)
% Reference: Ling, Y. T. T., Pease, M. E., Jefferys, J. L., Kimball, E. C., Quigley, H. A., 
% & Nguyen, T. D. (2020). Pressure-Induced Changes in Astrocyte GFAP, Actin, and Nuclear 
% Morphology in Mouse Optic Nerve. Investigative Ophthalmology & Visual Science, 61(11), 14-14.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [~,~,branch_lb]=MergeNodes(SkelI,bini);
    branch_lb1=branch_lb;
    no_branch=max(max(branch_lb)); % number of identified branches
    %pre-allocate output variables
    blstats = regionprops(branch_lb,'Perimeter','Orientation');
    bOrient=[blstats.Orientation];
       
       if showfig==1
            for l=1:no_branch
               branch_lb1(branch_lb==l)=bOrient(l);
            end
            branch_lb1=imdilate(branch_lb1,strel('disk',bini));
            dilatedbin=branch_lb1.*BWI;
            dilatedbin(BWI==0)=NaN;
            figure
            imagesc(dilatedbin)
            colormap jet
            colorbar
            title('colour map of beam orientation')
            figure
            histogram(bOrient)
            title('historgram of beam orientation')
       end
       branch_distance=[blstats.Perimeter]./2;
       % average orientation weighted by beam length
       %mean_orient=circ_mean(bOrient,branch_distance,2);
       [BeamKap,BeamOrient]=fitvonmises_cdf(bOrient,branch_distance,1);

end
       