function [beamNodeI,beamI,beamlabel]=MergeNodes(SkelI,bini)
%     This function MergeNodes merges nodes that are too close to each other 
%     after skelotonization procedure. 
%     i.e. distance between centroids of two nodes should be larger than 
%     initial estimate of beam width  
%
%     [beamNodeI,beamI]=MergeNodes(SkelI,bini)
% inputs,
%     skelI:        The 2D skeleton image from the network, where pixels in
%                   region of interest = 1 and pixels in the backgroun =0
% 
%     bini:         approximate average beam width(in pixels)
% outputs,
%     beamNodeI:    2D binary image of beam skeleton plus dilated node map
% 
%     beamI:        2D binary image of beam skeleton after removing dilated
%                   nodes, suitatble for measuring beam width, length,
%                   tortuosity etc.
%
% example, 
%     I=imread('samplegreen.tif'); % greyscale image
%     BWI=imbinrize(I,graythresh(I)); % binarize the image with beams
%     load('samplemask.mat');
%     SkelI = bwmorph(lastBW.*samplemask,'thin',Inf); % Skeletonize the beams
%     [beamNodeI,beamI,beamlabel]=MergeNodes(SkelI,8);
%     or refer to runAnalysis.m


% Function is written by YikTungTracy Ling, Johns Hopkins University (July 2019)
% Reference: Ling, Y. T. T., Pease, M. E., Jefferys, J. L., Kimball, E. C., Quigley, H. A., 
% & Nguyen, T. D. (2020). Pressure-Induced Changes in Astrocyte GFAP, Actin, and Nuclear 
% Morphology in Mouse Optic Nerve. Investigative Ophthalmology & Visual Science, 61(11), 14-14.
 
%% start function
    nodes = bwmorph(SkelI,'branchpoints');
    no_nodes=sum(sum(nodes));
    nse = strel('disk',round(bini/2)); 
    node_dil = imdilate(nodes,nse);
   
    %%%%% combine nodes to find the new centroid %%%%%
    %     node_dil1 = bwlabel(node_dil);
    %     stats = regionprops(node_dil1,'Centroid','PixelList');
    %     no_node=size(stats,1);
    %     for id=1:no_node
    %     cen=round(getfield(stats,{id},'Centroid'));
    %     nodes_positions_round(id,1)= cen(1,2);
    %     nodes_positions_round(id,2)= cen(1,1);
    %     end
    %%%%%
    
    % beam skeleton with dilated node map
    beamNodeI=logical(SkelI) |  logical(node_dil);

    % beam skeleton after removing dilated nodes
    beamI=SkelI.*(~node_dil);
    beamlabel=bwlabel(beamI);

end