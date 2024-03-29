function [connectivity]=connectiv(node_dil,no_nodes,beamI,beamlabel,showfig)
%     This function connectiv counts number of beams connecting to each dialated node. 
%     i.e. connectivity
%
%     [connectivity]=connectiv(node_dil,no_nodes,beamI,beamlabel)
% inputs,
%     node_dil:      The 2D binary image of nodes that were dilated
%     no_nodes:      number of nodes or write as max(max(bwlabel(node_dil))
%     beamI:         2D binary image of beam skeleton after removing dilated
%                   nodes, suitatble for measuring beam width, length,
%                   tortuosity etc.
%     beamlabel:    labeled image of beamI, or write as bwlabel(beamI)
%     showfig:      Optional, indicator, 1= show figures when running the 
%                   program, 0= don't show the figures, default =0
% outputs,
%     connectivity: 2D binary image of beam skeleton plus dilated node map
%       
% example, 
%     refer to MergeNodes.m


% Function is written by YikTungTracy Ling, Johns Hopkins University (July 2019)
% Reference: Ling, Y. T. T., Pease, M. E., Jefferys, J. L., Kimball, E. C., Quigley, H. A., 
% & Nguyen, T. D. (2020). Pressure-Induced Changes in Astrocyte GFAP, Actin, and Nuclear 
% Morphology in Mouse Optic Nerve. Investigative Ophthalmology & Visual Science, 61(11), 14-14.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 4
    error('not enough input variables')
elseif nargin == 4
    showfig=0;    
end
%% start function
    nodelbl=bwlabel(node_dil);
    stats = regionprops(nodelbl,'PixelList');
    connectivity=nan(no_nodes,1);
    nodetobranch=[];
        for k=1:no_nodes
            dots=getfield(stats,{k},'PixelList'); % list in [x y] format
            dots_p=size(dots,1); % number of pixels for this dilated node
            for dd=1:dots_p % for each pixel, find surrouding 8 squares for connecting node
                for searchy=(dots(dd,2)-1):(dots(dd,2)+1)
                    for searchx=(dots(dd,1)-1):(dots(dd,1)+1)
                        if beamI(searchy,searchx)~=0
                            nodetobranch=[nodetobranch,beamlabel(searchy,searchx)]; % then record branch number
                        end
                    end
                end   
            end
            nodetobranch=unique(nodetobranch); % delete repeated branch number
            no_br=numel(nodetobranch); % number of branches from node k
        % create connectivity table --------------------------------------------
            %conntable(k,1:no_br)=nodetobranch;
            connectivity(k,1)=no_br; % number of beams connecting to this node
            nodetobranch=[];
        %
        end
        % visualize
        if showfig==1
            connvis=nan(size(beamlabel));
            for k=1:no_nodes
                connvis(nodelbl==k)=connectivity(k,1);
            end
            figure
            imagesc(connvis)
            colormap jet
            colorbar
            caxis([1 5])
        end   
   
end