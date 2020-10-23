function [thickImask1,thickIden,histthickI,outline,thickstat]=process_thick(thickI,showfig,res,options,color)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This function process_thick processes images with thick beam network. It
%  binarizes grey scale images, smooths them and calculates network
%  properties. Process_thick also separates the area of interest into
%  central, peipheral and rim regions, refer to process_thin for an example
%  without separating regions
% 
% [thickImask1,thickIden,outline,thickstat]=process_thick(thickI,showfig,res,options,color)
% inputs,
%   I:          The 2D greyscale or binary image to be converted to colour
%               image if the image is binary, input should be I*255;
%   color:      must be 'r','g' or 'b' ,text for indicating whether output
%               image should be red green or blue in colour
%   showfig:    indicator, 1= show figures when running the program, 0=
%               don't show the figures 
% inputs,
%   thickI :        The original greyscale image with thick network
%   maskI :         The 2D logical/binary image indicating area of interest, 
%                   where pixels in region of interest = 1 and pixels in 
%                   the background =0
%   res:            resolution of the image, um/pixel
%   showfig:        Optional, indicator, 1= show figures when running the 
%                   program, 0= don't show the figures
%   color:          must be 'r','g' or 'b' ,text for indicating whether output
%                   image should be red green or blue in colour
%   options:        Stuct with input options:
%       .Lwin:      Large window size for more global contrast enhancement
%       .Swin:      Small window size for local contrast enhancement
%       .threshpix: radius of disk to be used for morphological closing for
%                   autoboundary detection of the area of interest, set to
%                   0 if not needed
%       .sthick:    radius of disk to be used for morphological opening of
%                   the branches, , set to 0 if not needed
%       .skelthr:   threshold for removing spurious small branches, larger number
%                   remove longer branch that are sticking out from the main
%                   branch, , set to 0 if not needed
%       .bini:      intial average beam width in pixels to merge nodes that are
%                   closed to each other for skeletonization
% outputs,
%    outline:      2D image with final identified network outlined on the
%                   median filtered image
%    thickIden:     Final 2D Binarized image after smoothing  
%                   for calculation of thin network properties
%    histthickI:    Final 2D greyscale image just before
%                   binarization, used for visualization
%    thickstat:     struct with measured properties 
%     .mask_area                         % total area of ROI (pixels^2)
%     .mask_AREA_um                      % total area of ROI (um^2)   
%     .mask_ar                           % aspect ratio of ROI    
%     .thickarea                         % total area of thick signals in ROI (pixels^2)
%     .thickIarea_um                     % total area of thick signals in ROI (um^2)
%     .thickdensity                      % area fraction (thickarea/mask_area) 
%     .thickcenden                       % area fraction in central zone
%     .thickperiden                      % area fraction in peripheral+rim zone
%     .thickperi90den                    % area fraction in rim zone
%     .thickrimden                       % area fraction in peripheral zone
%     .bwidth                            % average beam width 
%     .bwistd                            % standard deviation of beam width
%     .bwimax                            % maximum beam width 
%     .bwicen                            % average beam width in central zone
%     .bwicenstd                         % standard deviation of beam width in central zone
%     .bwiperi                           % average beam width in peripheral+rim zone
%     .bwiperistd                        % standard deviation of beam width in peripheral +rim zone    
%     .bwiperi90                         % average beam width in peripheral zone
%     .bwiperi90std                      % standard deviation of beam width in peripheral zone
%     .bwirim                            % average beam width in rim zone
%     .bwirimstd                         % standard deviation of beam width in rim zone
%     .thickkfit                         % network anisotropy (fitted k from von-mises distribution)                            
%     .thickori                          % circular mean orientation of the network
%
% example, 
%   I=imread('samplegreen.tif'); % can be greyscale image
%   options = struct('Lwin',500 , 'Swin', 30, 'threshpix', 300, 'sthick', 3, 'skelthr', 20, 'bini',8);
%   [thickImask1,thickIden,outline,thickstat]=process_thick(I,1,1,options,'g');
%   figure, imshow(outline);
%   or refer to Ex1_ImAnalysis.m

% Function is written by YikTungTracy Ling, Johns Hopkins University (July 2019)
% Reference: Ling, Y. T. T., Pease, M. E., Jefferys, J. L., Kimball, E. C., Quigley, H. A., 
% & Nguyen, T. D. (2020). Pressure-Induced Changes in Astrocyte GFAP, Actin, and Nuclear 
% Morphology in Mouse Optic Nerve. Investigative Ophthalmology & Visual Science, 61(11), 14-14.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    defaultoptions = struct('Lwin',500 , 'Swin', 30, 'threshpix', 300, 'sthick', 3, 'skelthr', 20, 'bini',8);

    % Process inputs
    if(~exist('options','var'))
        options=defaultoptions; 
    else
        tags = fieldnames(defaultoptions);
        for i=1:length(tags)
             if (~isfield(options,tags{i}))
                 options.(tags{i})=defaultoptions.(tags{i}); 
             end
        end
        if  (length(tags)~=length(fieldnames(options)))
            error('process_thick: unknownoption','unknown options found');
        end
    end
    % %%%%%%%%%%%%%%%%% preprocessing and binarization %%%%%%%%%%%%%%%%%
    [histthickI,B1]=PreProcess_Binarize(thickI,options.Lwin,options.Swin,color,showfig);
    
    %%%%%%%%%%%%%%%%% auto detect boundary & %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%% seperate central pheroheral and rim zone %%%%%%%%%%%%%%%%%%%%
    [thickImask1,outline,thickIcen,thickIperi,thickIperi90,thickIrim,mask_area,mask_ar,mcen_area,mperi_area,mperi90_area,mrim_area]=AutoBoundary(histthickI,B1,options.threshpix,showfig);
% %     if showfig ==1
% %         msgfig = msgbox('Auto Boundary detection and regional mask generation','Step2','modal');
% %         uiwait(msgfig)
% %         disp('Continue Image processing.');
% %         close all
% %     end
    
    %%%%%%%%%%% morphological smoothing %%%%%%%%%%%%%%%%%%%%
    B1=thickImask1.*B1;
    % figure
    % imshow(B1)
    thickI1= bwmorph(B1,'clean');
    thickI1 = bwmorph(thickI1,'majority');
    thickI1=bwmorph(thickI1,'bridge');
    se=strel('disk',options.sthick);
    thickI1=imopen(thickI1,se);
    
    bound=bwmorph(thickI1,'remove');
    thicksurf=sum(sum(bound)); % perimeter of processes
    thickarea=sum(sum(thickI1)); % 
    thickcenarea=sum(sum(thickI1.*thickIcen)); % beam area in central zone
    thickperiarea=sum(sum(thickI1.*thickIperi));% beam area in peri +rim zone
    thickperi90area=sum(sum(thickI1.*thickIperi90));% beam area in peri zone
    thickrimarea=sum(sum(thickI1.*thickIrim));% beam area in rim zone
    

    %%%%%%%%%%%%%%%%% figure: overlay BW tracing on hist%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %     hisfilterthickI1=histthickI;
    %     hisfilterthickI1(bound==1)=255;
    %     hisfilterthickI2=histthickI;
    %     hisfilterthickI2(bound==1)=0;
    %     thickIden(:,:,1)=hisfilterthickI1;
    %     thickIden(:,:,2)= hisfilterthickI2;
    %     thickIden(:,:,3)=hisfilterthickI2;
    thickIden=thickI1;
    if showfig==1
        figure,imshow(thickIden)
        title('final binarized image after morphological smoothing')
    end
    
    
    clear B1 medfilterthickI1 medfilterthickI2 bound medfilterthickI medfilterthickI1 medfilterthickI2 bound
   
    %%%%%%%%%%%%%%%%%%%%%% skeletonization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    lastBW=thickI1;
    Skel2 = bwmorph(lastBW,'thin',Inf);
    Skel2=bwmorph(Skel2,'spur',options.skelthr);
    Skel2=bwmorph(Skel2,'clean');
    
    % Identify nodes from beam network and merge nodes that are close together   
    [beamNodeI,beamI,beamlabel]=MergeNodes(Skel2,options.bini);
    
    % beamI can be used to measure beam width, length, tortuosity etc
    %%%%%%%%%%%%%%%%%%%%%% calculate beam width %%%%%%%%%%%%%%%%%%%%%%
    dilsize=options.bini;
    %calculate average beam width in each zone
    %all
    [bwidth,bwistd,bwimax]=MeasWidthAll(beamI,lastBW,res,lastBW,dilsize,showfig);
    %zones
    [bwicen,bwicenstd,~]=MeasWidthAll(beamI,lastBW,res,thickIcen,dilsize,showfig);
    [bwiperi,bwiperistd,~]=MeasWidthAll(beamI,lastBW,res,thickIperi,dilsize,showfig);
    [bwiperi90,bwiperi90std,~]=MeasWidthAll(beamI,lastBW,res,thickIperi90,dilsize,showfig);
    [bwirim,bwirimstd,~]=MeasWidthAll(beamI,lastBW,res,thickIrim,dilsize,showfig);
    % % % %% optional: calculate average beam width by fitting gaussian distribution %%%
    % % % [muall,sigall]=fitgaus(distskel, Maxwid); % output distskel from
    % % % function MeasWidth
    
    % figure
    % plot([1 2 3],[bwicen bwiperi90 bwirim],'o')   
    %%%%%%%% calculate orientation anisotropy %%%%%%%%
    [kfit,mean_orient]=network_kap(lastBW,10,floor((min(size(lastBW))-1)/2),showfig);
    mean_orient=rad2deg(mean_orient);
    %%%%%%%%%%% record thickI features %%%%%%%%%%%%%%%%%%%%
    %record lamina area 
    thickstat.mask_area=mask_area;
    thickstat.mask_AREA_um=mask_area*res^2;
    %record aspect ratio of lamina
    thickstat.mask_ar=mask_ar;
    %thickstat.thicksurf=thicksurf;
    thickstat.thickarea=thickarea;
    thickstat.thickIarea_um=thickarea*res^2;
    thickstat.thickdensity=thickarea/mask_area;
    thickstat.thickcenden=thickcenarea/mcen_area;
    thickstat.thickperiden=thickperiarea/mperi_area;
    thickstat.thickperi90den=thickperi90area/mperi90_area;
    thickstat.thickrimden=thickrimarea/mrim_area;
    thickstat.bwidth=bwidth; 
    thickstat.bwistd=bwistd; 
    thickstat.bwimax=bwimax;
    thickstat.bwicen=bwicen; 
    thickstat.bwicenstd=bwicenstd;
    thickstat.bwiperi=bwiperi;
    thickstat.bwiperistd=bwiperistd;
    thickstat.bwiperi90=bwiperi90;
    thickstat.bwiperi90std=bwiperi90std;
    thickstat.bwirim=bwirim;
    thickstat.bwirimstd=bwirimstd;
    thickstat.thickkfit=kfit;
    thickstat.thickori=mean_orient;
    clear  beamNodeI beamI Skel2 lastBW
end