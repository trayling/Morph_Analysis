function [noutline,nuclayer,histnuc,nuc1_calc,nucstat]=process_nuc(nucI,maskI,showfig,res,options,color)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This function process_nuc processes images stained for cell nucleus. It
%  binarizes grey scale images, smooths them and calculates nuclear
%  properties
% 
% [noutline,nuclayer,nuc1_calc,nucstat]=process_nuc(nucI,maskI,showfig,res,options,color)
% 
% inputs,
%   nucI :          The original greyscale image with cell nuclei
%   maskI :         The 2D logical/binary image indicating area of interest, 
%                   where pixels in region of interest = 1 and pixels in 
%                   the background =0
%   res:            resolution of the image, um/pixel
%   showfig:        Optional, indicator, 1= show figures when running the 
%                   program, 0= don't show the figures
%   color:          must be 'r','g' or 'b' ,text for indicating whether output
%                   image should be red green or blue in colour
%   options:        Stuct with input options:
%       .Lwin:      Large window size for more global contrast enhancement,
%                   base on size of area of interest
%       .snuc:      size of diamond structure for morphological opening
%       .threshpix: Threshold for bwareaopen, remove area under this number
%                   of pixels
% outputs,
%    noutline:      2D image with identified cell nuclei outlined on the
%                   median filtered image
%    nuclayer:      2D Binarized nuclear image for removal of negative
%                   space from thin network Image
%    histnuc:       Final 2D greyscale image just before
%                   binarization, used for visualization
%    nuc1_calc:     Final 2D Binarized image after smoothing and bwareaopen 
%                   for calculation of nuclear properties
%    nucstat:       struct with measured nuclear properties
%     .Nucarea_total:   total area of all identified nuclei (pixel^2)
%     .Nucarea_med:     median area of all nuclei (pixel^2)
%     .Nucarea_mean:    mean area of all nuclei (pixel^2)
%     .Nucarea_std:     std of area of all nuclei (pixel^2)
%     .NucAR_med:       median aspect ratio of all nuclei
%     .NucAR_mean:      mean aspect ratio of all nuclei
%     .NucAR_std:       standard deviation of aspect ratio of all nuclei
%     .nuccount:        number of nuclei 
%     .Nucden:          number of nuclei/ mask area
%     .NucAREAden:      Nucarea_total/ mask area
%
% example, 
%   I=imread('sampleblue.tif');    %  greyscale image
%   load('samplemask.mat');         % import  ROI mask
%   options = struct('Lwin',500 , 'snuc',5,'threshpix',300);
%   [noutline,nuclayer,nuc1_calc,nucstat]=process_nuc(I,samplemask,1,1,options,'b');
%   figure, imshow(noutline);
%   or refer to Ex1_ImAnalysis.m

% Function is written by YikTungTracy Ling, Johns Hopkins University (July 2019)
% Reference: Ling, Y. T. T., Pease, M. E., Jefferys, J. L., Kimball, E. C., Quigley, H. A., 
% & Nguyen, T. D. (2020). Pressure-Induced Changes in Astrocyte GFAP, Actin, and Nuclear 
% Morphology in Mouse Optic Nerve. Investigative Ophthalmology & Visual Science, 61(11), 14-14.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
defaultoptions = struct('Lwin',500 , 'snuc',4,'threshpix',300);

if (~exist('options','var')) 
    options=defaultoptions; 
else
    tags = fieldnames(defaultoptions);
    for i=1:length(tags)
         if (~isfield(options,tags{i}))
             options.(tags{i})=defaultoptions.(tags{i}); 
         end
    end
    if  (length(tags)~=length(fieldnames(options)))
        error('process_nuc: unknown option','unknown options found');
    end
end
%%
    %%%%%%%%%%%%%%%%%%% median filter %%%%%%%%%%%%%%%%%%%%%%%%%%
    medfiltnuc=medfilt2(nucI);
    %%%%%%%%%%%%%%% contrast enhancement %%%%%%%%%%%%%%%%%%%%%
    histnuc=adapthisteq(medfiltnuc,'NumTiles',[options.Lwin/2 options.Lwin/2]);
    showsingle(histnuc,color,showfig);
    %%%%%%%%%%%%%%%%%% otsu binarization %%%%%%%%%%%%%%%%%%%%
    threshnuc=graythresh(histnuc);
    nucB1= imbinarize(histnuc,threshnuc);
    nucB1=maskI.*nucB1;
    % figure 
    % imshow(nucB1)
    % title('nucB1')
    %%%%%%%%%%% morphological smoothing %%%%%%%%%%%%%%%%%%%%
    nucB1 = bwmorph(nucB1,'majority');
    nucB1 = bwmorph(nucB1,'clean');
    %nucB1 = imfill(nucB1,'holes');
    se4 = strel('diamond',options.snuc);
    nucB1open=imopen(nucB1,se4);
    nuclayer=imfill(nucB1open,'holes'); %for removal of nuclear pixels from actin pores
    % morphological opening remove region smaller than threshold number of
    % pixels
    nuc1_calc=bwareaopen(nuclayer,options.threshpix); 
    bound=bwmorph(nuc1_calc,'remove');

    clear nucB1 

    medfilter1=medfiltnuc;
    medfilter1(bound==1)=255;
    medfilter2=medfiltnuc;
    medfilter2(bound==1)=0;
    % clear nucden
    noutline(:,:,1)=medfilter1;
    noutline(:,:,2)= medfilter2;
    noutline(:,:,3)=medfilter2;

    NUClabel=bwlabel(nuc1_calc);
    if showfig==1
        figure
        imshow(nuc1_calc)
        figure
        imshow(noutline)
        labelrgb = label2rgb(NUClabel, 'jet', 'w', 'shuffle');
        figure
        imshow(labelrgb)
%         msgfig = msgbox('nuclear channel','Step3','modal');
%         uiwait(msgfig)
%         disp('Continue Image processing.');
%         close all
    end
    
    clear medfiltered1 medfiltered2 bound medfilternuc
    
    nucstats=regionprops(nuc1_calc,'Area','Perimeter','MajorAxisLength','MinorAxisLength');
    %nucstats2=regionprops(nucB1open,'PixelIdxList');

    nuccount=numel(nucstats);
    

    Nucarea = [nucstats.Area];
    Nucarea_total=sum(Nucarea);
    % figure
    % boxplot(Nucarea)
    %Nucperi = [nucstats.Perimeter];
    Nucmajor = [nucstats.MajorAxisLength];
    Nucminor = [nucstats.MinorAxisLength];
    nucstat.Nucarea_total=Nucarea_total;
    nucstat.Nucarea_med=median(Nucarea);
    nucstat.Nucarea_mean=mean(Nucarea);
    nucstat.Nucarea_std=std(Nucarea);

    clear nucstats2 nucstats NUClabel

    %circularity=Nucperi.^2./(4.*pi.*Nucarea);
    
    NucAR=(Nucmajor)./(Nucminor);
    nucstat.NucAR_med=median(NucAR);
    nucstat.NucAR_mean=mean(NucAR);
    nucstat.NucAR_std=std(NucAR);
    nucstat.nuccount=nuccount;
    mask_area=sum(sum(maskI));
    nucstat.Nucden=nuccount/mask_area;
    nucstat.NucAREAden=Nucarea_total/mask_area;

    clear NucAR Nucmajor Nucminor NUClabel1 NUClabel nucB4open_forcalc
end