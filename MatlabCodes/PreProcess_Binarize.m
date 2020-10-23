function [hist,B1]=PreProcess_Binarize(Im,Lwin,Swin,color,showfig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This function measures an average beam width within a masked area
% 
% [Avgwid,Stdwid,Maxwid]=MeasWidthAll(beamI,BWI,res,maskI,dilsize,showfig);
% 
% inputs,
%     Im:        The orginal 2D greyscale image 
%     showfig:      Optional, indicator, 1= show figures when running the 
%                   program, 0= don't show the figures
% outputs,
%   Avgwid:         Average beam width in the masked area
%   Stdwid:         Standard deviation of beam width in the masked area
%   Maxwid:         Maximum beam width
%
% example, 
%   Im=imread('samplegreen.tif'); % can be greyscale image
%   [hist,B1]=PreProcess_Binarize(Im,500,30,'g',1);
%   or refer to process_thin.m


% Function is written by YikTungTracy Ling, Johns Hopkins University (July 2019)
% Reference: Ling, Y. T. T., Pease, M. E., Jefferys, J. L., Kimball, E. C., Quigley, H. A., 
% & Nguyen, T. D. (2020). Pressure-Induced Changes in Astrocyte GFAP, Actin, and Nuclear 
% Morphology in Mouse Optic Nerve. Investigative Ophthalmology & Visual Science, 61(11), 14-14.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% start 
 %%%%%%%%%%%%%%%%%%% median filter %%%%%%%%%%%%%%%%%%%%%%%%%%
    medfiltI=medfilt2(Im);
    showsingle(medfiltI,color,showfig);
    %title('after median filter')
    %%%%%%%%%%%%%%% contrast enhancement %%%%%%%%%%%%%%%%%%%%%
    % overall contrast enhancement for lamina
    hist=adapthisteq(medfiltI,'NumTiles' ,[Lwin Lwin]);
    % local contrast enhancement
    hist=adapthisteq(hist,'NumTiles' ,[Swin Swin]);
    showsingle(hist,color,showfig);
    %title('after CLAHE')
    %%%%%%%%%%%%%%%%%% otsu binarization %%%%%%%%%%%%%%%%%%%%
    threshI=graythresh(hist);
    B1= imbinarize(hist,threshI);  
end