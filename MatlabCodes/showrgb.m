function rgbI=showrgb(r,g,b,show)
%  This function showrgb displays a coloured image by combining 3 greyscale
%  images representing the red, green and blue channels
% 
% rgbI=showrgb(r,g,b,show);
% 
% inputs,
%   r,g,b :     The 2D greyscale or binary image for red green and blue
%               labeling, if the image is binary, input should be r*255;
%   show:       indicator, 1= show figures when running the program, 0=
%               don't show the figures
%
% outputs,
%   rgbI :      The coloured image with 3 channels combined (uint8)
%
% example, 
%   r=imread('samplered.tif');
%   g=imread('samplegreen.tif');
%   b=adapthisteq(imread('sampleblue.tif'));
%   rgbI=showrgb(r,g,b,1);
% Function is written by YikTungTracy Ling, Johns Hopkins University (July 2019)
if nargin < 4
  show = 1;
end


rgbI(:,:,1)=r;
rgbI(:,:,2)=g;
rgbI(:,:,3)=b;
rgbI=uint8(rgbI);
if show==1
    figure
    imshow(rgbI)
end
end

