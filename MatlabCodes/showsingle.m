function colI=showsingle(I,color,show)

%  This function showrgb displays a coloured image by combining 3 greyscale
%  images representing the red, green and blue channels
% 
% rgbI=showrgb(r,g,b,show);
% 
% inputs,
%   I:          The 2D greyscale or binary image to be converted to colour
%               image if the image is binary, input should be I*255;
%   color:      must be 'r','g' or 'b' ,text for indicating whether output
%               image should be red green or blue in colour
%   show:       indicator, 1= show figures when running the program, 0=
%               don't show the figures
%
% outputs,
%   rgbI :      The coloured image with 3 channels combined (uint8)
%
% example, 
%   r=imread('samplered.tif');
%   redI=showsingle(r,'r',1);
% Function is written by YikTungTracy Ling, Johns Hopkins University (July 2019)
if nargin < 3
  show = 1;
end

[sizey,sizex]=size(I);
    I1=zeros(sizey,sizex);
    I2=I;
if color=='r'
    colI(:,:,1)=I2;
    colI(:,:,2)=I1;
    colI(:,:,3)=I1;
    colI=uint8(colI);
    if show==1
        figure
        imshow(colI)
    end

elseif color=='b'
    colI(:,:,1)=I1;
    colI(:,:,2)=I1;
    colI(:,:,3)=I2;
    colI=uint8(colI);
    if show==1
        figure
        imshow(colI)
    end
    
elseif color=='g'
    colI(:,:,1)=I1;
    colI(:,:,2)=I2;
    colI(:,:,3)=I1;
    colI=uint8(colI);
    if show==1
        figure
        imshow(colI)
    end
else
    error('input colour must be in text r,g or b ')

end

