function [kfit,mean_orient]=network_kap(BW,lowth,highth,showfig)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This function calculates an orientation dispersion parameter kfit by
%  fitting a the semicircular von Mises probability
%  density function to the histogram of orientation angles centered 
%  about the circular average orientation. Images should be larger than
%  100x 100 pixels
% 
% [kfit,mean_orient]=network_kap(BW,lowth,highth,showfig);
% 
% inputs,
%   BW :            The 2D logical/binary image, where pixels in
%                   region of interest = 1 and pixels in the backgroun =0,
%                   image should be larger than 100x100 pixels
%   lowth:          lower limit for band pass filter to remove large variations
%                   in background, to include everyting, this should be 0
%   highth:         higher limit for band pass filter to remove noise, must
%                   be smaller or equal to 0.5x size of the image to match
%                   the radius of the images, to include everything, this
%                   should be floor((min(size(BW))-1)/2)
%   showfig:        optional, indicator, 1= show figures when running the program, 0=
%                   don't show the figures
% outputs,
%   kfit:           fitted dispersion parameter from semi-circular
%                   von-mises distribution
%   mean_orient:    weighted circular mean of all orientations in radians
%
% additional inputs,
%   lowfreq:        lower limit for band pass filter to remove large variations
%                   in background
%   highfreq:       higher limit for band pass filter to remove noise, must
%                   be smaller or equal to 0.5x size of the image to match
%                   the radius of the images
%
% example 1, 
%   I=imread('samplegreen.tif'); % can be greyscale image
%   BW=imbinrize(I,graythresh(I)); % binarize the image
%   load('samplemask.mat');  % import  ROI mask
%   [kfit,mean_orient]=network_kap(BW.*samplemask,10,floor((min(size(BW))-1)/2),1);
%   or refer to runAnalysis.m

% example 2,
% % create image with horizontal or vertical lines
% data = zeros(100,100);
% data([5:5:end],:) = 1;
% % data(:,[5:5:end]) = 1;
% figure
% imshow(data)
% BW=data;
% [kfit,mean_orient]=network_kap(BW);

%  Function is written by YikTungTracy Ling, 
%  Johns Hopkins University (July 2019)
%  Reference: 
    % Band pass filter to eliminate high (noise) and low (large backgrounds) frequency components;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% parse input
    if nargin < 2
      showfig = 1;
      lowth=0;
      highth=floor((min(size(BW))-1)/2);
    end
    %% start function
    vert=size(BW,1);
    hori=size(BW,2);
    % make sure image is in odd number size 
    if mod(vert,2)==0
        vert=vert-1;
    end
    if mod(hori,2)==0
        hori=hori-1;
    end
    % crop image into a square, outer boundary is usually not within ROI
    if hori>vert
        startpt=(hori-vert)/2;
        cropBW=BW(1:vert,startpt+1:hori-startpt);
    elseif vert>hori
        startpt=(vert-hori)/2;
        cropBW=BW(startpt+1:vert-startpt,1:hori);
    else
        cropBW=BW(1:vert,1:vert);
    end

    %BWr = imrotate(BW,-90);
    BWf=fftshift(fft2(cropBW)); % fourier transform followed by shifting so that zero freq is at the center
    PS=abs(BWf).^2; % power spectrum

    if showfig==1
        figure
        imshow(cropBW)
        lnPS=log(1+(PS));
        figure
        imagesc((lnPS))
        title('power spectrum centered')
    end

    vert=size(cropBW,1);
    center=[((vert+1)/2),((vert+1)/2)];
    % conversion needed to be fliped up-down
    [xx,yy] = meshgrid(1:vert, 1:vert);

    % shift coordinates so that center is at the middle of the image
    coords = [yy(:)-center(1), xx(:)-center(2)];
    clear yy xx
    [theta, R] = cart2pol(coords(:,1),coords(:,2));
    clear coords
    % Convert to degrees and round them to the nearest degree
    degrees = mod(round(rad2deg(theta)), 360);
    dRange = 0:359;
    %close all

    % Average everything in same degrees (sum over PS and average by the number of pixels)
    for j = dRange
        PS_deg(j+1) = mean(PS(degrees == j & R > lowth & R < highth)); 
    end
        PS_deg(PS_deg==NaN)=0;
    if showfig==1
        % show 'histogram' of intensity in each degree
        figure
        polarplot(deg2rad(dRange),PS_deg)
        title('360 view')
    end

    % in our definition, 1 degree=181 degree, combine intensity in the same
    % degree
    PS_deg180=PS_deg(1:180)+PS_deg(181:360);
    clear PS_deg
    dRange = [0:179];
    dRange(dRange>90)=(dRange(dRange>90)-180);

    [orient1,index] = sort(dRange);
    PS_deg180=PS_deg180(index);

    orient1=deg2rad([-90,orient1]);
    PS_deg180=[PS_deg180(end),PS_deg180];

    if showfig==1
        % show 'histogram' of intensity in each degree from 0-179%
        figure
        polarplot(orient1,PS_deg180)
        title('180 view')
    end

    %mean_orient=sum(branch_orientation.*branch_distance)/sum(branch_distance);
    weight=PS_deg180';
    mean_orient=circ_mean(orient1',weight);
    % calculate orientation 


    binheight = cumsum(PS_deg180);
    cpdf=binheight./binheight(end);

    clear ps_integral180
    % fit to cumulative probability function of the semi-circular von-mises
    % distritbution
    fun=@(k)(cumsum(2*exp(k.*cos(2*orient1-mean_orient))./(2*pi*besseli(0,k)))./max(cumsum(2*exp(k.*cos(2*orient1))./(2*pi*besseli(0,k)))))-cpdf;
    %fun=@(k)(cumsum(2*exp(k.*cos(2*bin_mid)))./(2*pi*besseli(0,k))-vpdf);
    kfit = lsqnonlin(fun,1);

    if showfig==1
        cpdf_fit=cumsum(2*exp(kfit.*cos(2*orient1-mean_orient))./(2*pi*besseli(0,kfit)))./max(cumsum(2*exp(kfit.*cos(2*orient1))./(2*pi*besseli(0,kfit))));
        %pdf_fit=2*exp(kfit.*cos(2*orient2-mean_orient))./(2*pi*besseli(0,kfit));
        figure
        plot(orient1,cpdf)
        hold on
        plot(orient1,cpdf_fit)
        hold off
        legend('cpdf','fitted cpdf')
        set(gca,'FontSize',18)
        %title(['kfit=' num2str(abs(kfit))])
    end
    kfit = abs(kfit);
end



