function [kfit,mean_orient]=fitvonmises_cdf(branch_orientation,branch_distance,showfig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This function calculates an orientation dispersion parameter kfit by
%  fitting a the semicircular von Mises probability
%  density function to the histogram of orientation angles centered 
%  about the circular average orientation. This function is different from
%  network_kap because it used individual beam orientation instead of
%  fourier transform for historgram of beam orientations, used for a
%  different paper
% 
% [kfit,mean_orient]=fitvonmises_cdf(branch_orientation,branch_distance);
% 
% inputs,
%   branch_orientation:  branch_orientation of each labeled beams
%   branch_distance:     lenght of each labeled beams
% outputs,
%   kfit:           fitted dispersion parameter from semi-circular
%                   von-mises distribution, 
%                   kappa =0 for equally distributed orientations
%   mean_orient:    weighted circular mean of all orientations in radians

% example 1, 
%   data = zeros(100,100);
%   data([10:10:end],:) = 1;
%   figure
%   imshow(data)
%   BW=data;
%   SkelI = bwmorph(BW,'thin',Inf); % Skeletonize the beams
%   [beamNodeI,beamI,beamlabel]=MergeNodes(SkelI,8);
%   blstats = regionprops(beamlabel,'Perimeter','Orientation');
%   branch_orientation=[blstats.Orientation];
%   branch_distance==[blstats.Perimeter]./2;
%   [kfit,mean_orient]=fitvonmises_cdf(branch_orientation,branch_distance);
%   or refer to MeasOrient.m


% Function is written by YikTungTracy Ling, Johns Hopkins University (2018)
% Referece: Ling, Y.T.T., Shi, R., Midgett, D.E., Jefferys, J.L., Quigley, H.A. and Nguyen, T.D., 
% 2019. Characterizing the collagen network structure and pressure-induced strains of the human 
% lamina cribrosa. Investigative ophthalmology & visual science, 60(7), pp.2406-2422.

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if isempty(branch_orientation) == 1 || numel(branch_orientation)==1
    % return Nan if the region doesn't have branches
    mean_orient=NaN;
    kfit=NaN;
else
%mean_orient=sum(branch_orientation.*branch_distance)/sum(branch_distance);
mean_orient=circ_mean(branch_orientation,branch_distance,2);
% calculate orientation 
weight=branch_distance;
% sort in ascending order
[orient1,index] = sort(branch_orientation);
% rearrange weight to match the order
weight=weight(index);
% figure
% plot(orient1,weight)
% hold on
% % shift the distribution to central around average orientation
% orient1=orient1-mean_orient;
% plot(orient1,weight)
% hold off

bin_edge=[-pi/2:0.02:pi/2];
%shift the orient axis to mid point of the bin
bin_mid=[bin_edge+0.02/2];
bin_mid=bin_mid(1:numel(bin_edge)-1);
%count historgram
[N,~] = histcounts(orient1,bin_edge);
% sum the count to find the absolute start index for extractin the weight
weightextract=cumsum(N);

for e=1:numel(weightextract)
    if e==1
        if N(e)==0
            binheight(e)=0;
        else
            binheight(e)=sum(weight(1:weightextract(e)));
        end
    else
    if N(e)==0
        binheight(e)=binheight(e-1);
    else
        binheight(e)=sum(weight(1:weightextract(e)));
    end
    end
end
cpdf=binheight./binheight(end);
        
    

fun=@(k)(cumsum(2*exp(k.*cos(2*bin_mid-mean_orient))./(2*pi*besseli(0,k)))./max(cumsum(2*exp(k.*cos(2*bin_mid))./(2*pi*besseli(0,k)))))-cpdf;
%fun=@(k)(cumsum(2*exp(k.*cos(2*bin_mid)))./(2*pi*besseli(0,k))-vpdf);
kfit = lsqnonlin(fun,1);
%%%%%%%%%%%comment out when running in batch%%%%%%%%%%%%%%%%%%%%%%%
orient=[-pi/2:0.01:pi/2];
cpdf_fit=cumsum(2*exp(kfit.*cos(2*orient-mean_orient))./(2*pi*besseli(0,kfit)))./max(cumsum(2*exp(kfit.*cos(2*orient))./(2*pi*besseli(0,kfit))));
orient2=[-pi:0.01:pi];
pdf_fit=2*exp(kfit.*cos(2*orient2-mean_orient))./(2*pi*besseli(0,kfit));

if showfig==1
    figure
    % plot(bin_mid,vpdf)
    plot(bin_mid,cpdf)
    hold on
    plot(orient,cpdf_fit)
    title(['k=',num2str(kfit)])
    hold off

    % 
    figure
    polarplot(orient1,weight,'o')
    title(['k=',num2str(kfit)])
    figure
    polarplot(orient2,pdf_fit)
    title(['k=',num2str(kfit)])
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kfit = abs(kfit);
end
end
