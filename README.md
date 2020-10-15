# Morph_Analysis: Network_Fiber_Nucleus_Analysis

This code was originally developed to analyze astrocyte GFAP, actin and nucleus but can be generalized to conduct analysis on fiber network and cell nuclei.
This folder contains code and images for measuring network features such as fiber coverage, fiber width, fiber length, fiber tortuosity, fiber  orientation anisotropy , aspect ratio of nucleus, etc. It also has an auto boundary detection function for segmentation of region of interest(ROI). 

Ex1_main_ImAnalysis_3channels.m  , make sure the current directory is the folder containing matlab codes
% This code initiates the analyses for all images in a folder, assuming each image has 3 channels

% Thick/ Green channel is GFAP network of astrocyte cells

% Thin/Red channel is actin network of astrocyte cells

% Nucleus/Blue channel is obtained using DAPI to stain for cells nuclei

% This code saves the output images and output data 

Ex2_main_ImAnalysis_singlechannel.m

% This code initiates the analyses for all images in a folder, assuming each image has 1 channel only

% Reference: Ling, Y.T.T., Pease, M.E., Jefferys, J.L., Kimball, E.C., Quigley, H.A., & Nguyen, T.D. (2020). Pressure-Induced Changes in Astrocyte GFAP, Actin, and Nuclear  Morphology in Mouse Optic Nerve. Investigative Ophthalmology & Visual Science, 61(11), 14-14.

# Morph_Analysis
