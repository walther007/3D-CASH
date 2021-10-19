% Will perform Gerchberg-Saxton search for a hologramm that reconstructs
% a given set of target positions of a multi-spot focus pattern

% Contact: walther.akemann@ens.fr

% Following files are required in the active diretory:
% -------------------------------------------------------
% Gerchberg_Saxton_XY.m
% AOD_halo.m

% Following files required in the 'Functions' folder of the active directory:
% --------------------------------------------------------
% GerchbergSaxton.m
% GerchbergSaxton4.m
% phase2freq2.m
% --------------------------------------------------------

close all;
clear all;

% *************************************************************************

% INPUTS
% ------------- DEFINE TARGET PATTERN HERE ----------------

% Target pattern #1: 1x1 (single focus)
% Target pattern #2: 3x3 (8x8 um^2)
% Target pattern #3: 5x5 (16x16 um^2)
% Target pattern #4: 2x5 (24x24 um^2)

patternSELECT = 3;


% OUTPUTS
% --------------------------------------------------------
% Figure 1  1D hologram + RF for X-AOD
% Figure 2  1D hologram + RF for Y-AOD
% Figure 3  2D amplitude and phase hologram 
% Figure 4  Intensity distribution in object (Fourier) space
% Figure 5  Intesnity^2 distribution in object space
% Figure 6  3D Intensity distribution in object space
% --------------------------------------------------------
%
% *************************************************************************

Gerchberg_Saxton_XY;
AOD_holo;
