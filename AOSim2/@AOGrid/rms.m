function value = rms(a)

% RMS: The AOGrid 'rms' function.
%
% usage: RMS = rms(a)
%
% Written by: Johanan L. Codona, Steward Observatory: CAAO
% Dec 1, 2005
% 20090407: JLCodona New version for AOSim2.

value = std(a.grid_(:));
