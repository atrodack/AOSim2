function a = addNoise(a,sigma)
  
% ADDNOISE: Add complex Gaussian noise to an AOGrid. 
%
% Syntax: AOG = addnoise(AOG,sigma);
%
% Written by: Johanan L. Codona, Steward Observatory: CAAO
% July 11, 2002
% 20090408 JLCodona: AOSim2.
  
  % Just go with the MATLAB karma!
  dims = size(a);
  a.grid_ = a.grid_ + sigma*(randn(dims) + 1i*randn(dims));
  
