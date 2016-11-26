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

  if(isempty(a.seed))
      a.grid_ = a.grid_ + (sigma/sqrt(2))*(randn(dims) + 1i*randn(dims));
  else
      rng(a.seed);
      a.grid_ = a.grid_ + (sigma/sqrt(2))*(randn(dims) + 1i*randn(dims));
      rng('default');
  end
