function F = propagate(F,dz)
  
% PROPAGATE: Propagate an AOField by the distance dz.
% 
% USAGE: F = propagate(F,dz)
% 
% INPUTS:
% F: The AOField to work on.
% dz: The distance to propagate (m).
%
% OUTPUTS:
% F: The modified AOField object.
% 
% This routine used the parabolic wave equation to perform the
% propagation.  Note that it is a free-space propagator.  To add
% medium effects, you must call interact(F,PS).
% The result is in the x domain.
%
% SEE ALSO:
% AOGrid, interact.
% 
% Written by: Johanan L. Codona, Steward Observatory: CAAO
% July 13, 2002

  global env;

  if(nargin < 2)
    error('requires 2 arguments.');
  end

  tic; F = transformK(F); kFT=toc;

  tic; [kx,ky] = coords(F); 
  [KX,KY] = meshgrid(kx,ky); Ksynth=toc;
  
%  dk = F.AOGrid.dk;
%  theta = kx/env.k*206265;
%  fprintf('THETA_X runs from %f to %f with step size %f (arcsec).\n', ...
%	  min(theta),max(theta),theta(2)-theta(1)); 
%  theta = ky/env.k*206265;
%  fprintf('THETA_Y runs from %f to %f with step size %f.\n', min(theta), ...
%	  max(theta),theta(2)-theta(1)); 
  
  tic; P = exp(-i*(KX.^2 + KY.^2)*dz/2/env.k); Psynth=toc;
  
  tic; F.AOGrid.grid = P .* F.AOGrid.grid; ftime=toc;
  
  tic; F = transformX(F); xFT=toc;

  F.z = F.z - dz;
  
  fprintf('propagate:times: k-FT=%f Ksynth=%f Psynth=%f filter=%f x-FT=%f\n', kFT,Ksynth,Psynth,ftime,xFT);
  
