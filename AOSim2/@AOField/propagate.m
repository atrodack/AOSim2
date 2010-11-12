function F = propagate(F,dz,REGULARIZE)
  
% PROPAGATE: Propagate an AOField by the distance dz.
% 
% USAGE: F = propagate(F,dz,REGULARIZE)
% 
% INPUTS:
% F: The AOField to work on.
% dz: The distance to propagate (m).
% REGULARIZE: and optional number to filter high angles.
% 
% OUTPUTS:
% F: The modified AOField object.
% 
% This routine used the parabolic wave equation to perform the
% propagation.  Note that it is a free-space propagator.  To add
% medium effects, you must call interact(F,PS).
% The result is in the x domain.
%
% For more information on propagating paraxial waves through phase screens,
% see my dissertation, chap 1, intro theory.  Or Tatarskii, a lot of other
% places.  My thesis is online:
% http://leitzel.as.arizona.edu/thesis
% 
% SEE ALSO:
% AOGrid, interact.
% 
% Written by: Johanan L. Codona, Steward Observatory: CAAO
% July 13, 2002
% Nov 12, 2010 JLCodona.  New version where I keep F in x-space.



  if(nargin < 2)
    error('requires at least 2 arguments: propagate(AOField,distance).');
  end

  
  if(nargin<3)
      REGULARIZE = 0.01;
  end
  
  
  SZ = F.size;
  dK = F.dk;
  
  Rf = sqrt(dz*F.lambda);
%   fprintf('DEBUG: The Fresnel scale of this jump is %g m.\n',Rf);
  kx1 = mkXvec(SZ(1),dK(1));
  kx2 = mkXvec(SZ(2),dK(2));
  [KX2,KX1] = meshgrid(kx2,kx1);
 
  F.grid_ = ifft2(exp((-(REGULARIZE+1i)*dz/2/F.k)*(KX1.^2+KX2.^2))...
      .*fft2(F.grid)); 
  F.z = F.z - dz;
  
