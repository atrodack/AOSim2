function F = propagate(F,dz,REGULARIZE,PADDED)
  
% PROPAGATE: Propagate an AOField by the distance dz.
% 
% USAGE: F = propagate(F,dz,REGULARIZE)
% 
% INPUTS:
% F: The AOField to work on.
% dz: The distance to propagate (m).
% REGULARIZE: an optional number to filter high angles.
% PADDING: PAD The array to PADDED before propagating.
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
% http://books.google.com/books?id=59seqhaUvPwC
% 
% SEE ALSO:
% AOGrid, interact.
% 
% Written by: Johanan L. Codona, Steward Observatory: CAAO
% July 13, 2002
% Nov 12, 2010 JLCodona.  New version where I keep F in x-space.
% Nov 15, 2010 JLCodona. I added padding to the propagation.  Next is grid
% resampling.

  NFRESNEL = 25;

  if(nargin < 2)
    error('requires at least 2 arguments: propagate(AOField,distance).');
  end

  if(nargin<3)
      REGULARIZE = 0.01;
  end
  
  if(nargin<4)
      PADDED = 0;
  end
  
  SZ = F.size;
  DX = F.spacing;
  %dK = F.dk;
  
  field = F.grid;
  
  % pad the field by some number of Fresnel scales.
  Rf = sqrt(abs(dz)*F.lambda);
  %fprintf('DEBUG: The Fresnel scale of this jump is %g m.\n',Rf);
  
  if(nargin>4 | PADDED>F.nx)
      NPAD = round(NFRESNEL*Rf./DX);
      fprintf('DEBUG: PADDING the field array by [%d,%d] pixels.\n',NPAD);
      %field = padarray(field,NPAD,'post');
      %field = padarray(field,NPAD,F.subGrid(1,1),'post'); % pad with what is at the edge.
      field = [field,fliplr(field(:,end-NPAD+1:end))];
      field = [field;flipud(field(end-NPAD+1:end,:))];
  end
  
  dK = 2*pi./(size(field) .* DX);
  
  kx1 = mkXvec(size(field,1),dK(1));
  kx2 = mkXvec(size(field,2),dK(2));
  [KX2,KX1] = meshgrid(kx2,kx1);
  % make sure the regularization damps the high angles even in reverse.
  PROPAGATOR = exp((-(abs(REGULARIZE*dz)+1i*dz)/2/F.k)*(KX1.^2+KX2.^2));
  field = ifft2(PROPAGATOR.*fft2(field)); 
  F.grid_ = field(1:SZ(1),1:SZ(2)); 
  F.z = F.z - dz;
  
