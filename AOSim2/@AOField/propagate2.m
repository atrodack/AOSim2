function F = propagate2(F,Z,MaxTheta)
  
% AOField.propagate2(dz,MaxTheta)
% 
% dz: propagation distance in m.
% MaxTheta: max included angles in arcsecs.
% 
% 20150424 JLCodona

  NFRESNEL = 10;
  LATERAL = 1/4;

  if(nargin<3)
      MaxTheta = sqrt(F.lambda/Z) * NFRESNEL;
      MaxTheta = min(MaxTheta,min(F.extent)*LATERAL/Z)*206265;
  end
  
  if(nargin < 2)
    error('requires at least 2 arguments: propagate(AOField,distance).');
  end
  
  % pad the field by some number of Fresnel scales.
  Rf = sqrt(abs(Z)*F.lambda);
  %fprintf('DEBUG: The Fresnel scale of this jump is %g m.\n',Rf);
  
  SZ = F.size;
  
  n1 = 1:SZ(1);
  n1 = fftshift(n1);
  n1 = n1-n1(1);

  n2 = 1:SZ(2);
  n2 = fftshift(n2);
  n2 = n2-n2(1);
  
  DK = 2*pi./(F.size .* F.spacing);
  
  [K1,K2] = meshgrid(DK(1)*n1,DK(2)*n2);
  KR2 = K1.^2 + K2.^2;

  
  %[KX2,KX1] = F.KCOORDS;
  %KR2 = KX1.^2 + KX2.^2;
  % Note that this is corner-centered by construction.
  PROPAGATOR = exp(-1i*Z/2/F.k*KR2);
  
  MaxKappa = MaxTheta/206265*F.k;
%  MaxKappa = MaxKappa/2; % Be conservative.
  
  PROPAGATOR = PROPAGATOR .* smoothUP(1-KR2/MaxKappa^2 ,1);
  
  F.grid(ifftshift(ifft2(PROPAGATOR.*fft2(fftshift(F.grid))))); 
  F.z = F.z - Z;
  
