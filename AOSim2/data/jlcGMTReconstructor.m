function ActToSlopes = jlcGMTReconstructor(D,Act,NsubAps)

% ActToSlopes = jlcGMTReconstructor(D,Act,NsubAps)
%
% Hoo boy.  JLC 20080624

Npixels = 16; % per subap;
CH = 1:Npixels;

Ninfl = NsubAps * Npixels; % even
dx = D/(Ninfl-1);
xx = -D/2 + (0:Ninfl-1)*dx;  % CENTERS of the pixels;
[Xi,Yi] = meshgrid(xx,xx);  % This is the influence grid.

nn = 1:NsubAps;
[Nsh,Msh] = meshgrid(nn,nn); % these are the x y indices of the subaps.

SH = [Nsh(:),Msh(:)];
clear Nsh Msh
Nsh = size(SH,1);
Nact = size(Act,1);
Nfft = Npixels*8;

POKE = zeros(size(Act,1),1);
INFL = griddata(Act(:,1),Act(:,2),POKE,Xi,Yi,'cubic');
MASK = double(~isnan(INFL));

INFL(isnan(INFL)) = 0;

INFL = sparse(INFL);

% Poke-to-Slopes

XSlopes = zeros(Nsh,Nact);
YSlopes = zeros(Nsh,Nact);

CENTROID0 = zeros(Nsh,2);

fprintf('Computing the reference psf locations...\n');
for nsh=1:Nsh
	fprintf('.');
	infl = INFL((SH(nsh,1)-1)*Npixels+CH,(SH(nsh,2)-1)*Npixels+CH);
	mask = MASK((SH(nsh,1)-1)*Npixels+CH,(SH(nsh,2)-1)*Npixels+CH);

	psf = abs(fftshift2d(fft2(mask.*exp(i*full(infl)),Nfft,Nfft))).^2;
	CENTROID0(nsh,:) = centroid(psf);
end
fprintf('\n');

fprintf('Computing the actuator responses...\n');
for nact=1:Nact
	fprintf('Poking actuator %d\n',nact);
	POKE = zeros(size(Act,1),1);
	POKE(nact) = 1;
	INFL = griddata(Act(:,1),Act(:,2),POKE,Xi,Yi,'cubic');
	INFL(isnan(INFL)) = 0;
	INFL = sparse(INFL);

	for nsh=1:Nsh
		infl = INFL((SH(nsh,1)-1)*Npixels+CH,(SH(nsh,2)-1)*Npixels+CH);
		mask = MASK((SH(nsh,1)-1)*Npixels+CH,(SH(nsh,2)-1)*Npixels+CH);

		if(mask(:)'*abs(infl(:)) ~= 0)
			fprintf('.');
			psf = abs(fftshift2d(fft2(mask.*exp(i*full(infl)),Nfft,Nfft))).^2;
			CENTROID = centroid(psf);
			SLOPE = CENTROID - CENTROID0(nsh);
			XSlopes(nsh,nact) = SLOPE(1);
			YSlopes(nsh,nact) = SLOPE(2);
		else
			XSlopes(nsh,nact) = 0;
			YSlopes(nsh,nact) = 0;
		end
	end
	fprintf('\n');

end

ActToSlopes = [XSlopes;YSlopes];

save JLCRECONWS Act ActToSlopes Xi Yi MASK SH

