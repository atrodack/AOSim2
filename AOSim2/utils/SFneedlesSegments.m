function [spacings,Dphi,rmsDphi] = SFneedlesSegments(PHI,PUPILS,dx,PISTONS,NPOINTS)

% function [spacings,Dphi,rmsDphi] = SFneedlesSegments(PHI,PUPILS,dx,PISTONS,NPOINTS)
%
% PUPILS(:,:,npupils)
% PISTONS(npupils): values to test by adding.
%
% JLC 20060118.

if(nargin<5)
	NPOINTS = 1000;
end

Npupils = size(PUPILS,3);
spacings = (0:floor(sqrt(2)*size(PHI,2)))*dx;
Dphi = zeros([length(spacings) Npupils Npupils]);
rmsDphi = zeros([length(spacings) Npupils Npupils]);

% figure(10);
% hold off;
% imagesc(PHI); sqar;
% drawnow;

points_name = sprintf('POINTS_%03d_%06d.mat',Npupils,NPOINTS);

if(exist(points_name) == 2)
	load(points_name);
else
	for n=1:Npupils
		POINTS_{:,:,n} = getValidRandomPixels(NPOINTS,PUPILS(:,:,n)); % pick valid pixels first
	end
	save(points_name,'POINTS_');
end

for segment1=1:Npupils
	P1 = POINTS_{:,:,segment1};
	for segment2=segment1:Npupils
			fprintf('Segments(%d,%d)\n',segment1,segment2);
		P2 = POINTS_{:,:,segment2};
		
		Nspacings = zeros(size(spacings));
		for n=1:size(P1,1)
			NX1 = P1(n,:);
			phi1 = PHI(NX1(1),NX1(2)) + PISTONS(segment1);
% 			plot(NX1(2),NX1(1),'k^','MarkerSize',15);
% 			fprintf('Segments(%d,%d): head %d of %d.\n',segment1,segment2,n,size(P1,1));
			for m=1:size(P2,1)
				NX2 = P2(m,:);
% 				plot(NX2(2),NX2(1),'yo'); drawnow;

				phi2 = PHI(NX2(1),NX2(2)) + PISTONS(segment2);
				dp2 = (phi1-phi2)^2;
				dpixels = fix(round(norm(NX2-NX1)));
				bin = dpixels+1;
				Nspacings(bin) = Nspacings(bin) + 1;
				Dphi(bin,segment1,segment2) = Dphi(bin,segment1,segment2) + dp2;
				rmsDphi(bin,segment1,segment2) = rmsDphi(bin,segment1,segment2) + dp2^2;
			end
		end

		goo = squeeze(Dphi(:,segment1,segment2))';
		Dphi(:,segment1,segment2) = goo./(Nspacings+eps);
		
		goo = squeeze(rmsDphi(:,segment1,segment2))';
		rmsDphi(:,segment1,segment2) = sqrt(goo./(Nspacings+eps) ...
			- squeeze(Dphi(:,segment1,segment2))'.^2);% sqrt(mean(x^2)-mean(x)^2)
	end
end


return
end

function PIXELS = filterPoints(PIXELS,PUPIL)
% figure(10);
%    clf;
% hold on;
% plot(PIXELS(:,2),PIXELS(:,1),'y.','MarkerSize',2); sqar;

OUTSIDE = false(size(PIXELS,1));

for n=1:size(PIXELS,1)
	pixel = PIXELS(n,:);
	OUTSIDE(n) = (PUPIL(pixel(1),pixel(2)) < 0.5);
end
PIXELS(OUTSIDE,:) = [];
% hold on;
% plot(PIXELS(:,2),PIXELS(:,1),'ro'); sqar;
% hold off;
% drawnow;
return
end


function PIXELS = getValidRandomPixels(N,PUPIL)
% n.b. column 1 is X and column 2 is Y.
PIXELS = 1+fix(floor([rand(N,1)*size(PUPIL,2),rand(N,1)*size(PUPIL,1)]));
PIXELS = filterPoints(PIXELS,PUPIL);

return
end

