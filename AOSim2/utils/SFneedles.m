function [spacings,Dphi,rmsDphi] = SFneedles(WAVEFRONT,APER,NPOINTS)

% function [spacings,Dphi,rmsDphi] = SFneedles(WAVEFRONT,Aperture,[NPOINTS])
% This method estimates the structure function by dropping random points
% within a pupil and then computing statistics on the phase value
% differences.
% 
% PHI is an AOATMO or an AOScreen.
% 
% NOTE: This computes the wavefront displacement structure function!
% 
% JLCodona: 20091118
% JLCodona: 20110805 Tailored for better use with AOSim2.

if(nargin<4)
	NPOINTS = 1000;
end

% dx = X(1,2)-X(1,1);
% [x,y] = coords(PHI);
[x,y] = APER.coords;
% dx = A.dx;
PHI = WAVEFRONT.grid;
PUPIL = APER.grid;

% spacings = (0:floor(sqrt(2)*size(PHI,2)))*dx;
if(numel(x) > numel(y))
    spacings = x-x(1);
else
    spacings = y-y(1);
end

Dphi = zeros(size(spacings));
rmsDphi = zeros(size(spacings));
Nspacings = zeros(size(spacings));

figure(10);
hold off;
imagesc(PHI); sqar;
colormap(gray);
drawnow;

POINTS_ = getValidRandomPixels(NPOINTS,PUPIL); % pick valid pixels first

% figure(2);
for n=1:length(POINTS_)
% for n=1:3
% 	figure(10);
% 	hold off;
% 	imagesc(PUPIL.*PHI); sqar;
% 	hold on;
% % 	plot(POINTS_(:,2),POINTS_(:,1),'ro'); sqar;
% 	drawnow;

	
	NX1 = POINTS_(n,:);
	phi1 = PHI(NX1(2),NX1(1));
% 	plot(NX1(2),NX1(1),'k^','MarkerSize',15);
% 	fprintf('set %d of %d.\r',n,length(POINTS_));
	for m=n+1:length(POINTS_)
		NX2 = POINTS_(m,:);
% 		plot(NX2(2),NX2(1),'yo'); drawnow;
		
		phi2 = PHI(NX2(2),NX2(1));
		dp2 = (phi1-phi2)^2;
		dpixels = fix(round(norm(NX2-NX1)));
		bin = dpixels+1;
		Nspacings(bin) = Nspacings(bin) + 1;
		Dphi(bin) = Dphi(bin) + dp2;
		rmsDphi(bin) = rmsDphi(bin) + dp2^2;
		
% 		if(dp2>1000)
% 			fprintf('%d,%d: (%d,%d) (%d,%d) dist=%d  phi1=%.1f phi2=%.1f dphase^2=%.1f\n',...
% 				n,m,NX1(2),NX1(1),NX2(2),NX2(1),dpixels,phi1,phi2,dp2);
% 			plot(NX2(2),NX2(1),'rx','MarkerSize',15); 
% 			plot(NX2(2),NX2(1),'rs','MarkerSize',15);
% 			drawnow;
% 		end
	end
	% 	bar(Nspacings);
% 	clf;
% 	
% 	plot(PIXELS(:,2),PIXELS(:,1),'ro'); sqar;
% 	drawnow;
end

Dphi = Dphi./(Nspacings+eps);
rmsDphi = sqrt(rmsDphi./(Nspacings+eps) - Dphi.^2);% sqrt(mean(x^2)-mean(x)^2)

return
end

function PIXELS = filterPoints(PIXELS,PUPIL)
%     figure(10);
    %clf;
%     hold on;
%     plot(PIXELS(:,2),PIXELS(:,1),'y.','MarkerSize',2); sqar;
    
    OUTSIDE = logical(zeros(size(PIXELS,1)));
    
    for n=1:size(PIXELS,1)
        pixel = PIXELS(n,:);
        OUTSIDE(n) = (PUPIL(pixel(2),pixel(1)) < 0.5);
    end
    PIXELS(OUTSIDE,:) = [];
%     hold on;
%     plot(PIXELS(:,2),PIXELS(:,1),'ro'); sqar;
%     hold off;
%     drawnow;
    return
end


function PIXELS = getValidRandomPixels(N,PUPIL)
  % n.b. column 1 is X and column 2 is Y.
  %PIXELS = 1+fix(floor([rand(N,1)*size(PUPIL,2),rand(N,1)*size(PUPIL,1)]));
  PIXELS = 1+fix(floor([rand(N,1)*size(PUPIL,1),rand(N,1)*size(PUPIL,2)]));
  PIXELS = filterPoints(PIXELS,PUPIL);

return
end

