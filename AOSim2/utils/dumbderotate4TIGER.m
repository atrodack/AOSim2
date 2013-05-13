function [PSF_,ETAS] = derotate(PSF,CENTER,SWING)

% PSF_ = derotate(PSF,CENTER,ExpStart,ExpDur,DEC,LAT)
% 
% PSF: Image of the PSF 
% CENTER: Pixel location of the center of the PSF. (Can be non-integer)
% ExpStart: in HOURS re meridian crossing
% ExpDur: in HOURS
% DEC: Observation Declination (degrees)
% LAT: Observer Latitude (degrees)
% 
% ETAS: Optional retval rotation angles.
% 
% 20110310 JLCodona.  Steward Observatory

DEBUG = 1;
SCALE = 1;

[X1,X2,R] = mkImageCoords(PSF,SCALE,CENTER);

PSF_ = 0;

% LAMBDAS = SED(1,1):.01:SED(end,1);
% BRIGHTS = interp1(SED(:,1),SED(:,2),LAMBDAS);

% for n=1:length(LAMBDAS)
%     L = LAMBDAS(n);
%     dPSF_ = interp2(X2,X1,PSF,L*X2,L*X1);
%     dPSF_(isnan(dPSF_)) = 0;
%     PSF_ = PSF_ + dPSF_;
%     if(DEBUG)
%         imagesc(log10(normalize(PSF_)),[-4 0]);
%         daspect([1 1 1]);
%         colorbar;
%         drawnow;
%     end
% end

% lat = LAT/180*pi;
% dec = DEC/180*pi;
% 
% TIMES = ExpStart+(0:.05:ExpDur);
% PA = TIMES*pi/12;
% 
% OBJECT = cos(dec)*[cos(PA);sin(PA)];
% OBJECT(3,:) = sin(dec);
% 
% ZENITH = [cos(pi/2-lat);0;sin(pi/2-lat)];
% ETAS = zeros(1,size(OBJECT,2));
% 
% for n=1:size(OBJECT,2)
%     OBJ = OBJECT(:,n);
%     dOBJ = OBJ - (ZENITH'*OBJ)*ZENITH;
%     plot(dOBJ(2),dOBJ(3),'x'); hold on;
%     ETAS(n) = atan2(dOBJ(3),dOBJ(2));
%     
% end
% hold off;
% plot(TIMES,ETAS*180/pi+90,'.-');grid;drawnow;

% PSF = PSF_;
PSF_ = 0;
START = -SWING/2;
STOP = SWING/2;

ETAS = (START:1:STOP)*pi/180;
for n=1:length(ETAS)
    eta = ETAS(n);
    dPSF_ = interp2(X2,X1,PSF,cos(eta)*X2-sin(eta)*X1,cos(eta)*X1+sin(eta)*X2);
    dPSF_(isnan(dPSF_)) = 0;
    PSF_ = PSF_ + dPSF_;
    if(DEBUG)
        imagesc(log10(normalize(PSF_)),[-4 0]);
        daspect([1 1 1]);
        colorbar;
        drawnow;
    end
end

