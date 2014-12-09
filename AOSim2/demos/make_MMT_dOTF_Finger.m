% I am assuming that the AOSim2 data directory is here.
% A symlink is good enough.

LAMBDA = AOField.JBAND;

make_the_MMT_AO_jlc

thld = LAMBDA/D*206265;
FOV = 25 * thld;
PLATE_SCALE = thld/3;


FINGER = AOSegment(A)
FINGER.name = 'The Finger!';

[X,Y] = FINGER.COORDS;

% Just put a mask into the grid directly...
% FINGER.grid(~(Y<-2.5 & abs(X-1)<0.25));
FINGER.grid(exp(1i*pi/2*double(Y<-2.5 & abs(X-1)<0.25)));

F = AOField(A);
F.lambda = LAMBDA;
F.FFTSize = 1024;

N1 = 3; 
N2 = 2;

subplot(N1,N2,1);
F.planewave*A;
F.show
PSF0 = F.mkPSF(FOV,PLATE_SCALE);
OTF0 = ifftshift(fft2(fftshift(PSF0)));

subplot(N1,N2,2);
F*FINGER;
F.show
PSF1 = F.mkPSF(FOV,PLATE_SCALE);
OTF1 = ifftshift(fft2(fftshift(PSF1)));

subplot(N1,N2,N2+1);
imagesc([PSF0 PSF1].^.25);
daspect([1 1 1]);
axis off;

subplot(N1,N2,N2+2);
imagesc(abs([OTF0 OTF1]));
daspect([1 1 1]);
axis xy;
axis off;

subplot(N1,N2,2*N2+1);
imagesc(abs(OTF0-OTF1));
daspect([1 1 1]);
axis xy;
axis off;

subplot(N1,N2,2*N2+2);
imagesc(angle(OTF0-OTF1));
daspect([1 1 1]);
axis xy;
axis off;

