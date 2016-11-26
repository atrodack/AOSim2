%% A model of the Quarter-Meter Telescope (Meade LX200)
% Defocus experiments.
% 
% JLCodona, 20150424.

LAMBDA = AOField.VBAND;

D = 10.0 * 0.0254;
dx = 0.001;
PUPIL = [ 0, 0, D,     1, dx,   0, 0, 0, 0, 0
    0, 0, 0.3*D, 0, dx/2, 0, 0, 0, 0, 0 ];
      
A = AOSegment;
A.name = 'QMT';
A.pupils = PUPIL;
A.spacing(0.001);
A.make;
A.show;
colormap(gray);

PS = AOScreen(2048);
PS.spacing(0.005);
% PS.setR0(0.15);
PS.setR0(0.03);
PS.make;

F = AOField(A);
F.lambda = AOField.VBAND;
F.FFTSize = 1024;
F.padBy(256);

F.planewave*PS*A;F.show

THld = F.lambda/D * 206265; % Lambda/D in arcsecs.
FOV = 25*THld; % arcsecs
PLATE_SCALE = THld/3;

R = 5000;

for Z=10:100:R
    F.sphericalWave(1,Z)*PS;
    F.propagate2(R-Z)*A;
    %plotComplex(F.grid,4);
    %drawnow;
    [PSF,thx,thy]=F.mkPSF(2*FOV,PLATE_SCALE);
    imagesc(thx,thy,PSF);sqar;
    title(sprintf('Screen Distance: %.0f',Z));
    drawnow;
end

Z = 1000;

for n=1:100
    PS.shiftPixels([0 1]);
    F.sphericalWave(1,Z)*PS;
    F.propagate2(R-Z)*A;
    %plotComplex(F.grid,4);
    %drawnow;
    [PSF,thx,thy]=F.mkPSF(2*FOV,PLATE_SCALE);
    imagesc(thx,thy,PSF);sqar;
    title(sprintf('Screen Distance: %.0f (%d)',Z,n));
    drawnow;
end

[x,y] = F.coords;
Z = 2000;
MAX_ANGLE = D/R*206265*1.5; % arcsecs

for Zback=[10:100:R,R]
    for n=1:4
        PS.shiftPixels([0 1]);
        F.sphericalWave(1,Z);
        F*PS;
        F.propagate2(R-Z,MAX_ANGLE)*A;
        %F.sphericalWave(1,Z)*PS*A;
        %F.propagate2(R-Z,MAX_ANGLE);
        
        F.propagate2(-Zback,MAX_ANGLE);
        
        imagesc(x,y,F.mag2);sqar;axis xy;
        title(sprintf('Focus Distance: %.0f (%.2f to screen, %.2f to src)',...
            Zback,(R-Zback)/Z,Zback/R));
        colorbar;
        drawnow;
    end
end
