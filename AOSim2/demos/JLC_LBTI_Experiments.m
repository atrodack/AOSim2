dTHETA = 0.1;
RING = 1.5;
AR = 0.25;
dAZ = pi/32;
FOV = 2;
% r0 = 0.2;

NXPIX = 3;
NYPIX = 2;

% make_the_LBTI_AO_jlc
% make_the_MMT_AO_jlc
% PS = AOScreen(round([10 64]/.04));
PS = AOScreen(round(64*[1 1]/.04));
PS.name = 'Just 1 Phase Screen';
% PS.setR0(200);
PS.setR0(r0);
PS.make;

% F = AOField(PS);

[x,y] = A.coords;

A1Center = [0 -7.2085];
A2Center = [0  7.2085];

CCD = 0;
ny = 0;
% for ny=-10:10
%     for nx=-10:10
%
%         THETA = [nx ny]*0.1;

THETAS = [0 0];
BRIGHTS = 1;

for AZ=0:dAZ:(2*pi-1e-3)
    THETA = [cos(AZ) AR*sin(AZ)]*RING;
    THETAS(end+1,:) = THETA;
    BRIGHTS(end+1) = 1e-3;
end

% for AZ=0:dAZ:(2*pi-1e-3)
%     THETA = [cos(AZ) AR*sin(AZ)]*RING;

for n=1:length(BRIGHTS)
    THETA = THETAS(n,:);
    AMP = sqrt(BRIGHTS(n));
    A.Offset = A1Center;
    F1 = AOField(A);
    F1.lambda = AOField.MBAND;
    % F1.planewave(1,THETA)*A;
    F1.planewave(AMP,THETA)*A*PS;
    [x1_1,x1_2] = F1.coords;
    [X1,Y1] = F1.COORDS;
    
    A.Offset = A2Center;
    F2 = AOField(A);
    F2.lambda = AOField.MBAND;
    % F2.planewave(1,THETA)*A;
    F2.planewave(AMP,THETA)*A*PS;
    [x2_1,x2_2] = F2.coords;
    [X2,Y2] = F2.COORDS;
    
    clf;
    subplot(NYPIX,NXPIX,1);
    % show(F1*PS);
    show(F1);
    title('Left');
    subplot(NYPIX,NXPIX,2);
    % show(F2*PS);
    show(F2);
    title('Right');
    %         drawnow;
    
    % subplot(2,2,[3 4]);
    subplot(NYPIX,NXPIX,3);
    % F.zero + F1 + F2;
    F2.grid(-F2.grid_);
    F1.Offset = [0 0];
    F2.Offset = [0 0];
    
    % F = AOField(PS);
    F = AOField(F1);
    F.lambda = AOField.MBAND;
    % F.FFTSize = 2048;
    F.FFTSize = 1024;
    
    F.zero + F1 + F2;
    F.show;
    title('Bracewell');
    subplot(NYPIX,NXPIX,6);
    imagesc(x,y,F.mag2,[0 4]*BRIGHTS(n));
    sqar;
    axis xy;
    title('Bracewell Intensity');
   
    [kx,ky] = F.kcoords;
    [PSF,thx,thy] = F.mkPSF(FOV,FOV/100);
    CCD = CCD + PSF;
    subplot(NYPIX,NXPIX,4);
    
    %imagesc(thx,thy,log10(normalize(PSF)),[-3 0]);
    imagesc(thx,thy,log10(normalize(CCD)),[-3 0]);
%     imagesc(thx,thy,CCD);
%     colorbar;
    sqar;axis xy;
    title('CCD');
    drawnow;
end
% end

HEADER = struct();
HEADER.r0 = r0;

fits_write_image(sprintf('LBTI_RING_r0_%g.fits',r0),CCD');
