%% Simple AOSim2 VLT model.
%  Johanan L. Codona, Steward Observatory, University of Arizona.
%  20120306

%% Definitions
SCIENCE_WAVELENGTH = 4.05e-6;
STREHL = 0.85;
D = 8.2;

PVLT = [ 0 0 8.2   1 0.05 0 0 0 0 0;
    0 0 1.066 0 0.02 0 0 0 0 0];
PVLT(:,3) = D*[1 0.14];

APPFILE = '/data/NewMath/VLT/bestOf_vlt_tests_20090915_203705_rmDD_APPS.fits'
APPMODES = 7;

%%

HEADER = struct;
HEADER.RUN = 'VLT Simulation'
HEADER.GOP = GOPdate;
HEADER.Strehl = STREHL;
HEADER.APPFILE = APPFILE;
HEADER.APPNUM = APPMODES;
HEADER.lambda = SCIENCE_WAVELENGTH;
%%

A = AOSegment;
A.name = 'VLT';
A.pupils = PVLT;
A.make.show;
colormap(gray);
drawnow;

APP = AOScreen(1);
APP.lambdaRef = SCIENCE_WAVELENGTH;

% APP.importAPP('data/MMT_20090102_130151_rmDD_APPS.fits',91);
% APP.importAPP('data/MMT_20090102_130151_rmDD_APPS_scaled.fits',91);

% APP.importAPP('/data/NewMath/BAPP/GMT/20110117_000443_GMT_TIGER_2_10ld_rmDD_APPS.fits',91);
APP.importAPP(APPFILE,APPMODES);
APP.spacing(.12*8.2/8.0);
APPMASK = AOSegment;
APPMASK.name = 'VLT APP MASK';
PVLTAPP = PVLT;
PVLTAPP(:,3) = D*[1 0.2];


APPMASK.pupils = PVLTAPP;
APPMASK.make;

F = AOField(A);
F.lambda = SCIENCE_WAVELENGTH;
F.FFTSize = 2048;

F.planewave*A*APP*APPMASK;F.show;

PSF = F.mkPSF(2,.01);

F = AOField(A);
F.lambda = SCIENCE_WAVELENGTH;
F.FFTSize = 2048;

RMSdZ = sqrt(-log(STREHL))*SCIENCE_WAVELENGTH/2/pi*1e9

PS = AOScreen(2048);
PS.spacing(APP.dx/2)
PS.make.show;
drawnow;
Z = PS.grid;

LPFx = 0.8:.4:2.5;
WFE = zeros(size(LPFx));

for n=1:length(LPFx)
    PS.grid(Z);
    PS.grid(Z-PS.LPF(LPFx(n)));
    WFE(n) = PS.rms*1e9;
    plot(LPFx,WFE,'o-'); grid; drawnow;
end

FILTER_SCALE = findZeros(LPFx,WFE-RMSdZ);
FILTER_SCALE = FILTER_SCALE(1)
% interp(LPFx,WFE-RMSdZ)

PS.grid(Z);
PS.grid(Z-PS.LPF(FILTER_SCALE));
PS.show;
drawnow;

clear Z;
SR = exp(-(PS.rms/SCIENCE_WAVELENGTH*2*pi)^2)

XTENT = PS.dx*PS.nx
dLOC = D*0.85;
NLOCS = round((XTENT-D)/dLOC);
LOCS = demean((1:NLOCS)*dLOC)

NCUBE = NLOCS^2

PSF = F.mkPSF(2,.02);
CUBE = single(zeros([size(PSF) NCUBE]));

say('Starting the direct imaging simulation.');

ncube = 0;
for n=1:NLOCS
    for m=1:NLOCS
        PS.Offset = [LOCS(n) LOCS(m)];
        F.planewave*PS*A;
        %F.show;drawnow;
        PSF = F.mkPSF(2,.02);
        imagesc(log10(normalize(PSF)),[-4 0]);sqar;drawnow;
        ncube = ncube + 1;
        CUBE(:,:,ncube) = PSF;
    end
    say(sprintf('%d percent.',round(100*n/NLOCS)));
end

say('Done. Saving the results');

fits_write_image(sprintf('VLTSim_direct_%s.fits',HEADER.GOP) ,CUBE,HEADER);
MEAN = mean(CUBE,3);
VAR = var(CUBE,0,3);
save(sprintf('VLTSim_direct_STATS_%s.mat',HEADER.GOP),'HEADER','MEAN','VAR');

say('Moving onto the A P P simulation.');

ncube = 0;
for n=1:NLOCS
    for m=1:NLOCS
        PS.Offset = [LOCS(n) LOCS(m)];
        F.planewave*PS*A*APPMASK*APP;
        %F.show;drawnow;
        PSF = F.mkPSF(2,.02);
        imagesc(log10(normalize(PSF)),[-4 0]);sqar;drawnow;
        ncube = ncube + 1;
        CUBE(:,:,ncube) = PSF;
    end
     say(sprintf('%d percent.',round(100*n/NLOCS)));
end

say('Done. Saving the results');

fits_write_image(sprintf('VLTSim_APP_%s.fits',HEADER.GOP) ,CUBE,HEADER);
MEAN = mean(CUBE,3);
VAR = var(CUBE,0,3);
save(sprintf('VLTSim_APP_STATS_%s.mat',HEADER.GOP),'HEADER','MEAN','VAR');

% LIMIT_APP=2.5*log10(sqrt(VAR(CH,CH))/max(MEAN(:))/2)-5.2;

% imagesc(2.5*log10(sqrt(VAR(CH,CH))/max(MEAN(:))/2)-5.2,[-14 -10]);colorbar
% bigtitle('VLT NACO APP Sensitivity in 1 hour 85% Strehl')
% saveJPEG('sensitivityVLTSim_1hour_SR85pct.jpg',150)
