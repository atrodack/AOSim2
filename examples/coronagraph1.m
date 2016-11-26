% A demonstration of using AOSim2 to watch the evolution of various
% statistics beyond a Kolmogorov phase screen.
% Use the Kuiper 61" to view the field.
% 
% 20150221 JLCodona

lambda = AOField.RBAND; % Red light.
r0 = 0.15; % r0 is 15 cm at 500 nm.
PHOTONS_PER_EXPOSURE = 1e4; % Just a test.

% WIND_SHIFT = [1 5]; % wind motion in phase screen pixels per exposure.
WIND_SHIFT = [0 2]; % wind motion in phase screen pixels per exposure.

D = 1.54;
secondary = 14.5/100;

SPACING = 0.01;            % fine spacing makes nice pupil images but is really overkill.
aa = SPACING;              % for antialiasing.
% aa = 0.04;
spider = 0.0254;
% spider = 0.01;

% Scales
THld = lambda/D * 206265; % Lambda/D in arcsecs.
FOV = 5; % arcsecs
PLATE_SCALE = THld/3;

PUPIL_DEFN = [
   0 0 D         1 aa 0 0 0 0 0
   0 0 secondary 0 aa/2 0 0 0 0 0
   0 0 spider   -2 aa 4 0 D/1.9 0 0
   ];

% Since this demo only uses one AOSegment, I will not use the AOAperture wrapper.  
% I will only use AOSegment.  This is fine for simple pupils.

A = AOSegment;
A.spacing(SPACING);
A.name = 'Kuiper 61inch Primary';
A.pupils = PUPIL_DEFN;
A.make;

clf;
% colormap(gray);
colormap(jet);

% A.show;
% input 'Press ENTER to Continue...'
% pause(3);

%% Make a Kolmogorov phase screen.
TURBULENCE = AOScreen(2048); % Make it big so we get good low-frequency behavior.
%TURBULENCE.lambdaRef = AOField.VBAND; %This is the default.

TURBULENCE.spacing(.02);
TURBULENCE.setR0(r0); 
TURBULENCE.make;

%% Make an AOField object.

F = AOField(A);
F.FFTSize = 1024; % Used to compute PSFs, etc.
F.resize(F.FFTSize);
F.lambda = lambda;

F.planewave*A;

%% Walk through the coronagraph steps and make crude masks.

F.planewave*A; % Just go through the pupil.
F.grid(F.fft/F.nx); % Go to the focal plane.

FPMASK = AOSegment(F);
FPMASK.grid(exp(-normalize(F.mag2)/0.1) ); % This is pretty ad hoc.

F.grid(F.fft/F.nx); % Go to the Lyot pupil plane.

LYOT = AOSegment(F);
% LYOT.grid(F.real); % This is a good way to bootstrap a Lyot Stop.

% This is a better way...
LYOTSTOP_DEFN = [
   0 0 (D*0.8)         1 aa 0 0 0 0 0  % undersize the Lyot stop
   0 0 (secondary*1.1) 0 aa/2 0 0 0 0 0 % oversize the secondary
   0 0 spider         -2 aa 4 0 D/1.9 0 0
   ];

LYOT.pupils = LYOTSTOP_DEFN;
LYOT.make;

[PSF0,thx,thy] = F.mkPSF(FOV,PLATE_SCALE); % This is the reference PSF.
PSFmax = max(PSF0(:)); % Save for normalizing.

PSF0 = PSF0/PSFmax; % make the brightest value =1.

%% Now do it again with the masks to test the coronagraph...

% Do an on-axis star
%  F.planewave(1,[0 1.5])*A; % Just go through the pupil. With Source offset
F.planewave*A; % Just go through the pupil.
F.grid(F.fft/F.nx); % Go to the focal plane.
F*FPMASK; % Pass through the focal plane mask.
F.grid(F.fft/F.nx); % Go to the Lyot pupil plane.
F*LYOT; % Pass through the Lyot Stop.

[PSF,thx,thy] = F.mkPSF(FOV,PLATE_SCALE); % This is the reference PSF.
PSFstar = PSF/PSFmax;  % reference to the non-coronagraph PSF peak.

%% Do a planet offset by some angle

F.planewave(1,[1 0]*THld*4)*A; % Just go through the pupil. With Source offset
% F.planewave*A; % Just go through the pupil.
F.grid(F.fft/F.nx); % Go to the focal plane.
F*FPMASK; % Pass through the focal plane mask.
F.grid(F.fft/F.nx); % Go to the Lyot pupil plane.
F*LYOT; % Pass through the Lyot Stop.

[PSF,thx,thy] = F.mkPSF(FOV,PLATE_SCALE); % This is the reference PSF.
PSFplanet = PSF/PSFmax;  % reference to the non-coronagraph PSF peak.

colormap(jet);
imagesc(log10([PSF0 PSFstar PSFplanet]),[-6 0]);sqar;colorbar;



%% Now add turbulence (and FakeAO (tm))

TURBULENCE.setR0(0.15).make;
wavefront=TURBULENCE.grid();

PLANET_CONTRAST = 1e-2; 
PLANET_OFFSET = [1 1]/sqrt(2)*THld*4;

% colormap(flipud(gray));
colormap((hot));
CCD = 0;
IMG_LYOT = 0;  % for studying the average irradiance in the Lyot plane.

for n=1:1000
    modprint(n,20);
    
    TURBULENCE.grid(wavefront-circshift(wavefront,WIND_SHIFT));
    wavefront = circshift(wavefront,WIND_SHIFT);
    
    %% Now do it again with the masks to test the coronagraph...
    
    % Do an on-axis star
    %  F.planewave(1,[0 1.5])*A; % Just go through the pupil. With Source offset
    F.planewave*TURBULENCE*A; % Just go through the pupil.
    F.grid(F.fft/F.nx); % Go to the focal plane.
    F*FPMASK; % Pass through the focal plane mask.
    F.grid(F.fft/F.nx); % Go to the Lyot pupil plane.
    IMG_LYOT = IMG_LYOT + F.mag2;
    
    F*LYOT; % Pass through the Lyot Stop.
    
    [PSF,thx,thy] = F.mkPSF(FOV,PLATE_SCALE); % This is the reference PSF.
    PSFstar = PSF/PSFmax;  % reference to the non-coronagraph PSF peak.
    
    %% Do a planet offset by some angle
    
    F.planewave(1,PLANET_OFFSET)*TURBULENCE*A; % Just go through the pupil. With Source offset
    % F.planewave*A; % Just go through the pupil.
    F.grid(F.fft/F.nx); % Go to the focal plane.
    F*FPMASK; % Pass through the focal plane mask.
    F.grid(F.fft/F.nx); % Go to the Lyot pupil plane.
    F*LYOT; % Pass through the Lyot Stop.
    
    [PSF,thx,thy] = F.mkPSF(FOV,PLATE_SCALE); % This is the reference PSF.
    PSFplanet = PSF/PSFmax;  % reference to the non-coronagraph PSF peak.
    
    subplot(2,1,1);
    imagesc(log10([PSFstar PSFplanet]),[-6 0]);sqar;colorbar;
    axis off;
    
    subplot(2,1,2);
    
    %CCD = CCD + (PSFstar + PLANET_CONTRAST*PSFplanet);
    CCD = CCD + photonz(PSFstar + PLANET_CONTRAST*PSFplanet,PHOTONS_PER_EXPOSURE);
    imagesc(log10(normalize(CCD)),[-4 0]);sqar;colorbar;
    axis off;
    
    drawnow;
    
end

fitswrite(CCD,'coronagraph_image.fits');






















N1=2; N2=2; % subplots


