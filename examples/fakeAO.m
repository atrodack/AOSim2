% A demonstration of using AOSim2 to watch the evolution of various
% statistics beyond a Kolmogorov phase screen.
% Use the Kuiper 61" to view the field.
% 
% 20150221 JLCodona

lambda = AOField.RBAND; % Red light.
r0 = 0.15; % r0 is 15 cm at 500 nm.
PHOTONS_PER_EXPOSURE = 1e4; % Just a test.

WIND_SHIFT = [1 5]; % wind motion in phase screen pixels per exposure.

D = 1.54;
secondary = 14.5/100;

SPACING = 0.01;            % fine spacing makes nice pupil images but is really overkill.
aa = SPACING;              % for antialiasing.
% aa = 0.04;
spider = 0.0254;
% spider = 0.01;

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
% colormap(jet);
colormap(hot);

% A.show;
% input 'Press ENTER to Continue...'
% pause(3);

%% Make a Kolmogorov phase screen.
TURBULENCE = AOScreen(2048); % Make it big so we get good low-frequency behavior.
%TURBULENCE.lambdaRef = AOField.VBAND; %This is the default.
TURBULENCE.name = 'Turbulent Layer';

TURBULENCE.spacing(.02);
TURBULENCE.setR0(r0); 
TURBULENCE.make;
[xx,yy] = TURBULENCE.coords;

CORRECTOR = AOScreen(TURBULENCE);
CORRECTOR.grid(-TURBULENCE.grid);  % This will be what we use to make the opposite effect.

%% Make an AOField object.

F = AOField(A);
F.resize(256); % do this to give a buffer around pupil for propagation.
F.FFTSize = 1024; % Used to compute PSFs, etc.
F.lambda = lambda;

[x,y] = F.coords;

F.planewave*A;

% Scales
THld = F.lambda/D * 206265; % Lambda/D in arcsecs.
FOV = 5; % arcsecs
PLATE_SCALE = THld/3;

% 
F.planewave*A; % Just go through the pupil.
[PSF0,thx,thy] = F.mkPSF(FOV,PLATE_SCALE);
PSFmax = max(PSF0(:)); % Save for normalizing.

PSF0 = PSF0/PSFmax; % make the brightest value =1.

N1=2; N2=2; % subplots

%% Plan of attack:
% We start with a plane wave, 
% go through a phase screen, 
% propagate by z, 
% pass through pupil, 
% correct the wavefront phase,
% compute the PSF.
% Add a number of them to simulate a time exposure.

% Do this at these ranges
% RANGES = [100 1e3 5e3 10e3 20e3 ];
RANGES = [500 1e3 5e3 10e3 20e3 50e3 ];

for z=RANGES
    fprintf('Range: %.1f m\n',z);
    % Assume the phase screen is in the z= plane.    

    TURBULENCE.make; % Make a new turbulent layer realization.
    
    clf; % Start with a fresh figure.

    subplot(N1,N2,1);
    imagesc(thx,thy,log10(PSF0),[-4 0]);
    daspect([1 1 1]);
    axis xy;
    colorbar;
    if(z<1000)
        title(sprintf('Ideal PSF z=%.1f m',z));
    else
        title(sprintf('Ideal PSF z=%.2f km',z/1000));
    end
    
    %
    CCD_noAO = 0; % Start the exposure.
    CCD_AO = 0; % Start the exposure.
    
    for t=1:1000  % t is just a counter here, not real time.
    %for t=1:100  % t is just a counter here, not real time.
        modprint(t,20);
        
        CORRECTOR.grid(-TURBULENCE.grid); % we want the corrector to be from the past by some lag.
        
        TURBULENCE.shiftPixels([0 1]); % simulate wind.
        
        F.planewave*TURBULENCE;
        F.propagate(z)*A;
    
        PSF = F.mkPSF(FOV,PLATE_SCALE);
        %CCD_noAO = CCD_noAO + PSF; % Add to the exposure.
        CCD_noAO = CCD_noAO + photonz(PSF,PHOTONS_PER_EXPOSURE); % Add photons to the exposure.
        
        F*CORRECTOR; % Fix up the phase (but not scintillation)
        PSF = F.mkPSF(FOV,PLATE_SCALE);
        %CCD_AO = CCD_AO + PSF; % Add to the exposure.
        CCD_AO = CCD_AO + photonz(PSF,PHOTONS_PER_EXPOSURE); % Add photons to the exposure.

        
        subplot(N1,N2,2);
        %imagesc(thx,thy,log10(PSF/PSFmax),[-4 0]);
        %daspect([1 1 1]);
        %axis xy;
        %colorbar;
        %title('Instantaneous PSF');
        %F.show;
        %title('Post-AO Field');
        wf = TURBULENCE.interpGrid(A);
        imagesc((x+WIND_SHIFT(1)*TURBULENCE.dx)/z*206265,(y+WIND_SHIFT(1)*TURBULENCE.dx)/z*206265,wf.*A.grid);sqar;
        colorbar;
        %setFoV(FOV);
        title('WF in angle');

        subplot(N1,N2,3);
        imagesc(thx,thy,log10(CCD_noAO/max(CCD_noAO(:))),[-4 0]);
        daspect([1 1 1]);
        axis xy;
        colorbar;
        
        title('No AO Correction');

        subplot(N1,N2,4);
        imagesc(thx,thy,log10(CCD_AO/max(CCD_AO(:))),[-4 0]);
        daspect([1 1 1]);
        axis xy;
        colorbar;
        
        title('With Fake AO (tm)');

        drawnow;
        
    end

    %     CUBE = [];
    %     CUBE(:,:,1) = CCD_noAO;
    %     CUBE(:,:,2) = CCD_AO;
    %     fitswrite(CUBE,sprintf('fake_AO_z%.0fm.fits',z));
    
    fitswrite([CCD_noAO CCD_AO],sprintf('fake_AO_z%.0fm.fits',z));
    
    %input 'Press ENTER to go on...'
    
end

