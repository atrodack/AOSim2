%% Kuiper 61inch example for AOSim2.
% Johanan L. Codona, Steward Observatory, University of Arizona.
% 20150217 JLCodona: Based this off of make_the_Omega61 in demos.

%% Some parameters for easy setting

D = 1.54;
secondary = 14.5/100;

SPACING = 0.01;
aa = SPACING;
% aa = 0.04;
% spider = 0.0254;
spider = 0.01;

PUPIL_DEFN = [
   0 0 D         1 aa 0 0 0 0 0
   0 0 secondary 0 aa/2 0 0 0 0 0
   %0 0 spider   -2 aa 4 0 D/1.9 0 0
   ];

% Since this demo only uses one AOSegment, I will not use the AOAperture wrapper.  
% I will only use AOSegment.  This is fine for simple pupils.

A = AOSegment;
A.spacing(SPACING);
A.name = 'Kuiper 61inch Primary';
A.pupils = PUPIL_DEFN;
A.make;

clf;
colormap(gray);

A.show;
input 'Continue...'

%% Make a lens to play with.

LENS = AOScreen(A);
LENS.name = 'Lens';

F = AOField(A);
F.lambda = AOField.VBAND;

LENS.zero.addZernike(2,0,-F.lambda/8,D);

[x,y] = F.coords;

F.planewave*A*LENS;
F.show;

input 'Continue...'

% This adds a reference wave to the field and computes the intensity.
% imagesc(x,y,F.interferometer(1),[0 3]);
% sqar;
% axis xy;
% drawnow;
% input 'Continue...'

F.FFTSize = 1024; % Used to compute PSFs, etc.
THld = F.lambda/D * 206265; % Lambda/D in arcsecs.

% 
F.planewave*A; % Just go through the pupil.
[PSF,thx,thy] = F.mkPSF(5,THld/4);
PSFmax = max(PSF(:)); % Save for normalizing.

PSF = PSF/PSFmax; % make the brightest value =1.

imagesc(thx,thy,log10(PSF),[-4 0]); 
axis square;
axis xy;
colorbar;

input 'Continue...'


%% Now pass through the lens and defocus...

for DEFOCUS=-2:.05:2
    
    % Set the lens defocus.  Use F.lambda as a reference.
    LENS.zero.addZernike(2,0,-F.lambda*DEFOCUS,D);
    LENS.addZernike(3,3,F.lambda,D);
    
    F.planewave*A*LENS;

    [PSF,thx,thy] = F.mkPSF(5,THld/4);
    PSF = PSF/PSFmax; % make the brightest value =1.

    imagesc(thx,thy,log10(PSF),[-4 0]);
    axis square;
    axis xy;
    colorbar;
    title(sprintf('%.2f lambda',DEFOCUS));
    
    drawnow; 
end


%% Now make a Kolmogorov phase screen and run through focus using it.
TURBULENCE = AOScreen(2048); % Make it big so we get good low-frequency behavior.
%TURBULENCE.lambdaRef = AOField.VBAND; %This is the default.

TURBULENCE.setR0(0.15); % A typical value on top of a mountain looking up.
TURBULENCE.make;

TURBULENCE.show;

input 'Continue...'

%% Now pass through the lens and defocus...

for DEFOCUS=-2:.1:2
    
    % Set the lens defocus.  Use F.lambda as a reference.
    LENS.zero.addZernike(2,0,-F.lambda*DEFOCUS,D);
    
    F.planewave*TURBULENCE*A*LENS;

    [PSF,thx,thy] = F.mkPSF(5,THld/4);
    PSF = PSF/PSFmax; % make the brightest value =1.

    imagesc(thx,thy,log10(PSF),[-4 0]);
    axis square;
    axis xy;
    colorbar;
    title(sprintf('%.2f lambda',DEFOCUS));
    
    drawnow; 
end

input 'Continue...'

%% Do it again with weaker turbulence...

TURBULENCE.setR0(0.5); 
TURBULENCE.make;

for DEFOCUS=-2:.1:2
    
    % Set the lens defocus.  Use F.lambda as a reference.
    LENS.zero.addZernike(2,0,-F.lambda*DEFOCUS,D);
    
    F.planewave*TURBULENCE*A*LENS;

    [PSF,thx,thy] = F.mkPSF(5,THld/4);
    PSF = PSF/PSFmax; % make the brightest value =1.

    imagesc(thx,thy,log10(PSF),[-4 0]);
    axis square;
    axis xy;
    colorbar;
    title(sprintf('%.2f lambda',DEFOCUS));
    
    drawnow; 
end



