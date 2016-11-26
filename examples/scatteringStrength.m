% A demo of what happens to the PSF seen through a telescope as the
% strength of a phase screen is increased.  
% This demo just reuses the same phase screen but with an increasing rms
% wavefront error.
% 
% 20150225 JLCodona

lambda = AOField.RBAND; % Red light.
r0 = 0.15; % r0 is 15 cm at 500 nm.

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
colormap(gray);

% A.show;
% input 'Press ENTER to Continue...'
% pause(3);

%% Make a Kolmogorov phase screen.
TURBULENCE = AOScreen(2048); % Make it big so we get good low-frequency behavior.
%TURBULENCE.lambdaRef = AOField.VBAND; %This is the default.

TURBULENCE.spacing(.01); % This is ridiculously small for demo purposes.
TURBULENCE.setR0(r0); 
TURBULENCE.make;

OPL = TURBULENCE.grid; % Keep this for scaling.

%% Make an AOField object.

F = AOField(A);
F.resize(1024); % make it big to study the field before the pupil.
F.FFTSize = 1024; % Used to compute PSFs, etc.
F.lambda = lambda;

F.planewave*A;
% F.show;

% input 'Continue...'

% This adds a reference wave to the field and computes the intensity.
% imagesc(x,y,F.interferometer(1),[0 3]);
% sqar;
% axis xy;
% drawnow;
% input 'Continue...'

THld = F.lambda/D * 206265; % Lambda/D in arcsecs.

F.planewave*A;
[PSF,thx,thy] = F.mkPSF(3,THld/5);

PSFmax = max(PSF(:));

for q=0.01:0.01:3

    r0q = TURBULENCE.r0/q^(6/5); % This is found by scaling the Kolmogorov structure function.
    fprintf('r0 = %.3f m\n',r0q);

    TURBULENCE.grid(OPL*q);
    
    F.planewave*TURBULENCE*A;
    
    %subplot(N1,N2,4);
    [PSF,thx,thy] = F.mkPSF(3,THld/5);
    imagesc(thx,thy,log10(PSF/PSFmax),[-4 0]);
    daspect([1 1 1]);
    axis xy;
    colorbar;

    title(sprintf('r_0=%.3f m,  D/r_0=%.1f\n',r0q,D/r0q));
    
    
    drawnow; 
end
