%% Steward Observatory White 21in Simulator. (Note: This telescope does not have adaptive optics.)  
% Johanan L. Codona, Steward Observatory, University of Arizona.
% 20090825 JLCodona: First-light version.

WFS_FPS = 1000;
GAMMA = 2;  % This is the gamma correction for the PSF image.
% SCIENCE_WAVELENGTH = AOField.RBAND;
SCIENCE_WAVELENGTH = 400e-9;

D = 21*0.0254;
lambda0 = 500e-9;

FOV = 2.5;
dFOV = lambda0/D*206265/2;

%% Define the Atmosphere model and winds aloft.
ATMO = AOAtmo(A);

WFlow = AOScreen(1024,0.10,500e-9);
WFlow.name = 'Lower altitude turbulence';
WFhigh = AOScreen(2048,0.15,500e-9);
WFhigh.name = 'High altitude turbulence';

ATMO.addLayer(WFlow,1000);
ATMO.addLayer(WFhigh,8000);

ATMO.layers{1}.Wind = [3 1];
ATMO.layers{2}.Wind = [1 -1]*20;

r0 = ATMO.totalFriedScale;
th_scat = AOField.VBAND/r0*206265;

fprintf('The total r0 is %.1f cm.\n',100*ATMO.totalFriedScale);
fprintf('The V-band seeing is %.2f arcsecs.\n',th_scat);

% Turning this off is like using dynamic refocus.
ATMO.GEOMETRY = false;

%% Guide star selection
SODIUM_LAYER = 90e3;

LGS_BEACON0 = [0 0 1] * SODIUM_LAYER;
LGS_BEACON = [0 1/206265 1] * SODIUM_LAYER;  % Offset by 1 arcsec.
STAR = [0 0 1e10];
LEO_TARGET = [0 0 400e3];
GEOSYNC_TARGET = [0 0 42e6];

% NGS CASE
GUIDE_STAR = STAR; % pick one.  
SCIENCE_OBJECT = STAR; % pick one.  

% Na LGS CASE
% GUIDE_STAR = LGS_BEACON; 
% SCIENCE_OBJECT = STAR; 

% Looking at GeoSynchronous satellites using LGS
% GUIDE_STAR = LGS_BEACON; 
% SCIENCE_OBJECT = GEOSYNC_TARGET; 

ATMO.BEACON = GUIDE_STAR; % Set this so ATMO knows how to compute the wavefront.

%% Create Science AOField objects...
F = AOField(A);
F.lambda = SCIENCE_WAVELENGTH;  % Pick your favorite Science Wavelength at the top.
F.FFTSize = 256*[1 1]; % This should be thought about for best performance.
PSF = F.mkPSF(FOV,dFOV); %TODO: This is a bug workaround.  FIXME!

%% Set your Primary and adaptive secondary initial conditions...
A.trueUp;
touch(A);

clf;
colormap(gray);  % Looks more official.

% This is to select the points that will be included in the phase
% histogram.
mask = (A.grid>0.5);

% This is the brightest pixel seen to date.
Ipeak = 0;  
HISTMAX = 2000;
[x,y]=coords(F);

N1=1;N2=2;  % Selects display geometry.

CCD = 0;
%% "The clock has started..."
for n=1:2000
	% 1kHz Frame rate.  Divide by 2000 to match CoDR.
    t = n/WFS_FPS;
    ATMO.time = t-1.5;
	
    %% Meanwhile, in the Science Band...
    ATMO.BEACON = SCIENCE_OBJECT;
	F.planewave*ATMO*A;
	
    %% Plot some interesting pictures...

    clf; % Don't screw around.  Just clear it.
    subplot(N1,N2,1);
    
	RNG = FOV * [-1 1];
	PSF = F.mkPSF(FOV,dFOV);
    
    CCD = CCD + PSF;
    if(mod(n,50)==0)
        fits_write_image('/tmp/SOWhite21inch_exposure.fits',CCD/max(CCD(:)));
    end
        
	Ipeak = max(Ipeak,max(PSF(:)));
	imagesc(RNG,RNG,(PSF/Ipeak).^(1/GAMMA));
    daspect([1 1 1]);
    title(sprintf('PSF (\\lambda=%.2g microns, \\gamma=%g) t=%.3f',...
        F.lambda*1e6,GAMMA,t));
    xlabel('arcsecs');
    ylabel('arcsecs');
    axis xy;
    
    subplot(N1,N2,2);
	
    imagesc(x,y,(F.interferometer(1.0)),[0 3]);
    title('Science Band Interferometer');
    daspect([1 1 1]);
    axis xy;
	
%     subplot(N1,N2,3);
%     xScale = linspace(-pi,pi,64);
%    
%     g=F.grid_;
% 	binDat = histc(angle(g(mask)),xScale);
% 	[vals,indx] = max(binDat);
% 	phase0 = xScale(indx);
% 	binDat = histc(angle(exp(-1i*phase0)*g(mask)),xScale);
% 	%bar(xScale,binDat);
% 	plot(xScale,binDat,'k.');
%     HISTMAX = max(HISTMAX,max(binDat));
%     ylim([0 1.1*HISTMAX]);  % Tweak this for your situation...
% 	title(sprintf('Phase Histogram'));
% 	xlabel('pupil phase');
% 	ylabel('frequency');
% 
    %% This saves the current picture as a JPEG.
    filename = sprintf('/tmp/FRAME_WHITE21in_%04d.jpg',n);
    rez = 120;

    resolution = sprintf('-r%d',rez);
    print(resolution,'-djpeg',filename);

    
    
end

