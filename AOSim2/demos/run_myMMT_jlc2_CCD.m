%% MMTSim.  
% Johanan L. Codona, Steward Observatory, University of Arizona.
% 20090426 JLCodona: First-light version.

% % Load in my pre-built MMT model.
% load data/MMTAO_Model_JLC_20090426
% load data/MMTAO_Model_withBadActs_JLC_20090427

WFS_FPS = 550;
EXPOSURE = 1/15;

% gain=1; % gain>2 is asking for trouble!
gain = 0.3; % gain>2 is asking for trouble!
GAMMA = 2;  % This is the gamma correction for the PSF image.
SCIENCE_WAVELENGTH = AOField.MBAND;

FOV_START = 1.5;
FOV_AO_ON = 1.5;

ZOOM_STARTTIME = 0.;
ZOOM_ENDTIME = 0.001;

AO_STARTTIME = 0.001;

%% Define the Atmosphere model and winds aloft.
ATMO = AOAtmo(A);

WFlow = AOScreen(1024,0.17,500e-9);
WFlow.name = 'Lower altitude turbulence';
WFhigh = AOScreen(2048,0.20,500e-9);
WFhigh.name = 'High altitude turbulence';

ATMO.addLayer(WFlow,1000);
ATMO.addLayer(WFhigh,8000);

ATMO.layers{1}.Wind = [3 1];
ATMO.layers{2}.Wind = [1 -1]*20;

r0 = ATMO.totalFriedScale;
th_scat = AOField.VBAND/r0*206265;

fprintf('The total r0 is %f cm.\n',100*ATMO.totalFriedScale);
fprintf('The seeing is %.2f arcsecs.\n',th_scat);

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

%% Create the WFS and Science AOField objects...
Fwfs = AOField(A);
Fwfs.lambda = RECON.lambda;  % The Reconstructor was calibrated at a certain wavelength.

F = AOField(A);
F.lambda = SCIENCE_WAVELENGTH;  % Pick your favorite Science Wavelength at the top.
F.FFTSize = 512*[1 1]; 
PSF = F.mkPSF(.2,.2/100); %TODO: This is a bug workaround.  FIXME!

%% Set your Primary and adaptive secondary initial conditions...
A.trueUp;
DM.setActs(0);
% DM.addRippleActs(.3*[1 1],500e-9,0);
% DM.addRippleActs(.5*[-1.2 .4],300e-9,pi/2);
% for n=1:5
%     DM.addRippleActs(.5*randn(1,2),600e-9*rand,2*pi*rand);
% end

% These touches should no longer be needed, but they don't hurt.
touch(DM);
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

N1=2;N2=3;  % Selects display geometry.

CCD = 0;
FRAMES_PER_EXP = round(WFS_FPS*EXPOSURE);
nframe = 1;

%% "The clock has started..."
for n=1:2000
	% 1kHz Frame rate.  Divide by 2000 to match CoDR.
    t = n/WFS_FPS;
    ATMO.time = t-1.5;
	
    %% This is the guts of the AO closed-loop integrating servo....
	
    ATMO.BEACON = GUIDE_STAR;
    WFS.sense(Fwfs.planewave*ATMO*A*DM);
    
    if(t>AO_STARTTIME)  % Suffer with seeing limit for 50ms.
        DM.bumpActs(-gain*RECON.RECONSTRUCTOR * WFS.slopes);
        DM.removeMean;
	end
    % That was it!
    	
    %% Meanwhile, in the Science Band...
    ATMO.BEACON = SCIENCE_OBJECT;
	F.planewave*ATMO*A*DM;
	
    %% Plot some interesting pictures...

    clf; % Don't screw around.  Just clear it.
    subplot(N1,N2,1);
    
	if(t < ZOOM_STARTTIME)
		FOV = FOV_START;
	elseif(t < ZOOM_ENDTIME)
		FOV = (FOV_AO_ON-FOV_START)*...
			(t-ZOOM_STARTTIME)/(ZOOM_ENDTIME-ZOOM_STARTTIME)...
			+FOV_START;  % In arcsecs.
	else
		    FOV = FOV_AO_ON;  % In arcsecs.
	end

	RNG = FOV * [-1 1];
	PSF = F.mkPSF(FOV,FOV/100);
	Ipeak = max(Ipeak,max(PSF(:)));
	
	CCD = CCD + PSF;
	
	%imagesc(RNG,RNG,(PSF/Ipeak).^(1/GAMMA));
	imagesc(RNG,RNG,CCD.^(1/GAMMA));
	
	
    daspect([1 1 1]);
    title(sprintf('PSF (\\lambda=%.2g microns, \\gamma=%g) t=%.3f',...
        F.lambda*1e6,GAMMA,t));
    xlabel('arcsecs');
    ylabel('arcsecs');
    
	subplot(N1,N2,2);
	A.show;
    colorbar off;
	WFS.quiver(1);
    %title('WFS Slopes (autoscaled)');
    title('WFS Slopes');
	
    subplot(N1,N2,3);
	
    imagesc(x,y,(F.interferometer(.75)),[0 3]);
    title('Science Band Interferometer');
	DM.plotActuators;
	
    subplot(N1,N2,4);
    xScale = linspace(-pi,pi,64);
   
    g=F.grid_;
	binDat = histc(angle(g(mask)),xScale);
	[vals,indx] = max(binDat);
	phase0 = xScale(indx);
	binDat = histc(angle(exp(-1i*phase0)*g(mask)),xScale);
	%bar(xScale,binDat);
	plot(xScale,binDat,'k.');
    HISTMAX = max(HISTMAX,max(binDat));
    ylim([0 1.1*HISTMAX]);  % Tweak this for your situation...
	title(sprintf('Phase Histogram: Correcting %d modes, gain=%.2f.',...
		RECON.Nmodes,gain));
	xlabel('pupil phase');
	ylabel('frequency');

	subplot(N1,N2,5:6);
	surf(x,y,DM.grid,A.grid,'LineStyle','none');
	zlim([-1 1]*3e-6);
	daspect([1 1 5e-6]);
	lt=light();
	set(lt,'Position',[-1 0 1]);
    % You may need this if you aren't saving the frames.
    drawnow;

    %% This saves the current picture as a JPEG.
	if(mod(n,FRAMES_PER_EXP)==1)
		
        CUBE(:,:,nframe) = CCD;
        fits_write_image('/tmp/MMT_Simjlc2.fits',CUBE(:,:,1:nframe));
		filename = sprintf('/tmp/FRAME_%04d.jpg',nframe); 
        nframe=nframe+1;
		rez = 160;
		
		resolution = sprintf('-r%d',rez);
		print(resolution,'-djpeg',filename);
		
        CCD = 0;
	end
	
	if(ATMO.time>2.0)
		break;
	end
end

%% Movie creation...
% Run this command after it is done to create the movie...
% mencoder "mf:///tmp/FRAME_*.jpg" -mf fps=10 -o MOVIE.avi -ovc lavc -lavcopts vcodec=wmv1
system('mencoder "mf:///tmp/FRAME_*.jpg" -mf fps=30 -o MOVIE_automake.avi -ovc lavc -lavcopts vcodec=wmv1');
