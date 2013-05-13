%% GMTSim2.  
% Johanan L. Codona, Steward Observatory, University of Arizona.
% 20090426 JLCodona: First-light version.

% % Load in my pre-built MMT model.
% load data/MMTAO_Model_JLC_20090426
% load data/MMTAO_Model_withBadActs_JLC_20090427
% load data/MMTAO_Model_working

%WFS = WFS17; % local alias.

WFS_FPS = 1000;

% gain=1; % gain>2 is asking for trouble!
% gain = 0.3; % gain>2 is asking for trouble!
pgain = 0.1; % pgain is the piston gain.

GAMMA = 2;  % This is the gamma correction for the PSF image.
SCIENCE_WAVELENGTH = AOField.KBAND;

FOV_START = 0.2;
FOV_AO_ON = 0.2;

ZOOM_STARTTIME = 0.01;
ZOOM_ENDTIME = 0.02;

AO_STARTTIME = 0.005;

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
F.FFTSize = 2048*[1 1]; % This needs to be HUGE for the GMT.
PSF = F.mkPSF(.2,.2/100); %TODO: This is a bug workaround.  FIXME!

%% Set your Primary and adaptive secondary initial conditions...

A.lambdaRef = F.lambda; % for plotting only.

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

TIMES = (1:2000)/WFS_FPS;

ALL_PISTONS = zeros(7,length(TIMES));

% Strehl plot setup
mask = (A.grid>0.5);
STREHL = zeros(size(TIMES));
maxStrehl = 0.3;
N1=2;N2=3;  % Selects display geometry.

CCD = 0;

%% "The clock has started..."
for n=1:2000
    t = TIMES(n);
    ATMO.time = t-TIMES(1000);
	
    %% This is the guts of the AO closed-loop integrating servo....
	
    ATMO.BEACON = GUIDE_STAR;
    WFS.sense(Fwfs.planewave*ATMO*A*DM);
    
    if(t>AO_STARTTIME)  % Suffer with seeing limit for 50ms.
        DM.bumpActs(-gain*RECON.RECONSTRUCTOR * WFS.slopes);
        DM.removeMean;
		
% 		% This is NEW!
		PISTONS = WFS.magicPistonSensor(Fwfs,A);
		A.bumpPistons(- pgain*PISTONS);
		
	end
    % That was it!
    	
    %% Meanwhile, in the Science Band...
    ATMO.BEACON = SCIENCE_OBJECT;
	F.planewave*ATMO*A*DM;
	
	if(t>AO_STARTTIME)  % Suffer with seeing limit for 50ms.
		% This is NEW!
		PISTONS = WFS.magicPistonSensor(F,A);
		A.bumpPistons(- pgain*PISTONS);
		
		ALL_PISTONS(:,n) = PISTONS;
	end
	
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
	
	if(t>AO_STARTTIME) 
		CCD = CCD + PSF;
	end
	
	Ipeak = max(Ipeak,max(PSF(:)));
	imagesc(RNG,RNG,(PSF/Ipeak).^(1/GAMMA));
	axis xy;
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
	axis xy;
	
    subplot(N1,N2,3);
	
	[x,y]=F.coords;
    imagesc(x,y,(F.interferometer(.75)),[0 3]);
    title('Science Band Interferometer');
	daspect([1 1 1]);
	axis xy;
	% DM.plotActuators;
	
    %%
    subplot(N1,N2,4);
	g = F.grid_;
    STREHL(n) = abs(mean(g(mask)))^2;
	maxStrehl = min(1,max(maxStrehl,STREHL(n)*1.2));
    plot(TIMES(1:n),STREHL(1:n),'k-');
    ylim([0 maxStrehl]);
    xlim([0 TIMES(end)]);
    %xlim([-0.25 0]+t);
    ylabel('Strehl');
    xlabel('t (secs)');
	title(sprintf('Strehl History: Correcting %d modes, gain=%.2f.',...
		RECON.Nmodes,gain));
	clear g

% 	subplot(N1,N2,4);
%     xScale = linspace(-pi,pi,64);
%    
%     g=F.grid_;
% 	binDat = histc(angle(g(mask)),xScale);
% 	[vals,indx] = max(binDat);
% 	phase0 = xScale(indx);
% 	binDat = histc(angle(exp(-1i*phase0)*g(mask)),xScale);
% 	bar(xScale,binDat);
% 	plot(xScale,binDat,'k.');
%     HISTMAX = max(HISTMAX,max(binDat));
%     ylim([0 1.1*HISTMAX]);  % Tweak this for your situation...
% 	title(sprintf('Phase Histogram: Correcting %d modes, gain=%.2f.',...
% 		RECON.Nmodes,gain));
% 	xlabel('pupil phase');
% 	ylabel('frequency');

	subplot(N1,N2,5:6);
% 	surf(x,y,DM.grid,A.grid,'LineStyle','none');
% 	surf(x,y,DM.grid,'LineStyle','none');
% 	zlim([-1 1]*3e-6);
% 	daspect([1 1 5e-6]);
% 	lt=light();
% 	set(lt,'Position',[-1 0 1]);
    % You may need this if you aren't saving the frames.
    DM.show; colorbar;
    drawnow;

    %% This saves the current picture as a JPEG.
    filename = sprintf('/tmp/FRAME_%04d.jpg',n);
    rez = 160;

    resolution = sprintf('-r%d',rez);
    print(resolution,'-djpeg',filename);
	
	if(ATMO.time>2.0)
		break;
	end
end

%% Movie creation...
% Run this command after it is done to create the movie...
% mencoder "mf:///tmp/FRAME_*.jpg" -mf fps=10 -o MOVIE.avi -ovc lavc -lavcopts vcodec=wmv1
% system('mencoder "mf:///tmp/FRAME_*.jpg" -mf fps=30 -o MOVIE_automake.avi -ovc lavc -lavcopts vcodec=wmv1');
