%% Canned GMTSim2 Demo routine.  
% This script does NOT rebuild the RECONSTRUCTOR or any of the parts.
%  For that, see "make_the_GMT_AO.m"
% 
% Johanan L. Codona, Steward Observatory, University of Arizona.
% 20090424 JLCodona: First-light version.

% Clean up from earlier.  May not be needed if you aren't editing the
% classes.

% fprintf('Comment out these lines if you want to use your own model...\n');
% close;
% clear classes;
% 
% % Load in my pre-built GMT model.
% load AOSim2/data/JLC_GMTAO_Model_AOSim2.mat
% 
% fprintf('Comment out to here...\n');



if(RECON.Nmodes > 200)
    beep;
    fprintf('WARNING: Your reconstructor probably has too many modes enabled.');
    fprintf('Run the command: "RECON.rebuild(100)" to remake it with 100 modes.');
    fprintf('Too many modes cause edge effects and sensitivity to noise, etc.');
    fprintf('You have been WARNED.');
    beep;
    % RECON.rebuild(100);
end

% I called my d/17 WFS WFS17.  If you did not, change this line....
WFS=WFS17;  % Because I use "handles", this is an alias, not a copy.

gain=1; % gain>2 is asking for trouble!
GAMMA = 2;  % This is the gamma correction for the PSF image.
SCIENCE_WAVELENGTH = AOField.MBAND;

FOV_START = 2;
FOV_AO_ON = 0.2;

ZOOM_STARTTIME = 0.1;
ZOOM_ENDTIME = 0.25;

AO_STARTTIME = 0.4;

%% Define the Atmosphere model and winds aloft.
ATMO = AOAtmo(A);

WFlow = AOScreen(1024,0.17,500e-9);
WFlow.name = 'Lower altitude turbulence';
WFhigh = AOScreen(2048,0.20,500e-9);
WFhigh.name = 'High altitude turbulence';

ATMO.addLayer(WFlow,1000);
ATMO.addLayer(WFhigh,8000);

ATMO.layers{1}.Wind = [5 0];
ATMO.layers{2}.Wind = [1 -1]*30;

% Turning this off is like using dynamic refocus.
ATMO.GEOMETRY = false;

%% Playing with boundary conditions...
r0 = ATMO.totalFriedScale
MaxActRadius = max(normRows(DM.actuators(:,1:2)));

DM.defineBC(MaxActRadius+3*r0,32);

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
HISTMAX = 0;

N1=2;N2=2;  % Selects display geometry.

%% "The clock has started..."
for n=1:2000
	% 1kHz Frame rate.  Divide by 2000 to match CoDR.
    ATMO.time = n/1000;
	
    %% This is the guts of the AO closed-loop integrating servo....
	
    ATMO.BEACON = GUIDE_STAR;
    WFS.sense(Fwfs.planewave*ATMO*A*DM);
    
    if(ATMO.time>AO_STARTTIME)  % Suffer with seeing limit for 50ms.
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
    
	if(ATMO.time < ZOOM_STARTTIME)
		FOV = FOV_START;
	elseif(ATMO.time < ZOOM_ENDTIME)
		FOV = (FOV_AO_ON-FOV_START)*...
			(ATMO.time-ZOOM_STARTTIME)/(ZOOM_ENDTIME-ZOOM_STARTTIME)...
			+FOV_START;  % In arcsecs.
	else
		    FOV = FOV_AO_ON;  % In arcsecs.
	end

	RNG = FOV * [-1 1];
	PSF = F.mkPSF(FOV,FOV/100);
	Ipeak = max(Ipeak,max(PSF(:)));
	imagesc(RNG,RNG,(PSF/Ipeak).^(1/GAMMA));
    daspect([1 1 1]);
    title(sprintf('PSF (\\lambda=%.2g microns, \\gamma=%g) t=%.3f',...
        F.lambda*1e6,GAMMA,ATMO.time));
    xlabel('arcsecs');
    ylabel('arcsecs');
    
	subplot(N1,N2,2);
	A.show;
    colorbar off;
	WFS.quiver(1);
    %title('WFS Slopes (autoscaled)');
    title('WFS Slopes');
	
    subplot(N1,N2,3);
    imagesc(sqrt(F.interferometer(1)),[0 2.5]);sqar;axis off;
    title('Science Band Interferometer');

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
	title('Science wavelength phase histogram');
	xlabel('pupil phase');
	ylabel('frequency');

    % You may need this if you aren't saving the frames.
    % drawnow;

    %% This saves the current picture as a JPEG.
    filename = sprintf('/tmp/FRAME_%04d.jpg',n);
    rez = 160;

    resolution = sprintf('-r%d',rez);
    print(resolution,'-djpeg',filename);
end

%% Movie creation...
% Run this command after it is done to create the movie...
% mencoder "mf:///tmp/FRAME_*.jpg" -mf fps=10 -o MOVIE.avi -ovc lavc -lavcopts vcodec=wmv1
system('mencoder "mf:///tmp/FRAME_*.jpg" -mf fps=30 -o MOVIE_MBAND.avi -ovc lavc -lavcopts vcodec=wmv1');
