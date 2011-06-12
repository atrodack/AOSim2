%% MMTSim.  
% Johanan L. Codona, Steward Observatory, University of Arizona.
% 20090426 JLCodona: First-light version.

% % Load in my pre-built MMT model.
% load data/MMTAO_Model_JLC_20090426
% load data/MMTAO_Model_withBadActs_JLC_20090427
% load data/MMTAO_Model_working

STROKE = 5.5e-6;

WFS_FPS = 500.;
WFS.qscale = 1; % smaller values make the WFS arrows bigger in the quiver plot.

gain=1; % gain>2 is asking for trouble!
% gain = 0.3; % gain>2 is asking for trouble!
GAMMA = 4;  % This is the gamma correction for the PSF image.
SCIENCE_WAVELENGTH = AOField.KBAND;

FOV_START = 1.5;
FOV_AO_ON = 1.5;

ZOOM_STARTTIME = 0.1;
ZOOM_ENDTIME = 0.2;

AO_STARTTIME = 0.1;
AO_STARTTIME = 0.010;

%% Define the Atmosphere model and winds aloft.
ATMO = AOAtmo(A);

WFlow = AOScreen(1024,0.15,500e-9);
WFlow.name = 'Lower altitude turbulence';
WFhigh = AOScreen(2048,0.17,500e-9);
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
F.FFTSize = 512*[1 1]; % This needs to be HUGE for the GMT.
PSF = F.mkPSF(2,.1); %TODO: This is a bug workaround.  FIXME!

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
HISTMAX = 200;
[x,y]=coords(F);

CCD = 0;

N1=2;N2=3;  % Selects display geometry.

%% "The clock has started..."
for n=1:1000
    if(mod(n,10)==0)
        fprintf('\n');
    end
    
    t = n/WFS_FPS;
    ATMO.time = t-1.5;
	
    fprintf('%d ',n);
    
    %% This is the guts of the AO closed-loop integrating servo....
	
    ATMO.BEACON = GUIDE_STAR;
    WFS.sense(Fwfs.planewave*ATMO*A*DM);
    
    if(t>AO_STARTTIME)  % Suffer with seeing limit for 50ms.
        DM.bumpActs(-gain*RECON.RECONSTRUCTOR * WFS.slopes);
        DM.removeMean.clip(STROKE);
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
% <<<<<<< .mine
    CCD = CCD + PSF;
	Ipeak = max(Ipeak,max(PSF(:)));
	imagesc(RNG,RNG,(PSF/Ipeak).^(1/GAMMA));
    daspect([1 1 1]);
    title(sprintf('PSF (\\lambda=%.2g microns, \\gamma=%g) t=%.3f',...
        F.lambda*1e6,GAMMA,t));
    xlabel('arcsecs');
    ylabel('arcsecs');
    axis xy;
%     
% 	subplot(N1,N2,2);
% 	A.show;
%     colorbar off;
% 	WFS.quiver(1);
%     %title('WFS Slopes (autoscaled)');
%     title('WFS Slopes');
% 	
%     subplot(N1,N2,3);
% 	
%     imagesc(x,y,(F.interferometer(.75)),[0 3]);
%     title('Science Band Interferometer');
% 	DM.plotActuators;
% 	
%     subplot(N1,N2,4);
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
% 	title(sprintf('Phase Histogram: Correcting %d modes, gain=%.2f.',...
% 		RECON.Nmodes,gain));
% 	xlabel('pupil phase');
% 	ylabel('frequency');
% 
% 	subplot(N1,N2,5:6);
% 	surf(x,y,DM.grid,A.grid,'LineStyle','none');
% 	zlim([-1 1]*3e-6);
% 	daspect([1 1 5e-6]);
% 	lt=light();
% 	set(lt,'Position',[-1 0 1]);
%     % You may need this if you aren't saving the frames.
%     drawnow;

    if(mod(n,20)==0)
        HEADER = struct;
        HEADER.DATA = GOPdate;
        HEADER.NFRAMES = n;
        
        fits_write_image('/tmp/CCD_Dump.fits',normalize(CCD),HEADER);
    end
% =======
	if(t>ZOOM_ENDTIME && t>AO_STARTTIME)
		CCD = CCD + PSF;
		if(mod(n,50)==0)
			fits_write_image('/tmp/PSF_Exposure.fits',CCD/max(CCD(:)));
		end
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
    if(isnan(RECON.Nmodes))
        title(sprintf('Phase Histogram: JLC ad hoc reconstructor, gain=%.2f.',gain));
    else
        title(sprintf('Phase Histogram: Correcting %d modes, gain=%.2f.',...
            RECON.Nmodes,gain));
    end
    
    xlabel('pupil phase');
	ylabel('frequency');

	subplot(N1,N2,5:6);
	surf(x,y,DM.grid,A.grid,'LineStyle','none');
	zlim([-1 1]*3e-6);
	%daspect([1 1 30e-6]);
	daspect([1 1 10e-6]);
	lt=light();
	set(lt,'Position',[-1 0 1]);
    % You may need this if you aren't saving the frames.
    % drawnow;

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
