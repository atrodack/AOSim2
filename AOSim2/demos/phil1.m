% iff needed...
% close ; clear classes
% load AOSim2/data/JLC_GMTAO_Model_AOSim2.mat
% RECON.rebuild(100)

GAIN = 1;
TTGAIN = 1;

% I called my d/17 WFS WFS17.  If you did not, change this line....
WFS=WFS17;  % Because I use "handles", this is an alias, not a copy.

gain=1; % gain>2 is asking for trouble!
GAMMA = 2;  % This is the gamma correction for the PSF image.
SCIENCE_WAVELENGTH = AOField.MBAND;

FOV_START = 2;
FOV_AO_ON = 0.2;

ZOOM_STARTTIME = 0.1;
ZOOM_ENDTIME = 0.15;

AO_STARTTIME = 0.2;

%% Create the WFS and Science AOField objects...
Fwfs = AOField(A);
Fwfs.lambda = RECON.lambda;  % The Reconstructor was calibrated at a certain wavelength.

F = AOField(A);
F.lambda = SCIENCE_WAVELENGTH;  % Pick your favorite Science Wavelength at the top.
F.FFTSize = 2048*[1 1]; % This needs to be HUGE for the GMT.
PSF = F.mkPSF(.2,.2/100); %TODO: This is a bug workaround.  FIXME!

A.trueUp;
DM.setActs(0).touch;

%put in 10 Hz wobble on segment one amplitude+/-0.1"
% n=2;
% A.setTipTilt(1,[5e-7*(sin(n/100*2*3.14)),0]).show;
A.setTipTilts(randn(7,2)*5e-7);

%% This is the guts of the AO closed-loop integrating servo....

for nn=1:100
	fprintf('%d\n',nn);
	
	if(nn<10)
		gain = 0;
	else
		gain = GAIN;
	end
	
	clf;
	colormap(gray);
	subplot(2,2,1);
	WFS.sense(Fwfs.planewave*A*DM).quiver;
	
	
	% if(ATMO.time>AO_STARTTIME)  % Suffer with seeing limit for 50ms.
	DM.bumpActs(-gain*RECON.RECONSTRUCTOR * WFS.slopes);
	DM.removeMean;
	
	% This is essentially a Tip-Tilt offloader from the DM to A...
	A.bumpTipTilts(TTGAIN*DM.estimateSegmentTT(1:7));
% 	DMTT = DM.estimateSegmentTT(1:7);
% 	A.bumpTipTilts(-0.05*DMTT);
	
	F.planewave*A*DM;
	
	subplot(2,2,2);
	F.show;
	
	subplot(2,2,3);
	DM.show;
	
	subplot(2,2,4);
	% 	imagesc(F.phase,[-1 1]); sqar; axis xy; colorbar;
	A.show;
	
	drawnow;
end
