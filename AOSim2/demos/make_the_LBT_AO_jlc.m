%% A prototype AOSim2 LBT AO model.
% This model uses two separate telescopes running independent AO systems.
% The design goal was to be efficient, so I actually use one AOAperture,
% one AOWFS, and one AODM.  The usage model is a bit different than we have
% seen in previous AO models, and all of the common elements get appied and
% operate at the common center.  This file just sets things up.  See, for
% example, JLC_LBT_Fizeau for an example of usage.  
% 
% This is an experiment. YMMV.
% 
% Johanan L. Codona, Steward Observatory, University of Arizona.
% 20101011 JLCodona: First version.

%% Start clean...
% close
% clear classes

%% Load in some definitions...

MAX_MODES = 56;

D = 8.4;
load data/PLBT_oneEye.mat
load data/LBT_actuators_plusCenter.mat
LBT_Actuators(end,:) = [];
LBT_Actuators = LBT_Actuators/12.84*D/2; % convert weird units to m.

ALCenter = [0 -7.2085];
ARCenter = [0  7.2085];
Center = [0 0];

% Start building the telescope...
Seg = AOSegment;
Seg.name = 'LBT Segment';
Seg.pupils = PLBT;
Seg.make;

clf;
% Seg.touch.make.show;
A = AOAperture;
A.name = 'LBT';
A.addSegment(Seg);
A.show;
colormap(gray);

%% Make our working fields.  These will be translated a lot using Offset.

Fwfs = AOField(A); % start centered on the origin.
Fwfs.name = 'WFS Field';
Fwfs.lambda = AOField.RBAND;

Fscience = AOField(A); % start centered on the origin.
Fscience.name = 'Science Field';
Fscience.lambda = AOField.MBAND;

Fcombined = AOField([20 20]./A.dx); % start centered on the origin.
Fcombined.name = 'LBT Combined Beam';
Fcombined + Fscience; Fcombined.show;

%% This model uses a common DM.  We will keep the actuators outside in a vector.
% Make a DM with an OPD grid matched to the Aperture A...
DM = AODM(A);
DM.name = 'LBT DM672';

%% Add in the actuators.
DM.addActs(LBT_Actuators,1,A);

% Specify the boundary conditions... 
DM.defineBC(8,16); % A circle of 16 null points at 8m radius.
DM.plotRegions; daspect([1 1 1]); drawnow;

%% Build the Shack-Hartmann WFS.
WFS = AOWFS(A,D/12);
WFS.name = 'LBT Shack-Hartmann WFS (common multiplexed)';
A.show; WFS.quiver(1); drawnow; % Show them.

%% Now for some real work.  Building the RECONSTRUCTOR...
RECON = AOReconstructor(A,DM,WFS);

%% Now program this crazy thing.

% RECON.adhocProgram(1);
% RECON.show;

% RECON.adhocProgram(1);
% RECON.show;

% OWD = sqrt(MAX_MODES/pi);
% RECON.program(D,6*sqrt(2)); % Use Fourier modes. OWD is ~6 lambda/D for programming.
% RECON.zprogram(D,12);  % program using Zernikes.
% RECON.rebuild(56).show;

