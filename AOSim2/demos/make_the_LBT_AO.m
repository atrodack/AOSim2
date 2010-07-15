%% How to build and program the AOSim2 LBT AO model.
% Johanan L. Codona, Steward Observatory, University of Arizona.
% 20090426 JLCodona: MMT and GMT First-light versions.
% 20090430 JLCodona: Started a special LBT design.  The LBT is different
% because it is like two completely separate and independent AO telescopes
% that have a beam combiner at the end.  This presents some unique
% challenges to model cleanly that differ from the multi-segment single
% logical pupil systems like the GMT.

%% Define any parameters here...

SUBAPS_ACROSS_SEGMENT = 20;

%% Start clean...
% close
% clear classes

%% Load in some definitions...

load data/PLBT_oneEye.mat

Seg = AOSegment;
Seg.name = 'LBT Primary';
Seg.pupils = PLBT;
Seg.make.show;
colormap(gray);drawnow;

% Seg.touch.make.show;
AL = AOAperture;
AL.name = 'LBT Left';
AL.addSegment(Seg);
AL.Offset = [0 -7.2085];
AL.show;

AR = AOAperture;
AR.name = 'LBT Right';
AR.addSegment(Seg);
AR.Offset = [0 7.2085];
AR.show;


%% Grab the actuator coordinates from my hexapolar DM design.  There are
% pix in the data/pix directory.  I have higher and lower actuator density
% options already designed as well.
load data/LBT_actuators_plusCenter.mat
LBT_Actuators(673,:) = [];
Rmax = sqrt(max(sum(LBT_Actuators.^2,2)));
D = PLBT(1,3);
Rseg = PLBT(1,3)/2;
ACT = LBT_Actuators/Rmax*(Rseg-0.01); % 1 cm in from edge. (TODO: Check schematics!)
% plot(LBT_Actuators(:,1),LBT_Actuators(:,2),'.'); 
% daspect([1 1 1]); 
% axis xy;

%% Make a DM with an OPD grid matched to the Aperture AL...
DML = AODM(AL);
DML.name = 'DM672L';
DML.centerOn(AL);

%% Add in the actuators.
DML.addActs(ACT,1,AL);

% Specify the boundary conditions... 
DML.defineBC(5,16); % A circle of 8 null points at 5m radius.

clf;
AL.show;
DML.plotRegions; 
daspect([1 1 1]); 
xlim([-1 1]*15);
ylim([-1 1]*6);
colorbar off;
drawnow;

%% Make a DM with an OPD grid matched to the Aperture AR...
DMR = AODM(AR);
DMR.name = 'DM672R';
DMR.centerOn(AR);

%% Add in the actuators.
DMR.addActs(ACT,1,AR);

% Specify the boundary conditions... 
DMR.defineBC(5,16); % A circle of 8 null points at 5m radius.

% clf;
hold on;
AR.show;
DMR.plotRegions; 
daspect([1 1 1]); 
xlim([-1 1]*15);
ylim([-1 1]*6);
colorbar off;
hold off;
drawnow;

%% Build the LEFT Shack-Hartmann WFS.
WFSL = AOWFS(AL,D/SUBAPS_ACROSS_SEGMENT);
WFSL.name = 'Left LBT SHWFS';
AL.show; WFSL.quiver(1); drawnow; % Show them.

%% Build the RIGHT Shack-Hartmann WFS.
WFSR = AOWFS(AR,D/SUBAPS_ACROSS_SEGMENT);
WFSR.name = 'Right LBT SHWFS';
AR.show; WFSR.quiver(1); drawnow; % Show them.

error('');

%% Now for some real work.  Building the RECONSTRUCTOR...
RECON = AOReconstructor(AL,DML,WFSL); % Both reconstructors will be the same.

% Now program this crazy thing.
% We can look at the singular values and choose things manually later.  
% For now, we will make default assumptions.  You can rebuild it quickly
% later.  The MMT currently runs with 56 modes corrected.

RECON.program(D,6*sqrt(2)); % OWD is ~8.5 lambda/D for programming.
semilogy(RECON.s/RECON.s(1));

FL = AOField(AL);
FL.lambda = RECON.lambda;
[x,y] = coords(FL);
% Look at the reconstructor for good "gut" feelings...
for n=1:size(RECON.RECONSTRUCTOR,2)
    DML.setActs(500*RECON.RECONSTRUCTOR(:,n)).touch.render;
    FL.planewave*AL*DML;
    imagesc(x,y,FL.interferometer(1),[0 3]);
    sqar;
    axis xy;
    drawnow;
end

RECON.rebuild(56);
RECON.Nmodes
DML.nActs

save data/JLC_LBTAO_Model AL AR DML DMR WFSL WFSR RECON D % Paranoid.

% Now run something like the canned_GMT script, but don't load anything in
% first.  
