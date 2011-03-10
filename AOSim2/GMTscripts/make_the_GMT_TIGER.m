%% How to build and program the AOSim2 GMT AO model.
% Johanan L. Codona, Steward Observatory, University of Arizona.
% 20090424 JLCodona: First-light version.

%% Start clean...
% close
% clear classes

%% Parameters

NMODES = 500; % this can be reprogrammed later.

%% Load in some definitions...
load data/NewGMTPupil.mat

% This was old.  Nuke it.
delete(DM)
clear DM

%% Grab the actuator coordinates from my hexapolar DM design.  There are
% pix in the data/pix directory.  I have higher and lower actuator density
% options already designed as well.
load data/GMT_Actuators_JLC631_justCoords.mat Act*
colormap(gray);
A.show;

%% Make a DM with an OPD grid matched to the Aperture A...
DM = AODM(A);
DM.name = 'GMT Hexapolar Adaptive Secondary';

%% Add in the actuators and tag them with their segment positions.  
% This loops over th outer segments...
for n=1:6
    DM.addActs(Act,n,A);
end

% This adds the actuators for the center segment.
DM.addActs(Act7,7,A);

% Specify the boundary conditions... 
DM.defineBC(15,5*6); % A circle of 5*6=30 null points at 15m radius.
DM.plotRegions; daspect([1 1 1]); drawnow;


%% Build the Shack-Hartmann WFS.
% The CoDR says to make the subaps 8.4m/17 ~ 0.5m.  
% We do this and fill the whole pupil bounding box with subaps.

d = 8.4;
BB=A.BBox;
% D = mean(BB(2,:)-BB(1,:));
D = max(BB(2,:)-BB(1,:)); % this is ensure that the training sets cover the entire pupil.

WFS = AOWFS(A,d/17);
WFS.name = 'GMT Shack-Hartmann WFS';
A.show; WFS.quiver(1); drawnow; % Show them.


%% Now for some real work.  Building the RECONSTRUCTOR...
RECON = AOReconstructor(A,DM,WFS);

% The reconstructor has all of its physical parts now, but is is unprogrammed.
% I am paranoid, so I'll save it here first...
save TIGER_GMT_Model A DM WFS RECON d D % Paranoid.

% Now program the crazy thing.
% We can look at the singular values and choose things manually later.  
% For now, we will make default assumptions.  You can rebuild it quickly
% later.

% RECON = program(RECON,D,OWD,step,lambda)
% RECON.program(D,16,0.5,AOField.RBAND);

% RECON.OWD = 12;
% RECON.program(D,16);
 %RECON.program(D,10);
% RECON.OWD
% RECON.program(D,3); % OWD is 3 lambda/D which is VERY SMALL.  This is for debugging only.

OWD = sqrt(NMODES/4/pi)

% RECON.program(D,OWD); % OWD is 2.5 lambda/D which is VERY SMALL.  This is for debugging only.
% RECON.zprogram(D,10); % OWD is 2.5 lambda/D which is VERY SMALL.  This is for debugging only.
% RECON.zprogram(1.05*D,9); % Trying to avoid the Zernike edge effects.
% RECON.rebuild(NMODES);
% semilogy(RECON.s/RECON.s(1));

save GMTAO_DemoModel A DM WFS RECON d D % Paranoid.

F = AOField(A);
F.lambda = RECON.lambda;
[x,y] = coords(F);
% % Look at the reconstructor for good "gut" feelings...
% for n=1:size(RECON.RECONSTRUCTOR,2)
%     DM.setActs(500*RECON.RECONSTRUCTOR(:,n)).touch.render;
%     F.planewave*A*DM;
%     imagesc(x,y,F.interferometer(1),[0 3]);
%     sqar;
%     axis xy;
%     drawnow;
% end

RECON.Nmodes
DM.nActs

% Now run something like the canned_GMT script, but don't load anything in
% first.  
