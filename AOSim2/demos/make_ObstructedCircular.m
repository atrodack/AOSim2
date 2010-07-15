%% How to build and program the AOSim2 MMT AO model.
% Johanan L. Codona, Steward Observatory, University of Arizona.
% 20090426 JLCodona: First-light version.
% 20090504 JLCodona: Had to add a Seg.make line.  This is a bug that needs
% to be found.  Sorry.
% 20090911 JLCodona: This is a generic pupil for the purposes of
% bootstrapping APP designs.

%% Start clean...
% close
% clear classes

%% Load in some definitions...

% This limits the reconstructor probe spatial frequencies during
% programming.  This is to compensate for my different mode ordering so
% that when we say "56 modes" we get them all, instead of tippish and
% tiltish being way up like modes 83 and 84 and therefore not getting
% included at all.  This is a tricky point and has to have been dealt with
% in the MMT reconstructor design.  All roads lead back to GUIDO BRUSA!
MAX_MODES = 56;

D = 8;
% obstruction = 0.2;
if(obstruction>0)
	PUPIL = [ 0 0 D 1 0.05 0 0 0 0 0
		0 0 obstruction*D 0 0.05 0 0 0 0 0 ];
else
	  PUPIL = [ 0 0 D 1 0.05 0 0 0 0 0 ];
end

% load data/PMMT.mat
% load data/MMT_DM336_Actuators.mat

Seg = AOSegment;
Seg.name = sprintf('Circular D=%g obs=%g',D,obstruction);
Seg.pupils = PUPIL;
Seg.make;

% Seg.touch.make.show;
A = AOAperture;
A.name = Seg.name;
A.addSegment(Seg);
A.show;
colormap(gray);

%% Grab the actuator coordinates from my hexapolar DM design.  There are
% pix in the data/pix directory.  I have higher and lower actuator density
% options already designed as well.
% load data/MMT_DM336_Actuators.mat ACT BAD

%% Make a DM with an OPD grid matched to the Aperture A...
% DM = AODM(A);
% DM.name = 'MMT DM336';

%% Add in the actuators.
% DM.addActs(ACT*6.5/2,1,A);

% Mark BAD actuators
% DM.disableActuators(BAD);
% DM.disableActuators(MMT_BADACTS_NO_CURRENT);
% DM.disableActuators(MMT_BADACTS_NO_POSITION);

% Specify the boundary conditions... 
% DM.defineBC(5,8); % A circle of 8 null points at 5m radius.
% DM.plotRegions; daspect([1 1 1]); drawnow;

%% Build the Shack-Hartmann WFS.
% WFS = AOWFS(A,D/12);
% WFS.name = 'MMT Shack-Hartmann WFS';
% A.show; WFS.quiver(1); drawnow; % Show them.

%% Now for some real work.  Building the RECONSTRUCTOR...
% RECON = AOReconstructor(A,DM,WFS);

% Now program this crazy thing.
% We can look at the singular values and choose things manually later.  
% For now, we will make default assumptions.  You can rebuild it quickly
% later.  The MMT currently runs with 56 modes corrected.

% OWD = sqrt(MAX_MODES/pi);
% RECON.program(D,6*sqrt(2)); % OWD is ~8.5 lambda/D for programming.
% RECON.program(D,OWD); % OWD is ~8.5 lambda/D for programming.
% semilogy(RECON.s/RECON.s(1));

% F = AOField(A);
% F.lambda = RECON.lambda;
% [x,y] = coords(F);
% % Look at the reconstructor for good "gut" feelings...
% for n=1:size(RECON.RECONSTRUCTOR,2)
%     DM.setActs(500*RECON.RECONSTRUCTOR(:,n)).touch.render;
%     F.planewave*A*DM;
%     imagesc(x,y,F.interferometer(1),[0 3]);
%     sqar;
%     axis xy;
%     drawnow;
% end
% 
% RECON.rebuild(56);
% RECON.Nmodes
% DM.nActs

save data/ObsCirc_Model_working A D  % Paranoid.

% Now run something like the canned_GMT script, but don't load anything in
% first.  
