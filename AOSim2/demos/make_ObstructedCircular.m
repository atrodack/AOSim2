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
obstruction = 0.2;
% dx = 0.02;
dx = 1/30;

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
A = AOAperture();
A.spacing(dx);
A.name = Seg.name;
A.addSegment(Seg);
A.show;
colormap(gray);

% save data/ObsCirc_Model_working A D  % Paranoid.

% Now run something like the canned_GMT script, but don't load anything in
% first.  
