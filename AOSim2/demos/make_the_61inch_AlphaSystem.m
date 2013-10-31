%% How to build and program the AOSim2 MMT AO model.
% Johanan L. Codona, Steward Observatory, University of Arizona.
% 20090426 JLCodona: First-light version.
% 20090504 JLCodona: Had to add a Seg.make line.  This is a bug that needs
% to be found.  Sorry.

%% Start clean...
% close
% clear classes

%% Load in some definitions...

%     * Elevation: 2510 m = 8235 ft.
%     * Latitude: +32° 24' 59.3" N
%     * Longitude: 110° 44' 04.3" W
%             These coordinates are accurate to about 1".
%     * Time Zone: +7 hours 
% 
%     * Primary Mirror Diameter: 1.54 m = 61 inches
%     * Primary Focal Ratio: f/4 
% 
%     * f/13.5 Cassegrain focus
%           o Plate Scale: 100 microns/arcsec = 10.0 arcsec/mm (nominal)
%           o Useful Field of View: >435 arcsec diameter
%           o Secondary Diameter: 40.96 cm 
%     * f/45 Cassegrain focus
%           o Plate Scale: 351 microns/arcsec = 2.85 arcsec/mm (nominal)
%           o Useful Field of View: >325 arcsec diameter
%           o Secondary Diameter: 14.5 cm 
% 
%     * Typical Seeing: 1-2" 

D = 1.54;
% secondary = 40.96/100;
secondary = 14.5/100;
aa = 0.04;
spider = aa;


% PUPIL_DEFN = [
%    0 0 D         1 aa 0 0 0 0 0
%    0 0 secondary 0 aa 0 0 0 0 0
%    0 0 spider   -2 aa 4 0 0 0 0];

PUPIL_DEFN = [
   0 0 D         1 aa 0 0 0 0 0
   0 0 secondary 0 aa 0 0 0 0 0 ];

Seg = AOSegment;
Seg.name = 'Kuiper 61inch Primary';
Seg.spacing(0.01);
Seg.pupils = PUPIL_DEFN;
Seg.make;

clf;
% Seg.touch.make.show;
A = AOAperture;
A.name = 'Kuiper 61 inch';
A.addSegment(Seg);
A.show;
colormap(gray);

%% Make a DM with an OPD grid matched to the Aperture A...
DM = AODM(A);
DM.name = 'BMC 140 element MEMS';

xx = (1:12)*D/12;
xx = xx - mean(xx);
[X,Y] = meshgrid(xx,xx);
ACTS = [X(:) Y(:)];

%% Add in the actuators.
DM.addActs(ACTS,1,A);

% Mark BAD actuators
%DM.disableActuators(BAD);
%DM.disableActuators(MMT_BADACTS_NO_CURRENT);
%DM.disableActuators(MMT_BADACTS_NO_POSITION);

% Specify the boundary conditions... 
% DM.defineBC(5,8); % A circle of 8 null points at 5m radius.
DM.defineBC(D,8); % A circle of 8 null points at D radius.
DM.plotRegions; daspect([1 1 1]); drawnow;

%% Build the Shack-Hartmann WFS.

WFS = AOWFS(A,D/12);
WFS.name = 'Omega SHWFS';
A.show; WFS.quiver(1); drawnow; % Show them.

%% Now for some real work.  Building the RECONSTRUCTOR...
RECON = AOReconstructor(A,DM,WFS);

% Now program this crazy thing.
% We can look at the singular values and choose things manually later.  
% For now, we will make default assumptions.  You can rebuild it quickly
% later.  The MMT currently runs with 56 modes corrected.

RECON.adhocProgram(D/12*3);
RECON.show;

% OWD = sqrt(MAX_MODES/pi);
% % RECON.program(D,6*sqrt(2)); % Use Fourier modes. OWD is ~6 lambda/D for programming.
% RECON.zprogram(D,12);  % program using Zernikes.
% RECON.rebuild(56).show;

% % RECON.dhprogram(D,11); % program using disk harmonics.
% semilogy(RECON.s/RECON.s(1));
% 
% F = AOField(A);
% F.lambda = RECON.lambda;
% [x,y] = coords(F);
% % Look at the reconstructor for good "gut" feelings...
% lim0 = min(RECON.RECONSTRUCTOR(:));
% lim1 = max(RECON.RECONSTRUCTOR(:));
% for n=1:size(RECON.RECONSTRUCTOR,2)
%     DM.setActs(5*RECON.RECONSTRUCTOR(:,n)).touch.render;
%     %F.planewave*A*DM;
%     %     imagesc(x,y,F.interferometer(1),[0 3]);
%     %     sqar;
%     %     axis xy;
%     imagesc(x,y,DM.grid .* A.grid,[lim0 lim1]);colorbar;
%     daspect([1 1 1]);
%     drawnow;
% end
% 
% RECON.rebuild(500);
% RECON.Nmodes
% DM.nActs
% 
% % Mark BAD actuators
% % DM.disableActuators(BAD);
% % DM.disableActuators(MMT_BADACTS_NO_CURRENT);
% % DM.disableActuators(MMT_BADACTS_NO_POSITION);
% 
% save MMTAO_Model_jlc A DM WFS RECON D % Paranoid.
% 
% % Now run something like the canned_GMT script, but don't load anything in
% % first.  
