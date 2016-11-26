% A demo of the AOAtmo class.
% 
% 20150225 JLCodona
% 20151204 JLCodona... play version.

THICKNESS = 1000;
% THICKNESS = 3000
% NSCREENS = 10;
NSCREENS = 5

% Cn2 =  ATMO.setThickness(THICKNESS).setR0(r0).Cn2
Cn2 = 3e-13

WIND = [0 5]; % All winds the same.

lambda = AOField.RBAND; % Red light.
r0 = 0.15; % r0 is 15 cm at 500 nm.

% D = 1.54;  % Kuiper
D = 0.254;  % QMT
secondary = D/3;

SPACING = 0.005;            % fine spacing makes nice pupil images but is really overkill.
aa = SPACING;              % for antialiasing.
% aa = 0.04;
spider = 0.0254;
% spider = 0.01;

PUPIL_DEFN = [
   0 0 D         1 aa 0 0 0 0 0
   0 0 secondary 0 aa/2 0 0 0 0 0
%    0 0 spider   -2 aa 4 0 D/1.9 0 0
   ];

% Since this demo only uses one AOSegment, I will not use the AOAperture wrapper.  
% I will only use AOSegment.  This is fine for simple pupils.

A = AOSegment;
A.spacing(SPACING);
A.name = 'QMT';
A.pupils = PUPIL_DEFN;
A.make;


clf;
colormap(gray);

% A.show;
% input 'Press ENTER to Continue...'
% pause(3);

%% Make a multi-layer atmosphere.

% AOAtmo is a class that holds multiple AOScreens and uses them t create a
% signle object that acts like a fancyAOScreen.

% Make an AOAtmo with 5 layers all of the same strength.

%ATMO = AOAtmo(A);
ATMO = AOAtmo(A);
ATMO.name = 'Layered Atmosphere';

% for n=1:NSCREENS
%     ps = AOScreen(2*1024);
%     ps.name = sprintf('Layer %d',n);
%     ps.spacing(0.02);
%     ps.setCn2(1e-13,THICKNESS/NSCREENS);
%     ATMO.addLayer(ps,THICKNESS*n/NSCREENS);
%     %ATMO.layers{n}.Wind = randn([1 2])*5; % random wind layers.
%     ATMO.layers{n}.Wind = WIND; % Uniform wind.
% end    

PSTHICK = THICKNESS/NSCREENS

for n=1:NSCREENS
    ps = AOScreen(1024);
    ps.name = sprintf('Layer %d',n);
    ps.spacing(0.02);
    ps.setCn2(Cn2,PSTHICK);
    ps.interpolate_method = 'cubic';
    %ps.interpolate_method = 'linear';
    
    ATMO.addLayer(ps,PSTHICK*n);
    
    %ATMO.layers{n}.Wind = randn([1 2])*5; % random wind layers.
    ATMO.layers{n}.Wind = WIND; % Uniform wind.
end    

ATMO.show;

r00 = ATMO.totalFriedScaleStar
NAPER = round(D/r00);
ATMO.spacing(r00);
ATMO.resize(NAPER);

ATMO.useGeometry(false);

% 
% % Define some beacons from which to calculate ATMO OPLs...
% %% Guide star selection
% 
% % LGS_BEACON0 = [0 0 100];
% % LGS_BEACON = [0 1/206265 1] * SODIUM_LAYER;  % Offset by 1 arcsec.
% STAR = [0 0 1e10];
% LEO_TARGET = [0 0 400e3];
% GEOSYNC_TARGET = [0 0 42e6];
% 
% % NGS CASE
% GUIDE_STAR = STAR; % pick one.  
% SCIENCE_OBJECT = STAR; % pick one.  
% 
% % Na LGS CASE
% % GUIDE_STAR = LGS_BEACON; 
% % SCIENCE_OBJECT = STAR; 
% 
% % Looking at GeoSynchronous satellites using LGS
% % GUIDE_STAR = LGS_BEACON; 
% % SCIENCE_OBJECT = GEOSYNC_TARGET; 
% 
% % ATMO.BEACON = GUIDE_STAR; % Set this so ATMO knows how to compute the wavefront.
% 
% ATMO.make;
% 
% N1=2; N2=2;
% 
% %% Make an AOField object.
% 
% F = AOField(A);
% F.name = 'Field';
% % F.resize(1024); % make it big to study the field before the pupil.
% F.FFTSize = 1024; % Used to compute PSFs, etc.
% F.lambda = lambda;
% 
% F.planewave*A;
% F.show;
% 
% [x,y] = F.coords;
% 
% % input 'Continue...'
% 
% % This adds a reference wave to the field and computes the intensity.
% % imagesc(x,y,F.interferometer(1),[0 3]);
% % sqar;
% % axis xy;
% % drawnow;
% % input 'Continue...'
% 
% THld = F.lambda/D * 206265; % Lambda/D in arcsecs.
% FOV = 4;
% PLATE_SCALE = THld/5;
% 
% F.planewave*A;
% [PSF,thx,thy] = F.mkPSF(FOV,PLATE_SCALE);
% PSFmax = max(PSF(:));
% 
% fprintf('Use light from a finite-distance beacon.\n')
% 
% ATMO.setBeacon([0 0 150000]);
% 
% % This includes the geometry as well as the OPD from the layers...
% ATMO.useGeometry(false);
% 
% t = 0;
% for RANGE=25:25:2*THICKNESS
%     ATMO.setBeacon([0 0 RANGE]);
%     ATMO.setObsTime(t);
%     F.planewave*ATMO*A;
%     
%     subplot(N1,N2,1);
%     ATMO.show;
%     title(sprintf('wavefront:time=%.3fs',t));
%     
%     subplot(N1,N2,2);
%     F.show;
%     title('Field');
% 
%     
%     subplot(N1,N2,3);
%     [PSF,thx,thy] = F.mkPSF(FOV,PLATE_SCALE);
%     imagesc(thx,thy,log10(PSF/PSFmax),[-4 0]);
%     daspect([1 1 1]);
%     axis xy;
%     colorbar;
%     title('PSF');
%     
%     subplot(N1,N2,4);
%     imagesc(x,y,F.interferometer(1),[0 4]); sqar;
%     axis xy;
%     title('interferometer');
%     
%     
%     drawnow;
% end
% 
% % for t=0:.01:1
% %     ATMO.setObsTime(t);
% %     F.planewave*ATMO*A;
% %     
% %     subplot(N1,N2,1);
% %     ATMO.show;
% %     title(sprintf('wavefront:time=%.3fs',t));
% %     
% %     subplot(N1,N2,2);
% %     F.show;
% %     title('Field');
% % 
% %     
% %     subplot(N1,N2,3);
% %     [PSF,thx,thy] = F.mkPSF(FOV,PLATE_SCALE);
% %     imagesc(thx,thy,log10(PSF/PSFmax),[-4 0]);
% %     daspect([1 1 1]);
% %     axis xy;
% %     colorbar;
% %     title('PSF');
% %     
% %     subplot(N1,N2,4);
% %     imagesc(x,y,F.interferometer(1),[0 4]); sqar;
% %     axis xy;
% %     title('interferometer');
% %     
% %     
% %     drawnow;
% % end
% 
% % This turns off the geometry and is like focusing on the beacon.
% 
% % % fprintf('Focus on the beacon (only OPL variations).\n')
% % % 
% % % ATMO.useGeometry(false);
% % % 
% % % for t=0:.01:1
% % %     ATMO.setObsTime(t);
% % %     F.planewave*ATMO*A;
% % %     
% % %     subplot(N1,N2,1);
% % %     ATMO.show;
% % %     title(sprintf('time=%fs',t));
% % %     
% % %     subplot(N1,N2,2);
% % %     F.show;
% % %     title('Field');
% % % 
% % %     
% % %     subplot(N1,N2,3);
% % %     [PSF,thx,thy] = F.mkPSF(FOV,PLATE_SCALE);
% % %     imagesc(thx,thy,log10(PSF/PSFmax),[-4 0]);
% % %     daspect([1 1 1]);
% % %     axis xy;
% % %     colorbar;
% % %     title('PSF');
% % %     
% % %     subplot(N1,N2,4);
% % %     imagesc(x,y,F.interferometer(1),[0 4]); sqar;
% % %     axis xy;
% % %     title('interferometer');
% % %     
% % %     
% % %     drawnow;
% % % end
% % % 
% % % % This computes the wavefronts across the pupil for an array of point
% % % % sources, all at the same time.
% % % 
% % % clf;
% N = 5;
% 
% 
% fprintf('Beacon array Phase Maps focusing on beacons.\n')
% 
% % BEACONX = linspace(-10,10,N);
% % % ATMO.useGeometry(true);
% % ATMO.useGeometry(false);
% % ATMO.setObsTime(0);
% % 
% % for n=1:N
% %     for m=1:N
% %         subplot(N,N,(n-1)*N+m);
% %         ATMO.setBeacon([BEACONX(n) BEACONX(m) 150000]);
% %         imagesc(ATMO.grid);
% %         axis xy;
% %         axis off;
% %         daspect([1 1 1]);
% %         
% %     end
% %     drawnow;
% % end
% % 
% % 
% % n=1; m=1;
% % ATMO.setBeacon([BEACONX(n) BEACONX(m) 150000]);
% % 
% % fprintf('Beacon array Phase Maps using [1,1] as a reference.\n')
% % 
% % reference = ATMO.grid; % Look at the residual from this reference...
% % 
% % for n=1:N
% %     for m=1:N
% %         subplot(N,N,(n-1)*N+m);
% %         ATMO.setBeacon([BEACONX(n) BEACONX(m) 150000]);
% %         imagesc(ATMO.grid-reference,[-1 1]*1e-6);
% %         axis xy;
% %         axis off;
% %         daspect([1 1 1]);
% %         %colorbar;
% %     end
% %     drawnow;
% % end
