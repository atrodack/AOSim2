% Define an Aperture.  
% You can use AOSegment if it is going to just be one piece.  
% Using AOAperture allows you to build a more complicated pupil out of
% separately-actuated AOSegments.

A = AOSegment;
A.name = 'Pupil';

D = 10*0.0254; % 10 inch primary.
secondary = D/4;

SPACING = D/100;
aa = SPACING/3; % This is smoothing or "antialiasing".
% aa = 0.04;
% spider = 0.0254;
spider = 0.01;

% This is the rather primative way I define pupils.

% Supported Apodizations: apod_type = pupils(6)
%
% 0: Cosine (arg 5 holds the width)
% 1: Sonine (nu is arg 7)  (beware of definitions!)
% 2: Elliptical: arg7 is Dy.
% 3: Angel-Cheng-Woolf
% 4: Angel2
% 5: Spergel2D: (arg 5 is gaussian equivalent diameter).
% 6: Woolf slab: (arg 5 is the height)
% 7: Specified circular mask. (set shadingMask and shadingRadii
%     mask) arg 7 is an overall complex phasor.)
% 8: ellipse: (use AOAperture/addEllipse)
%
%
% "Transmission" types:
% 0: circular hole
% 1: mirror
% -1: wedge
% -2: spider vanes (3: width, 6: nvanes, 7: starting theta)

PUPIL_DEFN = [
   0 0 D         1 aa 0 0 0 0 0
   0 0 secondary 0 aa/2 0 0 0 0 0
   %0 0 spider   -2 aa 4 0 D/1.9 0 0
   ];

A.spacing(SPACING);
A.pupils = PUPIL_DEFN;
A.make; % This builds the grid from the definitions in A.pupils.
colormap(gray(256));
A.show;

input 'Press a key to continue...'

% This is how to load in an image as a pupil map...

img = imread('doubleSlit.png');
% AOGrid will only use the first plane if it has more dims.

whos img
size(img)

A.grid(img).show;

% Note that if you say A.make it will overwrite your manual grid with a
% rendered one based on the pupils definitions.  It will not do this unless
% you tell it to.

D = max(A.extent)

%% Now to make an AOField.

F = AOField(A);  % Use A as the template.
F.name = 'Field';

% This is pretty small, so make it bigger.

F.resize(256);
%F.resize(64);

F.lambda = 1e-6;
F.FFTSize = 1024;

THld = F.lambda/D * 206265;

F.planewave*A;

subplot(1,2,1);
F.show;

[PSF,thx,thy] = F.mkPSF(25*THld,THld/4);
PSFmax = max(PSF(:));

subplot(1,2,2);
% imagesc(thx,thy,log10(PSF/PSFmax),[-4 0]);
imagesc(thx,thy,(PSF/PSFmax));
daspect([1 1 1]);

%% Make a phase screen

PS = AOScreen(2048);
PS.setR0(0.15);
PS.make.show;

clf;

for n=1:100
    % Move the phase screen manually...
    PS.grid(circshift(PS.grid,[0 1]));
    
    F.planewave*PS*A;
    
    subplot(1,2,1);
    F.show;
    
    [PSF,thx,thy] = F.mkPSF(10*THld,THld/4);
    PSFmax = max(PSF(:));
    
    subplot(1,2,2);
    imagesc(thx,thy,PSF/PSFmax);
    daspect([1 1 1]);
    
    drawnow;
end

% Again with weaker turbulence...

PS.setR0(0.5);
PS.make.show;

clf;

for n=1:100
    % Move the phase screen manually...
    PS.grid(circshift(PS.grid,[0 1]));
    
    F.planewave*PS*A;
    
    subplot(1,2,1);
    F.show;
    
    [PSF,thx,thy] = F.mkPSF(10*THld,THld/4);
    
    subplot(1,2,2);
    imagesc(thx,thy,PSF);
    daspect([1 1 1]);
    
    drawnow;
end


% Do it again with two holes...

%img = imread('doubleHoles.png');
A.grid(imread('doubleHoles.png')).show;

clf;

PS.setR0(0.15);
PS.make;

CCD = 0; % Add a long exposure.
for n=1:200
    % Move the phase screen manually...
    PS.grid(circshift(PS.grid,[1 2]));
    
    F.planewave*PS*A;
    
    subplot(2,2,1);
    F.show;
    
    [PSF,thx,thy] = F.mkPSF(10*THld,THld/4);
    CCD = CCD + PSF;
    
    subplot(2,2,2);
    imagesc(thx,thy,PSF/PSFmax);
    %imagesc(thx,thy,CCD);
    daspect([1 1 1]);
    
    subplot(2,2,4);
    %imagesc(thx,thy,PSF/PSFmax);
    imagesc(thx,thy,CCD);
    daspect([1 1 1]);

    drawnow;
end


%% Now do that again with our pupil...

A = AOSegment;
A.name = 'Pupil';

D = 10*0.0254; % 10 inch primary.
secondary = D/4;

SPACING = D/100;
aa = SPACING/3; % This is smoothing or "antialiasing".
% aa = 0.04;
% spider = 0.0254;
spider = 0.01;

PUPIL_DEFN = [
   0 0 D         1 aa 0 0 0 0 0
   0 0 secondary 0 aa/2 0 0 0 0 0
   %0 0 spider   -2 aa 4 0 D/1.9 0 0
   ];

A.spacing(SPACING);
A.pupils = PUPIL_DEFN;
A.make; % This builds the grid from the definitions in A.pupils.
A.show;

F = AOField(A);
F.lambda = 1e-6;
F.FFTSize = 1024;

THld = F.lambda/D * 206265;

PS.spacing(0.01);
PS.setR0(0.08);
PS.make;

CCD = 0; % Add a long exposure.
for n=1:200
    % Move the phase screen manually...
    PS.grid(circshift(PS.grid,[1 2]));
    
    F.planewave*PS*A;
    
    subplot(2,2,1);
    F.show;
    title('pupil field');
    
    [PSF,thx,thy] = F.mkPSF(10*THld,THld/4);
    CCD = CCD + PSF;
    
    subplot(2,2,2);
    imagesc(thx,thy,log10(PSF/max(PSF(:))),[-3 0]);
    %imagesc(thx,thy,PSF/PSFmax);
    %imagesc(thx,thy,CCD);
    daspect([1 1 1]);
    title('I(t)');
    
    subplot(2,2,4);
    %imagesc(thx,thy,PSF/PSFmax);
    imagesc(thx,thy,log10(CCD/max(CCD(:))),[-3 0]);
    daspect([1 1 1]);
    title('log exposure');
    
    drawnow;
end

