%% Clear Workspace
clear all;
clc;
close all;

%% Some Starting Comments
% This code makes a model of an IrisAO PTT111 type mirror.  It does not
% include the limits on stroke or tip/tilt angle that a real mirror has.
% The size of the segments is accurate to within 3 microns of a real
% mirror, and the piston, tip, and tilt functionality is correct (to the
% best of my knowledge right now, after lots of debugging).
%
% The scripts that are used in creating the mirror model and using it are
% all in the folder you found this one.  The ones that are important to
% look at should you want to adapt this script for your own purposes are:
%******************************
% makeIrisAODM.m
% IrisAOComputeZernPositions.m
%******************************
% These functions are what are used to make, and then update the shape of
% the mirror. IrisAOComputeZernPositions.m is not required to be used, but
% can be helpful if you want to apply zernike polynomials.
%
% Please also note that this REQUIRES AOSim2 to be on your current pathway.
% The model is completely built in, and depends on the functionality of
% AOSim2.
%
% If you find a bug, or that it is not working accurately, please let me
% know what the problem you encountered is and I will work to fix it.  
%
% ATR (email: atrodack@email.arizona.edu)



%% Set Initial Parameters and Flags
%Parameters (currently set to Specs from UAWFC Testbed)
segpitch = 606e-6;
magnification = 1;
lambda = AOField.VBAND;

%Flags
verbose_makeDM = false; %turns on/off plotting the mirror as it is constructed

Scalloped_Field = true; %turns on/off returning an AOField Object that 
%encodes the actual surface shape of the segments.

%% Make the Mirror
if Scalloped_Field == true
    [DM,F_scal] = makeIrisAODM(magnification,verbose_makeDM,Scalloped_Field);
else
    DM = makeIrisAODM(magnification,verbose_makeDM,Scalloped_Field);
end
DM.lambdaRef = lambda;

clc; %clears the command window of the text generated from adding segments to the AOAperture object
fprintf('Mirror Constructed\n');

%Save the Coordinate vectors of the DM
[x,y] = DM.coords;

%Compute the extent of the Segment's Grid to Construct Aperture
extentx = abs(min(x)) + abs(max(x));
extenty = abs(min(y)) + abs(max(y));
fprintf('X Width = %0.4f mm\n',extentx * 10^3);
fprintf('Y Width = %0.4f mm\n',extenty * 10^3);

%% Make an Aperture and a Field
A = AOSegment(DM);
D = (extentx);
dx = (magnification *segpitch) / 100;
PNECO = [0 0 D 1 1.5*dx 0 0 0 0 0];
A.pupils = PNECO;
A.make;

DM.trueUp;
A.centerOn(DM);

if Scalloped_Field == false
    F = AOField(A);
    F.name = 'IrisAO DM';
    F.FFTSize = [2048 2048];
    F.lambda = lambda;
end

THld = lambda/D * 206265; % Lambda/D in arcsecs.
FOV = 50*THld; % FOV for PSF computation
PLATE_SCALE = THld/5; % Pixel Size for PSF computation

%% Get Something to Put on Mirror
%These are all options of what can be sent to the mirror, in the syntax an
%actual IrisAO Mirror would like them, with a minor difference. The
%actaul mirror wants its inputs in microns and milliradians, but this
%wants meters and radians.

% Flatten the Mirror
% PTTpos = zeros(37,3);

%Apply only a Tip
% PTTpos = horzcat(horzcat(zeros(37,1)*10^-6,ones(37,1)*10^-3),zeros(37,1)*10^-3);

% Apply a Tip and a Tilt
% PTTpos = horzcat(horzcat(zeros(37,1)*10^-6,ones(37,1)*10^-3),ones(37,1)*10^-3);

% Set a Random Mirror Shape
% PTTpos = horzcat(horzcat(randn(37,1)*10^-6,randn(37,1)*10^-3),randn(37,1)*10^-3);

% Set Random Piston
% PTTpos = horzcat(horzcat(randn(37,1)*10^-6,zeros(37,1)*10^-3),zeros(37,1)*10^-3);

% Apply Random Tip and a Tilt
% PTTpos = horzcat(horzcat(zeros(37,1)*10^-6,randn(37,1)*10^-3),randn(37,1)*10^-3);

% Set a Zernike Polynomial
Zernike_Number = 5;
Zernike_Coefficient_waves = 1;
PTTpos = IrisAOComputeZernPositions( lambda, Zernike_Number, Zernike_Coefficient_waves);


%% Map the PTTpos Matrix
load('IrisAO_SegMap.mat');
PTT = zeros(37,3);

for ii = 1:37
    mapped_segment = IrisAO_SegMap(ii);
    PTT(ii,1:3) = PTTpos(mapped_segment,:);
end

DM.PTT(PTT);
DM.touch;
DM.render;
figure(1);
DM.show;

if Scalloped_Field == true
    F = F_scal.copy;
    F * DM * A;
else
    F.planewave * DM * A;
end

figure(2);
[PSF,thx,thy] = F.mkPSF(FOV,PLATE_SCALE);
PSFmax = max(PSF(:));
% imagesc(thx,thy,PSF/PSFmax);
imagesc(thx,thy,log10(PSF/PSFmax),[-3,0]);
colormap(gray);
axis xy;
sqar;
title(sprintf('PSF\n'));
