%%Script estimating the effect of the ASM print through on the GMT PSF
%PMH 20090901
%Script used to show that the quilting error does not affect the PSF
%even for HCAO.

clear;
%%Set parameters likely to be fiddled with
GAMMA = 2;  % This is the gamma term for the PSF grayscale stretch.
SCIENCE_WAVELENGTH = AOField.JBAND;
%% Create Aperture
    CENTERX= [8.70  4.3500 -4.3500 -8.70 -4.3500 4.3500 0];
    CENTERY= [0.00 -7.5344 -7.5344  0.00  7.5344 7.5344 0];
    CENTERS=[CENTERX' CENTERY'];
    SegmentPupil=[0 0 8.4 1 0.001 0 0 0 0 0];
    SEG=AOSegment;
    SEG.pupils=SegmentPupil;
    SEG.make;
    A=AOAperture;
    A.spacing(0.03);
    for n=1:7;
        A.addSegment(SEG,[CENTERX(n) CENTERY(n)]);
    end
    A.show;
    % Make a DM with an OPD grid matched to the Aperture A...
    DM = AODM(A);
    DM.name = 'GMT Hexapolar Adaptive Secondary';
    touch(DM);
    %% Add in the actuators and tag them with their segment positions.  
    % Import GMT actuator positions for a single segment.
    load data/GMT_DM_673_actuators_hexapolar.mat;
    LBT_Actuators=LBT_Actuators/3.15;
    % This loops over all segments...
    for n=1:7
        DM.addActs(LBT_Actuators,n,A);
    end 
numact=length(DM.actuators);
fprintf('Number of actuators is %d \n',numact);  

F = AOField(A);
F.lambda = SCIENCE_WAVELENGTH;  % Science Wavelength
F.FFTSize = 2048*[1 1]; % How should this be scaled??

%%Set up plotting parameters
% This is the brightest pixel seen to date.
Ipeak = 0;  
[x,y]=coords(F);
%First value is Peak to value in meters, second is ~width of bumps
DM.addQuilt(30e-9,0.1);



    %% Now calculate image quality for the science object
	F.planewave*A*DM;  
    F.show;
    drawnow;
    %% Plot some interesting pictures...
    %clf; % Clear the Figure each time step
    FOV = 3;  % In arcsecs.
	RNG = FOV * [-1 1];
	PSF = F.mkPSF(FOV,FOV/1000);
    mask = (A.grid>0.5);    
    g = F.grid_;
    STREHL = abs(mean(g(mask)))^2;     
	Ipeak = max(Ipeak,max(PSF(:)));
    MAG=log10(PSF/Ipeak);
    imagesc(RNG,RNG,MAG)
	%imagesc(RNG,RNG,(PSF/Ipeak).^(1/GAMMA));
	axis xy;
    daspect([1 1 1]);
    title(sprintf('PSF (\\lambda=%.2g microns, Strehl=%.3g)',...
        F.lambda*1e6,STREHL));
    xlabel('arcsecs');    ylabel('arcsecs');
