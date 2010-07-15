%%Script to create a time series of actuator positions using AOSim2
%PMH 20090802
%This script was used in August 2009 to generate phase screens for
%the Microgate/ADS estimation of 
%The script has a bug due to the edge of the screen being reached.
%This could be fixed, but is left as written for documentation purposes.

clear;
% Load in aperture definitions
load data/NewGMTPupil.mat
% Import GMT segment coordinates.
load data/GMT_672act_segment_PMH.mat;
% Make a DM with an OPD grid matched to the Aperture A...
DM = AODM(A);
DM.name = 'GMT Hexapolar Adaptive Secondary';
%% Add in the actuators and tag them with their segment positions.  
% This loops over th outer segments...
for n=1:6
    DM.addActs(GMT_Actuators,n,A);
end
% This adds the actuators for the center segment.
DM.addActs(GMT_Actuators,7,A);
ActPositions=DM.actuators(:,1:2);
%Save actuator positions
save GMT_Actuator_Positions.mat ActPositions


%% Now Generate the Phase screens and loop through in time

%Will carry out numsim simulations, each numframes frames long
numframes=1000; %numfer of frames in each simulation
numsim=10;      %number of simulations
WFS_FPS=1000;   %running at 1 kHz
numact=length(DM.actuators);
fprintf('Number of actuators is %d \n',numact);
Position=zeros(numact,numframes,numsim);   %units is m
Phaseabs=zeros(numact,numframes,numsim);   %radians at V band
for a=1:numsim    %S
    % Define the Atmosphere model and winds aloft.
    ATMO = AOAtmo(A);
    WFlow = AOScreen(1024,0.21,500e-9);     %50th percentile
    WFhigh = AOScreen(2048,0.30,500e-9);    %50th percentile
    %WFlow = AOScreen(2048,0.13,500e-9);     %90th percentile
    %WFhigh = AOScreen(4096,0.20,500e-9);    %90th percentile
    ATMO.addLayer(WFlow,1000);
    ATMO.addLayer(WFhigh,8000);
    ATMO.layers{1}.Wind = [20 0];           %50th percentile
    ATMO.layers{2}.Wind = [0 45];           %50thpercentile
    %ATMO.layers{1}.Wind = [30 0];           %90th percentile
    %ATMO.layers{2}.Wind = [0 60];           %90thpercentile
    WFlow.name = 'Lower altitude turbulence';
    WFhigh.name = 'High altitude turbulence';
    r0 = ATMO.totalFriedScale;
    th_scat = AOField.VBAND/r0*206265;
    fprintf('The total r0 is %.2f cm.\n',100*ATMO.totalFriedScale);
    fprintf('The seeing is %.2f arcsecs.\n',th_scat);
    beaconheight=1e7;
    STAR = [0 0 beaconheight];
    ATMO.BEACON = STAR; % Set this so ATMO knows how to compute the wavefront.
    % Turning this off is like using dynamic refocus.
    ATMO.GEOMETRY = false;
    for n=1:numframes
        t = n/WFS_FPS;
        ATMO.time = t;    
        clf; % Clear figure
        %[x,y]=ATMO.coords;
        %colormap(gray);
        %phase=ATMO.grid-beaconheight;
        DM.actuators(:,3) = -ATMO.interpGrid(DM.actuators(:,1),DM.actuators(:,2));
        DM.removeMean;
        acts=DM.actuators;
        %imagesc(x,y,phase);
        %daspect([1 1 1]);
        %axis xy;
        %drawnow;
        Position(:,n,a)=DM.actuators(:,3);
        %Below used for estimating t_0
        %Phaseabs(:,n,a)=abs((Position(:,n,a)-Position(:,1,a))*2*3.14159/AOField.VBAND);
        %meanphase=mean(Phaseabs(:,n,a));
        %fprintf('Iteration %d,  Time=%.0f ms   Phase %.3f radians\n',a,(t*1000)-1,meanphase);
        fprintf('Iteration %d,  Time=%.0f ms \n',a,(t*1000)-1);
    
    end
end
%save Positions value
save GMT_ActuatorvsTime_50thSeeing.mat Position