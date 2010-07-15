%%Script for initial testing of AOSim2 results
%This is a modified version of the run_myGMT scripts by JLC in demos
%PMH 20091215

%Script used in "Preliminary Report" in December 2009

clear;
numframes=100;     %number of frames to generate
WFS_FPS=1000;       %running at  "WFS_FPS" Hz
AO_STARTTIME=0.002; %Starting time of AO correction
PHASE_START=0.05;  %Starting time of phase correction
PHASE_END=numframes/WFS_FPS;
FOV = 0.1;  % In arcsecs.
numtrials=3;

STREHL = zeros(numframes,numtrials);
%STREHLp= zeros(numframes,1);


for i=1:numtrials  %Run complete program multiple times
    %%Set parameters likely to be fiddled with
    
    if (i==1)  %Run each iteration with exact same random number seed
        defaultStream = RandStream.getDefaultStream;
        savedState = defaultStream.State;
    else
        defaultStream.State=savedState
    end
    

    %Gain on the error feedback in the AO loop
    gain = 1; % Seems to be optimum at ~0.3
    if(i==1)
        pgain = 0.2; % This only works at low gain 
    else
        pgain=0.2;
    end

    GAMMA = 2;  % This is the gamma term for the PSF grayscale stretch.
    SCIENCE_WAVELENGTH = AOField.KBAND;

    %%First decide whether to load in an existing reconstructor or build one.
    REBUILD=0;

    if(REBUILD) 
        %This is derived from make_the_GMT_AO script in demos by JLC
        NMODES=700;
        % Load in aperture definitions
        load data/NewGMTPupil.mat
        %Defines:   A, AOAperture object, defining the pupil
        %           CENTERS, position of the center of each GMT segment
        %           DM, AODM object, defining the actuator attributes
        %           Seg1, AOSegment object for outer segments
        %           seg7, AOSegment object for center segment        

        % Make a DM with an OPD grid matched to the Aperture A...
        DM = AODM(A);
        DM.name = 'GMT Hexapolar Adaptive Secondary';
        %% Add in the actuators and tag them with their segment positions.  
        % Import GMT actuator positions for a single segment.
        load GMT_672act_segment_PMH.mat;
        % This loops over all segments...
        for n=1:7
            DM.addActs(GMT_Actuators,n,A);
        end    
        % Specify the boundary conditions... 
        DM.defineBC(50,5*6); % A circle of 5*6=30 null points at 30m radius.
        DM.plotRegions; daspect([1 1 1]); drawnow;

        %% Build the Shack-Hartmann WFS.
        % The CoDR says to make the subaps 8.4m/17 ~ 0.5m.  
        % We do this and fill the whole pupil bounding box with subaps.
        d = 8.417;
        BB=A.BBox;
        D= max(BB(2,:)-BB(1,:));
        WFS = AOWFS(A,d/17);
        WFS.name = 'GMT Shack-Hartmann WFS';
        A.show; WFS.quiver(1); drawnow; % Show them.

        fprintf('D is %g\n',D);
        %Define and build the RECONSTRUCTOR...
        RECON = AOReconstructor(A,DM,WFS);
        %Set Scutoff high to make the cutoff be whatever the modes explored is.
        RECON.Scutoff=30000;
        ORDER = sqrt(NMODES*2);
        RECON.zprogram(D,ORDER);
        %Use rebuild to restrict reconstructor to exact number of modes
        %RECON.rebuild(NMODES);  
        %semilogy(RECON.s/RECON.s(1));
        save GMTAO_DemoModel A DM WFS RECON d D % 
        F = AOField(A);
        F.lambda = RECON.lambda;
        [x,y] = coords(F);
    else  %Load in pre-made reconstructors
        %load GMT66z_PMH
        %load GMT231z_PMH
        if(i==1)
            load GMTAO548 
        elseif(i==2)
            %load GMTAO548; 
            load GMTAO548;
        else
            %load GMTAO548; 
            load GMTAO548;
        end
        %load GMTAO757
        %load GMTAO1068
        
    end
    DM.setActs(0);
    A.trueUp;
    numact=length(DM.actuators);
    fprintf('Number of actuators is %d \n',numact);

    %% Now Generate the Phase screens 
    %Define the Atmosphere model and winds aloft.
    ATMO = AOAtmo(A);
    WFlow = AOScreen(1024,0.21,500e-9);     %50th percentile
    WFhigh = AOScreen(2048,0.30,500e-9);    %50th percentile
    %WFlow = AOScreen(2048,0.13,500e-9);     %90th percentile
    %WFhigh = AOScreen(4096,0.20,500e-9);    %90th percentile
    ATMO.addLayer(WFlow,1000);
    ATMO.addLayer(WFhigh,8000);
    ATMO.layers{1}.Wind = [20 0]*1;           %50th percentile
    ATMO.layers{2}.Wind = [0 45]*1;           %50thpercentile
    %ATMO.layers{1}.Wind = [30 0];           %90th percentile
    %ATMO.layers{2}.Wind = [0 60];           %90thpercentile
    WFlow.name = 'Lower altitude turbulence';
    WFhigh.name = 'High altitude turbulence';
    r0 = ATMO.totalFriedScale;
    th_scat = AOField.VBAND/r0*206265;
    fprintf('The total r0 is %.2f cm.\n',100*ATMO.totalFriedScale);
    fprintf('The seeing is %.2f arcsecs.\n',th_scat);
    %TODO: Calculate t_0 from phase screen

    %The ATMO.BEACON is defined in meters in x and y, and z from the primary
    %Examples
    beaconheight=1e7;
    NAHEIGHT=9e4;
    STAR = [0 0 beaconheight];
    LGS = [0 1/206265 1]*NAHEIGHT;  %Sodium guide star 1 arcsec off-axis
    ATMO.BEACON = STAR; % Set this so ATMO knows how to compute the wavefront.
    ATMO.GEOMETRY = false;  %false=dynamic focus of BEACON

    GUIDE_STAR=STAR;        %Position of beacon used for AO correction
    SCIENCE_OBJECT=STAR;    %Position of object used for image quality estimate

    %%Now Set up the field objects
    % Create the WFS and Science AOField objects...
    Fwfs = AOField(A);
    % The Reconstructor was calibrated at a certain wavelength.
    Fwfs.lambda = RECON.lambda;  

    F = AOField(A);
    F.lambda = SCIENCE_WAVELENGTH;  % Science Wavelength
    F.FFTSize = 2048*[1 1]; % How should this be scaled??
    PSF = F.mkPSF(.2,.2/100); %TODO: This is a bug workaround.  FIXME!

    %%Set up plotting parameters
    % This is the brightest pixel seen to date.
    Ipeak = 0;  
    Cpeak=0;
    [x,y]=coords(F);
    TIMES = (1:numframes)/WFS_FPS;

    % Strehl plot setup
    mask = (A.grid>0.5);
    
    xang = zeros(numframes,1);
    yang = zeros(numframes,1);
    maxStrehl = 0.3;
    CCD = 0;

    Dap=25.4;
    sigmafit=AOField.VBAND/AOField.KBAND.*sqrt(0.2944*(RECON.Nmodes^(-sqrt(3)/2))*(Dap/r0)^(5/3));
    Strehlfit=exp(-(sigmafit^2));
    fprintf('Expected Fitting error %6.2f radians, Strehl =%2.2f\n',sigmafit,Strehlfit);
    fprintf('Dap=%f   Num. modes=%d WFS wavelength=%6.2e\n', ...
             Dap,RECON.Nmodes,RECON.lambda);

    touch(DM);     
    %% Start the time loop
    for n=1:numframes
        t = n/WFS_FPS;
        %go from -t/2 to t/2 to center phase screens 
        ATMO.time = t-numframes/WFS_FPS/2.0;      
        %% This is the guts of the AO closed-loop integrating servo....
        ATMO.BEACON = GUIDE_STAR;
        %fprintf('Sensing WF\n');
        WFS.sense(Fwfs.planewave*ATMO*A*DM);
        if(t>AO_STARTTIME)  % Suffer with seeing limit until AO_STARTTIME.
            DM.bumpActs(-gain*RECON.RECONSTRUCTOR * WFS.slopes);
            DM.removeMean;
            if ((t>PHASE_START)&&(t<PHASE_END))
            %Sense and correct segment pistons
            %PISTONS = WFS.magicPistonSensor(Fwfs,A);
            %A.bumpPistons(- pgain*PISTONS);	
                PISTONS = WFS.magicPistonSensor(F,A);
                A.bumpPistons(-pgain*PISTONS);
            end %End of Phase servo    
        end  %End of AO servo

        %% Now calculate image quality for the science object
        %fprintf('Calculating science image\n');
        ATMO.BEACON = SCIENCE_OBJECT;
        F.planewave*ATMO*A*DM;  

        g = F.grid_;
        STREHL(n,i) = abs(mean(g(mask)))^2;

     %% Plot some interesting pictures...
        clf; % Clear the Figure each time step
        subplot(2,3,1);
        imagesc(x,y,DM.grid .* mask); daspect([1 1 1]);
        axis xy;
        
        
        subplot(2,3,2);
        %Calculate centroid
        PSF = F.mkPSF(FOV,FOV/100);
        sizepsf=floor(length(PSF)/2);
        xa=-sizepsf:sizepsf;
        PSF(isnan(PSF))=0;
        norm = sum(PSF(:));
        if(t>(AO_STARTTIME+0.002))
            xang(n) = mean(xa*PSF)/norm*FOV;
            yang(n) = mean(PSF*xa')/norm*FOV;
        end
        xlim([-FOV FOV]);
        ylim([-FOV FOV]);
        ylabel('xcentroid');
        xlabel('ycentroid');
        title(sprintf('Centroid Variation'));
       
        plot(xang(1:n),yang(1:n),'k-'); 
        %A.show;
        %colorbar off;
        %WFS.quiver(1);
        %title('WFS Slopes');
        %imagesc(x,y,(F.interferometer(.75)),[0 3]);
        %title('Science Band Interferometer');
        daspect([1 1 1]);
        %axis xy;

        subplot(2,3,3);
        [x,y]=F.coords;
        imagesc(x,y,angle(g));  %angle is atan2(imag(g)/real(z))
        %imagesc(x,y,(F.interferometer(.75)),[0 3]);
        %title('Science Band Interferometer');
        
        daspect([1 1 1]);
        axis xy;

        subplot(2,3,4);
        RNG = FOV * [-1 1];
        
        Ipeak = max(Ipeak,max(PSF(:)));
        imagesc(RNG,RNG,(PSF/Ipeak).^(1/2));
        axis xy;
        daspect([1 1 1]);
        title(sprintf('PSF Strehl=%.3f t=%.3f',...
            STREHL(n,i),t));
        xlabel('arcsecs');
        ylabel('arcsecs');
        
        
        subplot(2,3,5);
         if(t>AO_STARTTIME) 
            CCD = CCD + PSF;
         end
         Cpeak = max(Cpeak,max(CCD(:)));
        imagesc(RNG,RNG,(CCD/Cpeak).^(1/2));
        axis xy;
        daspect([1 1 1]);
        title(sprintf('Average \\lambda=%.2g microns',...
            F.lambda*1e6));
        xlabel('arcsecs');
        ylabel('arcsecs');
        
        %%
        subplot(2,3,6);
        maxStrehl = min(1,max(maxStrehl,STREHL(n,i)*1.2));
        plot(TIMES(1:n),STREHL(1:n,i),'b-',...
            TIMES(1:numframes),STREHL(1:numframes,2),'r-',...
            TIMES(1:numframes),STREHL(1:numframes,1),'k-');
        ylim([0 maxStrehl]);
        xlim([0 TIMES(end)]);
        %xlim([-0.25 0]+t);
        ylabel('Strehl');
        xlabel('t (secs)');
        title(sprintf('Strehl History:  gain=%.2f.',...
            gain));
        legend1=legend ('1068 modes', '757 modes','548 modes');
set(legend1,'Location','SouthEast');
        
        
        drawnow;
        clear g;
        %if(n==10);
         fprintf('Time: %3.0f ms, Strehl= %0.3f xpos=%2.4f ypos=%2.4f \n',t*1000,STREHL(n,i),STREHL(n,1),yang(n));
        %end

    end
    STREHLp=STREHL;
end

