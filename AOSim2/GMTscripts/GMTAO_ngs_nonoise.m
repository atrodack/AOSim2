%%Script to simulate initial NGS comparison
% Simulation parameters are from specification defined by van Dam
%PMH 20100315

clear;
numframes=200;     %number of frames to generate
WFS_FPS=1000;       %running at  "WFS_FPS" Hz
AO_STARTTIME=0.001; %Starting time of AO correction
PHASE_START=0.001;  %Starting time of phase correction
PHASE_END=numframes/WFS_FPS;
hFOV = 0.05;  % In arcsecs. Half-width
numtrials=1;



STREHL = zeros(numframes,10);
%STREHLp= zeros(numframes,1);


for i=1:numtrials  %Run complete program multiple times
    %%Set parameters likely to be fiddled with
    
    if (i==1)  %Run each iteration with exact same random number seed
        defaultStream = RandStream.getDefaultStream;
        savedState = defaultStream.State;
    else
        defaultStream.State=savedState
    end
    
    load binaries/GMTAO_dh1246
    

    %Gain on the error feedback in the AO loop
    gain=1;
    %phase error gain
    pgain=0.1;

    GAMMA = 2;  % This is the gamma term for the PSF grayscale stretch.
    SCIENCE_WAVELENGTH = AOField.KBAND;

        
        
    DM.setActs(0);
    A.trueUp;
    numact=length(DM.actuators);
    fprintf('Number of actuators is %d \n',numact);
    
    
    
    

    %% Now Generate the Phase screens 
    %Define the Atmosphere model and winds aloft.
    %Atmospheric parameters.
    
    %From van Dam specifications
    fractionlow=0.7;
    fractionhigh=0.3;
    R_0=0.15;
    rlow=R_0*(1+fractionhigh/fractionlow)^(3/5);
    rhigh=rlow*(fractionlow/fractionhigh)^(3/5);
    
    ATMO = AOAtmo(A);
    WFlow = AOScreen(1024,rlow,500e-9);     %
    WFhigh = AOScreen(2048,rhigh,500e-9);    %
    ATMO.addLayer(WFlow,300);
    ATMO.addLayer(WFhigh,5000);
    ATMO.layers{1}.Wind = [10 0]*1;           %
    ATMO.layers{2}.Wind = [3.5 6.06]*1;           % 7 m/s at 60 deg. angle
    WFlow.name = 'Lower altitude turbulence';
    WFhigh.name = 'High altitude turbulence';
    r0 = ATMO.totalFriedScale;
    th_scat = AOField.VBAND/r0*206265;
    fprintf('The total r0 is %.2f cm.\n',100*ATMO.totalFriedScale);
    fprintf('The seeing is %.2f arcsecs.\n',th_scat);
    %TODO: Calculate t_0 from phase screen

    %The ATMO.BEACON is defined in meters in x and y, and z from the
    %primary
    beaconheight=1e7;
    starangle=30*(3-i)*0;  %angle in arcsec
    xstar=starangle/206265*beaconheight;
    STAR = [xstar 0 beaconheight];
    SCIENCE=[0 0 beaconheight];
    ATMO.BEACON = STAR; % Set this so ATMO knows how to compute the wavefront.
    ATMO.GEOMETRY = false;  %false=dynamic focus of BEACON

    GUIDE_STAR=STAR;        %Position of beacon used for AO correction
    SCIENCE_OBJECT=SCIENCE;    %Position of object used for image quality estimate

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
            %Segment corrections
            %DM.bumpActs(-sgain*RECONseg.RECONSTRUCTOR * WFS.slopes);
            %Global correction
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
        
        %Calculate centroid and FWHM from moments
        pix=200;
        angleperpix=2*hFOV/(pix);
        PSF = F.mkPSF(hFOV,angleperpix);
        xa=1:pix+1;
        PSF(isnan(PSF))=0;
        M00 = sum(PSF(:));
        M10=sum(PSF*xa')/M00;
        M01=sum(xa*PSF)/M00;
        xa2=(xa'-M10).^2;
        ya2=(xa-M01).^2;
        M20 = sum(PSF*xa2)/M00;
        M02 = sum(ya2*PSF)/M00;
        FWHM=sqrt((M20+M02)/2)*2.35482*angleperpix;
        

     %% Plot some interesting pictures...
        clf; % Clear the Figure each time step
        
        %Plot DM shape
        subplot(2,3,1);
        imagesc(x,y,DM.grid .* mask); daspect([1 1 1]);
        axis xy;
        title('DM Shape');
        
        %Plot residual phase at science wavelength
        subplot(2,3,2);
        [x,y]=F.coords;
        imagesc(x,y,angle(g));  %angle is atan2(imag(g)/real(g))
        title('Residual Phase');
        
        daspect([1 1 1]);
        axis xy;
        
        
        %Plot Strehl Variation
        subplot(2,3,3);
        maxStrehl = min(1,max(maxStrehl,STREHL(n,i)*1.2));
        plot(TIMES(1:n),STREHL(1:n,i),'b-',...
            TIMES(1:numframes),STREHL(1:numframes,2),'r-',...
            TIMES(1:numframes),STREHL(1:numframes,1),'k-');
        ylim([0 maxStrehl]);
        xlim([0 TIMES(end)]);
        ylabel('Strehl');
        xlabel('t (secs)');
        title(sprintf('Strehl History'));
        %legend1=legend ('1068 modes', '757 modes','548 modes');
        %set(legend1,'Location','SouthEast');
        
        
        %Instantaneous PSF
        subplot(2,3,4);
        RNG = hFOV * [-1 1];  
        Ipeak = max(Ipeak,max(PSF(:)));
        imagesc(RNG,RNG,(PSF/Ipeak).^(1/2));
        axis xy;
        daspect([1 1 1]);
        title(sprintf('Instant. PSF Strehl=%.3f t=%.3f',...
            STREHL(n,i),t));
        xlabel('arcsecs');
        ylabel('arcsecs');
        
        %Plot cut through PSF
        subplot(2,3,5);
            center=floor(pix/2);
            cut=PSF(center,:);
            an=(xa-(pix/2))*2*hFOV/pix;
            plot(an,cut);
            xlim([-hFOV hFOV]);
            xlabel('arcsecs'); 
            title(sprintf('xcut through PSF'));
        
        
        %Average PSF
        subplot(2,3,6);
        if(t>AO_STARTTIME) 
            CCD = CCD + PSF;
        end
        Cpeak = max(Cpeak,max(CCD(:)));
        imagesc(RNG,RNG,(CCD/Cpeak).^(1/2));
        axis xy;
        daspect([1 1 1]);
        title(sprintf('Average PSF  \\lambda=%.2g microns',...
            F.lambda*1e6));
        xlabel('arcsecs');
        ylabel('arcsecs');        
        
           
        
        drawnow;
        
        clear g;
        %if(n==10);
         fprintf('Time: %3.0f ms, Strehl= %0.3f xpos=%2.4f ypos=%2.4f FWHM=%2.4f \n' ...
             ,t*1000,STREHL(n,i),xang(n),yang(n),FWHM);
        %end
        
           %% This saves the current picture as a JPEG.
     filename = sprintf('/tmp/FRAME_%04d.jpg',n);
     rez = 160;
 
     resolution = sprintf('-r%d',rez);
     print(resolution,'-djpeg',filename);
	

    end
    STREHLp=STREHL;
end

%% Movie creation...
% Run this command after it is done to create the movie...
% mencoder "mf:///tmp/FRAME_*.jpg" -mf fps=10 -o MOVIE.avi -ovc lavc -lavcopts vcodec=wmv1
% system('mencoder "mf:///tmp/FRAME_*.jpg" -mf fps=30 -o MOVIE_automake.avi
% -ovc lavc -lavcopts vcodec=wmv1');
