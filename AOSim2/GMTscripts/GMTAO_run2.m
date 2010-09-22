%%Script to simulate initial NGS comparison
% Simulation parameters are from specification defined by van Dam
%This script was used to generate NGS initial modeling results for GMT
%PMH 20100315

clear;
numframes=20;     %number of frames to generate
WFS_FPS=1000;       %running at  "WFS_FPS" Hz
AO_STARTTIME=0.001; %Starting time of AO correction
PHASE_START=0.002;  %Starting time of phase correction
%PHASE_END=0.1;
PHASE_END=numframes/WFS_FPS;
hFOV = 0.05;  % In arcsecs. Half-width
pix=201;       %pixels in PSF image
numtrials=1;
GRAPHICS=1;  



STREHL = zeros(numframes,10);
avSTREHL=0;
%STREHLp= zeros(numframes,1);


for i=1:numtrials  %Run complete program multiple times
    %%Set parameters likely to be fiddled with
    
    if (i==1)  %Run each iteration with exact same random number seed
        defaultStream = RandStream.getDefaultStream;
        savedState = defaultStream.State;
    else
        defaultStream.State=savedState
    end
    
    if (i==1)
        PHASE_START=0.011;
    else
        PHASE_START=0.4;
    end
    %if (i==1)
    %    load binaries/GMTAO_dh630;
    %elseif (i==2)
    %    load binaries/GMTAO_dh924;
    %else
        load binaries/GMTAO_dh1246;
    %end

    %Gain on the error feedback in the AO loop
    gain=1;
    %phase error gain
    pgain=0.2;

    GAMMA = 2;  % This is the gamma term for the PSF grayscale stretch.
    SCIENCE_WAVELENGTH = AOField.KBAND;
    
    %output file
    FITSfile='test.fits';

        
        
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
    PSF = F.mkPSF(hFOV,hFOV/pix*2); %TODO: This is a bug workaround.  FIXME!
    

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
    
    PLooptime=0;
    CLooptime=0;
    %% Start the time loop
    for n=1:numframes
        
        TLoop=tic; %time complete loop
        t = n/WFS_FPS;
        %go from -t/2 to t/2 to center phase screens 
        ATMO.time = t-numframes/WFS_FPS/2.0;      
        %% This is the guts of the AO closed-loop integrating servo....
        ATMO.BEACON = GUIDE_STAR;
        
        GLoop=tic; %time generating atmosphere 
        ftemp=Fwfs.planewave*ATMO*A*DM;
        GLooptime=toc(GLoop);
        SLoop=tic; %time sensing
        WFS.sense(ftemp);
        SLooptime=toc(SLoop);
        if(t>AO_STARTTIME)  % Suffer with seeing limit until AO_STARTTIME.
            %Segment corrections
            %DM.bumpActs(-sgain*RECONseg.RECONSTRUCTOR * WFS.slopes);
            %Global correction
            CLoop=tic; %time correction
            DM.bumpActs(-gain*RECON.RECONSTRUCTOR * WFS.slopes);
            DM.removeMean;
            CLooptime=toc(CLoop);
                    
            if ((t>PHASE_START)&&(t<PHASE_END))
            %Sense and correct segment pistons
            %PISTONS = WFS.magicPistonSensor(Fwfs,A);
            %A.bumpPistons(- pgain*PISTONS);	
                PLoop=tic;
                PISTONS = WFS.magicPistonSensor(F,A);
                A.bumpPistons(-pgain*PISTONS);
                PLooptime=toc(PLoop);
              
            end %End of Phase servo 
        
            
        end  %End of AO servo
        
        if (GRAPHICS)
            %% Now calculate image quality for the science object
            %fprintf('Calculating science image\n');
            ATMO.BEACON = SCIENCE_OBJECT;
            F.planewave*ATMO*A*DM;
            
            g = F.grid_;
            STREHL(n,i) = abs(mean(g(mask)))^2;
            
            if (t>AO_STARTTIME)
                avSTREHL=mean(STREHL(AO_STARTTIME*WFS_FPS+1:n,i));
            else
                avSTREHL=0;
            end
             
            
           
           
            %% Plot some interesting pictures...
            figure(1);
            clf; % Clear the Figure each time step
            
            %Plot DM shape
            subplot(2,3,1);
            imagesc(x,y,DM.grid .* mask); daspect([1 1 1]);
            axis xy;
            title(sprintf('t=%.3f \nDM Shape',t));
            
            %Instantaneous PSF
            
            xa=1:pix+2;
            angleperpix=2*hFOV/(pix);
            PSF = F.mkPSF(hFOV,angleperpix);
            
            subplot(2,3,2);
            RNG = hFOV * [-1 1];
            Ipeak = max(PSF(:));
            imagesc(RNG,RNG,(PSF/Ipeak).^(1/2));
            axis xy;
            daspect([1 1 1]);
            title(sprintf('Strehl=%.3f \nInstant. PSF',STREHL(n,i)));
            xlabel('arcsecs');
            ylabel('arcsecs');
            
            %Average PSF
            subplot(2,3,3);
            if(t>AO_STARTTIME)
                CCD = CCD + PSF;
            end
            Cpeak = max(CCD(:));
            imagesc(RNG,RNG,(CCD/Cpeak).^(1/2));
            axis xy;
            daspect([1 1 1]);
            title(sprintf(' avg. Strehl=%.3f \nAverage PSF  ',...
                avSTREHL));
            xlabel('arcsecs');
            ylabel('arcsecs');
            
            
            
            
            %Plot residual phase at science wavelength
            if ((t>PHASE_START)&&(t<PHASE_END))
                phaseStatus='on';
            else
                phaseStatus='off';
            end
            
            subplot(2,3,4);
            [x,y]=F.coords;
            imagesc(x,y,angle(g));  %angle is atan2(imag(g)/real(g))
            title(sprintf('Seg. phasing = %s \n Residual Phase',phaseStatus));
            
            daspect([1 1 1]);
            axis xy;
            
            
            
            %Plot cut through PSF
            subplot(2,3,5);
            center=floor(pix/2);
            cut=PSF(center,:)/Ipeak*STREHL(n,i);
            an=(xa-(pix/2))*2*hFOV/pix;
            plot(an,cut);
            xlim([-hFOV hFOV]);
            ylim([0 1]);
            xlabel('arcsecs');
            title(sprintf('\\lambda=%.2g \\mum \n x-cut through inst. PSF', ...
                F.lambda*1e6));
            
            
            %Plot cut through PSF
            subplot(2,3,6);
            if(t>AO_STARTTIME)
                center=floor(pix/2);
                cut=CCD(center,:)/Cpeak*avSTREHL;
                an=(xa-(pix/2))*2*hFOV/pix;
                plot(an,cut);
                xlim([-hFOV hFOV]);
                ylim([0 1]);
                xlabel('arcsecs');
                title(sprintf('x-cut through avg. PSF'));
            end
            
            
            
            
              
            %% This saves the current picture as a JPEG.
            filename = sprintf('/tmp/%1dFRAME_%04d.jpg',i,n);
            rez = 160;
            resolution = sprintf('-r%d',rez);
            print(resolution,'-djpeg',filename);
            
            
            
            figure(2);
            %Plot Strehl Variation
            
            maxStrehl = min(1,max(maxStrehl,STREHL(n,i)*1.2));
            plot(TIMES(1:n),STREHL(1:n,i),'b-',...
                TIMES(1:numframes),STREHL(1:numframes,2),'r-',...
                TIMES(1:numframes),STREHL(1:numframes,1),'k-');
            ylim([0 maxStrehl]);
            xlim([0 TIMES(end)]);
            ylabel('Strehl');
            xlabel('t (secs)');
            title(sprintf('Strehl History'));
            legend1=legend ('1246 modes', '924 modes','630 modes');
            set(legend1,'Location','SouthEast');
          
            drawnow;
            
        end

        clear g;
        %if(n==10);
        
    
        
        
        
        Looptime=toc(TLoop);  %Timing of sim. loop speed
        fprintf('T=%3.0f ms, S=%0.3f avS=%2.3f ' ...
             ,t*1000,STREHL(n,i),avSTREHL);
         fprintf('Gen=%2.2f Sense=%0.2f, Corr.= %0.3f, Phase=%0.2f s  Tot=%0.3f\n' ...
             ,GLooptime,SLooptime,CLooptime,PLooptime,Looptime);
        %end
        
   
	

    end
    STREHLp=STREHL;
    fitswritesimple(CCD,FITSfile);
end

%% Movie creation...
% Run this command after it is done to create the movie...
% mencoder "mf:///tmp/FRAME_*.jpg" -mf fps=10 -o MOVIE.avi -ovc lavc
% -lavcopts vcodec=msmpeg4v2
% system('mencoder "mf:///tmp/FRAME_*.jpg" -mf fps=30 -o MOVIE_automake.avi
% -ovc lavc -lavcopts vcodec=msmpeg4v2');

%encode in msmpeg4v2 rather than wmv1 to allow use by macs and windows.
