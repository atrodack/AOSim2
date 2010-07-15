%%Script demonstrating the setup and operation of the LBT AO simulation
%This is a modified version of the run_myGMT scripts by JLC in demos
%PMH 20100327
clear;
numframes=50;     %number of frames to generate
numtrials=3;  %number of time to do simulation with same random seed
STREHL = zeros(numframes,numtrials);

for i=1:numtrials %Run complete program multiple times
    if (i==1)  %Run each iteration with exaxt same random number seed
        defaultStream = RandStream.getDefaultStream;
        savedState = defaultStream.State;
    else
        defaultStream.State=savedState
    end
    %%Set parameters likely to be fiddled with
    
    WFS_FPS=1000;       %running at  "WFS_FPS" Hz
    AO_STARTTIME=0.001; %Starting time of AO correction
    
    gain = 0.4*i;         % Fraction of error to feedback
    SCIENCE_WAVELENGTH = AOField.KBAND;
    gsangle=0;           %Angle of guide star from axis in arcsec
    windfactor=1;       % for slowing wind down to check lag
    r0factor=1;         % for increasing r0 to check turbulence effects
    hFOV=0.2;            %half FOV of PSF in arcsec
    subapps=30;          %pyr WFS number of subapps
    REBUILD=0;          %Should reconstructor be generated first?
    ORDER = 2;
    fprintf('maximum modes correctable is ~%3.0f \n',subapps^2/4*pi);

   
    if(REBUILD) 
    disp('Rebuilding')
    NModes=round((ORDER*3.14159)*2+1)*ORDER;
        %This is based on the make_the_GMT_AO script in demos by JLC
        % Create aperture definitions
        load data/PLBT_oneEye.mat
        Seg = AOSegment;
        Seg.name = 'LBT Primary';
        Seg.pupils = PLBT;
        Seg.make;

        AL = AOAperture;
        AL.name = 'LBT Left';
        AL.addSegment(Seg);
        DM = AODM(AL);
        DM.name = 'LBT Hexapolar Adaptive Secondary';
        %% Use actuators from GMT positions  
        load GMT_672act_segment_PMH.mat;
        %Add them to the DM
        DM.addActs(GMT_Actuators)
        %Define a radius of 10 m and 8azimuthal points as boundary
        DM.defineBC(20,30);
        
        %% Build the WFS.
        d = 8.25;
        BB=AL.BBox;
        %D= mean(BB(2,:)-BB(1,:));
        WFS = AOWFS(AL,d/subapps);
        WFS.name = 'LBT Pyramid WFS';
        %% Now for some real work.  Building the RECONSTRUCTOR...
        RECON = AOReconstructor(AL,DM,WFS);
        %Set Scutoff high to make the cutoff be whatever the modes explored is.
        RECON.Scutoff=30000;
        show=false;
        %RECON.dhprogram(d,ORDER,show);
        order=12;
        RECON.zprogram(d,order,show);
        matfile = sprintf(['LBTAOzorder%03.0f'],order);
        save (matfile, 'AL', 'DM', 'WFS', 'RECON', 'd') % 
    else  %Load in pre-made reconstructor
       %load binaries/LBTAOdh104
       load LBTAOzorder012
    end

    %% Now Generate the Phase screens 
    %Define the Atmosphere model and winds aloft.
    %From van Dam specifications
    fractionlow=0.7;
    fractionhigh=0.3;
    R_0=0.15;
    rlow=R_0*(1+fractionhigh/fractionlow)^(3/5);
    rhigh=rlow*(fractionlow/fractionhigh)^(3/5);
    
    ATMO = AOAtmo(AL);
    WFlow = AOScreen(1024,rlow*r0factor,550e-9);     %
    WFhigh = AOScreen(2048,rhigh*r0factor,550e-9);    %
    ATMO.addLayer(WFlow,300);
    ATMO.addLayer(WFhigh,5000);
    ATMO.layers{1}.Wind = [10 0]*windfactor;           
    ATMO.layers{2}.Wind = [3.5 6.06]*windfactor;     %7 m/s at 60 deg.angle
    WFlow.name = 'Lower altitude turbulence';
    WFhigh.name = 'High altitude turbulence';
    r0 = ATMO.totalFriedScale;
    th_scat = AOField.VBAND/r0*206265;
    fprintf('The total r0 is %.2f cm.\n',100*ATMO.totalFriedScale);
    fprintf('The seeing is %.2f arcsecs.\n',th_scat);

    %The ATMO.BEACON is defined in m in x and y, and z from the primary
    %Examples
    beaconheight=1e7;
    xstar=beaconheight*gsangle/206265;
    STAR = [xstar 0 beaconheight];
    SC_OBJECT=[0 0 beaconheight];
    ATMO.BEACON = STAR; % Set this so ATMO knows how to compute the wavefront.
    ATMO.GEOMETRY = false;  %false=dynamic focus of BEACON
    GUIDE_STAR=STAR;        %Position of beacon used for AO correction
    SCIENCE_OBJECT=SC_OBJECT;    %Position of object used for image quality estimate

    %%Now Set up the field objects
    % Create the WFS and Science AOField objects...
    Fwfs = AOField(AL);
    % The Reconstructor was calibrated at a certain wavelength.
    Fwfs.lambda = RECON.lambda;  

    F = AOField(AL);
    F.lambda = SCIENCE_WAVELENGTH;  % Science Wavelength
    F.FFTSize = 1024*[1 1]; 
    %PSF = F.mkPSF(.2,.2/100); %TODO: This is a bug workaround.  FIXME!

    %%Set up plotting parameters
    % This is the brightest pixel seen to date.
    Ipeak = 0;  
    Cpeak=0;
    [x,y]=coords(F);
    xvals=padarray(x,length(y)-1,'replicate','post');
    yvals=padarray(y',[0,length(x)-1],'replicate','post');
    
    TIMES = (1:numframes)/WFS_FPS;

    % Strehl plot setup
    mask = (AL.grid>0.5);
    xang = zeros(numframes,1);
    yang = zeros(numframes,1);
    maxStrehl = 0.3;
    CCD = 0;
    Dap = 8.4;
    sigmafit=AOField.VBAND/SCIENCE_WAVELENGTH.*sqrt(0.2944*(RECON.Nmodes^(-sqrt(3)/2))*(Dap/r0)^(5/3));
    sigmatime=AOField.VBAND/SCIENCE_WAVELENGTH.*sqrt(0.3*(30/r0/WFS_FPS)^(5/3));
    Strehlfit=exp(-(sigmafit^2));
    Strehltime=exp(-(sigmatime^2));
    fprintf('Expected Fitting error %6.2f radians, Strehl =%2.2f\n',sigmafit,Strehlfit);
    fprintf('Expected Temporal error %6.2f radians, Strehl =%2.2f\n',sigmatime,Strehltime);
    fprintf('D=%f   Num. modes=%d WFS wavelength=%6.2e\n', ...
             Dap,RECON.Nmodes,RECON.lambda);

    DM.zero;     

    %% Start the time loop
    for n=1:numframes
    % for n = 1:1
        t = n/WFS_FPS;
        %go from -t/2 to t/2 to center phase screens 
        ATMO.time = t-numframes/WFS_FPS/2.0;      
        %% This is the guts of the AO closed-loop integrating servo....
        ATMO.BEACON = GUIDE_STAR;
        %fprintf('Sensing WF\n');
        %WFS.sense(Fwfs.planewave*ATMO*AL*DM);
        WFS.sensePyramid(Fwfs.planewave*ATMO*AL*DM);
        if(t>AO_STARTTIME)  % Suffer with seeing limit until AO_STARTTIME.
            if ((t-AO_STARTTIME)<0.005)
                ggain=gain*10;
            else
                ggain=gain;
            end
                
            DM.bumpActs(-ggain*RECON.RECONSTRUCTOR * WFS.slopes);
            DM.removeMean;
        end  %End of AO servo

        %% Now calculate image quality for the science object
        ATMO.BEACON = SCIENCE_OBJECT;
        F.planewave*ATMO*AL*DM;  
        g = F.grid_;
        STREHL(n,i) = abs(mean(g(mask)))^2;
        
        %xslope=angle(g)./xvals;
        %yslope=angle(g)./yvals;
        %xslope(isinf(xslope))=0;
        %yslope(isinf(yslope))=0;
        %xslope(isnan(xslope))=0;
        %yslope(isnan(yslope))=0;
   
        %xang=mean(mean(xslope));
        %yang=mean(mean(yslope));


%% Now calculate image quality for the science object
        %fprintf('Calculating science image\n');
        ATMO.BEACON = SCIENCE_OBJECT;
        F.planewave*ATMO*AL*DM;  

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
     figure(1);
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
        plot(TIMES(1:n),STREHL(1:n,i),'k-',...
            TIMES(1:numframes),STREHL(1:numframes,2),'r-',...
            TIMES(1:numframes),STREHL(1:numframes,1),'b-');
        ylim([0 maxStrehl]);
        xlim([0 TIMES(end)]);
        ylabel('Strehl');
        xlabel('t (secs)');
        title(sprintf('Strehl History:  gain=%.2f. iter=%1.1f',...
            gain,i));
        %legend1=legend ('1068 modes', '757 modes','548 modes');
        %set(legend1,'Location','SouthEast');
        
        
        %Instantaneous PSF
        subplot(2,3,4);
        RNG = hFOV * [-1 1];  
        Ipeak = max(Ipeak,max(PSF(:)));
        imagesc(RNG,RNG,(PSF/Ipeak).^(1/2));
        axis xy;
        daspect([1 1 1]);
        title(sprintf('PSF Strehl=%.3f t=%.3f',...
            STREHL(n,i),t));
        xlabel('arcsecs');
        ylabel('arcsecs');
        
        %Average PSF
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
        
        %Plot cut through PSF
        subplot(2,3,6);
        if(t>AO_STARTTIME)
            center=floor(pix/2);
            cut=CCD(center,:);
            an=(xa-(pix/2))*2*hFOV/pix;
            plot(an,cut);
            xlim([-hFOV hFOV]);
            xlabel('arcsecs');   
        end    
        
        drawnow;
        
        clear g;
        %if(n==10);
         fprintf('Time: %3.0f ms, Strehl= %0.3f xpos=%2.4f ypos=%2.4f \n' ...
             ,t*1000,STREHL(n,i),xang(n),yang(n));
        %end

    end
    STREHLp=STREHL;       

end