%%Script to simulate basic phase offsets of the GMT with a pyramid sensor 
%PMH 100921. Basec on pyramid code in @AOWFS.m by VB


clear;
numframes=40;     %number of frames to generate
WFS_FPS=1000;       %running at  "WFS_FPS" Hz
hFOV = 0.05;  % In arcsecs. Half-width
pix=201;       %pixels in PSF image
 
%Load in a GMT package of files.  We only use the AOaperture, A.
load binaries/GMTAO_dh1246;
    
%%Set up a field
F = AOField(A);
F.lambda = AOField.KBAND;  % Science Wavelength
F.FFTSize = 2048*[1 1]; % How should this be scaled??
%Create an initial PSF with the field (the first one is junk)
PSF = F.mkPSF(hFOV,hFOV/pix*2); %TODO: This is a bug workaround.  FIXME! 

%%Set up plotting parameters
% This is the brightest pixel seen to date.
Ipeak = 0;  
Cpeak=0;
CCD=0;
[x,y]=coords(F);
%% Start the time loop
for n=1:numframes
    t = n/WFS_FPS;     
    %Note transpose to make column vector
    f=25;  %freqency of phase variation
    a=F.lambda/4/2;   %F.lambda corresponds to half wave (??) scaling is not right
    %put in a sinusoidal phase error on one aperture only
    phase=a*sin(2*3.14*t*f);
    PISTONS = [phase 0 0 0 0 0 0]';
    A.setPistons(PISTONS);

    %% Now calculate image quality for the science object
    % Strehl plot setup
    mask = (A.grid>0.5);
    F.planewave*A;
    g = F.grid_;
    STREHL(n) = abs(mean(g(mask)))^2;
    avSTREHL=mean(STREHL(1:n));
 
%%%%%%%%Pyramid calculation from VB in AOWFS%%%%%%%%%%%%%%
    field=F.grid;
    field(1,:)=[]; %get rid of first row to make number of rows even %%%HACK!
    fieldPad=padarray(field,size(field));
    if(mod(size(fieldPad,1), 2)==0 || mod(size(fieldPad,2),2)==0)
        display(' ')
        display('Field has an even number of elements.  FFT2 may produce an error.')
        display(' ')
        %pause
    end
    psfAtTipPad = fftshift(fft2(fieldPad));
    
    Nelq = (size(psfAtTipPad)+1)/2;

    %Now, for the wobbling!
    ww =0;
    w(1,:) = [ ww  ww  -ww -ww];
    w(2,:) = [ ww -ww  -ww  ww];

    ee = size(psfAtTipPad);

    V1 = 1 : (ww + 1);
    V2 = (ee(1) - ww) : ee(1) ;

    for i=1:length(w)

        clear psfAtTip
        psfAtTip = circshift( psfAtTipPad , w(:,i) );

        for j = 1:2
            if w(j,i) < 0
                dd(j,:) = V1; %This was V2.  I don't understand this PMH
            else
                dd(j,:) = V1;
            end
        end

        psfAtTip(dd(1,:),:)=0;
        psfAtTip(:,dd(2,:))=0;

        %Zoom in on quadrant "psfs"
        %figure(2);    
        %subplot(1,2,1);
        %[sizex,sizey]=size(psfAtTipPad);
        %hf=sizex/100;
        %imagesc(abs(psfAtTipPad(sizex/2-hf:sizex/2+hf,sizey/2-hf:sizey/2+hf)));
        %daspect([1 1 1]); title('orig');
        %subplot(1,2,2);
        %[sizex,sizey]=size(psfAtTip);
        %imagesc(abs(psfAtTip(sizex/2-hf:sizex/2+hf,sizey/2-hf:sizey/2+hf)));
        %daspect([1 1 1]); title(sprintf('shifted %d',i));
        %pause

        %split PSF into quadrants
        %make sure each quadrant has an odd number of elements for FFT2?

        clear psfA psfB psfC psfD

        psfA = psfAtTip( 1:Nelq(1)  , 1:Nelq(2) );
        psfB = psfAtTip( 1:Nelq(1)  , Nelq(2):end );
        psfC = psfAtTip( Nelq(1):end , 1:Nelq(2) );
        psfD = psfAtTip( Nelq(1):end , Nelq(2):end );

        if(mod(size(psfA,1), 2)==0 || mod(size(psfA,2),2)==0)
            display(' ')
            display('Quadrant of field has an even number of elements.  FFT2 may produce an error.')
            display(' ')
            %pause;
        end

        
        %psfpyr = [psfA,psfB;psfC,psfD];
        %figure(2)
        %imagesc(log(abs(psfpyr))); daspect([1 1 1]);
        %drawnow;

        %transform back to pupil plane of WFS and downsample 
        %to correct number of subapps
        clear pupilA pupilB pupilC pupilD

        pupilAo = abs(ifft2(psfA)).^2 ; 
        pupilBo = abs(ifft2(psfB)).^2 ;
        pupilCo = abs(ifft2(psfC)).^2 ;
        pupilDo = abs(ifft2(psfD)).^2 ;

        %chop off extra pupil created during initial padding
        %chop off slightly less than 1/3 on each side b/c
        %mask has dark border - use scaling factor to adjust
        %sc = 31/29;  
        nn= ceil(size(pupilAo)/3-.01*size(pupilAo,1));
        %n= floor(0.5 * ( size(pupilAo) - ceil( sc*size(pupilAo)/3 ) ));
        pupilA = pupilAo(nn(1):end-nn(1)+1 , nn(2):end-nn(2)+1 );
        pupilB = pupilBo(nn(1):end-nn(1)+1 , nn(2):end-nn(2)+1 );
        pupilC = pupilCo(nn(1):end-nn(1)+1 , nn(2):end-nn(2)+1 );
        pupilD = pupilDo(nn(1):end-nn(1)+1 , nn(2):end-nn(2)+1 );


        %pupilpyr = [pupilA,pupilB;pupilC,pupilD];
        %figure(2)
        %imagesc(pupilpyr); daspect([1 1 1]);
        %drawnow;
        
        %PMH calculate slopes and display at nominal resolution
        %Sx=(pupilA-pupilB+pupilC-pupilD);
        %Sy=(pupilA-pupilC+pupilB-pupilD);

        %pupilslopes = [Sx,Sy];
        %figure(2);
        %imagesc(pupilslopes); daspect([1 1 1]);
        %title(sprintf('time = %2.0f ms, phase=%.3f waves',t*1e3,phase*1e6*0.5));
        %drawnow;
        %pause;
        
        %interpolate to final WFS sampling
        xp=linspace(1,size(pupilA,2),length(WFS.masked));
        yp=linspace(1,size(pupilA,1),length(WFS.masked));   
        %xp=linspace(1,size(pupilA,1),length(WFS.masked)-1);
        %yp=linspace(1,size(pupilA,2),length(WFS.masked)-1);
        [xxp,yyp]=meshgrid(xp,yp);
        %%%%%%%%%%%%%%%%%%%%%
        %FIX ME!!!
        %%%%%%%%%%%%%%%%%%%%

        %Don't use interp to downsample.  It will just find best fit point at the location you
        %request, not best fit value to area between points Need to upsample to
        %an integer ratio of the final size, then bin.  
        %Ex: Final Nsubap=25.  PupilA=91x91.  Upsample to 100x100 and do 4x4 binning (summing).  


        pupilAi(:,:,i) = interp2(pupilA, xxp,yyp);
        pupilBi(:,:,i) = interp2(pupilB, xxp,yyp);
        pupilCi(:,:,i) = interp2(pupilC, xxp,yyp);
        pupilDi(:,:,i) = interp2(pupilD, xxp,yyp);


    end

    pupilAint = sum(pupilAi,3);
    pupilBint = sum(pupilBi,3);
    pupilCint = sum(pupilCi,3);
    pupilDint = sum(pupilDi,3);

    figure(2);
    pupilint = [ pupilAint , pupilBint ; pupilCint, pupilDint ];
    imagesc(pupilint); daspect([1 1 1]);
    drawnow;


    %Then normalize the intensity in each pupil by the ratio of the 
    % (sum of all 4 pupils) / sum of original field intensity.
    %Do this so that photon noise is correct.

    fieldSum = sum(sum(abs(fieldPad)));
    pupilSumI = sum(sum(pupilAint+ pupilBint + pupilCint + pupilDint));
    pupilNorm = fieldSum / pupilSumI;

    
    %pupilAint = pupilAint*pupilNorm;
    %pupilBint = pupilBint*pupilNorm;
    %pupilCint = pupilCint*pupilNorm;
    %pupilDint = pupilDint*pupilNorm;
    
    
    %PMH calculate slopes at nominal resolution
    Sx=(pupilAint-pupilBint+pupilCint-pupilDint);
    Sy=(pupilAint-pupilCint+pupilBint-pupilDint);


    %pupilslopes = [Sx,Sy];
    %figure(2);
    %imagesc(pupilslopes, [-60 60]); daspect([1 1 1]);
    
    %title(sprintf('time = %2.0f ms, phase=%.3f waves',t*1e3,phase/F.lambda*4)); %why *4?
    %colorbar;
    %drawnow;
% This saves the current picture as a JPEG.
    filename = sprintf('/tmp/FRAME_%04d.jpg',n);
    rez = 160;
    resolution = sprintf('-r%d',rez);
    print(resolution,'-djpeg',filename);


%%%End of pyramid calculation from VB


    
    

    %% Plot some interesting pictures...
    figure(1);
    clf; % Clear the Figure each time step

    %Instantaneous PSF
    xa=1:pix+2;
    angleperpix=2*hFOV/(pix);
    PSF = F.mkPSF(hFOV,angleperpix);
    subplot(1,2,1);
    RNG = hFOV * [-1 1];
    Ipeak = max(PSF(:));
    imagesc(RNG,RNG,(PSF/Ipeak).^(1/2));
    axis xy;
    daspect([1 1 1]);
    title(sprintf('Strehl=%.3f \nInstant. PSF',STREHL(n)));
    xlabel('arcsecs');
    ylabel('arcsecs');

  
    %Plot residual phase at science wavelength
    subplot(1,2,2);
    [x,y]=F.coords;
    imagesc(x,y,angle(g));  %angle is atan2(imag(g)/real(g))
    title(sprintf('time = %2.0f ms, phase=%.3f um',t*1e3,phase*1e6*4));  %*4???
    daspect([1 1 1]);
    axis xy;

    drawnow; 
    clear g;
end
%    fitswritesimple(CCD,FITSfile);

%% Movie creation...
% Run this command after it is done to create the movie...
% mencoder "mf:///tmp/FRAME_*.jpg" -mf fps=10 -o MOVIE.avi -ovc lavc -lavcopts vcodec=msmpeg4v2
% system('mencoder "mf:///tmp/FRAME_*.jpg" -mf fps=30 -o MOVIE_automake.avi
% -ovc lavc -lavcopts vcodec=msmpeg4v2');

%encode in msmpeg4v2 rather than wmv1 to allow use by macs and windows.
