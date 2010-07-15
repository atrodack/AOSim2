close;clf;clear classes
colormap(gray);

load data/TestPupil.mat
SEG0 = AOSegment;
SEG0.name = 'Center Segment';
SEG0.pupils = PAEOS;
SEG0.make;

SEG = AOSegment;
SEG.name = 'Outer Segment';
SEG.pupils = PAEOS(1,:);
SEG.make;

% show(SEG)

A = AOAperture;
A.name = 'Toy GMT';
A.addSegment(SEG0)
for n=1:6
    A.addSegment(SEG,[cos(2*pi*n/6) sin(2*pi*n/6)]*3.7);
end

show(A);

% CH=256+(-31:32);
CH=256+(-63:64);
HALO=A.fft(512);
imagesc(log10(normalize(abs(interp2(HALO(CH,CH),1)).^2)),[-6 0]);

TT = cumsum(randn([7,2,100]),1)*1e-8;
PISTONS = cumsum(randn([7,100]),2)*1e-8;

maxHalo = max(HALO(:));
HALO=A.fft(512);
imagesc(log10(abs(interp2(HALO(CH,CH)/maxHalo,1).^2)),[-6 0]);
sqar;
% 
for nt=1:100
	%     for s=1:7
	%         A.segList{s}.tiptilt = TT(nt,:,s);
	%     end
	% touch(A);
    
	A.setTipTilts(TT(:,:,nt));
	A.setPistons(PISTONS(:,nt));
	
	subplot(1,2,1);
	show(A); 
%     sqar; 
	
	subplot(1,2,2);
	HALO = A.fft(512);
	imagesc(log10(abs(interp2(HALO(CH,CH)/maxHalo,1).^2)),[-4 0]);
	sqar;

	drawnow;
end

lambda = AOField.HBAND;

PS = AOScreen(512,1.5);
PS.make;

F = AOField(A);

% CCD = 0;
% maxHalo = 0;
% for nt=1:5
% 	%     for s=1:7
% 	%         A.segList{s}.tiptilt = TT(nt,:,s);
% 	%     end
% 	% touch(A);
%     
% 	A.setTipTilts(TT(:,:,nt));
% 	A.setPistons(PISTONS(:,nt));
% 	
% 	F.planewave*A*PS;
% 	
% % 	subplot(1,2,1);
% % 	show(F); sqar; 
% 	
% % 	subplot(1,2,2);
% 	HALO = F.fft(512);
% 	HALO = HALO(CH,CH);
% % 	maxHalo = max(maxHalo,max(HALO(:)));
% % 	imagesc(log10(abs(interp2(HALO/maxHalo,1).^2)),[-4 0]);
% % 	sqar;
% 
% 	CCD = CCD + abs(HALO).^2;
% 	fprintf('.');
% 	drawnow;
% end
% fprintf('\n');
% 
% imagesc(log10(normalize(CCD)),[-4 0]);sqar;


for nt=1:size(TT,3)
    A.setTipTilts(TT(:,:,nt));
    A.setPistons(PISTONS(:,nt));
    
    subplot(1,2,1);
    F.lambda = AOField.VBAND;
    imagesc(mag2(F.planewave*A+0.5),[0 2.5]);sqar;axis off;
    
    subplot(1,2,2);
    F.lambda = AOField.HBAND;
    F.planewave*A;HALO=F.fft(512);imagesc(log10(normalize(abs(interp2(HALO(CH,CH),1).^2))),[-6 0]);sqar;drawnow;
end


clf;
% colormap(jet);

A.trueUp;


ATMO = AOAtmo(A);
% WF1.altitude = 1000;
WFlow = AOScreen(1024,0.15,500e-9);
WFhigh = AOScreen(2048,0.20,500e-9);
% WF1.altitude = 8000;
ATMO.addLayer(WFlow,1000);
ATMO.addLayer(WFhigh,8000);
% ATMO.addLayer(WFlow,1000).addLayer(WFhigh,8000);
% ATMO.addLayer(AOScreen(512,0.12,500e-9));
% ATMO.addLayer(AOScreen(512,0.15,500e-9));

ATMO.GEOMETRY = false;

BEACON = [0 1/206265 1] * 90e3;
ATMO.BEACON = BEACON;

DM = AODM(A);
DM.defineBC(30,8);

% x = -6:2/3:6;
x = -6:0.5:6;
[X,Y] = meshgrid(x,x);
Xact = X(:); Yact = Y(:);
Aact = A.interpGrid(Xact,Yact);
SEL = Aact>0.05;
Xact = Xact(SEL); Yact = Yact(SEL);
DM.addActs([Xact Yact]);
clear Xact Yact Aact SEL

A.show;
DM.plotActuators;


% DM.setDM(ATMO);
DM.setActs(ATMO);
DM.render;
% DM.addFocus(-400e3); % focus on LEO. (TODO: Do this on the actuators so this is automatic.)

F.lambda = AOField.HBAND;
(F.planewave*A*DM)*ATMO;
HALO=F.fft;
imagesc(abs(interp2(HALO,1)),[0 5]);
axis square;
colorbar;axis off;
drawnow;

CH=256+(-63:64);

CCD = 0;
% for x = -500:150:500
% for x = -5:0.5:5
% % 	for y = -5:0.5:5
% % 	for y = -200:50:200
% 	for y=x 
% 		ATMO.BEACON = [x,y,400e3];
% 		ATMO.GEOMETRY = true;
% 		(F.planewave*ATMO*A*DM);
% % 		(F.planewave*A)*ATMO;
% 		HALO=F.fft;
% % 		CCD = CCD + abs(interp2(HALO,1)).^2;
% % 		CCD = CCD + abs(HALO(CH,CH)).^2;
% % 		CCD = CCD + abs(interp2(HALO(CH,CH),1)).^2;
% % 		CCD = CCD + abs(HALO).^2;
% 		CCD = CCD + abs(interp2(HALO,1)).^2;
% % 		imagesc(abs(interp2(HALO,1)),[0 5]);
% 		imagesc(CCD);
% 		axis square;
% 		colorbar;
% % 		axis off;
% 		drawnow;
% 	end
% end

x = -5;
y = -5;

ATMO.GEOMETRY = false;

% for HEIGHT=logspace(1,5.5,30)
for HEIGHT=logspace(1,5.5,30)
	ATMO.BEACON = [x,y,HEIGHT];
	show(ATMO);
	bigtitle(sprintf('Height: %f m',HEIGHT));
	drawnow;
end

% HEIGHT = 60000
% for x = -5:0.5:5
	for y = -5:1:5
		ATMO.BEACON = [x,y,HEIGHT];
		show(ATMO);drawnow;
	end
% end

% WFS = AOWFS(A,.25);
WFS = AOWFS(A,1/3);
imagesc(WFS.masked);sqar;
drawnow;

show(F);

WFS.sense(F).quiver;

RECON = AOReconstructor(A,DM,WFS);
RECON.program;

% CCD = 0;
% for nt=1:100
%     for s=1:7
%         A.segList{s}.tiptilt = TT(nt,:,s);
%     end
%     touch(A);
%     HALO=A.fft(512);
% %     imagesc(log10(abs(interp2(HALO(CH,CH)/maxHalo,1).^2)),[-6 0]);
%     CCD = CCD + abs(HALO).^2;
% 	imagesc(log10(normalize(CCD(CH,CH))),[-4 0]);
% 	drawnow;
% end

% imagesc(log10(normalize(abs(interp2(CCD(CH,CH),1)).^2)),[-6 0]);
% imagesc(log10(normalize(CCD(CH,CH))),[-6 0]);
% 
% TT(1,:,:)
% 
% squeeze(TT(1,:,:))

% % for nt=1:100
% %     for s=1:7
% %         A.segList{s}.tiptilt = TT(nt,:,s);
% %     end
% %     touch(A)
% %     HALO=A.fft(512);imagesc(log10(abs(interp2(HALO(CH,CH)/maxHalo,1).^2)),[-6 0]);drawnow;
% % end
% % for nt=1:100
% %     for s=1:7
% %         A.segList{s}.tiptilt = TT(nt,:,s);
% %     end
% %     touch(A)
% %     plotCAmpl(A.grid);drawnow;
% % end
% F = AOField(A);
% F
% F.size
% F = AOField(A);
% obj.name
% obj.defaultSize
% obj.FFTSize
% size(obj.grid_)
% nxy.axis_
% obj.AXIS_PIXEL
% F = AOField(A);
% nxy.axisPixel
% show(F)
% plotCAmpl(F.fft(512),1/8)
% plotCAmpl(F.fft(1024),1/8)
% imagesc(log10(normalize(F.fft(1024))),[-6 0]);
% imagesc(log10(normalize(abs(F.fft(1024)).^2)),[-6 0]);
% imagesc(log10(normalize(abs(F.fft(1024)).^2)),[-4 0]);
