% IrisAO pitch is 606 microns.

MAGNIFICATION = 1.75;
ACTUAL_D = 3.4e-3;

SegPitch = MAGNIFICATION * 606e-6;
dx = SegPitch/100;
GAP = SegPitch*0.02;

Seg = AOSegment(ceil((SegPitch+2*GAP)/dx*1.25));
Seg.spacing(dx);

SCALLOPING = 100e-9;

[x,y] = Seg.coords;
[X,Y] = Seg.COORDS;

%Q = sqrt(3)*max(x)*exp(2*pi*1i*(0:6)'/6);
Q = (SegPitch/2)*exp(2*pi*1i*(0:6)'/6);
V = [real(Q) imag(Q)];
dV = diff(V,1,1);

SEG = 1;

% for n=1:6
%     SEG=SEG.*smoothedge(dV(n,1)*(X-V(n,1))+dV(n,2)*(Y-V(n,2)),dx/40);
%     imagesc(x,y,SEG);sqar;axis xy;
%     colormap(gray);
%     drawnow;
%     input ' por que? '
% end
for n=1:6
    th = 60/180*pi*n;
    %SEG=SEG.*((X*cos(th)-Y*sin(th))<SegPitch/2);
    SEG=SEG.*smoothedge(SegPitch/2- (X*cos(th)-Y*sin(th)),dx);
    %     imagesc(x,y,SEG);sqar;axis xy;
    %     colormap(gray);
    %     drawnow;
    %     input ' por que? '
end

Seg.grid(SEG);
Seg.name = 'IrisAO DM Segment';

DM = AOAperture();
Seg.name = 'IrisAO DM';
DM.spacing(dx);
DM.show; drawnow;
% BB = Seg.BBox
% a2 = max(abs(BB(:,2))) + GAP
% a1 = max(abs(BB(:,1))) + GAP
a2 = (SegPitch/2+GAP)*1.333*1.0;
a1 = a2*sqrt(3)/2;

ua = a1 * [ 0  1.7 ]
th1 = 60/180*pi;
c = cos(th1); s = sin(th1);
ub = ([c -s; s c] * ua')'
ub_ = ([c s; -s c] * ua')'

for nrow=-3:0
    NumCol = 7-abs(nrow)
    
    START = -4*ua - nrow*ub;
    
    for ncol=1:NumCol
        DM.addSegment(Seg,START+ncol*ua);
    end
end

for nrow=1:3
    NumCol = 7-abs(nrow)
    START = -4*ua + nrow*ub_;
    
    for ncol=1:NumCol
        DM.addSegment(Seg,START+ncol*ua);
    end
end


SEGCOORDS = zeros(37,2);
for n=1:37
    SEGCOORDS(n,:) = DM.segList{n}.Offset;
end

imagesc(DM.grid);sqar;
DM.show; drawnow;

A = AOSegment(DM);
% D=5.8*SegPitch;
D = MAGNIFICATION *  ACTUAL_D;
PNECO = [0 0 D 1 1.5*dx 0 0 0 0 0];
A.pupils = PNECO;
A.make;

[x,y] = A.coords;
% [X,Y] = Seg.COORDS;

F = AOField(DM);
F.name = 'Science Field';
F.lambda = AOField.HeNe_Laser;
F.FFTSize = 1024*2;

DM.lambdaRef = F.lambda;

% DM.segList{10}.piston = 150e-9
% DM.segList{20}.tiptilt = [1 2]/2000/NN;

DM.touch;
F.planewave*DM*A;
F.show;

win = chebwin(9,20);
WIN = win*win';

PS = AOScreen(DM);
[X,Y] = PS.COORDS;

PS.zero;
w = SegPitch/2.2;

for n=1:length(DM.segList)
    LOC = DM.segList{n}.Offset;
    R = sqrt((X-LOC(2)).^2+(Y-LOC(1)).^2);
    PS.grid(PS.grid+exp(-(R/w).^2));
end

PS.grid(SCALLOPING*normalize(PS.grid-min(PS.grid_(:))));

FOV=30*F.lambda/D*206265;
dFOV=0.333*F.lambda/D*206265;

% FOV=1
dFOV = FOV/100
PSF = F.mkPSF(FOV,dFOV);
SCALE = .2

DM.trueUp.touch.make;

PS.show;

F.planewave*DM.touch.make*A*PS;
F.show;

PSF = F.mkPSF(FOV,dFOV);
% PSF_ = PSF_ + PSF;
