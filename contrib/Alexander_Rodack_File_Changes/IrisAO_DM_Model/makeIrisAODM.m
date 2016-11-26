function [varargout] = makeIrisAODM(MAGNIFICATION, verbose,Scalloped_Field)
% [DM] = makeIrisAODM(MAGNIFICATION, verbose, Scalloped_Field)
% Makes a 37 Segment IrisAO DM.
%
%Inputs:
% MAGNIFICATION scales the DM Segment Size. SegPitch of 606 microns is
% hardcoded because it is what standard IrisAO Mirrors use. If your mirror
% is larger, MAGNIFICATION will scale that SegPitch accordingly.
%
% verbose turns on/off plotting various things as the construction takes
% place. You should set it to true once....it is really cool.
%
% Scalloped_Field turns on returning an AOField Object that encodes the
% actual surface shape of the segments.
%
%Output:
% DM is an AOAperture class object that has 37 segments stored within it.
% This is the IrisAO Model
%

SegPitch = 606e-6;

if nargin == 0
    MAGNIFICATION = 1;
    SegPitch = MAGNIFICATION * 606e-6;
    verbose = true;
    Scalloped_Field = false;
elseif nargin == 1
    SegPitch = MAGNIFICATION * SegPitch;
    Scalloped_Field = false;
    verbose = false;
elseif nargin == 2
    SegPitch = MAGNIFICATION * SegPitch;
    Scalloped_Field = false;
elseif nargin == 3
    SegPitch = MAGNIFICATION * SegPitch;
else
    error('Inputs are Incorrect');
end

% set the grid spacing
dx = SegPitch/100;
% set the gap width
GAP = SegPitch*0.02;

% set up a segment
Seg = AOSegment(ceil((SegPitch+2*GAP)/dx*1.25));
Seg.spacing(dx);

SCALLOPING = 100e-9;

% Grab the coordinate systems
[x,y] = Seg.coords;
[X,Y] = Seg.COORDS;

SEG = 1;

% Make a Hexagon
for n=1:6
    th = 60/180*pi*n;
    SEG=SEG.*smoothedge(SegPitch/2- (X*cos(th)-Y*sin(th)),dx);
    if verbose == true
        imagesc(x,y,SEG);sqar;axis xy;
        colormap(gray);
        drawnow;
        pause(0.2)
    end
end

%Set the Segment Grid to the hexagon
Seg.grid(SEG);
Seg.name = 'IrisAO DM Segment';

%Initialize the Model
DM = AOAperture();
Seg.name = 'IrisAO DM';
DM.spacing(dx);

% Some important numbers, and some unit vectors to point to where the
% segments need to be placed as they are constructed
a2 = (SegPitch/2+GAP)*1.333*1.0;
a1 = a2*sqrt(3)/2;

ua = a1 * [ 0  1.7 ];
th1 = 60/180*pi;
c = cos(th1); s = sin(th1);
ub = ([c -s; s c] * ua')';
ub_ = ([c s; -s c] * ua')';

% Build the DM, start at bottom row, and add in segments sweeping from left
% to right. The axis of symmetry for the rows is row #0. Below this, ub is
% used, and above this ub_ is used.  These unit vectors account for the
% necessary sign change in the rotation now that the segments are in
% positive space vertically.
for nrow=-3:0
    NumCol = 7-abs(nrow);
    
    START = -4*ua - nrow*ub;
    
    for ncol=1:NumCol
        DM.addSegment(Seg,START+ncol*ua);
        if verbose == true
            DM.show;
            drawnow;
        end
    end
end

for nrow=1:3
    NumCol = 7-abs(nrow);
    START = -4*ua + nrow*ub_;
    
    for ncol=1:NumCol
        DM.addSegment(Seg,START+ncol*ua);
        if verbose == true
            DM.show;
            drawnow;
        end
    end
end

% Store the location of the Center of each Segment in the global coordinate
% system
SEGCOORDS = zeros(37,2);
for n=1:37
    SEGCOORDS(n,:) = DM.segList{n}.Offset;
end

%return the DM Model
varargout{1} = DM;


%% Scalloping on Field
if Scalloped_Field == true;
    F = AOField(DM);
    F.name = 'Scalloped Field';
    F.lambda = AOField.VBAND;
    F.FFTSize = 1024*2;
    
    DM.lambdaRef = F.lambda;
    
    DM.touch;
    
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
    
    F.planewave*DM*PS;
    if verbose == true;
        F.show;
    end
    varargout{2} = F;
end
end
