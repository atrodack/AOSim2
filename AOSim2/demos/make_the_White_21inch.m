%% How to build and program the Steward Observatory 61 inch model.
% Johanan L. Codona, Steward Observatory, University of Arizona.
% 20090825 JLCodona: First-light version. NOTE: This telescope does not
% have AO.

%% Start clean...
% close
% clear classes

D = 21.0 * 0.0254;
dx = 0.01;
PUPIL = [
            0            0          D            1         dx            0            0            0            0            0
            0            0         0.2*D       0         dx/2            0            0            0            0            0
%             0            0        0.012           -2     dx            4            0            0            0            0
];

Seg = AOSegment;
Seg.name = 'SO White Primary';
Seg.pupils = PUPIL;
Seg.spacing(0.01);
Seg.make;

clf;
% Seg.touch.make.show;
A = AOAperture;
A.spacing(0.01);
A.name = 'Steward White Telescope (21 inch)';
A.addSegment(Seg);
A.show;
colormap(gray);
