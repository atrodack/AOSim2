%% How to build JLC's 10" LX200 f/6.3 model.
% Johanan L. Codona, Steward Observatory, University of Arizona.
% 20090827 JLCodona: First-light version. NOTE: This telescope does not
% have AO.

%% Start clean...
% close
% clear classes

D = 10.0 * 0.0254;
d = 3.75 * 0.0254;
dx = 1/30;
PUPIL = [
            0            0          D            1         dx            0            0            0            0            0
            0            0          d            0         dx/2            0            0            0            0            0
%             0            0        0.012           -2     dx            4            0            0            0            0
];

Seg = AOSegment;
Seg.name = 'LX200 Primary';
Seg.pupils = PUPIL;
Seg.spacing(dx);
Seg.make;

clf;
% Seg.touch.make.show;
A = AOAperture;
A.spacing(dx);
A.name = 'LX200';
A.addSegment(Seg);
A.show;
colormap(gray);
