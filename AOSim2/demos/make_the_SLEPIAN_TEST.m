D = 6.5;
% OBS = 0.1;

PDEFN = [ 0 0 D        1 0.05 0 0 0 0 0
          0 0 (OBS*D)  0 0.01 0 0 0 0 0 ];

Seg = AOSegment;
Seg.name = 'Slepian Mirror';
Seg.pupils = PDEFN;
Seg.make;

A = AOAperture;
A.name = 'Slepian Test Mirror';
A.addSegment(Seg);
