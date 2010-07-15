%# Non-Redundant Mask Configurations for the MMT
%#
%# Josh Eisner
%# 5/4/2009
%#
%# Parameters: D=6.35 m == 4.44 mm; D_Obscuration=1.016 m
%# Parameters: 

D = 6.35;
Dmask = 4.44e-3;
D_Obscuration = 1.016;

%#
%#
EISNER_MASKS = {};

%#1.89-m sub-apertures: Mask throughput = 27.3%
mask = struct;
mask.diam = 1.89;
mask.centers = [
-2.18847      1.27437      1.01494
-0.194592     -1.81706      1.96874 ]';

P = zeros(size(mask.centers,1),9);
P(:,3) = mask.diam;
P(:,1:2) = mask.centers;
P(:,4) = 1;
P(:,5) = 0.04;

mask.pupils = P;

EISNER_MASKS{end+1} = mask;


%#
%#1.23-m sub-apertures: Mask throughput = 15.4%
mask = struct;
mask.diam = 1.23;
%#1.23-m sub-apertures: Mask throughput = 15.4%
mask.centers = [
2.12178     -1.36475     -2.25535    -0.676999
-1.39517     -2.09474     0.200582      2.46633]';

P = zeros(size(mask.centers,1),9);
P(:,3) = mask.diam;
P(:,1:2) = mask.centers;
P(:,4) = 1;
P(:,5) = 0.04;

mask.pupils = P;

EISNER_MASKS{end+1} = mask;

%#
%#0.94-m sub-apertures: Mask throughput = 11.2%
mask = struct;
mask.diam = 0.94;
mask.centers = [
1.25355     -2.45217      1.35699      2.64936    -0.729136
2.22212    -0.895969     -2.12183    -0.534693      2.52970]';

P = zeros(size(mask.centers,1),9);
P(:,3) = mask.diam;
P(:,1:2) = mask.centers;
P(:,4) = 1;
P(:,5) = 0.04;

mask.pupils = P;

EISNER_MASKS{end+1} = mask;

%#
%#0.80-m sub-apertures: Mask throughput =  9.8%
mask = struct;
mask.diam = 0.80;
mask.centers = [
-2.69742      2.21345    -0.635066     -1.39067      1.17673      1.81940
0.500169     0.453061     -1.88000      1.63989     -2.45520      2.08526]';

P = zeros(size(mask.centers,1),9);
P(:,3) = mask.diam;
P(:,1:2) = mask.centers;
P(:,4) = 1;
P(:,5) = 0.04;

mask.pupils = P;

EISNER_MASKS{end+1} = mask;

%#
%#0.62-m sub-apertures: Mask throughput =  6.8%
mask = struct;
mask.diam = 0.62;
mask.centers = [
-2.54513     0.657616      1.93198     -1.66641      2.30691     0.575806      1.01183
0.713915     -1.32296     -1.45635     -1.90264      1.10164      2.63698     -2.60825]';

P = zeros(size(mask.centers,1),9);
P(:,3) = mask.diam;
P(:,1:2) = mask.centers;
P(:,4) = 1;
P(:,5) = 0.04;

mask.pupils = P;

EISNER_MASKS{end+1} = mask;

%#
%#0.50-m sub-apertures: Mask throughput =  5.1%
mask = struct;
mask.diam = 0.50;
mask.centers = [
1.99167      2.65663     0.144214     -2.52787      1.37193     -2.69511    -0.699613 -0.950475
-1.36485     0.377195     -2.15341    -0.491367      1.80742     0.636230     -2.79628  2.51453]';

P = zeros(size(mask.centers,1),9);
P(:,3) = mask.diam;
P(:,1:2) = mask.centers;
P(:,4) = 1;
P(:,5) = 0.04;

mask.pupils = P;

EISNER_MASKS{end+1} = mask;




