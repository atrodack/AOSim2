% function d = dh_dhfast(dhm, l, a, r, theta)
%
% dh_dhfast: disk harmonic function evaluation (fast)
%
% dhm   - Bessel order number
% l     - disk harmonic spatial frequency
% a     - disk harmonic normalization constant
% r     - radial coordinate
% theta - azimuthal angle coordinate
%
% returns: value of disk harmonic function at specified point
%          (real part)
%
% Norman Mark Milton                             August 25, 2005
%
function d = dh_dhfast(dhm, l, a, r, theta)
mabs = abs(dhm);
k    = 2 * pi * l;
d    = a * besselj(mabs, k * r);
if dhm > 0
    d = sqrt(2) * d .* sin(mabs * theta);
elseif dhm < 0
    d = sqrt(2) * d .* cos(mabs * theta);
end