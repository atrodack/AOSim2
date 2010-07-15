% function d = dh_dh(dhn, dhm, r, theta)
%
% dh_dh: disk harmonic function evaluation
%
% dhn   - Zero number
% dhm   - Bessel order number
% r     - radial coordinate
% theta - azimuthal angle coordinate
%
% returns: value of disk harmonic function at specified point
%          (real part)
%
% Norman Mark Milton                             May 6, 2008
%
function d = dh_dh(dhn, dhm, r, theta)
if dhn < 0
    error('dh_dh: invalid zero number');
end
d = dh_dhfast(dhm, dh_bjprime_zero(dhn, dhm), dh_norm(dhn, dhm), r, theta);