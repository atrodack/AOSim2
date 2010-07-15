% function a = dh_norm(dhn, dhm)
%
% dh_norm: disk harmonic normalization constant
%
% dhn - Zero number
% dhm - Bessel order number
%
% returns: normalization constant for disk harmonic function
%
% Norman Mark Milton                             August 25, 2005
%
function a = dh_norm(dhn, dhm)
if dhn < 0
    error('dh_norm: invalid zero number');
end
if (dhn == 0) & (dhm ~= 0)
    error('dh_norm: invalid Bessel number for order zero');
end
if dhn == 0
    a = 1;
else
    mabs = abs(dhm);
    l    = dh_bjprime_zero(dhn, dhm);
    k    = 2 * pi * l;
    a    = sqrt(1/((1 - (mabs / k)^2) * besselj(mabs, k)^2));
end
