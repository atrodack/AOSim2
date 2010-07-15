% function l = dh_bjprime_zero(dhn, dhm)
%
% dh_bjprime_zero: disk harmonic BesselJ' zero
%
% dhn - Zero number
% dhm - Bessel order number
%
% returns: dhn'th zero to first derivative of dhm'th order BesselJ
%          function (first kind)
%
% Norman Mark Milton                             August 25, 2005
%
function l = dh_bjprime_zero(dhn, dhm)
global dh_bjprime_zero_tab
if dhn < 0
    error('dh_bjprime_zero: invalid zero number');
end
mabs = abs(dhm);
if dhn == 0
    l = 0;
else
    [r, c] = size(dh_bjprime_zero_tab);
    if dhn > c
        error('dh_bjprime_zero: invalid zero number');
    elseif mabs + 1 > r
        error('dh_bjprime_zero: invalid order number');
    else
        l = dh_bjprime_zero_tab(mabs + 1, dhn);
    end
end