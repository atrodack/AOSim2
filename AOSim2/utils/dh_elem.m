% function n = dh_elem(o)
%
% dh_elem: disk harmonic elements
%
% o - disk harmonic order number
%
% returns: number of disk harmonic basis functions from
%          order 0 through order o
%
% Norman Mark Milton                             August 25, 2005
%
function n = dh_elem(o)
if o < 0,
    error('dh_elem: invalid order number');
else
    n = fix((o + 1) * (o + 2) / 2);
end