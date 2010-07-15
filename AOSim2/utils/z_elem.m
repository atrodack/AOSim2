% function n = z_elem(o)
%
% z_elem: Zernike elements
%
% o - Zernike order number
%
% returns: number of Zernike basis functions from
%          order 0 through order o
%
% Norman Mark Milton                             August 25, 2005
%
function n = z_elem(o)
if o < 0,
    error('z_elem: invalid order number');
else
    n = fix((o + 1) * (o + 2) / 2);
end