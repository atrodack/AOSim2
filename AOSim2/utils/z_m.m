% function m = z_m(n, k)
%
% z_m: Zernike function number m
%
% n - Zernike order number
% k - Zernike function of order m
%
% returns: function number M of Zernike basis function with
%          order N and function number K
%
% Norman Mark Milton                             August 25, 2005
%
function m = z_m(n, k)
if n < 0,
    error('z_m: invalid order number');
elseif (k < 0) | (k > n)
    error('z_m: invalid function number');
else
    m = n - (2 * k);
end