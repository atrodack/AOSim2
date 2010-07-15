% function [n, k] = z_num(J)
%
% z_num: Zernike function numbers
%
% J - Zernike function index number
%
% returns: function numbers Zernike basis function with
%          index J
%
% Norman Mark Milton                             August 25, 2005
%
function [n, k] = z_num(J)
if J < 1
    error('z_num: invalid function index');
end
n      = 0;
n_elem = 1;
while n_elem < J
    n      = n + 1;
    n_elem = n_elem + n + 1;
end
if n == 0,
    k = 0;
else
    k = J - z_elem(n - 1) - 1;
end