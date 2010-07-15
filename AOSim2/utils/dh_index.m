% function [i, dhn, dhm] = dh_index(n, k)
%
% dh_index: disk harmonic function index
%
% n - disk harmonic order number
% k - disk harmonic mode number
%
% returns: index of disk harmonic basis function with
%          order n and mode number k
%          optionally return DH indices (dhn, dhm)
%
% Norman Mark Milton                             August 25, 2005
%
function [i, dhn, dhm] = dh_index(n, k)
if (nargout ~= 1) & (nargout ~= 3)
    error('dh_index: invalid number of ouptut arguements');
end
if n < 0
    error('dh_index: invalid order number');
end
if (k < 0) | (k > n)
    error('dh_index: invalid mode number');
end
    
if n == 0
    i = 1;
else
    i = z_elem(n - 1) + k + 1;
end

if nargout == 3
    [dhn, dhm] = dh_dhindex(n, k);
end