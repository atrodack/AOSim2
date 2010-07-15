% function [n, k, dhn, dhm] = dh_num(i)
%
% dh_num: disk harmonic function numbers
%
% i - disk harmonic function index number
%
% returns: function numbers disk harmonic basis function with
%          index i
%          optionally return DH indices (dhn, dhm)
%
% Norman Mark Milton                             August 25, 2005
%
function [n, k, dhn, dhm] = dh_num(i)
if (nargout ~= 2) & (nargout ~= 4)
    error('dh_num: invalid number of ouptut arguements');
end
if i < 1
    error('dh_num: invalid function index');
end
n      = 0;
n_elem = 1;
while n_elem < i
    n      = n + 1;
    n_elem = n_elem + n + 1;
end
if n == 0,
    k = 0;
else
    k = i - dh_elem(n - 1) - 1;
end

if nargout == 4
    [dhn, dhm] = dh_dhindex(n, k);
end