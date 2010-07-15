% function [n, k, m] = z_noll_num(J)
%
% z_noll_num: Noll Zernike function numbers
%
% J - Zernike function index number
%
% returns: Noll function numbers Zernike basis function with
%          index J
%
% Norman Mark Milton                             August 25, 2005
%
function [n, k, m] = z_noll_num(J)
if J < 1,
    error('z_noll_num: invalid function index');
end
[n, k] = z_num(J);            % not Noll k, just correct n %
if n > 0,
    i = J - z_elem(n - 1) - 1;
else
    i = 0;
end
if even(n),
    if ~even(i),
        i = i + 1;
    end
else
    if even(i),
        i = i + 1;
    end
end
if even(J),
    m = -i;
else
    m = i;
end
k = fix((n - m) / 2);         % compute correct k %