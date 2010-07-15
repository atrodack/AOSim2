% function [dhn, dhm] = dh_dhindex(n, k)
%
% dh_dhindex: disk harmonic radial and azimuthal index
%
% n - disk harmonic order number
% k - disk harmonic mode number
%
% returns: return DH indices (dhn, dhm)
%
% Norman Mark Milton                             August 25, 2005
%
function [dhn, dhm] = dh_dhindex(n, k)
% assume n and k are valid
% if n < 0
%     error('dh_dhindex: invalid order number');
% end
% if (k < 0) | (k > n)
%     error('dh_dhindex: invalid mode number');
% end
    
dhm = n - (2 * k);
if dhm == 0
    dhn = fix(n / 2);
elseif dhm > 0
    dhn = k + 1;
else
    dhn = n - k + 1;
end