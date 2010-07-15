function s = ZernikeStringP(n,f)
%<PRE>
%ZernikeStringP  Generate a string matrix for a Zernike circle polynomial in 
%	polar coordinates.
%
% Use: s = ZernikeStringP(n,f)
%
% Arguments: n = radial "degree" or "order", [a positive integer]
%  			 f = azimuthal frequency, [a signed integer, magnitude <= n]
%
% Method:  use eqn. A2.6 from "Optical Shop Testing" (1978) by D. Malacara.
%
% Terms:	n = radial "degree" or "order" of the polynomial 
%			f = azimuthal frequency (Note: Malacara uses letter l, which 
%				invites confusion)
%			m = (n-|f|)/2 = convenient parameter used by Malacara
%
% Required subroutines: functions FACTORIAL and BINCOEFF by LNT.
%
% LNT 8-Mar-98

% Error checking on input arguments
if nargin < 2, 
    error('Requires two input arguments');
end
if prod(size(n)) ~= 1 
      error('First input argument must be a scalar.')
end
if (imag(n) ~= 0) | (n < 0)
      error('First input argument must be real and nonnegative.')
end
if n ~= round(n)
      error('First input argument must be an integer.')
end
if prod(size(f)) ~= 1 
      error('Second input argument must be a scalar.')
end
if (imag(f) ~= 0) 
      error('Second input argument must be real.')
end
if f ~= round(f)
      error('Second input argument must be an integer.')
end
if abs(f) > n
      error('Magnitude of second argument must not exceed first argument.')
end
if ( 2*round((n-f)/2) ~= (n-f) )
	  error('Arguments must both be even, or both odd.')
end

% Malacara's eqn. A2.6 (1978), or eqn. 13.20 (1992) says to form product 
% of radial and trig functions. To Get the radial polynomial, use 
% eqn. A2.4 (1978) or eqn. 13.15 (1992) using absolute value of frequency.

m = (n-abs(f))/2;		% eqn. A2.4 is only valid for f>0, but same result 
%							applies for f<0
s='';					% initialize the output string

for i = 0:m
	c = (-1)^i * factorial(n-i)/(factorial(i)*factorial(m-i)*factorial(n-m-i));	
	power = n-2*i;
	
	% prepare a string representation of this term in the summation
	if c>=0
		s = [s,' +',num2str(c)];
	end
			
	if c<0
		s = [s,' -',num2str(-c)];
	end

	% first do the radial term
	if power > 0
		if power == 1
			s = [s,'.*r'];
		else
			s = [s,'.*r.^',num2str(power)];	
		end
	end

	% now append the angular term
	if f>0
		s = [s,'.*sin(',num2str(f),'*t)'];	
	elseif f<0
		s = [s,'.*cos(',num2str(-f),'*t)'];
	end	
end
%</PRE>	
