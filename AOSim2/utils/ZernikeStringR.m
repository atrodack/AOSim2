function [s,dx,dy] = ZernikeStringR(n,f)
%<PRE>
%ZernikeStringR  Generate string definitions of Zernike's circle polynomials  
%  and also their partial derivatives in rectangular form.
%
% Use: [s,dx,dy] = ZernikeStringR(n,f)
%
% NOTE: UNLIKE MALACARA, WE ADOPT THE MATHEMATICAL CONVENTION OF X-AXIS AS 
% REFERENCE MERIDIAN FOR SPECIFYING ANGLES.  i.e. we adopt the convention
%  x=r*cos(theta), y=r*sin(theta). This is achieved by exchanging X & Y 
% in Malacara's eqn. A2.17 (1978) ar eqn. 13.40 (1992)
%
% Output matrices: 	s is the Zernike equation
%					dx is its partial derivatives in x, 
%					dy is its partial derivative in y.
%
% Input args: 	n = radial "degree" or "order", [a positive integer]
%  			 	f = azimuthal frequency, [a signed integer, magnitude <= n]
%
% Method:  implement eqn. A2.17 from "Optical Shop Testing" (1978) by Malacara.
%			which is meant to be the same as eqn. 13.40 in 2nd edition (1992).
% 	IMPORTANT NOTE: there are errors in both editions of Malacara's book!
% 	1978 edition has an error in Table A2.2, case n=4,m=-2 which is corrected 
%	  in 1992.
%	1978 edition has errors in the bottom row of table A2.3 which are corrected
%	  in 1992.
% 	1992 version of Table 13.3 introduced an error for the case of n=even, 
%		(i.e. sine term).  Refer to 1978 version for correct formula.
%	1992 eqn. 13.40 has an equals sign on the second line which should be
%		 multiplication sign. Refer to 1978 version for correct formula.
%	1992 edition contains a critical phrase on p. 471 that is missing from 1978.
%		"where m must be replaced by n-m when n-2m =<0".  This is the key to 
%		handling negative azimuthal frequencies.  However, the statement only 
%		refers to the limits of summation, and one binomial term, but not to 
%		the other terms in the formula. Also, the sign of frequency must be 
%		reversed when f<0.
%
% 	Two solutions are required: one for positive azimuthal frequency, another 
% 	for negative. Malacara's text is very confusing on this point. He says on
% 	p. 471 of 1992 edition to replace m with n-m in summation formula. By trial
% 	& error, I determined that this should only be done for summantion limits, 
% 	and for just one of the algebraic terms. I also figured the absolute value 
%	of frequency needed to be used to get binomial coeff. I presume these are  
%	all typographical errors in the printed book since Malacara's tables of
%	functions in polar and cartesian forms are in agreement with each other.
% 	This MATLAB routine generates formulas that match Malacara's published tables.
%
% Verification: compare output with Malacara's Table A2.2 (1978) and 13.2 (1992).
%	Note that this program sometimes produces terms of the same order which may  
%	be combined manually prior to comparing results with Malacara's tables.  
% 	See test routine MALACARA_VS_AUTO for the method used to compare output 
%	of this routine with published table of Malacara (1992).  Extrapolation 
%	to orders > 20 needs careful testing.  It isn't obvious that coefficients for 
%	higher orders will be correctly computed by routines FACTORIAL and BINCOEFF.
%
% Background:
%	In polar form, a wavefront is represented by a sum of Zernike modes 
%	of the form:
%	 R[superscript_f, subscript_n] * ( C[n,f]*cos(f*theta) + D[n,f]*sin(f*theta) )
%	(this is eqn. A2.8 of Malacara (1978) or eqn. 13.22 of Malacara (1992)).
%
%	Amazingly, this polar equation can be converted to rectangular form without 
%	recourse to brute-force transformation eqns. x=rcos(q), y=rsin(q) plus the
%	 repeated use of trig identities to unpack the harmonic trig terms.  
%	The key is eqn. A2.14 (1978) or eqn. 13.37 (1992) which gives power series 
%	representation of cos(n*theta) and sin(n*theta) in rectangular form.  
%	(NOTE: there is a typographcial error in both of these eqns.: in the 
%	binomial coefficient, 2i+p should be 2j+p.)  According to eqn. A2.17, 
%	the un-normalized Zernike circle polynomials in monomial (i.e.rectangular) 
%	form is given by a triple sum over indices i,j,k which determine the values 
%	of the coefficients of each term in the polynomial, as well as the powers 
%	of x,y.  (Note: the corresponding eqn. 13.40 (1992) has a typographical error:
%	the equals sign ("=") in the last line should be a multiplication sign ("*").)
%
% 	Implementation of Malacara's eqn. A2.17 requires careful attention to detail in 
% 	handling the 4 cases of even & odd order, sine and cosine terms.
%
% Terms:	n = radial "degree" or "order" of the polynomial 
%			f = azimuthal frequency (Note: Malacara uses letter l in his eqns., 
%				which invites confusion).
%			m = (n-f)/2 = a convenient combination of n,f used by Malacara
%
% Required subroutines: functions FACTORIAL and BINCOEFF by LNT.
%
% LNT 7-Mar-98.  Add derivatives 8-MAR-98.

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

m = (n-f)/2;					% a positive index parameter

% Set up parameters that depend on sign of azimuthal frequency
% Evaluate parameters p,q based on Malacara's revised table 13.4 in 1991 edition.
% WARNING!  Malacara's table 13.4 has an error in n-even, sine condition.

ff = abs(f);				% frequency variable for binomial coeff
if f>0 						% sine version
	B = m;					% a loop parameter
	if (2*round(n/2) == n)	% true if n is even
		p=1;
		q=round(f/2)-1;		% use 1978 version
	else					% true if n is odd
		p=1;
		q=round((f-1)/2);	% nervous about floating point operations
	end
else						% cosine version
	B = n-m; 					% change loop param & one term when f<0
	if (2*round(n/2) == n)	% true if n is even
		p=0;
		q=-round(f/2);		% change sign of f
	else					% true if n is odd
		p=0;
		q=-round((f+1)/2);	% change sign in 2 places
	end
end

% Implement Malacara's eqn. A2.17 as 3 nested loops

s='';						% initialize the output strings
dx='';
dy='';

for i = 0:q
	for j = 0:B
		for k = 0:B-j
			c = (-1)^(i+j) * bincoeff(ff,2*i+p) * bincoeff(B-j,k) * factorial(n-j);
			c =	c /( factorial(j) * factorial (m-j) * factorial (n-m-j) );
			y_power = 2*(i+k)+p;
			x_power = n-2*(i+j+k)-p;

			% prepare a string representation of this Z term

			if c>=0
				s = [s,' +',num2str(c)];
			end
			
			if c<0
				s = [s,' -',num2str(-c)];
			end
			
			if x_power > 0
				if x_power == 1
					s = [s,'.*x'];
				else
					s = [s,'.*x.^',num2str(x_power)];	
				end
			end
			
			if y_power > 0
				if y_power == 1
					s = [s,'.*y'];
				else
					s = [s,'.*y.^',num2str(y_power)];	
				end
			end	
			
			% differentiate w/r to x
			cx = c*x_power;
			dx_power = x_power - 1;
			
			% prepare a string representation of dZ/dx term
			if cx ~= 0
				if cx > 0
					dx = [dx,' +',num2str(cx)];
				end
			
				if cx < 0
					dx = [dx,' -',num2str(-cx)];
				end
			
				if dx_power > 0
					if dx_power == 1
						dx = [dx,'.*x'];
					else
						dx = [dx,'.*x.^',num2str(dx_power)];	
					end
				end
			
				if y_power > 0
					if y_power == 1
						dx = [dx,'.*y'];
					else
						dx = [dx,'.*y.^',num2str(y_power)];	
					end
				end	
			end
			
			% differentiate w/r to y
			cy = c*y_power;
			dy_power = y_power - 1;

			% prepare a string representation of dZ/dy term
			if cy ~= 0
				if cy>0
					dy = [dy,' +',num2str(cy)];
				end
			
				if cy<0
					dy = [dy,' -',num2str(-cy)];
				end
					
				if x_power > 0
					if x_power == 1
						dy = [dy,'.*x'];
					else
						dy = [dy,'.*x.^',num2str(x_power)];	
					end
				end
			
				if dy_power > 0
					if dy_power == 1
						dy = [dy,'.*y'];
					else
						dy = [dy,'.*y.^',num2str(dy_power)];	
					end
				end	
			end
			
		end	%end of k-loop
	end	%end of j-loop
end %end of i-loop
%</PRE>
