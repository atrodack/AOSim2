function r0 = scaleR0(r0,lambda0,lambda)
  
% scaleR0: Scale r0 from a reference wavelength to another wavelength.
%
% usage: r0 = scaleR0(r0,lambda0,lambda)
  
  
  alpha = 6/5;
  
  r0 = r0 .* (lambda./lambda0).^alpha;
  
  
