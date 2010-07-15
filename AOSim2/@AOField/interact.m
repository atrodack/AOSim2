function F = interact(F,medium)
  
% INTERACT: The AOField interaction routine.
%
% Written by: Johanan L. Codona, Steward Observatory: CAAO
% July 14, 2002
% March 5, 2003  Trimmed the checking somewhat.
  
  Fdom = F.AOGrid.domain;
  Mdom = medium.AOGrid.domain;
  
  if(~strcmp(Fdom,Mdom))
    error(sprintf('domain mismatch: F:%s med:%s',Fdom,Mdom));
  end

  % Center both arrays.
  F = center(F);
  medium = center(medium);
  
  if(F.AOGrid.origin == medium.AOGrid.origin & (F.AOGrid.size) == (medium.AOGrid.size)) % Just do it and go.
    if(isa(medium,'AOPhaseScreen'))
      medium = expi(medium);
    end

    F = F .* medium;
    return;
  end
  
  error('grid size mismatch??? Finish writing the dang routine!!!');
  
					     
