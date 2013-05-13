function setFoV(halfwidth)

% usage: setFoV(halfwidth)
if(halfwidth == 0)
	fprintf('WARNING: Tried to set a zero Field of View.\n');
	return;
end

halfwidth = abs(halfwidth);

xlim([-halfwidth halfwidth]);
ylim([-halfwidth halfwidth]);
