function a = subsasgn(a,index,value)

% SUBSASGN: Set AOField fields by named indexing.

global env;

switch index.type
	case '()'
		fprintf('  >>> type ()\n');
		disp(index.subs);
	case '.'
		switch index.subs
			case 'Cn2'
				a.Cn2 = value;
				a.r0 = (0.423*(2*pi/a.lambdaRef)^2*a.Cn2*a.thickness)^(-3/5); % Roddier
				
			case 'thickness' 
				a.thickness = value;
				% note: I compute r0 here rather than Cn2 because Cn2 is a 
				% material property and r0 is an effect.
				a.r0 = (0.423*(2*pi/a.lambdaRef)^2*a.Cn2*a.thickness)^(-3/5); % Roddier

			case 'r0'
				a.r0 = value;
				a.Cn2 = (value^(-5/3))/0.423/(2*pi/a.lambdaRef)^2/a.thickness;
				
			case 'altitude'
				a.altitude = value;
				
			case 'grid'
				a.AOGrid.grid = value;
				
			case 'phasor'
				a.type = 'phasor';
				a.AOGrid.grid = value;
				
			case 'phase'
				a.type = 'phase';
				a.AOGrid.grid = value;
				
			case 'type'
				a.type = value;
				
			case 'show'
				show@AOGrid(a);
				
			otherwise
				error(['unknown or illegal subscript: ' index.subs]);
		end
		
	otherwise
		fprintf('  >>> unsupported subscript type: %s\n',index.type);
		disp(index.subs);
end

