function M = Vprod(M,v)

% M = Vproduct(M,v)
%
% Do what it takes to do a vector equivalent of an element-by-element
% product of vector v with matrix M.
% 
% 20081228: JLCodona

if(dims(M) ~= 2)
	error('M is not a matrix');
end

if(~isvector(v))
	error('v is not a vector');
end

if(isscalar(v))
	M = M*v;
	return;
end

if(size(v,1)>1) % column vector
	if(size(v,1) ~= size(M,1))
		whos v M
		error('column vector does not match first M dimension');
	end
	
	M = M .* (v*ones(1,size(M,2)));
	
else % row vector
	if(size(v,2) ~= size(M,2))
		whos v M
		error('row vector does not match second M dimension');
	end
	
	M = M .* (ones(size(M,1),1)*v);
	
end

return;
