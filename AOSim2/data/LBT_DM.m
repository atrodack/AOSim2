% number of actuators, radius in degrees from center (seen from sphere
% center)

LBT_meta = [
	9 1.2615
	15 2.1498
	21 3.0382
	27 3.9266
	33 4.8150
	39 5.7034
	45 6.5918
	51 7.4802
	57 8.3686
	63 9.2570
	69 10.1454
	75 11.0338
	81 11.9222
	87 12.8105];

LBT_Actuators = [];

for arow=1:size(LBT_meta,1)
	r = LBT_meta(arow,2);
	N = LBT_meta(arow,1);
	for n=1:N
		th = (n-1)/N*2*pi;
		LBT_Actuators(end+1,:) = r*[cos(th) sin(th)];
	end

end

LBT_Actuators(end+1,:) = [0 0];
