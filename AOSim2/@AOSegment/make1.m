function [AOA] = make1(AOA)

% MAKE: Build the AOAperture grid.
% USAGE:
% AOA = make(AOA)
%
% Supported Apodizations: apod_type = pupils(6)
%
% 0: Cosine (arg 5 holds the width)
% 1: Sonine (nu is arg 7)  (beware of definitions!)
% 2: Elliptical: arg7 is Dy.
% 3: Angel-Cheng-Woolf
% 4: Angel2
% 5: Spergel2D: (arg 5 is gaussian equivalent diameter).
% 6: Woolf slab: (arg 5 is the height)
% 7: Specified circular mask. (set shadingMask and shadingRadii
%     mask) arg 7 is an overall complex phasor.)
% 8: ellipse: (use AOAperture/addEllipse)
%
%
% "Transmission" types:
% 0: circular hole
% 1: mirror
% -1: wedge
% -2: spider vanes (3: width, 6: nvanes, 7: starting theta)
%
% Johanan L. Codona, Steward Observatory, CAAO
% July 15, 2002     - First version
% July 2, 2002      - Modified to use the new sim infrastructure.
% Aug 20, 2002      - Added support for wedges and spiders.
% Nov. 22, 2002     - Added circular shaded masks.
% Dec 6, 2002       - Added support for multiple shaded masks.
% 20090406: JLC: make the grid bigger if it is too small to hold the pupil.

% global env;

% Prepare to build the pupils...
P = AOA.pupils;
S = size(P);
N = S(1);

% if(~AOA.touched)
% 	warning('AOSegment:UNTOUCHED','Not remaking an untouched pupil.');
% 	return;
% end
	
if(N==0)
	warning('AOSegment:UNDEF_PUPIL','The pupil has no entries.');
	return;
end

tr_glass   = 1;
tr_hole    = 0;
tr_wedge   = -1;
tr_spider  = -2;

AOA.tX;
AOA.center;

% maxX = max(P(:,1)+P(:,3)/2);
% maxY = max(P(:,2)+P(:,3)/2);
% minX = min(P(:,1)-P(:,3)/2);
% minY = min(P(:,2)-P(:,3)/2);
% 
% BBox = [minY minX; maxY maxX];
% BBox = AOA.BBox;
% BBoxN = ceil((BBox(2,:)-BBox(1,:))./AOA.spacing);
% BBoxN = 2*ceil(max(abs(BBox(:)))/min(AOA.spacing));
% if(sum(BBoxN > AOA.size)>0)
% 	%AOA.resize(BBoxN(2)+2,BBoxN(1)+2);
% 	AOA.resize(BBoxN+2);
% end
AOA.setBBox(AOA.BBox,AOA.dx);


% [x,y] = coords(AOA,true);
% % [X,Y] = meshgrid(x,y);
% [X,Y] = COORDS(AOA,true);
[x,y] = coords(AOA);
% [X,Y] = meshgrid(x,y);
[X,Y] = COORDS(AOA);

% Start with a clean slate...
AOA.zero;
%AOA = constant(AOA,1);

% decode the smoothing setting.

if(AOA.smooth<0)
	smooth = abs(AOA.smooth) * AOA.dx;
else
	smooth = abs(AOA.smooth);
end

for p=1:N
	pupil = P(p,:);

	x0 = pupil(1);
	y0 = pupil(2);

	D  = pupil(3);
	d = D/2;

	tr = pupil(4);

	apod = pupil(5);
	apod_type = pupil(6);

	if(tr == tr_glass)			% aperture
		switch apod_type
			case 0			% normal round
				fprintf('make: circular COSINE pupil(%d)\n',p);

				R = sqrt((X-x0).^2 + (Y-y0).^2);
				%AOA.AOGrid.grid = AOA.AOGrid.grid + (R<=(d-apod)) + ...
				%    (R>(d-apod)).*(R<d).*cos(.5*pi*(R-d-apod)/apod);

				g = AOA.grid;

				INSIDE  = R<(d-apod);
				OUTSIDE = R>d;

				g = g + INSIDE + (~INSIDE & ~OUTSIDE) .* cos(.5*pi*(R-d+apod)/apod);

				%whos
				AOA.grid_ = g;

			case 100			% SIMPLE ROUND
				fprintf('make:  SIMPLE ROUND (%d)\n',p);

				mask = 1-smoothedge(sqrt((X-x0).^2 + (Y-y0).^2) - d/2, apod);
				AOA.grid_ = AOA.grid + mask;
				clear mask;

			case -1			% normal round
				fprintf('make: circular COSINE-SQUARED pupil(%d)\n',p);

				R = sqrt((X-x0).^2 + (Y-y0).^2);
				%AOA.AOGrid.grid = AOA.AOGrid.grid + (R<=(d-apod)) + ...
				%    (R>(d-apod)).*(R<d).*cos(.5*pi*(R-d-apod)/apod);

				g = AOA.grid;

				INSIDE  = R<(d-apod);
				OUTSIDE = R>d;

				g = g + INSIDE + (~INSIDE & ~OUTSIDE) .* cos(.5*pi*(R-d+apod)/apod).^2;

				%whos
				AOA.grid_ = g;

			case 1				% sonine
				fprintf('make: circular SONINE pupil(%d)\n',p);

				R = sqrt((X-x0).^2 + (Y-y0).^2);
				nu = pupil(7);
				AOA.grid_ = AOA.grid + (R<d).*((1-(R./d)).^nu);
				%AOA.AOGrid.grid = AOA.AOGrid.grid + (R<d).*(cos(pi.*R./d/2)).^4;

			case 2				% elliptical
				fprintf('make: ELLIPTICAL pupil(%d)\n',p);

				Dx = pupil(3)/2;
				Dy = pupil(7)/2;
				R = ((X-x0)./Dx).^2 + ((Y-y0)/Dy).^2;

				AOA.grid_ = 1 - (1-AOA.grid) .* (1-.5*erfc(50*(R-1)));

			case 3			% Angel-Woolf notch mask.
				fprintf('make: Angel-Cheng-Woolf notch mask. (%d)\n',p);

				R = sqrt((X-x0).^2 + (Y-y0).^2);
				AOA.AOGrid.grid = AOA.AOGrid.grid + ...
					1 ...			% start from the middle.
					- smoothedge(R/(0.663*d)-1,smooth) ...
					+ smoothedge(R/(0.900*d)-1,smooth) ...
					- smoothedge(R/d-1,smooth);

			case 4			% Angel2 notch mask.
				fprintf('make: Angel2 notch mask. (%d)\n',p);

				R = sqrt((X-x0).^2 + (Y-y0).^2);
				AOA.AOGrid.grid = AOA.AOGrid.grid + ...
					1 ...			% start from the middle.
					- smoothedge(R/(0.778*d)-1,smooth) ...
					+ smoothedge(R/(0.915*d)-1,smooth) ...
					- smoothedge(R/d-1,smooth);


			case 40			% JLC TWEAKED ACW mask.
				fprintf('make: JLC-tweaked ACW notch mask. (%d)\n',p);

				R = sqrt((X-x0).^2 + (Y-y0).^2);
				AOA.AOGrid.grid = AOA.AOGrid.grid + ...
					1 ...			% start from the middle.
					- smoothedge(R/(0.6494*d)-1,smooth) ...
					+ smoothedge(R/(0.9081*d)-1,smooth) ...
					- smoothedge(R/d-1,smooth);

				%	    - smoothedge(R/(0.6725330*d)-1,smooth) ...
				%	    + smoothedge(R/(0.9112400*d)-1,smooth) ...
				%	    - smoothedge(R/d-1,smooth);

				%	    - smoothedge(R/(0.649*d)-1,smooth) ...
				%	    + smoothedge(R/(0.908*d)-1,smooth) ...
				%	    - smoothedge(R/d-1,smooth);


			case 50			% Spergel-Ge-Debes
				fprintf('make: Spergel-Ge-Debes. (%d)\n',p);
				% based on Ge and Debes wacky equations.
				R = pupil(3);
				alpha = pupil(7);
				a = pupil(8);
				b = pupil(9);

				smth = pupil(5);

				xx = X-x0;
				yy = Y-y0;

				cutoff = exp(-alpha.^2);
				top =  R*(exp(-(alpha.*xx./R).^2)-cutoff);
				bottom = yy+b.*top;
				top = yy-a*top;

				mask = smoothedge(bottom,smth) .* (1-smoothedge(top,smth));

				AOA.AOGrid.grid = AOA.AOGrid.grid + mask;

				clear top bottom mask;

			case 51			% Angel2 notch mask.
				fprintf('make: Spergel2D. (%d)\n',p);

				R = sqrt((X-x0).^2 + (Y-y0).^2);
				Reff = pupil(5)/2;
				clip = exp(-(d/Reff)^2);
				norm = 1/(1-clip);
				AOA.AOGrid.grid = AOA.AOGrid.grid + ...
					norm*(exp(-(R.^2)/(Reff^2)) - clip) .* (R<d);


			case 6
				fprintf('make: Woolf slab. (%d)\n',p);

				h = pupil(5)/2;

				AOA.AOGrid.grid = AOA.AOGrid.grid + ...
					(1-smoothedge(abs(X-x0)-d,smooth)).*(1-smoothedge(abs(Y-y0)-h,smooth));

			case 7
				phasor = pupil(7);
				if(phasor == 0)
					phasor = 1;
				end

				fprintf('make: circularly symmetric APODIZED mask. (%d) |phasor|=%f phase=%fpi\n', p, abs(phasor), angle(phasor)/pi);

				R = sqrt((X-x0).^2 + (Y-y0).^2);

				tr = AOA.shadingValues;
				radii = AOA.shadingRadii;

				if(radii(1) ~= 0)
					radii = [0,radii(1)-eps,radii];
					tr = [0,0,tr];
				end

				tr(end+1) = 0;
				tr(end+1) = 0;
				radii(end+1) = 1000;

				%whos tr radii

				mask = interp1(radii,tr,R);

				AOA.AOGrid.grid = AOA.AOGrid.grid + phasor*mask;

				clear mask tr radii R phasor;



			case 8		% an ELLIPTICAL segment

				% pupil(:,1:2): center
				% pupil(:,4):   transparency
				% pupil(:,5):   anti-aliasing width (1.5 * spacing)
				% pupil(:,7):   x radius
				% pupil(:,8):   y radius
				% pupil(:,9):   rotation (radians)
				fprintf('make: elliptical pupil(%d)\n',p);

				%Xp = (X-x0)./pupil(7);
				%Yp = (Y-y0)./pupil(8);
				Xp = X-x0;
				Yp = Y-y0;

				th = pupil(9);

				Xrot = ( Xp*cos(th) + Yp*sin(th))./(pupil(8)./2);
				Yrot = (-Xp*sin(th) + Yp*cos(th))./(pupil(7)./2);

				R = sqrt(Xrot.^2 + Yrot.^2);

				g = AOA.AOGrid.grid;

				apod = pupil(5) / mean(pupil(7:8));

				d = 1;

				INSIDE  = R<(d-apod);
				OUTSIDE = R>d;

				g = g + INSIDE + (~INSIDE & ~OUTSIDE) .* cos(.5*pi*(R-d+apod)/apod);

				AOA.AOGrid.grid = g;

				clear Xp Yp Xrot Yrot g;

				% END of ellipses



			case 9			% an elliptical segment

				fprintf('make: shadow line segment(%d)\n',p);

				% pupil(:,1:2): endpoint1XY
				% pupil(:,3):   width
				% pupil(:,4):   transparency
				% pupil(:,5):   anti-aliasing width (1.5 * spacing)
				% pupil(:,6):   9
				% pupil(:,7:8): endpoint2XY
				% pupil(:,9):   NaN

				X1 = pupil(1);
				Y1 = pupil(2);

				Z1 = X1 + i*Y1;

				X2 = pupil(7);
				Y2 = pupil(8);

				Z2 = X2 + i*Y2;

				fprintf('LINE from (%f,%f) to (%f,%f)\n', X1,Y1,X2,Y2);

				len = abs(Z2-Z1);
				Zmid = (Z1+Z2)/2;
				thick = pupil(3)/2;
				theta = atan2(Y2-Y1,X2-X1);

				Z = X + Y*i;
				Z = (Z-Zmid) .* exp(-i*theta);

				smth = pupil(5)/2;

				mask = 1-smoothedge(abs(imag(Z))-thick/2,-smth).*smoothedge(abs(real(Z))-len/2,-smth);

				AOA.AOGrid.grid = AOA.AOGrid.grid .* mask;

				clear Z,mask;


				% END of line segment shadows.



			otherwise
				error(sprintf('unknown apod_type %d',apod_type));
		end				% switch



		%%% OBSTRUCTIONS

	elseif(tr == tr_hole)		% obstruction
		fprintf('make: circular HOLE.\n');

		R = sqrt((X-x0).^2 + (Y-y0).^2);

		AOA.grid_ = AOA.grid .* smoothedge(R-d,apod);




		% AOA.AOGrid.grid = AOA.AOGrid.grid .* (1-erfc((R-D/2)/smooth)/2);

	elseif(tr == tr_wedge)		% wedge obstruction

		theta0 = pi * pupil(7) / 180;
		dtheta = pi * pupil(5) / 180/2;

		R = sqrt((X-x0).^2 + (Y-y0).^2);
		ANG = atan2(Y-y0,X-x0);

		AOA.grid_ = AOA.grid .* (abs(ANG-theta0)>dtheta);


	elseif(tr == tr_spider)		% spider obstruction

		th0 = pi * pupil(7) / 180;
		nvanes = pupil(6);
		width = pupil(3);

		fprintf('SPIDER VANES: %d %f wide vanes at angle %f\n',nvanes, width,pupil(7));

		X0 = X-x0;
		Y0 = Y-y0;

		R = sqrt((X-x0).^2 + (Y-y0).^2);

		SPIDER = ones(size(X0));

		dth = 2*pi/nvanes;

		% 		w = max(AOA.dx,width/2);
		w = width;
		if(width < AOA.dx/2)
			peak = 2*width/AOA.dx;
		else
			peak = 1;
		end

		for n=0:(nvanes-1)
			th = th0+n*dth;

			Xrot =  X*cos(th) + Y*sin(th);
			Yrot = -X*sin(th) + Y*cos(th);

% 			SPIDER = SPIDER .* (1-peak*exp(-(Yrot/w).^2) .* double(Xrot>0));
			SPIDER = SPIDER .* double(~(abs(Yrot)<w/2 & Xrot>0));
        end

        if(pupil(8) > 0)
            SPIDER = 1 - (1-SPIDER).*(R<=pupil(8));
        end
		AOA.grid_ = AOA.grid .* SPIDER;
		clear SPIDER X0 Y0 Yrot;

	end
end

AOA.touched = false;

