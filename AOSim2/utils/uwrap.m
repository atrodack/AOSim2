function uwphase = uwrap(phase,algorithm)

% function uwphase = uwrap(phase,[algorithm])
% Unwrap phase using fancy external programs.
%
%

if(nargin>1)
    USE_FFT = true;
else
    USE_FFT = false;
end

[nx,ny] = size(phase);

if(USE_FFT)
    NX = 2^floor(log2(nx)+1)+1;
    NY = 2^floor(log2(ny)+1)+1;
    phase = padarray(phase,[(NX-nx) (NY-ny)],'post');
else
    NX = nx;
    NY = ny;
end

%TEMPDIR = '/mnt/ramdisk';
TEMPDIR = '/tmp';
RAND = floor(1000*rand);

RAW = sprintf('%s/_uwrap_%dx%d_%010d.phase',TEMPDIR,NX,NY,RAND);
RESULT = sprintf('%s/_uwrap_%dx%d_%010d.uwphase',TEMPDIR,NX,NY,RAND);

fid = fopen(RAW,'wb');
fwrite(fid, phase(:), 'float');
fclose(fid);

% fprintf('/home/jlc/LAB/gold -input pengaha1.phase -format float -xsize 513 -ysize 513 -dipole no -output pengahaout1')
% system('/home/jlc/LAB/gold -input pengaha1.phase -format float -xsize 513 -ysize 513 -dipole no -output pengahaout1')


if(USE_FFT)
    CMD = sprintf('/home/jlc/LAB/unwt -input "%s" -format float -xsize %d -ysize %d -output "%s"  > uwrap.out 2> uwrap.err',RAW,NX,NY,RESULT);
else
    CMD = sprintf('/home/jlc/LAB/gold -input "%s" -format float -xsize %d -ysize %d -output "%s" > uwrap.out 2> uwrap.err',RAW,NX,NY,RESULT);
end

% fprintf('CMD: %s\n',CMD);
ret = system(CMD);

fid = fopen(RESULT,'r');
uwphase = fread(fid,NX*NY,'float');
uwphase = reshape(uwphase,NX,NY);
uwphase = uwphase(1:nx,1:ny);
fclose(fid);

delete(RAW);
delete(RESULT);

return;
