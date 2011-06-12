function [pupilmask_trim,detector2d] = wrapmask(indpup)

%load list of illuminated pixels
%indpup = load('indpup20090925-113454');

%create a detector
detector = ones(6400,1);
for i=1:length(indpup)
  detector(indpup(i))=0;
end  

%illuminate pupils    
%for i=0:79
%detector2d(1:80,i+1)=detector(i*80+1:(i+1)*80);
%end

detector2d = reshape(detector,80,80);

pupilmask = detector2d(1:40,1:40);



%%

for j = 1:length(pupilmask)
   r0(j) = j*max(pupilmask(j,:)==0);
   c0(j) = j*max(pupilmask(:,j)==0);    
end

rr = find(r0~=0);
cc = find(c0~=0);

pupilmask_trim = pupilmask(min(rr):max(rr) , min(cc):max(cc));
pupilmask_trim = padarray(pupilmask_trim,[1,1],1);


subplot(1,3,1)
imagesc(detector2d); daspect([1 1 1]); axis xy;
subplot(1,3,2)
imagesc(pupilmask); daspect([1 1 1]); axis xy;
subplot(1,3,3)
imagesc(pupilmask_trim); daspect([1 1 1]); axis xy;


end