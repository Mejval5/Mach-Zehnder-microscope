clear a
for i=1:120
  a(i,:,:)=getsnapshot(vid);
end

a=squeeze(sum(a,1))/60;
imagesc(a)