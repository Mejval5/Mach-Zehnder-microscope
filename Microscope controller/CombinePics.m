mashedPic=zeros(512,608);
for i=1:512
    for ii=1:608
        picSelect=ceil(ii/(608/5));
        indexX=ceil(ii/5)+180;
        mashedPic(i,ii)=I(i,ii,picSelect);
    end
end
imagesc(mashedPic)