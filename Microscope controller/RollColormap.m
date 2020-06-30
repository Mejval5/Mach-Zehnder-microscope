
base = jet(256);
base2 = hsv(256);
base3 = winter(256);
base4 = hot(256);
base5 = colorcube(256);
base6 = parula(256);

base7 = pink(256);
base8 = hsv(256);
base9 = prism(256);
base10 = gray(256);
base11 = autumn(256);
base12 = lines(256);

base13 = cool(256);
base14 = spring(256);
base15 = bone(256);
base16 = copper(256);
base17 = jet(256).^2;
base18 = parula(256).^2;

holder = base;
filename = 'im';

for n = (128*6):(128*18)
    for k=1:256
        l = rem(k+n*2,256)
        if l == 0, l = 256; end
        choice = (k+2*n)/256;
        if choice < 2
            holder(k,:)=base(l,:);
        elseif choice<3
            holder(k,:)=base2(l,:);
        elseif choice<4
            holder(k,:)=base3(l,:);
%         elseif choice<5
%             holder(k,:)=base4(l,:);
%         elseif choice<6
%             holder(k,:)=base5(l,:);
%         elseif choice<7
%             holder(k,:)=base6(l,:);
%         elseif choice<8
%             holder(k,:)=base7(l,:);
%         elseif choice<9
%             holder(k,:)=base8(l,:);
%         elseif choice<10
%             holder(k,:)=base9(l,:);
%         elseif choice<11
%             holder(k,:)=base10(l,:);
%         elseif choice<12
%             holder(k,:)=base11(l,:);
%         elseif choice<13
%             holder(k,:)=base12(l,:);
%         elseif choice<14
%             holder(k,:)=base13(l,:);
%         elseif choice<15
%             holder(k,:)=base14(l,:);
%         elseif choice<16
%             holder(k,:)=base15(l,:);
%         elseif choice<17
%             holder(k,:)=base16(l,:);
%         elseif choice<18
%             holder(k,:)=base17(l,:);
%         elseif choice<19
%             holder(k,:)=base18(l,:);
        end
    end
    
    
    image = (data.phasePic+pi)/2/pi*256;
    
    name = strcat('../gifs/',filename,num2str(n),'.png');
    imwrite(image,holder,name);
    
    %     name = strcat('../gifs/',filename,'.gif');
    %       if n == 1
    %           imwrite(image,holder,name,'gif', 'Loopcount',inf,'DelayTime',0.01);
    %       else
    %           imwrite(image,holder,name,'gif','WriteMode','append');
    %       end
    
end