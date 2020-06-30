function [P2] = Phase_shift_5steps_2D(ref,circleSelect,modulationBool)


%% Read "5 Intensity images" and Stack in 3D array I
if nargin<3
    modulationBool = 0;
    if nargin<2
        circleSelect = 0;
        if nargin<1
            ref = 0;
        end
    end
end

I = load('PhaseImages','pic');
I = I.pic;
[m, n, ~]=size(I);
P = zeros(m, n); %extracted Phase array
C = zeros(m, n); %Complex number array(cos + i sin)$
if modulationBool
    Gamma = zeros(m, n);
end
if ref
    Iref = load('PhaseImages','picRef');
    Iref = Iref.picRef;
    Pref = zeros(m, n); %extracted Phase array
    Cref = zeros(m, n); %Complex number array(cos + i sin)
    
    
    for i=1:m
        for ii=1:n
            SIN_phase=2*((Iref(i,ii,4)-Iref(i,ii,2)));       %sin(theta)
            COS_phase=Iref(i,ii,1)+Iref(i,ii,5)-2*Iref(i,ii,3); %cos(theta)
            Pref(i,ii) = atan2(SIN_phase,COS_phase); % 5 step PSI phase calculation
            Cref(i,ii) = complex(COS_phase, SIN_phase); % Store complex number
        end
    end
    
end

%% Generation Phase array and Complex array(cosine_phase + i-sine_phase)
ifFive = size(I);
ifFive = ifFive(3)==5;

% Calculate each pixel of the array one by one
for i=1:m
    for ii=1:n
        if ifFive
            SIN_phase=2*((I(i,ii,4)-I(i,ii,2)));       %sin(theta)
            COS_phase=I(i,ii,1)+I(i,ii,5)-2*I(i,ii,3); %cos(theta)
            if modulationBool
                Gamma(i,ii) = 2/((I(i,ii,1)+I(i,ii,5))/2+I(i,ii,2)+I(i,ii,3)+I(i,ii,4))*sqrt((I(i,ii,2)-I(i,ii,4))^2+((I(i,ii,1)+I(i,ii,5))/2-I(i,ii,3))^2);
            end
        else
            SIN_phase=mean([I(i,ii,4), I(i,ii,8)])-mean([I(i,ii,2), I(i,ii,6)]);       %sin(theta)
            COS_phase=mean([I(i,ii,1), I(i,ii,5)])-mean([I(i,ii,3), I(i,ii,7)]); %cos(theta)
            if modulationBool
                Gamma(i,ii) = 2/(mean([I(i,ii,1), I(i,ii,5)])+mean([I(i,ii,2), I(i,ii,6)])+mean([I(i,ii,3), I(i,ii,7)])+mean([I(i,ii,4), I(i,ii,8)]))*sqrt((mean([I(i,ii,2), I(i,ii,6)])-mean([I(i,ii,4), I(i,ii,8)]))^2+(mean([I(i,ii,1), I(i,ii,5)])-mean([I(i,ii,3), I(i,ii,7)]))^2);
            end
        end
        P(i,ii) = atan2(SIN_phase,COS_phase); % 5 step PSI phase calculation
        C(i,ii) = complex(COS_phase, SIN_phase); % Store complex number
        
    end
end



%% Save Phase as "tiff" and Complex in MAT file
if modulationBool
    figure(55);
    imagesc(Gamma);
    title('Modulation')
    set(gca,'YDir','normal')
    axis xy
    axis equal
    xlabel('position [pixel]')
    ylabel('position [pixel]')
    colormap jet
    colorbar
end

figure(50);
C=angle(C);
imagesc(C);
title('Extracted Phase')
set(gca,'YDir','normal')
axis xy
axis equal
xlabel('position [pixel]')
ylabel('position [pixel]')
colormap jet
colorbar

if circleSelect
    
    h = drawcircle('Position',[1 100]);
    try
        [pos,rad] = customWait(h);
    catch
        P2 = P;
        return
    end
    
    [xgrid, ygrid] = meshgrid(1:size(P,2), 1:size(P,1));
    mask = ((xgrid-pos(1)).^2 + (ygrid-pos(2)).^2) <= rad.^2;
    P2 = P.*mask;
    
    square1x=max([int32(pos(1)-rad) 1]);
    square2x=min([int32(pos(1)+rad) size(P2,2)]);
    square1y=max([int32(pos(2)-rad) 1]);
    square2y=min([int32(pos(2)+rad) size(P2,1)]);
    P2 = P2(square1y:square2y,square1x:square2x);
    P2 (P2 == 0) = NaN;
    
    figure(50);
    imagesc(P2);
    title('Selected Phase');
    axis equal
    set(gca,'YDir','normal')
    colormap jet
    xlabel('position [pixel]')
    ylabel('position [pixel]')
    
else
    P2 = P;
end

if ref
    if circleSelect
        Pref2 = Pref.*mask;
        Pref2 = Pref2(square1y:square2y,square1x:square2x);
        Pref2 (Pref2 == 0) = NaN;
    else
        Pref2 = Pref;
    end
    
    figure(51);
    imagesc(Pref2);
    title('Selected Phase 2');
    axis equal
    colormap jet
    set(gca,'YDir','normal')
    
    
    save('ExtractPhase');
    disp('..............Press ENTER to continue..............')
    pause
    answer = questdlg('Would you like to continue?', 'Continue Menu', 'Yes','No','No');
    
    if strcmp(answer,'Yes')
        Goldstein_2D_unwrap(false,false,false,false,false,ref,circleSelect,false);
    end
else
    save('ExtractPhase');
end

%% Save calculated complex number data in mat file for 2D unwrapping

    function [pos,rad] = customWait(hROI)
        
        % Listen for mouse clicks on the ROI
        l = addlistener(hROI,'ROIClicked',@clickCallback);
        
        % Block program execution
        uiwait;
        
        % Remove listener
        delete(l);
        
        % Return the current position
        pos = hROI.Position;
        rad = hROI.Radius;
        
    end

    function clickCallback(~,evt)
        
        if strcmp(evt.SelectionType,'double')
            uiresume;
        end
        
    end

end