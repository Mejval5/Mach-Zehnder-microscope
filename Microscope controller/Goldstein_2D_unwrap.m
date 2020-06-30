function Goldstein_2D_unwrap(calcZernike, subCoefs, showCoefs, visualCoefs,subtractList,ref,circleSelect,linePlotBool)

% config
if nargin < 7
    linePlotBool = false;
    if nargin <6
        ref = 0;
        if nargin <5
            subtractList = [1,2,3,5];
            if nargin < 4
                visualCoefs = false;
                if nargin < 3
                    showCoefs = false;
                    if nargin < 2
                        subCoefs = false;
                        
                        if nargin < 1
                            calcZernike = false;
                        end
                    end
                end
            end
        end
    end
end


%% Load dependencies
addpath('Unwrap');

%% Load Complex image from .mat file

C = load('ExtractPhase.mat','P2');
C = C.P2;
C (isnan(C))  = 0;
sizeC = size(C);

%% complex field to amplitude and pahse
C_phase=C;                 % calculated phase from complex data C

%% Creat binary mask
% It is necessary in PhaseUnwrapFunction)
IM_mask=ones(size(C));

%%  Set parameters (default, MSK does not change anything here)
max_box_radius=4;                           %Maximum search box radius (pixels)

%% Unwrap
residue_charge=PhaseResidues(C_phase, IM_mask);                            %Calculate phase residues
branch_cuts=BranchCuts(residue_charge, max_box_radius, IM_mask);            %Place branch cuts
[C_unwrapped, ~, ~, coords]=FloodFill(C_phase, branch_cuts, IM_mask);   %Flood fill phase unwrapping
C_unwrapped=C_unwrapped(2:sizeC(1)-1,2:sizeC(2)-1); % remove zeros from edges
sizeC = size(C_unwrapped);

%% Display results

if ~ref
    
    [xgrid, ygrid] = meshgrid(1:sizeC(2), 1:sizeC(1));
    mask = ((xgrid-sizeC(2)/2).^2 + (ygrid-sizeC(1)/2).^2) <= (min([sizeC(1) sizeC(2)])/2)^2;
    C_plot = C_unwrapped.*mask;
    C_plot(C_plot == 0) = NaN;
    
    figure(2)
    mesh(C_plot),
    colormap(gray),
    title('Unwrapped phase');
    axis tight
    colormap jet
    colorbar;
    
    if linePlotBool
        figure(86);
        imagesc(C);
        title('Phase map');
        axis equal
        colormap jet
        set(gca,'YDir','normal')
        h = drawline('Position',[1,sizeC(1)/2;sizeC(2),sizeC(1)/2],'DrawingArea',[1,1,sizeC(2),sizeC(1)]);
        try
            [pos] = customWait(h);
        catch
            return
        end
        pos = floor(pos);
        
        Xmin=max(min([pos(1,1), pos(1,2)]),1);
        Xmax=min(max([pos(1,1), pos(1,2)]),sizeC(2));
        Ymin=max(min([pos(2,1), pos(2,2)]),1);
        Ymax=min(max([pos(2,1), pos(2,2)]),sizeC(1));
        deltaX=Xmax-Xmin;
        deltaY=Ymax-Ymin;
        len = max(deltaX,deltaY);
        slope = deltaY/deltaX;
        C_line = zeros(len);
        for i=1:len
            if deltaX > deltaY
                xx=Xmin+i;
                yy=Ymin+floor(i*slope);
                C_line(i)=C_unwrapped(xx,yy);
            else
                yy=Ymin+i;
                xx=Xmin+floor(i/slope);
                C_line(i)=C_unwrapped(xx,yy);
            end
        end
        plot(C_line)
        
    end
    
    
    
    if calcZernike
        
        coeffs = 1:15;
        coeffCalc = ZernikeCalc(coeffs, C_unwrapped, mask, 'standard');
        zernikeFit = ZernikeCalc(coeffs,coeffCalc,mask, 'standard');
        zernikeFitSummed = sum(zernikeFit,3);
        
        zernikeFitSummed = zernikeFitSummed.*mask;
        zernikeFitSummed (zernikeFitSummed == 0) = NaN;
        
        
        
        figure(3);
        mesh(zernikeFitSummed)
        colormap jet
        colorbar;
        axis tight
        title('Zernike fit');
        
        if showCoefs
            
            figure(4);
            title('Zernike fit coeff');
            plot(coeffs,coeffCalc,'*')
            ax = gca;
            ax.XGrid = 'on';
            ax.XTick = coeffs;
            bar(coeffCalc)
            
        end
        
        
        if visualCoefs
            
            figure(5);
            title('Zernike polynomials');
            
            
            subx = floor((size(coeffs,2)/4*3)^0.5);
            suby = ceil(subx*16/9);
            
            for i = coeffs
                subplot(subx,suby,i);
                zern = zernikeFit(:,:,i).*mask;
                zern(zern == 0) = NaN;
                mesh(zern)
                colormap jet
                colorbar;
                axis tight
                title(num2str(i-1));
            end
        end
        
        if subCoefs
            
            % C_subbed = C_unwrapped - zernikeFit(:,:,1) - zernikeFit(:,:,2) - zernikeFit(:,:,3) - zernikeFit(:,:,5);
            C_subbed = C_plot;
            for i=subtractList
                C_subbed = C_subbed - zernikeFit(:,:,i);
            end
            
            coeffCalcSub = ZernikeCalc(coeffs, C_subbed, mask, 'standard');
            zernikeFitSub = ZernikeCalc(coeffs,coeffCalcSub,mask, 'standard');
            % zernikeFitSummedSub = sum(zernikeFitSub,3);
            % zernikeFitSummedSub = zernikeFitSummedSub.*mask;
            % zernikeFitSummedSub (zernikeFitSummedSub == 0) = NaN;
            
            figure(6)
            mesh(C_subbed)
            colormap jet
            colorbar;
            axis tight
            title('Zernike fit Subbed');
            
            if showCoefs
                
                figure(4);
                title('Zernike subbed fit coeff');
                plot(coeffs,coeffCalcSub,'*')
                ax = gca;
                ax.XGrid = 'on';
                ax.XTick = coeffs;
                bar(coeffCalcSub)
                
            end
            
            
            if visualCoefs
                figure(7)
                
                for i = coeffs
                    subplot(subx,suby,i);
                    zern = zernikeFitSub(:,:,i).*mask;
                    zern(zern == 0) = NaN;
                    mesh(zern)
                    colormap jet
                    colorbar;
                    axis tight
                    title(num2str(i-1));
                end
                
                title('Zernike fit Subbed');
            end
        end
        
    end
    
    
else
    
    
    Pref2 = load('ExtractPhase.mat','Pref2');
    Pref2 = Pref2.Pref2;
    Pref2 (isnan(Pref2))  = 0;
    sizeCref = size(Pref2);
    C_phaseRef=Pref2;
    
    %% Creat binary mask
    % It is necessary in PhaseUnwrapFunction)
    IM_maskRef=ones(size(Pref2));
    
    %%  Set parameters (default, MSK does not change anything here)
    max_box_radiusRef=4;                           %Maximum search box radius (pixels)
    
    %% Unwrap
    residue_chargeRef=PhaseResidues(C_phaseRef, IM_maskRef);                            %Calculate phase residues
    branch_cutsRef=BranchCuts(residue_chargeRef, max_box_radiusRef, IM_maskRef);            %Place branch cuts
    [C_unwrappedRef, ~, ~]=FloodFill(C_phaseRef, branch_cutsRef, IM_maskRef, coords);   %Flood fill phase unwrapping
    C_unwrappedRef=C_unwrappedRef(2:sizeCref(1)-1,2:sizeCref(2)-1); % remove zeros from edges
    %     sizeCRef = size(C_unwrappedRef);
    
    if circleSelect
        [xgrid, ygrid] = meshgrid(1:sizeC(2), 1:sizeC(1));
        mask = ((xgrid-sizeC(2)/2).^2 + (ygrid-sizeC(1)/2).^2) <= (min([sizeC(1) sizeC(2)])/2)^2;
        save('test');
        coeffs = 1;
        coeffCalc = ZernikeCalc(coeffs, C_unwrapped, mask, 'standard');
        zernikeFit = ZernikeCalc(coeffs,coeffCalc,mask, 'standard');
        
        coeffCalcRef = ZernikeCalc(coeffs, C_unwrappedRef, mask, 'standard');
        zernikeFitRef = ZernikeCalc(coeffs,coeffCalcRef,mask, 'standard');
        
        C_orig = C_unwrapped.*mask;
        C_orig(C_unwrapped == 0) = NaN;
        C_ref = C_unwrappedRef.*mask;
        C_ref(C_unwrappedRef == 0) = NaN;
        
        
        C_orig = C_orig - zernikeFit;
        C_ref = C_ref - zernikeFitRef;
    else
        C_orig = C_unwrapped-mean(C_unwrapped,'all');
        C_ref = C_unwrappedRef-mean(C_unwrappedRef,'all');
    end
    
    C_diff = C_orig - C_ref;
    C_diff (isnan(C_diff))  = 0;
    RMS = 1/sqrt(2) * sqrt(sum(C_diff.^2,'all')/sizeC(1)/sizeC(2));
    C_diff(C_diff == 0) = NaN;
    C_diff(abs(C_diff) > RMS) = NaN;
    
    figure(8)
    mesh(C_diff)
    colormap jet
    colorbar;
    axis tight
    title(strcat('RMS: ',num2str(RMS)));
    
end


save('Unwrap');

    function [pos] = customWait(hROI)
        
        % Listen for mouse clicks on the ROI
        l = addlistener(hROI,'ROIClicked',@clickCallback);
        
        % Block program execution
        uiwait;
        
        % Remove listener
        delete(l);
        
        % Return the current position
        pos = hROI.Position;
        
    end

    function clickCallback(~,evt)
        
        if strcmp(evt.SelectionType,'double')
            uiresume;
        end
        
    end

end