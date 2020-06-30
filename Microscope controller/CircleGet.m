imagesc(P2)
axis equal
set(gca,'YDir','normal')
colormap jet
h = drawcircle('Position',[1 100]);
    try
        [pos,rad] = customWait(h);
    catch
        P2 = P;
        return
    end
    
    [xgrid, ygrid] = meshgrid(1:size(P,2), 1:size(P,1));
    mask = ((xgrid-pos(1)).^2 + (ygrid-pos(2)).^2) <= rad.^2;
    P2 = P2.*mask;
    
    square1x=max([int32(pos(1)-rad) 1]);
    square2x=min([int32(pos(1)+rad) size(P2,2)]);
    square1y=max([int32(pos(2)-rad) 1]);
    square2y=min([int32(pos(2)+rad) size(P2,1)]);
    P2 = P2(square1y:square2y,square1x:square2x);
    
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
