imagesc(P2)
axis equal
set(gca,'YDir','normal')
colormap jet
sizeC=size(P2);
h = drawline('Position',[1,sizeC(1)/2;sizeC(2),sizeC(1)/2],'DrawingArea',[1,1,sizeC(2),sizeC(1)]);
        try
            [pos] = customWait(h);
        catch
            return
        end
        pos = floor(pos);
        
        Xmin=max(min([pos(1,1), pos(2,1)]),1);
        Xmax=min(max([pos(1,1), pos(2,1)]),sizeC(2));
        Ymin=max(min([pos(1,2), pos(2,2)]),1);
        Ymax=min(max([pos(1,2), pos(2,2)]),sizeC(1));
        deltaX=Xmax-Xmin;
        deltaY=Ymax-Ymin;
        len = max(deltaX,deltaY);
        slope = (max(pos(1,2),1)-min(pos(2,2),sizeC(1)))/(min(pos(1,1),sizeC(2))-max(pos(2,1),1));
        C_line = zeros(len,1);
        for i=1:len
            if deltaX > deltaY
                xx=Xmin+i
                yy=Ymin+floor(i*slope)
                C_line(i)=P2(xx,yy);
            else
                yy=Ymin+i
                xx=Xmin+floor(i/slope)
                C_line(i)=P2(xx,yy);
            end
        end
        plot(C_line)
        
        
        
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