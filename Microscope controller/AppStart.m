function AppStart

load('config','config')
% s = serial('COM6');
% pause(1);
% SendCommand('RS');
% pause(1);
% SendCommand('1OR');

% Handles for the camera
imaqreset
% vid = videoinput('pointgrey', 1,'F7_Mono8_608x512_Mode5');
% vid2 = videoinput('winvideo', 1);
% src = getselectedsource(vid);
% 
% % Config for camera
src.Brightness = config.cameraConfig.Brightness;
src.Exposure = config.cameraConfig.Exposure;
src.ExposureMode = config.cameraConfig.ExposureMode;
src.FrameRate = config.cameraConfig.FrameRate;
src.FrameRateMode = config.cameraConfig.FrameRateMode;
src.Gain = config.cameraConfig.Gain;
src.GainMode = config.cameraConfig.GainMode;
src.Shutter = config.cameraConfig.Shutter;
src.ShutterMode = config.cameraConfig.ShutterMode;

% RampValueChanged(config.piezo.piezoRampValue);

% Data loading hacks
data.image = 0;
data.loadDataBool = 0;
data.phaseBool = 0;
data.FarFieldBool = 0;
data.NearFieldBool = 0;
updateWait = 0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Button functions to use camera.

% Value changed function: ExposureMode
    function ExposureModeValueChanged(app, ~)
        value = app.Value;
        updateWait  = 0.3;
        if value
            src.ExposureMode = 'Auto';
        else
            src.ExposureMode = 'Manual';
        end
        
    end

% Value changed function: GainMode
    function GainModeValueChanged(app, ~)
        value = app.Value;
        updateWait  = 0.3;
        if value
            src.GainMode = 'Auto';
        else
            src.GainMode = 'Manual';
        end
        
    end

% Value changed function: FrameRateMode
    function FrameRateModeValueChanged(app, ~)
        value = app.Value;
        updateWait  = 0.3;
        if value
            src.FrameRateMode = 'Auto';
        else
            src.FrameRateMode = 'Manual';
        end
        
    end

% Value changed function: ShutterMode
    function ShutterModeValueChanged(app, ~)
        value = app.Value;
        updateWait  = 0.3;
        if value
            src.ShutterMode = 'Auto';
        else
            src.ShutterMode = 'Manual';
        end
        
    end

% Value changed function: Brightness
    function BrightnessValueChanged(app, ~)
        updateWait  = 0.3;
        src.Brightness = BoundValue(app.Value,'Brightness');
        
    end

% Value changed function: Exposure
    function ExposureValueChanged(app, ~)
        updateWait  = 0.3;
        src.ExposureMode = 'Manual';
        src.Exposure = BoundValue(app.Value,'Exposure');
        ExposureMode.Value = false;
        
    end

% Value changed function: FrameRate
    function FrameRateValueChanged(app, ~)
        updateWait  = 0.3;
        src.FrameRateMode = 'Manual';
        src.FrameRate = BoundValue(app.Value,'FrameRate');
        FrameRateMode.Value = false;
        
    end

% Value changed function: Gain
    function GainValueChanged(app, ~)
        updateWait  = 0.3;
        src.GainMode = 'Manual';
        src.Gain = BoundValue(app.Value,'Gain');
        GainMode.Value = false;
        
    end

% Value changed function: Shutter
    function ShutterValueChanged(app, ~)
        updateWait  = 0.3;
        src.ShutterMode = 'Manual';
        src.Shutter = BoundValue(app.Value,'Shutter');
        ShutterMode.Value = false;
        
    end

% Control behaviour of option buttons for the measuring functions
    function calcZernike_Change(~, ~)
        if calcZernike.Value == false
            visualCoefs.Value = false;
            showCoefs.Value = false;
            subCoefs.Value = false;
        end
    end

    function coefs_Change(app, ~)
        if app.Value
            calcZernike.Value = true;
        end
    end

% SnapshotButton with subtracting background selection
    function a = SnapshotButtonPushed(~, ~)
        if load_BG.Value
            b = getsnapshot(vid)-config.BG;
            b(b<0) = 0;
            a = b;
        else
            a = getsnapshot(vid);
        end
    end

% Get background to use later
    function GetBG(~, ~)
        config.BG = getsnapshot(vid);
        save('config','config');
    end

    function GetFarField(~, ~)
        data.farField = getsnapshot(vid);
        data.FarFieldBool = 1;
        figure(99);
        imagesc(data.farField);
        title('Far field data');
        axis equal
        set(gca,'YDir','normal')
        colormap bone
        FarFieldButton.BackgroundColor = [.96 .51 .45];
    end

    function GetNearField(~, ~)
        data.nearField = getsnapshot(vid);
        data.NearFieldBool = 1;
        figure(100);
        imagesc(data.nearField);
        title('Near field data');
        axis equal
        set(gca,'YDir','normal')
        colormap bone
        NearFieldButton.BackgroundColor = [.96 .51 .45];
    end


% Livepreview for MEASURE camera
    function LivepreviewButtonPushed(~, ~)
        updateWait  = 0.3;
        preview(vid);
        
    end

% Livepreview for TIP camera
    function LivepreviewButtonPushed2(~, ~)
        %         preview(vid2);
    end

% Change camera mode
    function CameraModeDropDownValueChanged(app, ~)
        updateWait  = 0.3;
        value = app.Value;
        delete(vid);
        vid = videoinput('pointgrey', 1, value);
        src = getselectedsource(vid);
        
    end

% New data
    function NewDataPushed(~, ~)
        data.loadDataBool = 0;
        data.phaseBool = 0;
        CurrentDataText.Text = '*new data*';
        data.FarFieldBool = 0;
        data.NearFieldBool = 0;
        fiveStepButton.BackgroundColor = [.96 .96 .96];
        NearFieldButton.BackgroundColor = [.96 .96 .96];
        FarFieldButton.BackgroundColor = [.96 .96 .96];
    end

% Save data
    function SaveDataPushed(~, ~)
        if (data.phaseBool && sum(SaveDataList.Value == 1,'all') || data.NearFieldBool && sum(SaveDataList.Value == 2,'all') || data.FarFieldBool && sum(SaveDataList.Value == 3,'all'))
            [file,path] = uiputfile('../Measurement/filename');
            if file ~= 0
                a = strsplit(file,'.');
                dataName = cell2mat(a(1));
                mkdir(strcat(path,dataName));
                
                save(strcat(path,dataName,'\data.mat'),'data')
                
                if data.phaseBool && sum(SaveDataList.Value == 1,'all')
                    imwrite((data.phasePic+pi)/2/pi*256,jet(256),strcat(path,dataName,'\PhaseImage.png'));
                end
                
                if data.NearFieldBool && sum(SaveDataList.Value == 2,'all')
                    imwrite((data.nearField),bone(256),strcat(path,dataName,'\NearField.png'));
                end
                
                if data.FarFieldBool && sum(SaveDataList.Value == 3,'all')
                    imwrite((data.farField),bone(256),strcat(path,dataName,'\FarField.png'));
                end
                
                SavedDataText.Text = dataName;
                NewDataPushed
            end
        else
            disp('No DATA')
        end
    end

% Load data
    function LoadDataPushed(~, ~)
        path = uigetdir('../Measurement/','Select data folder');
        if path ~= 0
            b = strsplit(path,'\');
            dataName = cell2mat(b(end));
            filePath = strcat(path,'\data.mat');
            
            data=load(filePath,'data');
            UIFigure.Visible = 'off';
            UIFigure.Visible = 'on';
            data = data.data;
            data.loadDataBool = 1;
            CurrentDataText.Text = dataName;
            if data.phaseBool
                fiveStepButton.BackgroundColor = [.54 .81 .94];
            end
            if data.FarFieldBool
                NearFieldButton.BackgroundColor = [.54 .81 .94];
            end
            if data.NearFieldBool
                FarFieldButton.BackgroundColor = [.54 .81 .94];
            end
        end
    end

% Preview data
    function PreviewDataPushed(~, ~)
        if data.phaseBool
            figure(50);
            imagesc(data.phasePic);
            title('Preview phase data');
            axis equal
            set(gca,'YDir','normal')
            colormap jet
            colorbar
        end
        
        if data.FarFieldBool
            figure(99);
            imagesc(data.farField);
            title('Preview far field data');
            axis equal
            set(gca,'YDir','normal')
            colormap bone
        end
        
        if data.NearFieldBool
            figure(100);
            imagesc(data.nearField);
            title('Preview near field data');
            axis equal
            set(gca,'YDir','normal')
            colormap bone
        end
        
        if (~data.phaseBool && ~data.FarFieldBool && ~data.FarFieldBool)
            disp('No DATA')
        end
    end

% 5-step phase shift algorithm
    function fiveStepAlgorithm(~, ~)
        if ~data.loadDataBool || (~data.phaseBool && data.loadDataBool)
            ref = getError.Value;
            step = config.piezo.voltagePiezoFiveStep(2)-config.piezo.voltagePiezoFiveStep(1);
            rampValue = RampValue;
            testPic = SnapshotButtonPushed;
            pic = zeros(size(testPic,1),size(testPic,2),8);
            
            for i=1:8
                MoveTo(config.piezo.voltagePiezoFiveStep(i));
                pause(2*step/(rampValue*1000))
                pic(:,:,i) = SnapshotButtonPushed;
            end
            
            SendCommand('1PA0.00')
            
            if ref
                pause(refWait.Value);
                picRef = zeros(size(testPic,1),size(testPic,2),8);
                
                for i=1:8
                    MoveTo(config.piezo.voltagePiezoFiveStep(i));
                    pause(2*step/(rampValue*1000))
                    picRef(:,:,i) = SnapshotButtonPushed;
                end
                save('PhaseImages','pic','picRef')
            else
                save('PhaseImages','pic')
            end
            
            SendCommand('1PA0.00')
            pause(0.1);
        else
            ref = 0;
            pic = data.image;
            save('PhaseImages','pic')
        end
        
        data.image = pic;
        [data.phasePic] = Phase_shift_5steps_2D(ref,circleSelect.Value,Modulation.Value);
        data.phaseBool = 1;
        fiveStepButton.BackgroundColor = [.96 .51 .45];
    end

% Unwrap phase
    function UnwrapphaseButtonPushed(~, ~)
        if data.loadDataBool && data.phaseBool
            P2 = data.phasePic;
            save('ExtractPhase','P2');
        end
        Goldstein_2D_unwrap(calcZernike.Value,subCoefs.Value,showCoefs.Value,visualCoefs.Value,subtractList.Value,0,circleSelect.Value,linePlot.Value);
    end

% Calibrate piezo voltages for 5-step algorithm
    function CalibratePiezo(~, ~)
        ParamsBox = uifigure;
        ParamsBox.Position = [100 100 400 300];
        ParamsBox.Name = 'Calibration parameters';
        
        % Create HorizontalMode checkbox
        HorizontalMode = uicheckbox(ParamsBox);
        HorizontalMode.Text = 'Horizontal';
        HorizontalMode.Position = [220 205 100 20];
        HorizontalMode.Value = true;
        
        % Create Number of steps Label
        NumberLabel = uilabel(ParamsBox);
        NumberLabel.HorizontalAlignment = 'left';
        NumberLabel.Position = [30 205 150 22];
        NumberLabel.Text = 'Number of steps';
        
        % Create Number of steps box
        Number = uieditfield(ParamsBox, 'numeric');
        Number.Position = [150 200 50 30];
        Number.Value = 10;
        
        % Create VoltageRange Label
        VoltageRangeLabel = uilabel(ParamsBox);
        VoltageRangeLabel.HorizontalAlignment = 'left';
        VoltageRangeLabel.Position = [30 155 150 22];
        VoltageRangeLabel.Text = 'Voltage MIN/MAX';
        
        % Create Voltage Min box
        VoltageMin = uieditfield(ParamsBox, 'numeric');
        VoltageMin.Position = [150 150 50 30];
        VoltageMin.Value = 3;
        
        % Create Voltage Max box
        VoltageMax = uieditfield(ParamsBox, 'numeric');
        VoltageMax.Position = [210 150 50 30];
        VoltageMax.Value = 10;
        
        % Create CalibratePiezo Button
        CalibratePiezoParamsButton = uibutton(ParamsBox, 'push');
        CalibratePiezoParamsButton.ButtonPushedFcn = @(CalibratePiezoParamsButton,event) Calibrate5step(CalibratePiezoParamsButton,event,Number.Value,VoltageMin.Value,VoltageMax.Value,HorizontalMode.Value);
        CalibratePiezoParamsButton.Position = [145 50 110 44];
        CalibratePiezoParamsButton.Text = 'Calibrate';
    end

% Calibrating function
    function Calibrate5step(~, ~,N,Vstart,Vend,horizontal)
        
        Vrange = linspace(Vstart,Vend,N+1);
        Vstep = (max(Vrange)-min(Vrange))/N;
        sizeIm = size(getsnapshot(vid));
        imageHolder = zeros(sizeIm(1),sizeIm(2),N+1);
        rampValue = RampValue;
        MoveTo(Vstart)
        pause(2*Vstart/(rampValue*1000))
        for i=1:(N+1)
            imageHolder(:,:,i) = getsnapshot(vid);
            if i < N+1
                MoveTo(Vrange(i+1));
                pause(2*Vstep/(rampValue*1000))
            end
        end
        MoveTo(0.0);
        
        param = 2;
        if horizontal
            param = 1;
        end
        
        summedHolder = squeeze(sum(imageHolder,param));
        sizeSum = size(summedHolder);
        summedHolder = summedHolder(floor(sizeSum/3):floor(sizeSum/3*2),:);
        xx = size(summedHolder);
        xx = 1:xx(1);
        Pholder = zeros(1,(N+1));
        PholderRelative = zeros(1,(N+1));
        for i=1:(N+1)
            [SineP]=sineFit(xx,summedHolder(:,i)');
            Pholder(i) = SineP(1,4);
        end
        
        piHolder = 0;
        PholderRelative(1,1) = Pholder(1,1);
        for i=1:(N)
            if abs(Pholder(1,i)-Pholder(1,i+1))>pi
                piHolder = piHolder-sign(Pholder(1,i+1)-Pholder(1,i))*2*pi;
            end
            PholderRelative(1,i+1) = Pholder(1,i+1)+piHolder;
        end
        
        figure(20)
        p = polyfit(PholderRelative,Vrange,1);
        phaseNeeded = [0,pi/2,pi,pi*3/2,2*pi,5*pi/2,3*pi,7*pi/2];
        yyy = polyval(p,phaseNeeded);
        Ymin = min(yyy);
        yyy=yyy-Ymin+Vstart;
        plot(Vrange,PholderRelative-PholderRelative(1),yyy,phaseNeeded)
        config.piezo.voltagePiezoFiveStep = yyy;
        save('config','config');
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create UIFigure and components

% UIFigure
UIFigure = uifigure;
UIFigure.Position = [100 100 877 738];
UIFigure.Name = 'Über controller';

% Brightness label
BrightnessEditFieldLabel = uilabel(UIFigure);
BrightnessEditFieldLabel.HorizontalAlignment = 'left';
BrightnessEditFieldLabel.Position = [15 690 60 22];
BrightnessEditFieldLabel.Text = 'Brightness';

% Brightness slider
Brightness = uislider(UIFigure);
Brightness.ValueChangedFcn = @(Brightness,event) BrightnessValueChanged(Brightness,event);
Brightness.Position = [100 700 150 3];
Brightness.Value = src.Brightness;
% a = propinfo(src,'Brightness');
% Brightness.Limits = [a.ConstraintValue(1) a.ConstraintValue(2)];

% Exposure label
ExposureLabel = uilabel(UIFigure);
ExposureLabel.HorizontalAlignment = 'left';
ExposureLabel.Position = [15 640 60 22];
ExposureLabel.Text = 'Exposure';

% Exposure slider
Exposure = uislider(UIFigure);
Exposure.ValueChangedFcn = @(Exposure,event) ExposureValueChanged(Exposure,event);
Exposure.Position = [100 650 150 3];
Exposure.Value = src.Exposure;
% a = propinfo(src,'Exposure');
% Exposure.Limits = [a.ConstraintValue(1) a.ConstraintValue(2)];

% FrameRate Label
FrameRateLabel = uilabel(UIFigure);
FrameRateLabel.HorizontalAlignment = 'left';
FrameRateLabel.Position = [15 590 80 22];
FrameRateLabel.Text = 'Frame Rate';

% FrameRate slider
FrameRate = uislider(UIFigure);
FrameRate.ValueChangedFcn = @(FrameRate,event) FrameRateValueChanged(FrameRate,event);
FrameRate.Position = [100 600 150 3];
FrameRate.Value = src.FrameRate;
% a = propinfo(src,'FrameRate');
% FrameRate.Limits = [a.ConstraintValue(1) a.ConstraintValue(2)];

% Gain Label
GainLabel = uilabel(UIFigure);
GainLabel.HorizontalAlignment = 'left';
GainLabel.Position = [15 540 60 22];
GainLabel.Text = 'Gain';

% Create Gain slider
Gain = uislider(UIFigure);
Gain.ValueChangedFcn = @(Gain,event) GainValueChanged(Gain,event);
Gain.Position = [100 550 150 3];
Gain.Value = src.Gain;
% a = propinfo(src,'Gain');
% Gain.Limits = [a.ConstraintValue(1) a.ConstraintValue(2)];

% Shutter Label
ShutterLabel = uilabel(UIFigure);
ShutterLabel.HorizontalAlignment = 'left';
ShutterLabel.Position = [15 490 60 22];
ShutterLabel.Text = 'Shutter';

% Shutter slider
Shutter = uislider(UIFigure);
Shutter.ValueChangedFcn = @(Shutter,event) ShutterValueChanged(Shutter,event);
Shutter.Position = [100 500 150 3];
Shutter.Value = src.Shutter;
% a = propinfo(src,'Shutter');
% Shutter.Limits = [a.ConstraintValue(1) a.ConstraintValue(2)];

% ShutterPrecise Label
ShutterPreciseLabel = uilabel(UIFigure);
ShutterPreciseLabel.HorizontalAlignment = 'left';
ShutterPreciseLabel.Position = [15 435 150 22];
ShutterPreciseLabel.Text = 'Shutter precise';

% ShutterPrecise box
ShutterPrecise = uieditfield(UIFigure, 'numeric');
ShutterPrecise.ValueChangedFcn = @(Shutter,event) ShutterValueChanged(Shutter,event);
ShutterPrecise.Position = [150 430 100 30];
ShutterPrecise.Value = src.Shutter;

% Ramp Label
ShutterLabel = uilabel(UIFigure);
ShutterLabel.HorizontalAlignment = 'left';
ShutterLabel.Position = [15 395 100 22];
ShutterLabel.Text = 'Ramp [V/us]';

% Ramp box
Ramp = uieditfield(UIFigure, 'numeric');
Ramp.ValueChangedFcn = @(Ramp,event) RampValueChanged(Ramp.Value,event);
Ramp.Position = [150 390 100 30];
% Ramp.Value = RampValue;

% ExposureMode checkbox
ExposureMode = uicheckbox(UIFigure);
ExposureMode.ValueChangedFcn = @(ExposureMode,event) ExposureModeValueChanged(ExposureMode,event);
ExposureMode.Text = 'Auto';
ExposureMode.Position = [280 640 50 20];
ExposureMode.Value = AutoCheck(src.ExposureMode == 'Auto');

% FrameRateMode checkbox
FrameRateMode = uicheckbox(UIFigure);
FrameRateMode.ValueChangedFcn = @(FrameRateMode,event) FrameRateModeValueChanged(FrameRateMode,event);
FrameRateMode.Text = 'Auto';
FrameRateMode.Position = [280 590 50 20];
FrameRateMode.Value = AutoCheck(src.FrameRateMode == 'Auto');

% GainMode checkbox
GainMode = uicheckbox(UIFigure);
GainMode.ValueChangedFcn = @(GainMode,event) GainModeValueChanged(GainMode,event);
GainMode.Text = 'Auto';
GainMode.Position = [280 540 50 20];
GainMode.Value = AutoCheck(src.GainMode == 'Auto');

% ShutterMode checkbox
ShutterMode = uicheckbox(UIFigure);
ShutterMode.ValueChangedFcn = @(ShutterMode,event) ShutterModeValueChanged(ShutterMode,event);
ShutterMode.Text = 'Auto';
ShutterMode.Position = [280 490 50 20];
ShutterMode.Value = AutoCheck(src.ShutterMode == 'Auto');

% TemperatureGauge Label
TemperatureGaugeLabel = uilabel(UIFigure);
TemperatureGaugeLabel.HorizontalAlignment = 'center';
TemperatureGaugeLabel.Position = [743.5 468 73 22];
TemperatureGaugeLabel.Text = 'Temperature';

% TemperatureGauge
TemperatureGauge = uigauge(UIFigure, 'circular');
TemperatureGauge.Position = [720 500 120 120];
% TemperatureGauge.Value = src.Temperature-273;

% Get BG
GetBGButton = uibutton(UIFigure, 'push');
GetBGButton.ButtonPushedFcn = @(GetBGButton,event) GetBG(GetBGButton,event);
GetBGButton.Position = [460 506 100 44];
GetBGButton.Text = 'Get BG';

% Preview MEASURE camera
LivepreviewButton = uibutton(UIFigure, 'push');
LivepreviewButton.ButtonPushedFcn = @(LivepreviewButton,event) LivepreviewButtonPushed(LivepreviewButton,event);
LivepreviewButton.Position = [400 446 100 44];
LivepreviewButton.Text = 'Live preview';

% Preview TIP camera
LivepreviewButton2 = uibutton(UIFigure, 'push');
LivepreviewButton2.ButtonPushedFcn = @(LivepreviewButton2,event) LivepreviewButtonPushed2(LivepreviewButton2,event);
LivepreviewButton2.Position = [520 446 100 44];
LivepreviewButton2.Text = 'Preview tip';

% Calibrate Piezo Button
CalibratePiezoButton = uibutton(UIFigure, 'push');
CalibratePiezoButton.ButtonPushedFcn = @(CalibratePiezoButton,event) CalibratePiezo(CalibratePiezoButton,event);
CalibratePiezoButton.Position = [455 386 110 44];
CalibratePiezoButton.Text = 'Calibrate piezo';

% Camera MODES label
CameraModeDropDownLabel = uilabel(UIFigure);
CameraModeDropDownLabel.HorizontalAlignment = 'right';
CameraModeDropDownLabel.Position = [354 568 82 22];
CameraModeDropDownLabel.Text = 'Camera Mode';

% Camera MODES
CameraModeDropDown = uidropdown(UIFigure);
% CameraModeDropDown.Items = {'F7_Mono12_1224x1024_Mode1',	'F7_Mono12_1224x1024_Mode2',	'F7_Mono12_2448x2048_Mode0',	'F7_Mono12_608x512_Mode5',	'F7_Mono16_1224x1024_Mode1',	'F7_Mono16_1224x1024_Mode2',	'F7_Mono16_2448x2048_Mode0',	'F7_Mono16_608x512_Mode5',	'F7_Mono8_1224x1024_Mode1',	'F7_Mono8_1224x1024_Mode2',	'F7_Mono8_2448x2048_Mode0',	'F7_Mono8_608x512_Mode5',	'F7_Raw16_1224x1024_Mode1'	,'F7_Raw16_1224x1024_Mode2',	'F7_Raw16_2448x2048_Mode0',	'F7_Raw8_1224x1024_Mode1',	'F7_Raw8_1224x1024_Mode2',	'F7_Raw8_2448x2048_Mode0'};
CameraModeDropDown.Items = {'F7_Mono8_1224x1024_Mode1',	'F7_Mono8_1224x1024_Mode2',	'F7_Mono8_2448x2048_Mode0',	'F7_Mono8_608x512_Mode5'};
CameraModeDropDown.ValueChangedFcn = @(CameraModeDropDown,event) CameraModeDropDownValueChanged(CameraModeDropDown,event);
CameraModeDropDown.Position = [451 568 233 22];
% CameraModeDropDown.Value = vid.VideoFormat;

% Five step Button
fiveStepButton = uibutton(UIFigure, 'push');
fiveStepButton.ButtonPushedFcn = @(fiveStepButton,event) fiveStepAlgorithm(fiveStepButton,event);
fiveStepButton.Position = [116 318 110 44];
fiveStepButton.Text = '5-step algorithm';

% Unwrap phase Button
UnwrapphaseButton = uibutton(UIFigure, 'push');
UnwrapphaseButton.ButtonPushedFcn = @(UnwrapphaseButton,event) UnwrapphaseButtonPushed(UnwrapphaseButton,event);
UnwrapphaseButton.Position = [116 259 110 44];
UnwrapphaseButton.Text = 'Unwrap phase';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data management

% Current data label
CurrentDataLabel = uilabel(UIFigure);
CurrentDataLabel.HorizontalAlignment = 'right';
CurrentDataLabel.Position = [350 650 82 22];
CurrentDataLabel.Text = 'Current data: ';

% Current data text
CurrentDataText = uilabel(UIFigure);
CurrentDataText.HorizontalAlignment = 'left';
CurrentDataText.Position = [450 640 300 50];
CurrentDataText.Text = '*new data*';
CurrentDataText.FontWeight = 'bold';
CurrentDataText.FontSize = 22;

% Current data label
SavedDataLabel = uilabel(UIFigure);
SavedDataLabel.HorizontalAlignment = 'right';
SavedDataLabel.Position = [350 620 82 22];
SavedDataLabel.Text = 'Saved data: ';

% Saved data text
SavedDataText = uilabel(UIFigure);
SavedDataText.HorizontalAlignment = 'left';
SavedDataText.Position = [450 610 300 50];
SavedDataText.Text = '';
SavedDataText.FontSize = 18;

% New data
NewDataButton = uibutton(UIFigure, 'push');
NewDataButton.ButtonPushedFcn = @(NewDataButton,event) NewDataPushed(NewDataButton,event);
NewDataButton.Position = [116 180 110 50];
NewDataButton.Text = 'New data';

% Load data
LoadDataButton = uibutton(UIFigure, 'push');
LoadDataButton.ButtonPushedFcn = @(LoadDataButton,event) LoadDataPushed(LoadDataButton,event);
LoadDataButton.Position = [116 120 110 44];
LoadDataButton.Text = 'Load data';

% Preview data
PreviewDataButton = uibutton(UIFigure, 'push');
PreviewDataButton.ButtonPushedFcn = @(PreviewDataButton,event) PreviewDataPushed(PreviewDataButton,event);
PreviewDataButton.Position = [240 120 110 44];
PreviewDataButton.Text = 'Preview data';

% Near Field
NearFieldButton = uibutton(UIFigure, 'push');
NearFieldButton.ButtonPushedFcn = @(NearFieldButton,event) GetNearField(NearFieldButton,event);
NearFieldButton.Position = [240 180 110 50];
NearFieldButton.Text = 'Near Field';

% Far Field
FarFieldButton = uibutton(UIFigure, 'push');
FarFieldButton.ButtonPushedFcn = @(FarFieldButton,event) GetFarField(FarFieldButton,event);
FarFieldButton.Position = [364 180 110 50];
FarFieldButton.Text = 'Far Field';

% Save data
SaveDataButton= uibutton(UIFigure, 'push');
SaveDataButton.ButtonPushedFcn = @(SaveDataButton,event) SaveDataPushed(SaveDataButton,event);
SaveDataButton.Position = [116 60 110 44];
SaveDataButton.Text = 'Save data';

% Select which data to save
SaveDataList = uilistbox(UIFigure);
SaveDataList.Position = [250 45 100 57];
SaveDataList.Multiselect = 'on';
SaveDataList.Items = {'Phase Image', 'Near Field', 'Far Field'};
SaveDataList.ItemsData = [1, 2, 3];
SaveDataList.Value = [1, 2, 3];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Checkboxes for measuring functions' options

%%%% 5-set options
% load_BG checkbox
load_BG = uicheckbox(UIFigure);
load_BG.Text = 'Subtract saved background';
load_BG.Position = [250 345 200 20];
load_BG.Value = false;

% calculate error checkbox
getError = uicheckbox(UIFigure);
getError.Text = 'Find error with';
getError.Position = [450 345 200 20];
getError.Value = false;

% Camera Delay for refence label
CameraModeDropDownLabel = uilabel(UIFigure);
CameraModeDropDownLabel.HorizontalAlignment = 'right';
CameraModeDropDownLabel.Position = [580 345 50 20];
CameraModeDropDownLabel.Text = 'delay [s]';

% Camera Delay for refence numeric field
refWait = uieditfield(UIFigure, 'numeric');
refWait.Position = [550 348 30 20];
refWait.Value = 0;

% Circle Select checkbox
circleSelect = uicheckbox(UIFigure);
circleSelect.Text = 'Select circular area';
circleSelect.Position = [250 315 200 20];
circleSelect.Value = false;

% Modulation checkbox
Modulation = uicheckbox(UIFigure);
Modulation.Text = 'Show modulation';
Modulation.Position = [450 315 200 20];
Modulation.Value = false;

%%%% Unwrap options
% Calculate Zernike checkbox
calcZernike = uicheckbox(UIFigure);
calcZernike.ValueChangedFcn = @calcZernike_Change;
calcZernike.Text = 'Calculate Zernike coefficients';
calcZernike.Position = [250 285 250 20];
calcZernike.Value = false;

% Subtract piston, tilt and defocus coefficients checkbox
subCoefs = uicheckbox(UIFigure);
subCoefs.ValueChangedFcn = @coefs_Change;
subCoefs.Text = 'Subtract these terms:';
subCoefs.Position = [450 255 250 20];
subCoefs.Value = false;

% Show coefficients checkbox
showCoefs = uicheckbox(UIFigure);
showCoefs.ValueChangedFcn = @coefs_Change;
showCoefs.Text = 'Show coefficients';
showCoefs.Position = [250 255 200 20];
showCoefs.Value = false;

% Visual coefficients checkbox
visualCoefs = uicheckbox(UIFigure);
visualCoefs.ValueChangedFcn = @coefs_Change;
visualCoefs.Text = 'Show coefficients visually';
visualCoefs.Position = [450 285 200 20];
visualCoefs.Value = false;

% Line plot checkbox
linePlot = uicheckbox(UIFigure);
linePlot.Text = 'Line plot select';
linePlot.Position = [650 285 200 20];
linePlot.Value = false;

% Select which coeffs to subtract
subtractList = uilistbox(UIFigure);
subtractList.Position = [590 192 100 78];
subtractList.Multiselect = 'on';
subtractList.Items = {'Piston', 'TiltX', 'TiltY', 'Defocus'};
subtractList.ItemsData = [1, 2, 3, 5];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helpful functions

% Return bounded value clipped by camera constraints
    function y = BoundValue(x,name)
        updateWait  = 0.3;
        a = propinfo(src,name);
        y=min(max(x,a.ConstraintValue(1)),a.ConstraintValue(2));
        
    end

% convert number to boolean
    function y = AutoCheck(x)
        if x == 1
            y = true;
        else
            y = false;
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Piezo commands

% Send commands to piezo
    function SendCommand(string)
        fopen(s);
        fprintf(s,string);
        fclose(s);
    end

% Send callbacks to piezo
    function out = ReadCommand(string)
        fopen(s);
        fprintf(s,string);
        out = fscanf(s);
        fclose(s);
    end

% Move piezo to absolute voltage
    function MoveTo(voltage)
        value = num2str(voltage);
        command = strcat('1PA', value);
        SendCommand(command)
    end


% Value changed function: Ramp
    function RampValueChanged(value, ~)
        value = num2str(value);
        command = strcat('1VA', value);
        SendCommand(command);
        Ramp.Value = RampValue;
        config.piezo.piezoRampValue = RampValue;
    end


% Read value of ramp
    function out = RampValue(~, ~)
        read = ReadCommand('1VA?');
        out = sscanf(read, '1VA%f');
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Keep the values in program updated and handle app closing
    function ReloadValues()
        while 1
            pause(0.1);
            
            if ishghandle(UIFigure) == 0
                save('config','config');
                delete(vid);
                %                 delete(vid2);
                instrreset;
                break
            end
            if updateWait <= 0
                try
                    if updateWait > 0, continue; end
                    TemperatureGauge.Value = src.Temperature-273;
                    if updateWait > 0, continue; end
                    Brightness.Value = src.Brightness;
                    if updateWait > 0, continue; end
                    if contains(src.ExposureMode,'Auto'), Exposure.Value = src.Exposure; end
                    if updateWait > 0, continue; end
                    if contains(src.FrameRateMode,'Auto'), FrameRate.Value = src.FrameRate; end
                    if updateWait > 0, continue; end
                    if contains(src.GainMode,'Auto'), Gain.Value = src.Gain; end
                    if updateWait > 0, continue; end
                    if contains(src.ShutterMode,'Auto'), Shutter.Value = src.Shutter; end
                    if updateWait > 0, continue; end
                    if contains(src.Shutter,'Auto'), ShutterPrecise.Value = src.Shutter; end
                    if updateWait > 0, continue; end
                catch
                    continue
                end
            else
                updateWait = updateWait - 0.1;
            end
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run reloading function
ReloadValues();
end