function [SineParams]=sineFit(x,y)
%Purpose: Estimation of noisy sine curve  by peak and zero cross positions
% and non linear ftting.
%
% Syntax:
%       [SineParams]=sineFit(x,y)       
%       Input: x and y values, y=offs+amp+sin(2*pi*f*x+phi)+noise
%       Output: SineParams(1): offset (offs)
%               SineParams(2): amplitude (amp)
%               SineParams(3): frequency (f)
%               SineParams(4): phaseshift (phi)
%       yOut=offs+amp*sin(2*pi*f*x+phi)
%
%You may want to comment/uncomment the last statement (PlotAnalysis) in this function.
%
% Example:
% % generate f(x)
    % x=-4:5;
    % y=1+2*(sin(2*pi*0.1*x+2)+0.3*randn(size(x)));%Sine + noise
    % [SineP]=sineFit(x,y)
    % figure;
    % xx=x(1):(x(end)-x(1))/222:x(end);%better resolution
    % plot(x,y,xx,SineP(1)+SineP(2)*sin(2*pi*SineP(3)*xx+SineP(4)));
    % %uncomment following line if you want to save y=f(x) and run it sineFitDemo
    % % save('xy.mat','x','y');
%
%Author: Peter Seibold

if nargin~=2 
    help sineFit
    return %fail
end

offsetE=mean(y);% Offset Estimated
ysmooth=y-offsetE;
numely=numel(y);
%find zero cross positions
ZPos = find(diff(sign(ysmooth)));
%Adjust offsetE if only one peak and one zero cross
[~,indxMaxy]=max(ysmooth);
[~,indxMiny]=min(ysmooth);
if numel(ZPos)==1 && (indxMaxy==1 || indxMiny==1 || ...
     indxMaxy==numely || indxMiny==numely)   
   offsetE=(y(1)+y(end))/2;
   ysmooth=y-offsetE;
   ZPos = find(diff(sign(ysmooth)));
end 
%Interpolate zero cross position and
%remove duplicates (happens when y =0 exactly)
ZPosInterpol=unique(ZPos+ysmooth(ZPos)./(ysmooth(ZPos)-ysmooth(ZPos+1)))';
diffZPosLowSum2=1;%dummy
diffZPosHighSum=0;%dummy
%Prepare some variables
if numel(ZPosInterpol)>2
    diffZPosS=sort(diff(ZPosInterpol));
    diffZPosLowSum2=2*sum(diffZPosS(1:round(end*0.5-0.5)));
    diffZPosHighSum=sum(diffZPosS(round(end*0.5+0.5):end));
    temp=ceil(0.1*numel(diffZPosS));
    %underweight last zero cross diffs
    diffZPosHighSumE=sum(diffZPosS(round(end*0.5+1.5):end-temp))...
        +diffZPosS(end-temp)*temp;    
end
%find Peaks
PeakPos = find(diff(sign(diff(ysmooth))))+1;
[PeakInterpol,PeakPosInterpol]=PeakInterpolation(ysmooth,PeakPos);
PeakPosRaw=[];%make sure variable exists
if numel(PeakPos)>0
    %Peaks and zeros in one array
    ZPosInterpol2=[ZPosInterpol,zeros(numel(ZPosInterpol),1),...
        (1:numel(ZPosInterpol))',zeros(numel(ZPosInterpol),1)];
    PeakPosInterpol2=[PeakPosInterpol',PeakInterpol',(1:numel(PeakPosInterpol))',ones(numel(PeakPosInterpol),1)];
    %Zero&Peakpos&Peakks:ZPposP(position of PeakPos&ZeroPos, value Peak,index,  0=Z & 1=P)
    ZPposP=sortrows([ZPosInterpol2;PeakPosInterpol2]);
    %same only with raw peakpos and raw peaks
    PeakPosPraw=[PeakPos',ysmooth(PeakPos)',(1:numel(PeakPosInterpol))'];
    PeakPosRaw=EliminateMultiplePeaks(ZPos,PeakPos,ysmooth);
end

%Determine if peaks are due to noise or to low sample rate (high frequency)
LowSampleRate=true;%if true, all peaks are part of sinus, no smooth
LowSampleRateProb=zeros(7,1);%propability measure for low sample rate
if numel(PeakPosInterpol)>5
    %check for dominating adjacent peaks
    [~,PeakMaxIndx]=max(PeakInterpol(2:end-1));
    PeakMaxIndx=PeakMaxIndx+1;
    [~,PeakMinIndx]=min(PeakInterpol(2:end-1));
    PeakMinIndx=PeakMinIndx+1;
    %check if the 2 max peaks are single between Z
    [PeakMax,PeakMaxIndxZP]=max(ZPposP(2:end-1,2));
    PeakMaxIndxZP=PeakMaxIndxZP+1;
    [PeakMin,PeakMinIndxZP]=min(ZPposP(2:end-1,2));
    PeakMinIndxZP=PeakMinIndxZP+1;    
    %left zero cross next to the max peak
    if ZPposP(PeakMaxIndxZP-1,4)==0
       LowSampleRateProb(2)=1;
    end   
    %right zero cross next to the max peak
    if ZPposP(PeakMaxIndxZP+1,4)==0
       LowSampleRateProb(2)=max(LowSampleRateProb(2)+0.5,1);
    end
    %left zero cross next to the min peak
    if ZPposP(PeakMinIndxZP-1,4)==0
       LowSampleRateProb(3)=1;
    end
    %right zero cross next to the min peak
    if ZPposP(PeakMinIndxZP+1,4)==0
       LowSampleRateProb(3)=max(LowSampleRateProb(3)+0.5,1);%1
    end
    
    %big left min peak next to max peak
    if PeakInterpol(PeakMaxIndx-1)+0.7*PeakInterpol(PeakMaxIndx)<0
       LowSampleRateProb(4)=LowSampleRateProb(4)+1;
    end
    %big right min peak next to max peak
    if PeakInterpol(PeakMaxIndx+1)+0.7*PeakInterpol(PeakMaxIndx)<0
       LowSampleRateProb(4)=LowSampleRateProb(4)*2+1;
    end
    %big left max peak next to min peak
    if PeakInterpol(PeakMinIndx-1)+0.7*PeakInterpol(PeakMinIndx)>0
       LowSampleRateProb(5)=1;
    end
    %big right max peak next to min peak
    if PeakInterpol(PeakMinIndx+1)+0.7*PeakInterpol(PeakMinIndx)>0
       LowSampleRateProb(5)=LowSampleRateProb(5)*2+1;
    end  
    %check for many big peak diff
    %and mean ysmooth on left and right side nearly equal
    diffPeaksS=sort(abs(diff(PeakInterpol)));
    center=(numely+1)/2;%find center zero cross
    meanLeft=abs(mean(ysmooth(1:floor(center))));
    meanRight=abs(mean(ysmooth(ceil(center):end)));
    maxPeakDiff2=(PeakMax-PeakMin)*0.1;
    if mean(diffPeaksS(floor(0.5*end):end))>0.3*(PeakMax-PeakMin)...
        && meanLeft<maxPeakDiff2 &&  meanRight<maxPeakDiff2   
       LowSampleRateProb(7)=1;
    else 
       LowSampleRateProb(7)=-1; 
    end
    %find number of peaks between zero crosses
    countPbetweenMax=0;
    countPbetween3=0;
    distance3=0; 
    LowSampleRateProb(6)=0;
    ZPosInterpolE=[1;ZPosInterpol;numel(y)];
    for i=1:numel(ZPosInterpolE)-1
        countPbetweenTemp=numel(PeakPos(PeakPos>ZPosInterpolE(i) &...
           PeakPos<ZPosInterpolE(i+1)));
       if countPbetweenTemp>2
           countPbetween3=countPbetween3+1;
           distance3=distance3+ZPosInterpolE(i+1)-ZPosInterpolE(i);
       end
        countPbetweenMax=max(countPbetweenMax,countPbetweenTemp);
    end
    if countPbetweenMax>3 || countPbetween3>3
        %one time more than 3 peaks or more than 3 times more than 2 peaks
        %between adjacent zero crosses
        LowSampleRateProb(6)=-3;
    end
   if distance3>0.4*numely
        LowSampleRateProb(6)=LowSampleRateProb(6)-1;
   end 
   if diffZPosLowSum2<diffZPosHighSum
        LowSampleRateProb(6)=LowSampleRateProb(6)-1;
   end     
    %check if half periods of raw peaks equal and bigger than 2 
    %and seldom extra peaks between zero crosses
    %improvement sine-param results: ~0.1%
    LowSampleRateProb(1)=1;%set to no influence
    if numel(ZPosInterpol)>2
        PposRawDiffS=sort(diff(PeakPos));
        PposRawDiffSm=PposRawDiffS(floor(0.5*end)+1);
        if PposRawDiffSm>1 ...
            && PposRawDiffS(round(0.3*end))+1>=PposRawDiffS(round(0.8*end))...
            && (countPbetween3-1)*20<numel(PeakPos)...
            && diffZPosHighSumE*1.1<diffZPosLowSum2
            LowSampleRateProb(1)=9;%force no smooth
        end
    end
    %Evaluate if sine is high frequency
    if sum(LowSampleRateProb)>2 
        LowSampleRate=true;%do not smooth if more than 5 peaks!
    else
        LowSampleRate=false;%allow smooth
    end
end

 %Determine if sine is noisy
NoiseLevel=0;
if numel(PeakPos)>numel(ZPos)+1 
    NoiseLevel=numel(PeakPos)/(numel(ZPos)+1);
end
if numel(ZPos)>2 && diffZPosS(1)*1.5<diffZPosS(end)
    NoiseLevel=NoiseLevel+diffZPosS(end)/diffZPosS(1);
end

%define filter value
filterValues=ones(4,1);
%special for low (noisy) period number (0 to 3)
if numel(ZPos)==1 && numel(PeakPosRaw)==1 && numel(PeakPos)>1
    filterValues(1)=min(abs(ZPos-PeakPosRaw)/3,3);  
    NoiseLevel=max(NoiseLevel,1);
elseif numel(ZPos)==1 && numel(PeakPosRaw)==2 && numel(PeakPosRaw)< numel(PeakPos)   
   filterValues(1)=min((PeakPosRaw(2)-PeakPosRaw(1))/2.5,3); 
   NoiseLevel=max(NoiseLevel,1);     
elseif numel(ZPos)==2 && numel(PeakPosRaw)< numel(PeakPos)   
   filterValues(1)=min((ZPos(2)-ZPos(1))/3,3); 
   NoiseLevel=max(NoiseLevel,1);
%check for adjacent peaks around Z
elseif numel(ZPos)==3 && numel(PeakPos)==3 
    diffPeakPos=diff(PeakPos);
    if diffPeakPos(1)>3*diffPeakPos(2)  || diffPeakPos(1)< 3*diffPeakPos(2)
        [~,indxPeakMax]=max(abs(PeakInterpol));
        filterValues(1)=min((PeakPos(indxPeakMax)-ZPos(2))/1.5,3);
    end
elseif numel(ZPos)==4 && numel(PeakPosRaw)< numel(PeakPos)
    %noisy, filter a little bit
    filterValues(1)=1.2;
    NoiseLevel=max(NoiseLevel,1);
end

if numel(ZPos)>2 && ~LowSampleRate
    if numel(ZPos)>3 &&...
        diffZPosLowSum2<diffZPosHighSum 
        %many small zero cross intervals
        NoiseLevel=max(NoiseLevel,1);
        filterValues(1)=min(diffZPosS(end)/2,3);
    end
    if sum(diffZPosS)<0.3*(numely-1)
        %zero crossDiff covers small part, propably only one zero cross
        NoiseLevel=max(NoiseLevel,1);
        filterValues(2)=(numely-1)/5;
    end
    if diffZPosS(end)>0.3*(numely-1)
        %One dominating zero cross: propably only one zero cross
        NoiseLevel=max(NoiseLevel,1);
        filterValues(3)=(numely-1)/5;
    end    
    if countPbetweenMax>4  
        NoiseLevel=max(NoiseLevel,1);
       filterValues(4)=diffZPosS(end)/3;
    elseif countPbetween3/(numel(ZPosInterpol)-1)>0.3 && ...
           distance3>0.5*(numel(ZPosInterpol)-1)
        NoiseLevel=max(NoiseLevel,1);
       filterValues(4)=diffZPosS(max(end-countPbetween3,1))/1.5;        
    end    
end
filterValue=max(filterValues);
%redo find zero crosses and peaks if function y is filtered:
if filterValue>1
    ysmooth=smooth1D(y,filterValue)-offsetE;
    dysmooth=diff(ysmooth);
    ZPos = find(diff(sign(ysmooth)));
    %Interpolate zero cross position
    ZPosInterpol=unique(ZPos+ysmooth(ZPos)./(ysmooth(ZPos)-ysmooth(ZPos+1)))';
    %find peaks
    PeakPos = find(diff(sign(dysmooth)))+1;
end
    %Eliminate multiple peaks between zero crosses
PeakPos=EliminateMultiplePeaks(ZPos,PeakPos,ysmooth);

%parabola interpolation of peak position
[PeakInterpol,PeakPosInterpol]=PeakInterpolation(ysmooth,PeakPos);
xTOindx=(x(end)-x(1))/(numely-1);%x/indx

%=================================================================
%Estimate frequency (period) and amplitude
if NoiseLevel==0 && numel(PeakPosInterpol)>2
    fToleranceHigh=1.1;
    fToleranceLow=0.9;
else
    fToleranceHigh=1.3;
    fToleranceLow=0.6;
end

fLowEe=0;%frequency Low Estimate extra
fHighEe=0;
if numel(ZPosInterpol)==1 && numel(PeakPosInterpol)==0
    %no peaks found, guess frequency and amplitude
    fLowE=1/((x(end)-x(1))*4);
    fHighE=fLowE*3;
    amplitudeE=(max(ysmooth)-min(ysmooth))/1.8;
    fE=fLowE*1.5;%fE for 1st esimate and phase
elseif numel(ZPosInterpol)==1 && numel(PeakPosInterpol)==1 
    periodindx=abs(PeakPosInterpol(1)-ZPosInterpol(1))*4;
     if NoiseLevel==0
        fE=0.8/(periodindx*xTOindx);
        fHighE=fE*1.5;
        fLowE=fE*0.8;  
        amplitudeE=abs(PeakInterpol(1))*1.5;
     else
        fE=0.6/(periodindx*xTOindx);
        fHighE=fE*1.9;
        fLowE=fE*0.8;
        amplitudeE=abs(PeakInterpol(1))*1.3; 
     end
elseif numel(ZPosInterpol)==1 && numel(PeakPosInterpol)==2
    periodPindx=(PeakPosInterpol(2)-PeakPosInterpol(1))*2;
    if NoiseLevel==0 || filterValue>1
        fE=1/(periodPindx*xTOindx);
        fLowE=fE*fToleranceLow;
        fHighE=fE*fToleranceHigh;  
        amplitudeE=sum(abs(PeakInterpol))/2;
    else
        fE=1/(periodPindx*xTOindx)*0.6;
        fLowE=fE*fToleranceLow;
        fHighE=fE*1.5;             
        amplitudeE=sum(abs(PeakInterpol))/1.5;
    end
elseif numel(ZPosInterpol)==2 && numel(PeakPosInterpol)==1 
    periodZindx=(ZPosInterpol(2)-ZPosInterpol(1))*2;
    fHighE=1/(periodZindx*xTOindx);
    fLowE=fHighE/2.5; 
    fE=fHighE/2;
    amplitudeE=abs(PeakInterpol)*1.3;
elseif numel(ZPosInterpol)==2 && numel(PeakPosInterpol)==2
    periodZindx=(ZPosInterpol(2)-ZPosInterpol(1))*2;
    periodPindx=(PeakPosInterpol(2)-PeakPosInterpol(1))*2;
    pBig=max(periodZindx,periodPindx);
    pSmall=min(periodZindx,periodPindx);
    fE=1/(pBig*xTOindx);
    fLowE=fE*fToleranceLow;
    fHighE=1/(pSmall*xTOindx)*fToleranceHigh; 
    amplitudeE=sum(abs(PeakInterpol))/2;
elseif numel(ZPosInterpol)==2 && numel(PeakPosInterpol)==3
    periodZindx=(ZPosInterpol(2)-ZPosInterpol(1))*2;
    periodPindx=PeakPosInterpol(3)-PeakPosInterpol(1);
    pBig=max(periodZindx,periodPindx);
    pSmall=min(periodZindx,periodPindx);
    fLowE=1/(pBig*xTOindx)*fToleranceLow;
    fHighE=1/(pSmall*xTOindx)*fToleranceHigh;
    fE=(fLowE+fHighE)/2;
    amplitudeE=sum(abs(PeakInterpol))/3;
elseif numel(ZPosInterpol)==3 && numel(PeakPosInterpol)==2 
    %check if one Peak very small
    PeakMax=max(abs(PeakInterpol));
    PeakMin=min(abs(PeakInterpol));
    if PeakMin<0.3*PeakMax
       ZPosInterpol(2)=[];%delete center Z, for phase calculation
       fLowE=1/((ZPosInterpol(2)-ZPosInterpol(1))*2*xTOindx)*fToleranceLow;
       fHighE=2.5*fLowE;
       fE=(fLowE+fHighE)/2;
       amplitudeE=PeakMax*1.5;
    else
        periodBig=max([PeakPosInterpol(2)-PeakPosInterpol(1),diffZPosS(2)])*2;
        periodSmall=min([PeakPosInterpol(2)-PeakPosInterpol(1),diffZPosS(1)])*2;
       fLowE=1/(periodBig*xTOindx)*fToleranceLow; 
       fHighE=1/(periodSmall*xTOindx)*fToleranceHigh; 
       fE=(fLowE+fHighE)/2;
       amplitudeE=abs((PeakMax+PeakMin))*0.5;
    end
elseif numel(PeakPosInterpol)==4 ...
        && abs(PeakInterpol(1))<abs(PeakInterpol(2))+abs(PeakInterpol(3))...
        && abs(PeakInterpol(1))*0.7<abs(PeakInterpol(4))...
        && abs(PeakInterpol(4))*0.7<abs(PeakInterpol(1))
    %no small peaks between big border peaks
    fLowE=1/((PeakPosInterpol(3)-PeakPosInterpol(1))*xTOindx)*0.60;
    fHighE=fLowE*2;
    fE=fLowE*1.67;
    amplitudeE=abs(PeakInterpol(1)-PeakInterpol(4))/2; 
    %check if largest peaks at borders
    [~,PeakInterpolSindx]=sort(abs(PeakInterpol));
    if PeakInterpolSindx(1)>1 && PeakInterpolSindx(1)<4 && ...
       PeakInterpolSindx(2)>1 && PeakInterpolSindx(2)<4 
        periodPindx=(PeakPosInterpol(4)-PeakPosInterpol(1))*2;
        periodZindx=diffZPosS(end)*2;
        fLowEe=1/(max(periodPindx,periodZindx)*xTOindx)*0.9;
        fHighEe=1/(min(periodPindx,periodZindx)*xTOindx)*1.1;
    end
elseif numel(ZPosInterpol)>2
    %find zero cross intervals of biggest peaks
        %treat peaks at border with only one Z cross
    Z2=ZPosInterpol2(:,1);
    if PeakPosPraw(1)<Z2(1)
      Z2=[2*PeakPosPraw(1,1)-Z2(1);Z2];
    end
    if PeakPosPraw(end,1)>Z2(end)
      Z2=[Z2;2*PeakPosPraw(end,1)-Z2(end)];
    end
    
    PeakPosPrawS=sortrows(abs(PeakPosPraw),2);
    PeakPosPrawS=PeakPosPrawS(round(0.8*end):end,1);
    ZdiffAtBigP=zeros(numel(PeakPosPrawS),1);
    for i=1:numel(PeakPosPrawS);
       ZdiffAtBigP(i)=Z2(find(Z2>PeakPosPrawS(i),1))...
           -Z2(find(Z2<PeakPosPrawS(i),1,'last'));
    end
        %exclude values far z cross diff at max peak
    ZdiffAtBigPm=ZdiffAtBigP(end);
    ZdiffAtBig=ZdiffAtBigP(ZdiffAtBigP>0.7*ZdiffAtBigPm & ZdiffAtBigP<1.3*ZdiffAtBigPm);
    fLowEe=1/(max(ZdiffAtBig)*2*xTOindx);
    fHighEe=1/(min(ZdiffAtBig)*2*xTOindx)*1.3;
    [diffPeakS]=sort(diff(PeakPosInterpol));

   diffZS=sort(diff(ZPosInterpol));
    if NoiseLevel>0
        periodPindx=diffPeakS(round(end*0.9))*2;
        periodZindx=diffZS(round(end*.9))*2;
        if numel(diffPeakS)>4 && diffPeakS(ceil(0.2*end))*1.5>diffPeakS(round(0.85*end))
          %majority of similar peak diffs, eliminate AfLowE
          fLowEe=0;
          periodPindx=diffPeakS(floor(0.5*end)+1)*2;
          periodZindx=periodPindx;
        end        
    else
        periodPindx=diffPeakS(floor(0.5*end)+1)*2;
        periodZindx=diffZS(floor(0.5*end)+1)*2;
        if numel(diffPeakS)>4 && diffPeakS(2)*1.5>diffPeakS(end-1)
          %low variation of peakpos, eliminate AfLowE
          fLowEe=0;
          periodPindx=diffPeakS(floor(0.5*end)+1)*2;
          periodZindx=periodPindx;          
        end
    end
    diffZSm2=diffZS(floor(0.5*end)+1)*2;%median
    pBig=max([periodZindx,periodPindx,diffZSm2]);
    pSmall=min([periodZindx,periodPindx,diffZSm2]);
    fLowE=1/(pBig*xTOindx)*fToleranceLow;
    fHighE=1/(pSmall*xTOindx)*fToleranceHigh;
    fE=(1/(pBig*xTOindx)+1/(pSmall*xTOindx))/2;
    if fE*(x(end)-x(1))>3 %num periods
      PeakInterpolS=sort(abs(PeakInterpol));  
      amplitudeE=mean(PeakInterpolS(ceil(0.3*end):end)); 
    else
      diffPeakS=sort(abs(diff(PeakInterpol)));  
      amplitudeE=mean(diffPeakS(floor(0.5*end)+1:end))/2;
    end
end

%Estimate phase Phi0
phaseshiftE=0;
    %Only good if x=0 in observed range and low noise, 
    %does not improve results with nlinfit!
    %Sometimes good as a rough estimate. 
if x(1)<=0 && x(end)>=0
    %Translate ZInterpol to x-values
    xZInterpol=xTOindx*ZPosInterpol+x(1)-xTOindx;
    [~,Zindx]=min(abs(xZInterpol));
    phaseshiftE=-2*pi*fE*xZInterpol(Zindx);
    if ysmooth(floor(ZPosInterpol(Zindx))+1)<0
         %zero cross direction: \
        phaseshiftE=phaseshiftE+pi;
    end
    %  make positive Phi0 values
    if phaseshiftE<0
       phaseshiftE=phaseshiftE+2*pi;
    end 
end    

%Regression
paramsE=[offsetE,amplitudeE,0,phaseshiftE];%paramsE[3): frequency fToTest
    %since regression is very sensitive on frequency, 
    %try several frequencies around the two estimated frequency ranges.
maxf=1/(x(3)-x(1));%max frequency, min 2 samples per period.
if fHighE>maxf;fHighE=maxf;end;
step=0.015;
if fLowEe==0
    %Only one range since extra f not estimated
    fToTest=fLowE:step*fHighE:fHighE;
else
    if fHighEe>maxf;fHighEe=maxf;end;
    fmin=sort([fLowE,fLowEe]);
    fmax=sort([fHighE,fHighEe]);
    if fmin(2)<fmax(1);fmin(2)=fmax(1);end;%avoid overlaping of frequency
    fToTest=[fmin(1):step*fmax(1):fmax(1),fmin(2):step*fmax(2):fmax(2)];
end
%Determine sine parameters
SineParams=RUNnlinfit(paramsE,fToTest,x,y);

%uncomment following statement if you want to use lsqcurvefit and 
%comment statement above.
% SineParams=RUNlsqcurvefit(paramsE,fToTest,x,y);

%======The section below may be deleted
%Display possibly better result if amplitude is unexpected high
NP=SineParams(3)*(x(end)-x(1));%num periods
SpP=1/(SineParams(3)*xTOindx);%Samples per period
if SpP<=3
  amplitudeE=max(abs(ysmooth));
end
if (SineParams(2)>3*amplitudeE && NP>1.2 && SpP>3) ...
        || SineParams(2)>15*amplitudeE...
        || (SineParams(2)>3*amplitudeE && SpP<=3)
    if x(1)<=0 && x(end)>=0
      fprintf(['The ouput amplitude is much higher than estimated. The initial sine analysis might be better:\n'...
        'y=  %0.2g + %0.2g * sin(2*pi* %0.2g *x + %0.2g )'... 
        '   alternative frequency: %0.2g \n'],offsetE,amplitudeE,fE,phaseshiftE,SineParams(3))
    else
      fprintf(['The ouput amplitude is much higher than estimated. The initial sine analysis might be better:\n'...
        'y=  %0.2g + %0.2g * sin(2*pi* %0.2g *x + ?)'... 
        '   alternative frequency: %0.2g \n'],offsetE,amplitudeE,fE,SineParams(3))
    end
    %following statement improves a few results
%     SineParams=[roundsd(offsetE,2),roundsd(amplitudeE,2),SineParams(3),roundsd(SineParams(4),2)];
end
%Following statement may be deleted. It displays some analysis
% PlotAnalysis(x,y,ysmooth,SineParams,ZPosInterpol,...
%     PeakPosInterpol,PeakInterpol,offsetE,amplitudeE,fE,phaseshiftE,...
%     filterValue,SineParams)
%======Above section may be deleted

function [PeakInterpol,PeakPosInterpol]=PeakInterpolation(ysmooth,PeakPos)
    if numel(PeakPos)>0
        %Parabola interpolation of peak value and peak position
        y1=ysmooth(PeakPos-1);
        y2=ysmooth(PeakPos);
        y3=ysmooth(PeakPos+1);
        indexTemp=(y1-y3)./(2*y1-4*y2+2*y3);
        PeakPosInterpol=PeakPos+indexTemp;
        PeakInterpol=(y3+y1-2*y2)/2.*indexTemp.^2+(y3-y1)/2.*indexTemp+y2;
    else
        PeakInterpol=[];
        PeakPosInterpol=[];
    end

function PeakPos=EliminateMultiplePeaks(ZPos,PeakPos,ysmooth)  
%eliminate peaks if more than one peak between 2 zero crosses
if numel(PeakPos)>1
    for i=1:numel(ZPos)-1
        pos1=ZPos(i)+0.5;
        pos2=ZPos(i+1)+0.5;
        indxPeaks=find(PeakPos>pos1 & PeakPos<pos2);
        if numel(indxPeaks)>1
            %keep the largest peak
           [~,indxindxMaxPeak]=max(abs(ysmooth(PeakPos(indxPeaks))));
           indxMaxPeak=indxPeaks(indxindxMaxPeak);
           toDelete=indxPeaks(indxPeaks~=indxMaxPeak);
           PeakPos(toDelete)=[];
        end
    end
    %allow only one peak on extreme side of 1st and last zero cross
    indxPeaks=find(PeakPos<ZPos(1)+0.5);
    if numel(indxPeaks)>1
        %keep largest peak
       [~,indxindxMaxPeak]=max(abs(ysmooth(PeakPos(indxPeaks))));
       indxMaxPeak=indxPeaks(indxindxMaxPeak);
       toDelete=indxPeaks(indxPeaks~=indxMaxPeak);
       PeakPos(toDelete)=[]; 
    end
    indxPeaks=find(PeakPos>ZPos(end)+0.5);
    if numel(indxPeaks)>1
        %keep largest peak
       [~,indxindxMaxPeak]=max(abs(ysmooth(PeakPos(indxPeaks))));
       indxMaxPeak=indxPeaks(indxindxMaxPeak);
       toDelete=indxPeaks(indxPeaks~=indxMaxPeak);
       PeakPos(toDelete)=[]; 
    end
end


function ys=smooth1D(y,n)
%Moving smoothing filter for 1D data with floating point filter values .
%(The build-in function smooth.m from matlab operates only with odd integer
%nummbers.)
%Input:
%   y is the vector to be smoothed. 
%   n is the filter value with any number.
%     Filter values n below/equal 1 will return the original input.
%     Filter values larger than the number of elements in the y vector are reduced.
%     The filter is symmetric and has therefore an odd number of elements.
%     Filter examples:
%        n=1.5:      0.25, 1, 0.25
%        n=2.5:      0.75, 1, 0.75
%        n=4.0:    0.5, 1, 1, 1, 0.5
%Author: Peter Seibold

    if n<=1;ys=y;return;end;%no filter
    numely=numel(y);
    if n>numely-2 %exact formula: n>floor((numel(y)-1)/2)*2+1, slower
        n=floor((numely-1)/2)*2+1;%closest odd number
        warning('number of filter elements too large.\n%s',...
            ['Filter value reduced to ' num2str(n)]);
    end
    spanL=ceil((n-1)/2);%num filter coeff left or right of filter center
    borderValue=(n-spanL*2+1)/2;%e.g. borderValue=0.25 for n=1.5 
    ys=y;%Unchanged: ys(1)=y(1) and ys(end)=y(end)
    %Filter left and right border of y
    div=2;
     for i=2:spanL
        ys(i)=sum(y(1:div))/div;
        ys(end-i+1)=sum(y(end-div+1:end))/div;
        div=div+2;
     end
     %Filter remaining center part of y
     for i=spanL+1:numely-spanL
         ys(i)=y(i-spanL)*borderValue+sum(y(i-spanL+1:i+spanL-1))+y(i+spanL)*borderValue;
     end
     ys(spanL+1:numely-spanL)=ys(spanL+1:numely-spanL)/n;
     
     
%roundsd may be replaced for sineFit by round(xx,n,'significant') with higher matlab versions   
function y=roundsd(x,n,method)
%ROUNDSD Round with fixed significant digits
%	ROUNDSD(X,N) rounds the elements of X towards the nearest number with
%	N significant digits.
%
%	ROUNDSD(X,N,METHOD) uses following methods for rounding:
%		'round' - nearest (default)
%		'floor' - towards minus infinity
%		'ceil'  - towards infinity
%		'fix'   - towards zero
%
%	Examples:
%		roundsd(0.012345,3) returns 0.0123
%		roundsd(12345,2) returns 12000
%		roundsd(12.345,4,'ceil') returns 12.35
%
%	See also Matlab's functions ROUND, ROUND10, FLOOR, CEIL, FIX, and 
%	ROUNDN (Mapping Toolbox).
%
%	Author: Franois Beauducel <beauducel@ipgp.fr>
%	  Institut de Physique du Globe de Paris
%
%	Acknowledgments: Edward Zechmann, Daniel Armyr, Yuri Kotliarov
%
%	Created: 2009-01-16
%	Updated: 2015-04-03

%	Copyright (c) 2015, Franois Beauducel, covered by BSD License.
%	All rights reserved.
%
%	Redistribution and use in source and binary forms, with or without 
%	modification, are permitted provided that the following conditions are 
%	met:
%
%	   * Redistributions of source code must retain the above copyright 
%	     notice, this list of conditions and the following disclaimer.
%	   * Redistributions in binary form must reproduce the above copyright 
%	     notice, this list of conditions and the following disclaimer in 
%	     the documentation and/or other materials provided with the distribution
%	                           
%	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
%	AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
%	IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
%	ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
%	LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
%	CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
%	SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
%	INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
%	CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
%	ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
%	POSSIBILITY OF SUCH DAMAGE.

if nargin < 2
	error('Not enough input arguments.')
end

if nargin > 3
	error('Too many input arguments.')
end

if ~isnumeric(x)
		error('X argument must be numeric.')
end

if ~isnumeric(n) || ~isscalar(n) || n < 0 || mod(n,1) ~= 0
	error('N argument must be a scalar positive integer.')
end

opt = {'round','floor','ceil','fix'};

if nargin < 3
	method = opt{1};
else
	if ~ischar(method) || ~ismember(opt,method)
		error('METHOD argument is invalid.')
	end
end
e = floor(log10(abs(x)) - n + 1);
og = 10.^abs(e);
y = feval(method,x./og).*og;
k = find(e<0);
if ~isempty(k)
	y(k) = feval(method,x(k).*og(k))./og(k);
end	
y(x==0) = 0;    


function PlotAnalysis(x,y,ysmooth,SineParams,ZPosInterpol,...
    PeakPosInterpol,PeakInterpol,offsetE,amplitudeE,fE,phaseshiftE,...
    filterValue,BetaBest)
global figSineFit %used in sineFitDemo
numelx=numel(x);
xindx=1:numelx;
figSineFit=figure;
hold on;
p1=plot(xindx,y,'color',[0.5,0.5,0.5]);grid on;grid minor;

xstart=x(1);
xend=x(end);
x3=(xend-xstart)/(numelx-1)*(xindx-1)+xstart;
x4b=1:0.01:numel(x);
x4=(xend-xstart)/(numelx-1)*(x4b-1)+xstart;    
y3=SineParams(1)+SineParams(2)*sin(2*pi*SineParams(3)*x3+SineParams(4));
y4=SineParams(1)+SineParams(2)*sin(2*pi*SineParams(3)*x4+SineParams(4));
p3=plot(xindx,y3,'k.',x4b,y4,'r-','LineWidth',1);%result with dots
plot(ZPosInterpol,offsetE,'b+','LineWidth',0.8);%plot Z crosses
plot(PeakPosInterpol,PeakInterpol+offsetE,'bx', 'MarkerSize',7,'LineWidth',1.0);%plot peaks
ymax=max(y);
ymin=min(y);
diffMinMaxy=ymax-ymin;
diffMinMaxy01=diffMinMaxy*0.1;
%low y plot limit
if offsetE-SineParams(2)<ymin-diffMinMaxy
   newYlim(1)=ymin-diffMinMaxy01;
else
   newYlim(1)=min([ymin,offsetE-SineParams(2)])-diffMinMaxy01; 
end
%high y plot limit
if offsetE+SineParams(2)>ymax+diffMinMaxy
   newYlim(2)=ymax+diffMinMaxy01;
else
   newYlim(2)=max([ymax,offsetE+SineParams(2)])+diffMinMaxy01; 
end
set(gca,'ylim',newYlim);
if filterValue>1
    p2=plot(1:numelx,ysmooth+offsetE,'c');%ysmooth
    legend([p1,p2,p3(2)],'y in','y smoothed','result','Location','best');
else
    legend([p1,p3(2)],'y in','result','Location','best');
end

%second x-axis with x values:
set(gca,'xlim',[1,numelx]);
ax2_xt=1:ceil(numelx/5):numelx;
ax2_xL =cellfun(@(xx) num2str(xx,3),num2cell(x(ax2_xt)),'UniformOutput', false);
axes('XAxisLocation','top','Color','none','xlim',[1,numel(x)],...
    'XTick',ax2_xt,'XTickLabel',ax2_xL,'ylim',newYlim,'YTick',[],'YTickLabel',[]);
hold off;
opts = statset('nlinfit');opts.MaxIter=2000;
modelfun = @(paramc,x) paramc(1) + paramc(2) * sin(2*pi* paramc(3)*x+paramc(4) );    
warning('off','all');    
    [~,~,~,~,MSE] = nlinfit(x,y,modelfun,BetaBest,opts);% only for MSE
warning('on','all');  
disp(['1st estimate (NO result): y= ' num2str(offsetE,2) ' + ' num2str(amplitudeE,2) ... 
    ' * sin(2*pi*' num2str(fE,2) '+' num2str(phaseshiftE,2) ')']);  
disp(['Calculation: y= ' num2str(SineParams(1)) ' + ' num2str(SineParams(2)) ... 
' * sin(2*pi*' num2str(SineParams(3)) '+' num2str(SineParams(4)) ')  MSE: ' num2str(MSE)]);

function SineParams=RUNnlinfit(paramsE,fToTest,x,y)
%this function is obsolete if you use RUNlsqcurvefit
Beta=zeros(numel(fToTest),5);
opts = statset('nlinfit');opts.MaxIter=100;
% MaxIter=50 decreases number of good results by 0.05% compared to opts=100
% MaxIter=20 decreases number of good results by 0.1% compared to opts=100
% MaxIter=500 increases number of good results by 0.01% compared to opts=100
modelfun = @(paramc,x) paramc(1) + paramc(2) * sin(2*pi* paramc(3)*x+paramc(4));
%test different frequencies
warning('off','all');
for j=1:numel(fToTest)
    paramsE(3)=fToTest(j);  
    %statistics_toolbox necessary
    [BETA,~,~,~,MSE] = nlinfit(x,y,modelfun,paramsE,opts);
    Beta(j,:)=[BETA,MSE];
end
warning('on','all');

weight=[4;1.8;1.7;1.6;1.5;1.4;1.3;1.2;1.1;1.0;0.9;0.8];
BetaS=sortrows(Beta,5);
if BetaS(1,5)==inf
    %in very seldom cases no result found by nlinfit with 100 iterations
    %it  may happen with less than 5 samples in total
    BetaS(1,5)=9999;%set dummy value
end
MSEFactor=1.1;
SrIndM=min(size(BetaS,1),12);%Sorted Beta rounded max index
%check for negative frequencies (often bad estimate) or very big amplitudes.
if BetaS(1,3)<0 && BetaS(2,5)< 1.05*BetaS(1,5) || abs(BetaS(1,2))>15*paramsE(2)
    weight(1)=1.8;
    if abs(BetaS(1,2))>15*paramsE(2)
       weight(1)=1.7;
       MSEFactor=3;
       if abs(BetaS(2,2))>15*paramsE(2) 
        weight(2)=1.7; 
       end
    end
    i=1;
    while i<=size(BetaS,1) && abs(BetaS(i,2))>15*paramsE(2)
        weight(i)=1;
        i=i+1;
    end
    SrIndM=min(size(BetaS,1),i+15);    
    weight=ones(SrIndM,1);    
end
BetaSr=[roundsd(BetaS(1:SrIndM,3),2),roundsd(BetaS(1:SrIndM,5),2)...
    ,weight(1:SrIndM)];
lastValidIndx= find(BetaSr(:,2)<=BetaSr(1,2)*MSEFactor,1,'last');%last good MSE
BetaSr=abs(BetaSr(1:lastValidIndx,:));%eliminate frequencies with bad MSE, make all f positive
fMost = mode(BetaSr(:,1));%find frq with most occurrence
BetaSrSum1=sum(BetaSr(BetaSr(:,1)==BetaSr(1,1),3));
BetaSrSum2=sum(BetaSr(BetaSr(:,1)==fMost,3));
if BetaSrSum1>=BetaSrSum2
    BetaSrBestIndx=1;
else
    BetaSrBestIndx=find(BetaSr(:,1)==fMost,1,'first');
end
BetaBest=BetaS(BetaSrBestIndx,1:4);   
opts = statset('nlinfit');opts.MaxIter=2000;
% MaxIter=2000 improves amount of good results by 0.01% compared to MaxIter=1000
% MaxIter=500 as MaxIter 1000
%no further nlinfit: worse by 0.01% compared to MaxIter=1000
warning('off','all');    
    BETA = nlinfit(x,y,modelfun,BetaBest,opts);
warning('on','all');

%make frequency positive
if BETA(3)<0
   BETA(3)=-BETA(3);
   BETA(4)=pi-BETA(4);%sin(2*pi*-f-phi)=sin(2*pi*f+phi+pi)
end
%make amplitude positive
if BETA(2)<0 
    BETA(2)=-BETA(2);
    BETA(4)=BETA(4)+pi; 
end
%make phase smaller and positive
BETA(4)=rem(BETA(4),2*pi);
if BETA(4)<0
    BETA(4)=BETA(4)+2*pi;
end
SineParams=BETA;

function SineParams=RUNlsqcurvefit(paramsE,fToTest,x,y)
%this function is obsolete if you use nlinfit.m
%Results are not as good as with nlinfit.
Beta=zeros(numel(fToTest),5);
opts = optimoptions('lsqcurvefit','Display','off','MaxIterations',100);
modelfun = @(paramc,x) paramc(1) + paramc(2) * sin(2*pi* paramc(3)*x+paramc(4));
lb=[paramsE(1)/2,paramsE(2)/2,paramsE(3)/2,0];
ub=[paramsE(1)*2,paramsE(2)*4,paramsE(3)*2,2*pi];
warning('off','all');
for j=1:numel(fToTest)
    %try differnt frequencies
    lb(3)=fToTest(j)/2;
    ub(3)=fToTest(j)*2;
    paramsE(3)=fToTest(j);  
    %Optimization Toolbox necessary
    BETA = lsqcurvefit(modelfun,paramsE,x,y,lb,ub,opts);
    resnorm=sum((modelfun(BETA,x)-y).^2);
    Beta(j,:)=[BETA,resnorm];
end
warning('on','all');
weight=[4;1.8;1.7;1.6;1.5;1.4;1.3;1.2;1.1;1.0;0.9;0.8];
BetaS=sortrows(Beta,5);
amplitudeE=paramsE(2);
if BetaS(1,3)<0 && BetaS(2,5)< 1.05*BetaS(1,5) || abs(BetaS(1,2))>15*amplitudeE
    BetaS(1,5)=BetaS(2,5);%set like second best MSE
    weight(1)=1.8;
    if abs(BetaS(1,2))>15*amplitudeE 
       weight(1)=1.7; 
       if abs(BetaS(2,2))>15*amplitudeE 
        weight(2)=1.7; 
        BetaS(1,5)=BetaS(3,5);%set like 3rd best MSE
        BetaS(2,5)=BetaS(3,5);%set like 3rd best MSE
       end
    end
end
BetaSr=[roundsd(BetaS(1:min(12,end),3),2),roundsd(BetaS(1:min(12,end),5),2)...
    ,weight(1:min(12,size(BetaS,1)))];
lastValidIndx= find(BetaSr(:,2)<BetaSr(1,2)*1.5,1,'last');
BetaSr=abs(BetaSr(1:lastValidIndx,:));%eliminate frequencies with bad MSE, make all f positive
fMost = mode(BetaSr(:,1));%find frq with most occurrence
BetaSrSum1=sum(BetaSr(BetaSr(:,1)==BetaSr(1,1),3));
BetaSrSum2=sum(BetaSr(BetaSr(:,1)==fMost,3));
if BetaSrSum1>=BetaSrSum2
    BetaSrBestIndx=1;
else
    BetaSrBestIndx=find(BetaSr(:,1)==fMost,1,'first');
end
BetaBest=BetaS(BetaSrBestIndx,1:4);
opts = optimoptions('lsqcurvefit','Display','off','MaxIterations',2000);
warning('off','all');    
    lb(3)=BetaBest(3)/2;
    ub(3)=BetaBest(3)*2;
    BETA = lsqcurvefit(modelfun,BetaBest,x,y,lb,ub,opts);
warning('on','all');
%make frequency positive
if BETA(3)<0
   BETA(3)=-BETA(3);
   BETA(4)=pi-BETA(4);%sin(2*pi*-f-phi)=sin(2*pi*f+phi+pi)
end
%make amplitude positive
if BETA(2)<0 
    BETA(2)=-BETA(2);
    BETA(4)=BETA(4)+pi; 
end
%make phase smaller and positive
BETA(4)=rem(BETA(4),2*pi);
if BETA(4)<0
    BETA(4)=BETA(4)+2*pi;
end
SineParams=BETA;