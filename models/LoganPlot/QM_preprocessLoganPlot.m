function [results, plots] = QM_preprocessLoganPlot(Ct,Cref,times,gError,k2,tR,tC)
% QM_preprocessLoganPlot function
% function [results plots] = QM_preprocessLoganPlot(Ct,Cref,times,gError,k2,tR,tC)
%
% QM_preprocessLoganPlot is a function which fit the Logan Plot model to the 
% current TACs and return the fitted TAC and the result parameters.
% The input data come from the previous calculated TACs and the 
% parameters introduced by the user in the GUI. 
%
% Parameters:
% Ct -> Vector containing the TAC of the current ROI
% Cref -> Vector containing the TAC of the Reference region
% times -> timetable of current study
% gError -> Maximun error permited in the fitting 
% k2 -> a difussion constant
% tR -> restriction value times 
% tC -> boolean variable containing the values of the check box times
%---------------------------------------------------------------------------------
% QM_preprocessLoganPlot is part of QModeling.
%
% QModeling is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% QModeling is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with QModeling.  If not, see <http://www.gnu.org/licenses/>.
%---------------------------------------------------------------------------------
% Copyright (C) 2015, 2017 Francisco Javier Lopez Gonzalez, Jose Paredes Pacheco, 
% Karl-Khader Thurnhofer Hemsi, Daniel Toro Flores

global QMmaindata
global QMlastPreprocess


num_files=QMmaindata.num_files;

% Variables to the restrictions of t*
tRest=tR(1)*60;
tLwr=tR(2)*60;
tUpr=tR(3)*60;
tRC=tC(1);
tCh=tC(2);

%Ecuation

k2= k2/60;

pY=cumsum(diff(times,1,2).*Ct)./Ct;
numX=(cumsum(diff(times,1,2).*Cref) + Cref/k2);
pX=numX./Ct;
f=0;

% Correct cancel errors (0/0)
for i = 1:num_files
    if(isnan(pY(i)) || isinf(pY(i)))
        pY(i)=0;
    end
    if(isnan(pX(i)) || isinf(pX(i)))
        pX(i)=0;
    end
    if(pX(i)==0 && pY(i)==0)
        f=f+1;
    end
end

%VECTORs X and Y without zeros

pX=pX(f+1:end);
pY=pY(f+1:end);
     
%Error calc.
[Eq,S]=polyfit(pX,pY,1);
Z=polyval(Eq,pX);% Error of all the values of regression
uError = max(abs(pY-Z)./abs(Z))*100;
if(uError>=100)
    uError=100;
end


%Store times.
uTimes=times(f+1:end,1);
%Last time value.
uTimes((num_files+1)-f)=times(num_files,2);
%Middle times.
mTimes=(times(f+1:num_files,1)+times(f+1:num_files,2))/2;
nError=false;


%Restriction time is inactive

if(tCh==1 && tRC==0)
    i=1;
    if(uError>gError)
        while (uError >= gError && i<=length(pX)-2)
           
            [uEq,S]=polyfit(pX(i:end),pY(i:end),1);
            uZ=polyval(uEq,pX(i:end));
            uError = max(abs(pY(i:end)-uZ)./abs(uZ))*100;
            i=i+1;
            if(uError<gError)
                nError=true;
            end
        end
        i=i-1;
        Z=polyval(uEq,pX);
        BPnd=uEq(1)-1;
        intercept=uEq(2);
    end
    
    if(uError==gError)
        nError=true;
        BPnd=Eq(1)-1;
        intercept=Eq(2);

    end
end


if(tRC==1)
    t=0;
    i=1;
    u=1;
    
    if tLwr > mTimes(length(pX)-1)
    ME = MException('PreprocessLoganPlot:LowerTimeIncorrect','Lower time restriction only take the last frame');
    throw(ME)
    elseif(tLwr<=uTimes(1))
        i=1;
    else
        while(t == 0)
            i=i+1; 
            if(tLwr <= uTimes(i))
            t=1;
            end
        end
        i=i-1;
    end
    if(tLwr > mTimes(i))
        i=i+1;
    end
    t=0;
    if tUpr > mTimes(length(pX))
    u=length(pX);
    else
        while(t == 0)
            u=u+1; 
            if(tUpr <= uTimes(u))
            t=1;
            end
        end
        u=u-1;
    end
    if(tUpr >= mTimes(u))
        u=u+1;
    end
    if(uError>gError)    
        while (uError >= gError&& i<=u+1)
            [uEq,S]=polyfit(pX(i:end),pY(i:end),1);
            uZ=polyval(uEq,pX(i:end));
            uError = max(abs(pY(i:end)-uZ)./abs(uZ))*100;
            i=i+1;
            if(i>u+1)
                ME = MException('PreprocessLoganPlot:MaxErrCondition','Not full fill MaxErr condition');
                throw(ME) 
            end
            if(uError<gError)
                nError=true;
            end 
        end
         i=i-1;
         Z=polyval(uEq,pX);
         BPnd=uEq(1)-1;
         intercept=uEq(2);       
    end
end

%Restriction time is active
if(tCh==0)
    i=1;
    t=0;
    gError=100;
     if tRest > mTimes(length(pX)-1)
    ME = MException('PreprocessLoganPlot:NumberDataPointsFitCondition','The t* value only takes the last frame');
    throw(ME);
    end 
    if tRest >= uTimes(length(pX)-1)
        tRest=uTimes(length(pX)-1);
    end
    
    if(tRest>0)
     while(t == 0)
         if(tRest <= uTimes(i))
             t=1;
         end
     i=i+1;   
     end 
    i=i-2; 
    end
    
    if(i<length(pX)-1)
        if(tRest > mTimes(i))
            i=i+1;
        end
    end
    
    if(i>length(pX)-2)
        warndlg('The restrict time only takes the last two time points');
    end
 
    if(uError>gError)    
        while (uError >= gError && i<=length(pX))
            [uEq,S]=polyfit(pX(i:end),pY(i:end),1);
            uZ=polyval(uEq,pX(i:end));
            uError = max(abs(pY(i:end)-uZ)./abs(uZ))*100;
            i=i+1;
            if(uError<gError)
                nError=true;
            end
        end
         i=i-1;
         Z=polyval(uEq,pX);
         BPnd=uEq(1)-1;
         intercept=uEq(2);
        
    end
    if(uError==gError)
        
        [uEq,S]=polyfit(pX(i:end),pY(i:end),1);
        
        BPnd=uEq(1)-1;
        intercept=uEq(2);
        Z=polyval(uEq,pX);
        nError=true;
    end
end
 
if (nError == false)   
    ME = MException('PreprocessLoganPlot:MaxErrCondition','Not full fill MaxErr condition');
    throw(ME);
else   
    
    %Saving results
    pX=pX./60;
    %NMSE=goodnessOfFit(pY(i:end),Z(i:end),'NMSE');
    NMSE=1-((pY(i:end)-Z(i:end)).^2)./((pY(i:end)-mean(pY(i:end))).^2);
    CCoef=corrcoef(pY(i:end),Z(i:end));
%     results=[BPnd intercept uTimes(i)/60 pX(i) CCoef(1,2) NMSE];

    %Covariance matrix (see polyfit manual)
    CovM=(inv(S.R)*inv(S.R)')*(S.normr^2)/S.df;

    results={BPnd intercept uTimes(i)/60 CCoef(1,2) mean(NMSE) diag(CovM)};
    plots=[pX pY Z];
    QMlastPreprocess.LoganPlot.numX=numX;
    QMlastPreprocess.LoganPlot.l=i+f;
    
end

