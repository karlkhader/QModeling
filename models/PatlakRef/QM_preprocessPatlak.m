function [results, plots] = QM_preprocessPatlak(Ct,Cr,times,t,restrict,Max_err,waitbarhandle)
% QM_preprocessPatlak function
% function [results plots] = QM_preprocessPatlak(Ct,Cr,times,t,restrict,Max_err,waitbarhandle)
%
% QM_preprocessPatlak is a function which fit the Patlak model to the 
% current TACs and return the fitted TAC and the result parameters.
% The input data come from the previous calculated TACs and the 
% parameters introduced by the user in the GUI. 
%
% Parameters:
% Ct -> Vector containing the TAC of the current ROI
% Cr -> Vector containing the TAC of the Reference region
% times -> timetable of current study
% t -> Vector of two values containing the t* parameter
% restrict -> Vector with values to restrict the t*
% Max_err -> Maximun error permited in the fitting
% waitbarhandle -> handle of the current waitbar
%---------------------------------------------------------------------------------
% QM_preprocessPatlak is part of QModeling.
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
% Copyright (C) 2015, 2016 Francisco Javier Lopez Gonzalez, Jose Paredes Pacheco, 
% Karl-Khader Thurnhofer Hemsi


num_files = length(Ct);
tmin = ((times(1:num_files,1)+times(1:num_files,2))/2)/60;

% Copute the integral:
intvalue = cumsum(diff(times,1,2).*Cr);

%Obtengo los valores normalizados
X = intvalue./Cr;
Y = Ct./Cr;

% Correct cancel errors (0/0)
X(isnan(X)|isinf(X)) = 0;
Y(isnan(Y)|isinf(Y)) = 0;
% First frame to be used to fitting
indfr_min = find(X,1,'first');

found = false;
enoughpoints = true;

tchecked = t(1);
if tchecked == 1
    
    if (restrict(1) == 1)

        % Variables to the restrictions of t*
        rest_min = restrict(2);
        rest_max = restrict(3);
        
        rtm = find((times(:,1)/60) >= rest_min,1,'first');
        if rtm ~= 1 && rest_min <= tmin(rtm-1)
            rtm = rtm-1;
        end
        rtM = find((times(:,1)/60) < rest_max,1,'last');    %I must take a value if there is between the middle point of PET times
        if rtM ~= num_files && rest_max >= tmin(rtM)
            rtM = rtM+1;
        end

        if (rtm > num_files-1) || (rtm > rtM) % General case

%             i = 1;
            i = indfr_min;
            while (i <= num_files-2) && (found == false)      % At least 3 point to ajust
                x1 = X(i:end);
                y1 = Y(i:end);
                [p,S] = polyfit(x1,y1,1);

                z = polyval(p,x1);          % Error of all the values of regression
                err_act = max(abs(y1-z)./abs(z))*100;

                waitbar(i/num_files,waitbarhandle)

                if err_act < Max_err
                    Max_err = err_act;
                    new_p = p;
                    time_ind_st = i;
                    found = true;
                    waitbar(1,waitbarhandle)
                end
                i = i+1;
            end

        else

            if (rtm == num_files-1)

                x1 = X(rtm:end);
                y1 = Y(rtm:end);
                [new_p,S] = polyfit(x1,y1,1);
                new_z = polyval(new_p,x1);          % Error of all the values of regression
                Max_err = max(abs(y1-new_z)./abs(new_z))*100;
                time_ind_st = rtm;

            else

                if (rtM >= num_files-1)
                    rtM = num_files-2;
                end

%                 i = rtm;
                i = max(rtm,indfr_min);
                while (i <= rtM) && (found == false)      % At least 3 point to ajust
                    x1 = X(i:end);
                    y1 = Y(i:end);
                    [p,S] = polyfit(x1,y1,1);

                    z = polyval(p,x1);          % Error of all the values of regression
                    err_act = max(abs(y1-z)./abs(z))*100;

                    waitbar(i/num_files,waitbarhandle)

                    if err_act < Max_err
                        Max_err = err_act;
                        new_p = p;
                        time_ind_st = i;
                        found = true;
                        waitbar(1,waitbarhandle)
                    end
                    i = i+1;
                end
            end
        end           

    else

%         i = 1;
        i = indfr_min;
        while (i <= num_files-2) && (found == false)      % At least 3 point to ajust

            x1 = X(i:end);
            y1 = Y(i:end);
            [p,S] = polyfit(x1,y1,1);

            z = polyval(p,x1);                      % Error of all the values of regression
            err_act = max(abs(y1-z)./abs(z))*100;

            waitbar(i/num_files,waitbarhandle)

            if err_act < Max_err
                Max_err = err_act;
                new_p = p;
                time_ind_st = i;
                found = true;
                waitbar(1,waitbarhandle)
            end
            i = i+1;
            
        end

    end
    
else
    
    tmin_ind = find((times(:,1)/60) >= t(2),1,'first');
    if isempty(tmin_ind)
        enoughpoints = false;
    else
        if tmin_ind ~= 1 && t(2) <= tmin(tmin_ind-1)
            tmin_ind = tmin_ind-1;
        end

        %Nuevo
        tmin_ind = max(tmin_ind,indfr_min);

        if (tmin_ind == num_files)
            enoughpoints = false;
    %         ME = MException('PreprocessPatlak:NumberDataPointsFitCondition','The t* value only takes the last frame');
    %         throw(ME);
        else

            x1 = X(tmin_ind:end);
            y1 = Y(tmin_ind:end);
            [new_p,S] = polyfit(x1,y1,1);

            z = polyval(new_p,x1);                      % Error of all the values of regression
            Max_err = max(abs(y1-z)./abs(z))*100;

            time_ind_st = tmin_ind; %O indfr_min; ESTO VEO MEJOR DEJARLO COMO ESTÁ
            found = true;

            waitbar(1,waitbarhandle)
        end
    end
    
end


if (found == false && enoughpoints == true)   
    ME = MException('PreprocessPatlak:MaxErrCondition','Not full fill MaxErr condition');
    throw(ME);
elseif (enoughpoints == false)
    ME = MException('PreprocessPatlak:NumberDataPointsFitCondition','The t* value only takes the last frame');
    throw(ME);
else
  
    %Saving results
    T = times(time_ind_st)/60;
    K = new_p(1)*60;
    intercept = new_p(2);
    Start = X(time_ind_st)/60;

    %Saving plots
    Z = polyval(new_p,X);
    plots = [X Y Z];
    
    % Goodness of fit
    % MSE = mean((Y(time_ind_st:end)-Z(time_ind_st:end)).^2);
    NMSE=1-((Y(time_ind_st:end)-Z(time_ind_st:end)).^2)./((Y(time_ind_st:end)-mean(Y(time_ind_st:end))).^2);
    CCoef=corrcoef(Y(time_ind_st:end),Z(time_ind_st:end));
    
    %Covariance matrix (see polyfit manual)
    CovM=(inv(S.R)*inv(S.R)')*(S.normr^2)/S.df;
    
    % results = [T K intercept Start Max_err MSE CCoef(1,2)];
    results = {T K intercept Start Max_err mean(NMSE) CCoef(1,2) diag(CovM)};

end

