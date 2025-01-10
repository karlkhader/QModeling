function [results plots B theta_3] = QM_preprocessCarson(Ct,Cr,times,k2_p,resampling,basisF,k2a_max,k2a_min,waitbarhandle)
% QM_preprocessCarson function
% function [results plots B theta_3] = QM_preprocessCarson(Ct,Cr,times,k2_p,resampling,basisF,k2a_max,k2a_min,waitbarhandle)
%
% QM_preprocessCarson is a function which fit the SRTM2 model to the 
% current TACs and return the fitted TAC, the result parameters,
% the matrix with the values of basis functions and the theta 3
% parameter for calculating the parametric images.
% The input data come from the previous calculated TACs and the 
% parameters introduced by the user in the GUI. 
%
% Parameters:
% Ct -> Vector containing the TAC of the current ROI
% Cr -> Vector containing the TAC of the Reference region
% times -> timetable of current study
% k2_p -> k2' parameter set by the user
% resampling -> resampling to define the grid for convolution
% basisF -> number of basis functions
% k2amax -> k2a maximum defined by the user
% k2amin -> k2a minimum defined by the user
% waitbarhandle -> handle of the current waitbar
%---------------------------------------------------------------------------------
% QM_preprocessCarson is part of QModeling.
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

format long;

RSS_menor = Inf; 
i_min = 1;
                       
num_files=length(Ct);
w = ones(1,num_files); %vector of weights initializing by one
tseg = zeros(1,num_files);

for i=1:num_files 
    tseg(i)=(times(i,1)+times(i,2))/2;
end

W=eye(num_files); %weight diagonal matrix. Is the identity matrix (in fact, the weigths are not consider but we follow Gunn algorithm

tgrid = 0:resampling:times(num_files,2);
Cr_grid = interp1([0 tseg],[0;Cr],tgrid,'linear');
Cr_grid(isnan(Cr_grid)) = Cr_grid(find(isnan(Cr_grid),1,'first')-1); %all nonvalues of the vector Cr_grid are changed by the real number which appears in the vector before the first NaN

if (k2_p(1)==1) %is selected so we must calculate and estimate k2'
    theta_3 = exp(log(k2a_min):(log(k2a_max)-log(k2a_min))/(basisF-1):log(k2a_max)); %rango logaritmico de valores para el theta_3        
    theta_3 = theta_3/60;    

    B=zeros(length(theta_3),num_files);
    M=zeros(2*length(theta_3),num_files); %for each theta_3 value, fill in two rows
    
    for i=1:1:length(theta_3) %for each basis function

        conv_i = resampling*conv(Cr_grid,exp(-theta_3(i)*tgrid));    
        B(i,:) = interp1(0:resampling:resampling*(length(conv_i)-1),conv_i,tseg,'linear'); %basis function i calculated for each frame

        A = [Cr B(i,:)'];      
        [Q,R]=qr(W*A,0); %QR decomposition reduced
        M(2*i-1:2*i,:)=R\Q'; %(2 x num_frames)

        theta=M(2*i-1:2*i,:)*W*Ct; %(2 x 1)
               
        Ct_i=theta(1)*Cr + theta(2)*B(i,:)'; %Stimated concentration for the basis function B_i in each frame (num_frames x 1)

        L=(Ct-Ct_i).^2; %Matrix to the square error (num_frames x 1)

        RSS_aux = w*L;        
        
        if (RSS_aux < RSS_menor)  %Look for the basis function's index which produce the minimum error
                RSS_menor = RSS_aux;                      
                i_min = i;
        end  
        waitbar(i/length(theta_3),waitbarhandle);

    end

    clear A Q R;

    theta=(M(2*i_min-1:2*i_min,:)*W*Ct);    
    theta_1=theta(1);   
    theta_2=theta(2);    

    %Obtain the parameters
    theta_3_fin=theta_3(i_min);
    
    R1_fin=theta_1;    
    k2a_fin=theta_3_fin;    
    k2_fin=theta_1*k2a_fin+theta_2;    
    BPnd_fin=(k2_fin/k2a_fin)-1;    
    k2_P = (k2_fin/R1_fin);  

    Ct_estimate=theta(1)*Cr+theta(2)*B(i_min,:)' ; %stimated concentration for each frame
    
    % Goodness of fit
    NMSE =1-((Ct-Ct_estimate).^2)./((Ct-mean(Ct)).^2);
    CCoef=corrcoef(Ct,Ct_estimate);
    
    % Modelo: Ct = R1_fin*Cr+(R1_fin*(k2_P-k2a_fin))*B(i_min,:)';
    Ct_R1increased = R1_fin*1.01*Cr+((R1_fin*1.01)*(k2_P-k2a_fin))*B(i_min,:)';
    
    Ct_k2_Pincreased = R1_fin*Cr+(R1_fin*((k2_P*1.01)-k2a_fin))*B(i_min,:)';

    conv_i = resampling*conv(Cr_grid,exp(-(k2a_fin*1.01)*tgrid)); %Convolution
    B_increased = interp1(0:resampling:resampling*(length(conv_i)-1),conv_i,tseg,'linear'); %Interpolate
    Ct_k2aincreased =  R1_fin*Cr+(R1_fin*(k2_P-(k2a_fin*1.01)))*B_increased';

    % Sensitivity functions
    Sens=100*[(Ct_R1increased-Ct_estimate)./Ct_estimate,(Ct_k2_Pincreased-Ct_estimate)./Ct_estimate,...
       (Ct_k2aincreased-Ct_estimate)./Ct_estimate];
    Sens(1,:)=0; %Quitamos NaN

    for i=1:1:3
        for j=1:1:3
            Aux=cumsum(diff(times,1,2).*(Sens(:,i).*Sens(:,j)));
            SM(i,j)=Aux(end);
        end
    end
    % Inv_SM=inv(SM);
    CovarMatrix=SM\eye(3);
    Variance_Param=diag(CovarMatrix);
    CorrM_Param=CovarMatrix./(sqrt(Variance_Param*Variance_Param'));
    
    
    
    %Variance calculation as SRTM (to skip round errors)
%     Ct_R1increased = R1_fin*1.01*Cr+(k2_fin-((R1_fin*1.01*k2_fin)/(1+BPnd_fin)))*B(i_min,:)';
% 
%     conv_i = resampling*conv(Cr_grid,exp(-(((k2_fin)*1.01)/(1+BPnd_fin))*tgrid)); %Convolution
%     B_increased = interp1(0:resampling:resampling*(length(conv_i)-1),conv_i,tseg,'linear'); %Interpolate
%     Ct_k2increased = R1_fin*Cr+((k2_fin)*1.01-((R1_fin*(k2_fin)*1.01)/(1+BPnd_fin)))*B_increased';
% 
%     conv_i = resampling*conv(Cr_grid,exp(-((k2_fin)/(1+BPnd_fin*1.01))*tgrid)); %Convolution
%     B_increased = interp1(0:resampling:resampling*(length(conv_i)-1),conv_i,tseg,'linear'); %Interpolate
%     Ct_BPndincreased = R1_fin*Cr+((k2_fin)-((R1_fin*(k2_fin))/(1+BPnd_fin*1.01)))*B_increased';
% 
%     % Sensitivity functions
%     Sens=100*[(Ct_R1increased-Ct_estimate)./Ct_estimate,(Ct_k2increased-Ct_estimate)./Ct_estimate,...
%        (Ct_BPndincreased-Ct_estimate)./Ct_estimate];
%     Sens(1,:)=0; %Remove NaN
% 
%     for i=1:1:3
%         for j=1:1:3
%             Aux=cumsum(diff(times,1,2).*(Sens(:,i).*Sens(:,j)));
%             SM(i,j)=Aux(end);
%         end
%     end
% 
%     CovarMatrix=SM\eye(3);% inv(SM);
%     Variance_Param=diag(CovarMatrix);
    
%     dof=size(times,1)-3;% degree of freedom (#number of measurements-fitted parameters)
%     chisquared=RSS_menor/dof; % reduced chi squared
%     dk2=2*sqrt(Variance_Param(2)*chisquared); % 2*standard error (from chisquared formula)
%     dBP=2*sqrt(Variance_Param(3)*chisquared);
%     dk2a=abs(k2a_fin)*(dk2/abs(k2_fin)+dBP/abs(1+BPnd_fin)); % from propagation error (k2a=k2/(1+BP))
%     vk2a=((dk2a/2)^2)/chisquared; % from chisquared formula
% 
%     dR1=2*sqrt(Variance_Param(1)*chisquared);
%     dk2_P=abs(k2_P)*(dk2/abs(k2_fin)+dR1/abs(R1_fin));
%     vk2_P=((dk2_P/2)^2)/chisquared;
    

    %Variance of k2_P and k2a in 1/min
    Variance_Param=[Variance_Param(2)*3600, Variance_Param(1), Variance_Param(3)*3600];
    
    results = {k2_P*60, R1_fin, k2a_fin*60, BPnd_fin, mean(NMSE), CCoef(1,2), CorrM_Param, Variance_Param};
    
    plots = Ct_estimate';
    
    % obtain the basis function using the k2' calculated
    for i=1:1:length(theta_3) %for each basis function
        conv_i = resampling*conv(Cr_grid,exp(-theta_3(i)*tgrid));    
        B(i,:) = Cr' + (k2_P-theta_3(i))*interp1(0:resampling:resampling*(length(conv_i)-1),conv_i,tseg,'linear'); %basis function i calculated for each frame
    end
    
else %k2' is not selected, so we consider the given value
    theta_3 = exp(log(k2a_min):(log(k2a_max)-log(k2a_min))/(basisF-1):log(k2a_max)); %logarithmic range of theta_3's values
    theta_3 = theta_3/60;
    
    k2_P=k2_p(2)/60;    

    B=zeros(length(theta_3),num_files);
        
    for i=1:1:length(theta_3) %for each basis function

        conv_i = resampling*conv(Cr_grid,exp(-theta_3(i)*tgrid));    
        B(i,:) = Cr' + (k2_P-theta_3(i))*interp1(0:resampling:resampling*(length(conv_i)-1),conv_i,tseg,'linear'); %basis function i calculated for each frame

        RI = (B(i,:)*Ct)/sum(B(i,:).^2);

        L=(Ct-RI*B(i,:)').^2; %Matrix to the square error (num_frames x 1)

        RSS_aux = w*L;
        if RSS_aux < RSS_menor %Look for the basis function's index which produce the minimum error
            RSS_menor = RSS_aux;
            i_min = i;
        end 

        waitbar(i/length(theta_3),waitbarhandle);
    end

        %Obtain and save the stimated parameters
        theta_3_fin=theta_3(i_min); 
        
        R1_fin=(B(i_min,:)*Ct)/sum(B(i_min,:).^2);
        k2a_fin=theta_3_fin;
        BPnd_fin=R1_fin*(k2_P/k2a_fin)-1;

        Ct_estimate=R1_fin*B(i_min,:); %stimated concentration for each frame
              
        % Goodness of fit
%         MSE=mean((Ct-Ct_estimate').^2);
        NMSE =1-((Ct-Ct_estimate').^2)./((Ct-mean(Ct)).^2);
        CCoef=corrcoef(Ct,Ct_estimate');
        
        %Modelo: Ct = R1_fin*B(i_min,:)';
        Ct_R1increased = R1_fin*1.01*B(i_min,:);
    
        conv_i = resampling*conv(Cr_grid,exp(-k2a_fin*tgrid));    
        B_increased = Cr' + ((k2_P*1.01)-k2a_fin)*interp1(0:resampling:resampling*(length(conv_i)-1),conv_i,tseg,'linear'); 
        Ct_k2_Pincreased = R1_fin*B_increased;

        conv_i = resampling*conv(Cr_grid,exp(-(k2a_fin*1.01)*tgrid)); %Convolution
        B_increased = Cr' + (k2_P-(k2a_fin*1.01))*interp1(0:resampling:resampling*(length(conv_i)-1),conv_i,tseg,'linear');
        Ct_k2aincreased =  R1_fin*B_increased';

        % Sensitivity functions
        Sens=100*[(Ct_R1increased'-Ct_estimate')./Ct_estimate',(Ct_k2_Pincreased'-Ct_estimate')./Ct_estimate',...
           (Ct_k2aincreased-Ct_estimate')./Ct_estimate'];
        Sens(1,:)=0; %Quitamos NaN

        for i=1:1:3
            for j=1:1:3
                Aux=cumsum(diff(times,1,2).*(Sens(:,i).*Sens(:,j)));
                SM(i,j)=Aux(end);
            end
        end
        % Inv_SM=inv(SM);
        CovarMatrix=SM\eye(3);
        Variance_Param=diag(CovarMatrix);
        CorrM_Param=CovarMatrix./(sqrt(Variance_Param*Variance_Param'));
        
%         dR1=sqrt(Variance_Param(1)); 
%         dk2_P=sqrt(Variance_Param(2)); 
%         dk2a=sqrt(Variance_Param(3));
%         dBP=abs(BPnd_fin)*(dR1/abs(R1_fin)+dk2a/abs(k2a_fin));
%         Variance_Param(4)=dBP^2;        
        
        %Variance of k2_P and k2a in 1/min       
        Variance_Param=[0,Variance_Param(1),Variance_Param(3)*3600];
        results = {k2_P*60, R1_fin, k2a_fin*60, BPnd_fin, mean(NMSE), CCoef(1,2), CorrM_Param, Variance_Param};


        %Model: Ct = R1_fin*B(i_min,:)'; BPnd version
%         conv_i = resampling*conv(Cr_grid,exp(-k2a_fin*tgrid));    
%         B_increased = Cr' + ((k2a_fin*(1+BPnd_fin)/(R1_fin*1.01))-k2a_fin)*interp1(0:resampling:resampling*(length(conv_i)-1),conv_i,tseg,'linear'); 
%         Ct_R1increased = R1_fin*1.01*B_increased;
%          
%         conv_i = resampling*conv(Cr_grid,exp(-k2a_fin*tgrid));    
%         B_increased = Cr' + ((k2a_fin*(1+(BPnd_fin*1.01))/R1_fin)-k2a_fin)*interp1(0:resampling:resampling*(length(conv_i)-1),conv_i,tseg,'linear'); 
%         Ct_BPincreased = R1_fin*B_increased;
% 
%         conv_i = resampling*conv(Cr_grid,exp(-(k2a_fin*1.01)*tgrid)); %Convolution
%         B_increased = Cr' + (((k2a_fin*1.01)*(1+BPnd_fin)/R1_fin)-(k2a_fin*1.01))*interp1(0:resampling:resampling*(length(conv_i)-1),conv_i,tseg,'linear');
%         Ct_k2aincreased =  R1_fin*B_increased';
%         
%         % Sensitivity functions
%         Sens=100*[(Ct_R1increased'-Ct_estimate')./Ct_estimate',(Ct_BPincreased'-Ct_estimate')./Ct_estimate',...
%            (Ct_k2aincreased-Ct_estimate')./Ct_estimate'];
%         Sens(1,:)=0; %Quitamos NaN
% 
%         for i=1:1:3
%             for j=1:1:3
%                 Aux=cumsum(diff(times,1,2).*(Sens(:,i).*Sens(:,j)));
%                 SM(i,j)=Aux(end);
%             end
%         end
%         % Inv_SM=inv(SM);
%         CovarMatrix=SM\eye(3);
%         Variance_Param=diag(CovarMatrix);
%         CorrM_Param=CovarMatrix./(sqrt(Variance_Param*Variance_Param'));
%         
%         dR1=sqrt(Variance_Param(1)); 
%         dBP=sqrt(Variance_Param(2)); 
%         dk2a=sqrt(Variance_Param(3));
%         
%         Variance_Param(3)=Variance_Param(3)*3600;
%         
%         Variance_Param=[0,Variance_Param(1),Variance_Param(3),Variance_Param(2)];
%         results = {k2_P*60, R1_fin, k2a_fin*60, BPnd_fin, mean(NMSE), CCoef(1,2), CorrM_Param, Variance_Param};
%         
        plots = Ct_estimate;
end   
end