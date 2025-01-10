function [results plots] = QM_preprocessSRTM(Ct,Cr,times,k2a_min,k2a_max,num_BF,resampling,waitbarhandle)
% QM_preprocessSRTM function
% function [results plots] = QM_preprocessSRTM(Cpet,Ct,Cr,times,k2a_min,k2a_max,num_BF,resampling,waitbarhandle)
%
% QM_preprocessSRTM is a function which fit the SRTM model to the 
% current TACs and return the fitted TAC and the result parameters.
% The input data come from the previous calculated TACs and the 
% parameters introduced by the user in the GUI. 
%
% Parameters:
% Cpet -> Matrix with concentrations of PET study
% Ct -> Vector containing the TAC of the current ROI
% Cr -> Vector containing the TAC of the Reference region
% times -> timetable of current study
% k2amin -> k2a minimum defined by the user
% k2amax -> k2a maximum defined by the user
% num_BF -> number of basis functions
% resampling -> resampling to define the grid for convolution
% waitbarhandle -> handle of the current waitbar
%---------------------------------------------------------------------------------
% QM_preprocessSRTM is part of QModeling.
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

global QMlastPreprocess
global QMmaindata

format long;
%Obtain info
num_files = length(Ct);
t = ((times(1:num_files,1)+times(1:num_files,2))/2)/60;

%Logarithmic scale
theta_3 = exp(log(k2a_min):(log(k2a_max)-log(k2a_min))/(num_BF-1):log(k2a_max)); 
theta_3 = theta_3/60;

%Create auxiliary matrix
w = ones(1,num_files); 
B=zeros(length(theta_3),num_files);
M=zeros(2*length(theta_3),num_files); 
W=eye(num_files); 

tseg = t'*60;
RSS_menor = Inf; 
i_min = 1;

% Interpolate
tgrid = 0:resampling:times(num_files,2);
Cr_grid = interp1([0 tseg],[0;Cr],tgrid,'linear');
Cr_grid(isnan(Cr_grid)) = Cr_grid(find(isnan(Cr_grid),1,'first')-1);

%Calculations for the algorithm
for i=1:1:length(theta_3) 

    conv_i = resampling*conv(Cr_grid,exp(-theta_3(i)*tgrid)); %Convolution
    B(i,:) = interp1(0:resampling:resampling*(length(conv_i)-1),conv_i,tseg,'linear'); %Interpolate

    A = [Cr B(i,:)'];      
    [Q,R]=qr(W*A,0); % QR factoritation
    M(2*i-1:2*i,:)=R\Q'; 
    
    theta=M(2*i-1:2*i,:)*W*Ct; 
    
    Ct_i=theta(1)*Cr+theta(2)*B(i,:)'; 
    
    L=(Ct-Ct_i).^2; 
    
    RSS_aux = w*L;
    if RSS_aux < RSS_menor
        RSS_menor = RSS_aux;
        i_min = i;
    end
    
     waitbar(i/length(theta_3),waitbarhandle)
end


clear A Q R;

theta=(M(2*i_min-1:2*i_min,:)*W*Ct);
theta_1=theta(1);    
theta_2=theta(2);    

Ct_estimate=theta_1*Cr+theta_2*B(i_min,:)'; 

% Get the parameters
theta_3_fin=theta_3(i_min);
R1_fin=theta_1;
k2_fin=theta_1*theta_3_fin+theta_2;
k2a_fin=theta_3_fin;
BPnd_fin=(k2_fin/k2a_fin)-1;

% Goodness of fit
% MSE=mean((Ct-Ct_estimate).^2);
NMSE =1-((Ct-Ct_estimate).^2)./((Ct-mean(Ct)).^2);
CorrM_C=corrcoef(Ct,Ct_estimate);

% Model: Ct = R1_fin*Cr+(k2_fin/60-((R1_fin*k2_fin/60)/(1+BPnd_fin)))*B(i_min,:)';
Ct_R1increased = R1_fin*1.01*Cr+(k2_fin-((R1_fin*1.01*k2_fin)/(1+BPnd_fin)))*B(i_min,:)';

conv_i = resampling*conv(Cr_grid,exp(-(((k2_fin)*1.01)/(1+BPnd_fin))*tgrid)); %Convolution
B_increased = interp1(0:resampling:resampling*(length(conv_i)-1),conv_i,tseg,'linear'); %Interpolate
Ct_k2increased = R1_fin*Cr+((k2_fin)*1.01-((R1_fin*(k2_fin)*1.01)/(1+BPnd_fin)))*B_increased';

conv_i = resampling*conv(Cr_grid,exp(-((k2_fin)/(1+BPnd_fin*1.01))*tgrid)); %Convolution
B_increased = interp1(0:resampling:resampling*(length(conv_i)-1),conv_i,tseg,'linear'); %Interpolate
Ct_BPndincreased = R1_fin*Cr+((k2_fin)-((R1_fin*(k2_fin))/(1+BPnd_fin*1.01)))*B_increased';

% Sensitivity functions
Sens=100*[(Ct_R1increased-Ct_estimate)./Ct_estimate,(Ct_k2increased-Ct_estimate)./Ct_estimate,...
   (Ct_BPndincreased-Ct_estimate)./Ct_estimate];
Sens(1,:)=0; %Remove NaN

for i=1:1:3
    for j=1:1:3
        Aux=cumsum(diff(times,1,2).*(Sens(:,i).*Sens(:,j)));
        SM(i,j)=Aux(end);
    end
end

CovarMatrix=SM\eye(3);% inv(SM);
Variance_Param=diag(CovarMatrix);
CorrM_Param=CovarMatrix./(sqrt(Variance_Param*Variance_Param'));

% dof=size(times,1)-3;% degree of freedom (#number of measurements-fitted parameters)
% % chisquared=RSS_menor/dof; % reduced chi squared
% dk2=2*sqrt(Variance_Param(2)*chisquared); % 2*standard error (from chisquared formula)
% dBP=2*sqrt(Variance_Param(3)*chisquared);
% dk2a=abs(k2a_fin)*(dk2/abs(k2_fin)+dBP/abs(1+BPnd_fin)); % from propagation error (k2a=k2/(1+BP))
% Variance_Param(4)=((dk2a/2)^2)/chisquared; % from chisquared formula

dk2=sqrt(Variance_Param(2));
dBP=sqrt(Variance_Param(3));
% dk2a=abs(k2a_fin)*(dk2/abs(k2_fin)+dBP/abs(1+BPnd_fin));
% Variance_Param(4)=dk2a^2;

%Variance of k2 and k2a in 1/min
Variance_Param(2)=Variance_Param(2)*3600;
% Variance_Param(4)=Variance_Param(4)*3600;

% Results
% results = [R1_fin k2a_fin k2_fin BPnd_fin MSE CCoef(1,2)];
results = {R1_fin,k2_fin*60,BPnd_fin,k2a_fin*60,mean(NMSE),CorrM_C(1,2),CorrM_Param,Variance_Param};


% b = [R1_fin,k2_fin,BPnd_fin];
% % x=[Cr B(i_min,:)'];
% % 
% % modelfun=@(b,x) (b(1)*x(:,1)+(b(2)-((b(1)*b(2))/(1+b(3))))*x(:,2));
% x=[Cr times(:,2) tseg'];
% 
% modelfun=@(b,x) fSRTM(b,x(:,1),x(:,2),x(:,3));
% y=modelfun(b,x);
% 
% [beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(x,Ct,modelfun,b)

% Save preprocess
QMlastPreprocess.SRTM.M = M;
QMlastPreprocess.SRTM.B = B;
QMlastPreprocess.SRTM.numBF = num_BF;

plots(:,1)=Ct;
plots(:,2)=Ct_estimate;
plots(:,3)=Cr;
