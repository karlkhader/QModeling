function [results, plots] = QM_preprocessTwoTCM(Ct,Cp,times,vB,k4_checked,alpha1,alpha2,resampling,basisF,waitbarhandle)
% QM_preprocessTwoTCM function
% function [results plots] = QM_preprocessTwoTCM(Ct,Cp,times,vB,k4_checked,alpha1,alpha2,resampling,basisF,waitbarhandle)
%
% QM_preprocessTwoTCM is a function which fits the 2TCM model to the 
% current TACs and return the fitted TAC and the result parameters.
% The input data come from the previous calculated TACs and the 
% parameters introduced by the user in the GUI. 
%
% Parameters:
% Ct -> Vector containing the TAC from the current ROI
% Cp -> Vector containing the TAC from the plasma region
% times -> timetable of current study
% vB -> to check if the volume of blood parameter will be fitted or it is
% fixed to sme value specified by the user
% k4_checked -> to check if the k4 para meter will be fitted or it is fixed to 0
% alpha1 -> conditions and minor/mayor values of alpha 1
% alpha2 -> conditions and minor/mayor values of alpha 2
% resampling -> resampling to define the grid for convolution
% basisF -> number of basis functions
% waitbarhandle -> handle of the current waitbar
%---------------------------------------------------------------------------------
% QM_preprocessTwoTCM is part of QModeling.
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
% Copyright (C) 2015, 2018 Francisco Javier Lopez Gonzalez, Jose Paredes Pacheco, 
% Karl-Khader Thurnhofer Hemsi

global QMlastPreprocess

format long;

if isempty(QMlastPreprocess.TwoTCM.Cblood)
    Cb=Cp;
else
    Cb=QMlastPreprocess.TwoTCM.Cblood;
end

num_files = length(Ct);
t = ((times(1:num_files,1)+times(1:num_files,2))/2)/60;
tseg = t'*60;
RSS_menor = Inf;
i_min = 1;
j_min = 1;
k_min=1;
k=1;
% Interpolate
tgrid = 0:resampling:times(num_files,2);
Cp_grid = interp1([0 tseg],[0;Cp],tgrid,'linear');
Cp_grid(isnan(Cp_grid)) = Cp_grid(find(isnan(Cp_grid),1,'first')-1);

%Create auxiliary matrix
w = ones(1,num_files); 
W = eye(num_files);

if k4_checked %k4 will be fitted

    if alpha1(1) && alpha2(1) %alpha1 and alpha2 restricts are checked
        %Logarithmic scales for alpha1 y alpha2
        a1 = exp(log(alpha1(2)):(log(alpha1(3))-log(alpha1(2)))/(basisF-1):log(alpha1(3)));
        a1 = a1/60;       

        a2 = exp(log(alpha2(2)):(log(alpha2(3))-log(alpha2(2)))/(basisF-1):log(alpha2(3)));
        a2 = a2/60;     
        
    elseif alpha1(1) %only alpha1 restrict is checked
        %Logarithmic scales for alpha1
        a1 = exp(log(alpha1(2)):(log(alpha1(3))-log(alpha1(2)))/(basisF-1):log(alpha1(3)));
        a1 = a1/60;

        a2 = alpha2(2)/60;

    elseif alpha2(1) %only alpha2 restrict is checked
        %Logarithmic scales for alpha2
        a2 = exp(log(alpha2(2)):(log(alpha2(3))-log(alpha2(2)))/(basisF-1):log(alpha2(3)));
        a2 = a2/60;

        a1 = alpha1(2)/60;

    else %neither alpha1 nor alpha2 restricts are checked
        a1 = alpha1(2)/60;
        a2 = alpha2(2)/60;

    end 
    
else % k4 == 0
%      a1=zeros(1,basisF);
     a1=0;
     if alpha2(1)
        a2 = exp(log(alpha2(2)):(log(alpha2(3))-log(alpha2(2)))/(basisF-1):log(alpha2(3)));
        a2 = a2/60;
     else
         a2=alpha2(2)/60;
     end     
end

B1=zeros(length(a1),num_files);
B2=zeros(length(a2),num_files);

for i=1:1:length(a1)
    conv_i = resampling*conv(Cp_grid,exp(-a1(i)*tgrid)); %Convolution
    B1(i,:) = interp1(0:resampling:resampling*(length(conv_i)-1),conv_i,tseg,'linear'); %Interpolate
end
for j=1:1:length(a2)
    conv_j = resampling*conv(Cp_grid,exp(-a2(j)*tgrid)); %Convolution
    B2(j,:) = interp1(0:resampling:resampling*(length(conv_j)-1),conv_j,tseg,'linear'); %Interpolate
end


for i=1:length(a1)
    for j=1:length(a2)
        A = [B1(i,:)' B2(j,:)' Cb];
        [Q,R]=qr(W*A,0); % QR factoritation
        M(3*k-2:3*k,:)=R\Q';
%             A = [B1(i,:)' B2(j,:)' Cb];      
%             [U,S,V]=svd(W*A,0); %Singular Value Decomposition
%             M(3*k-2:3*k,:)=V*(S\U'); 
        sol=M(3*k-2:3*k,:)*W*Ct;
        
        if vB(1)
            Ct_ij=sol(1)*B1(i,:)'+sol(2)*B2(j,:)'+sol(3)*Cb;
        else
            Ct_ij=sol(1)*B1(i,:)'+sol(2)*B2(j,:)'+vB(2)*Cb;
        end
        
        L=(Ct-Ct_ij).^2;

        RSS_aux=w*L;
        if RSS_aux < RSS_menor
            RSS_menor = RSS_aux;
            i_min =i;
            j_min =j;
            k_min=k;
        end 
        k=k+1;                
    end
    waitbar(k/(length(a1)*length(a2)),waitbarhandle)
end
clear A Q R;

sol=M(3*k_min-2:3*k_min,:)*W*Ct;
theta_1 = sol(1);
theta_2 = sol(2);
if vB(1)
    vB_fin = sol(3);    
else
    vB_fin = vB(2);
end
            
Ct_estimate=theta_1*B1(i_min,:)'+theta_2*B2(j_min,:)'+vB_fin*Cb;

%Get the parameters
a1_fin=a1(i_min);
a2_fin=a2(j_min);

K1 = (theta_1 + theta_2)/(1-vB_fin);
% K1 = K1*(K1>0) + 0*(K1<=0);
k2 = (theta_1 * a1_fin + theta_2 * a2_fin)/(theta_1 + theta_2);
k4 = (a1_fin * a2_fin)/ k2;
k3 = a1_fin + a2_fin - k2 -k4;
% Ki = (K1*k3)/(k2+k3);

% Goodness of fit
NSME = 1- ((Ct-Ct_estimate).^2)./((Ct-mean(Ct)).^2);
CorrM_C = corrcoef(Ct,Ct_estimate);

Ct_vBincreased = ((1-vB_fin*1.01)*K1/(a2_fin-a1_fin))*((k3+k4-a1_fin)*B1(i_min,:)'+ (a2_fin-k3-k4)*B2(j_min,:)')+vB_fin*1.01*Cb;
Ct_K1increased = ((1-vB_fin)*K1*1.01/(a2_fin-a1_fin))*((k3+k4-a1_fin)*B1(i_min,:)'+ (a2_fin-k3-k4)*B2(j_min,:)')+vB_fin*Cb;


conv_i = resampling*conv(Cp_grid,exp(-((k2*1.01+k3+k4-sqrt((k2*1.01+k3*k4)^2-4*k2*1.01*k4))/2)*tgrid)); %Convolution
B1_increased = interp1(0:resampling:resampling*(length(conv_i)-1),conv_i,tseg,'linear'); %Interpolate
conv_j = resampling*conv(Cp_grid,exp(-((k2*1.01+k3+k4+sqrt((k2*1.01+k3*k4)^2-4*k2*1.01*k4))/2)*tgrid)); %Convolution
B2_increased = interp1(0:resampling:resampling*(length(conv_j)-1),conv_j,tseg,'linear'); %Interpolate

Ct_k2increased = ((1-vB_fin)*K1/(sqrt((k2*1.01+k3*k4)^2-4*k2*1.01*k4)))*((k3+k4-(k2*1.01+k3+k4-sqrt((k2*1.01+k3*k4)^2-4*k2*1.01*k4))/2)*B1_increased'+...
    ((k2*1.01+k3+k4+sqrt((k2*1.01+k3*k4)^2-4*k2*1.01*k4))/2-k3-k4)*B2_increased')+vB_fin*Cb;

conv_i = resampling*conv(Cp_grid,exp(-((k2+k3*1.01+k4-sqrt((k2+k3*1.01*k4)^2-4*k2*k4))/2)*tgrid)); %Convolution
B1_increased = interp1(0:resampling:resampling*(length(conv_i)-1),conv_i,tseg,'linear'); %Interpolate
conv_j = resampling*conv(Cp_grid,exp(-((k2+k3*1.01+k4+sqrt((k2+k3*1.01*k4)^2-4*k2*k4))/2)*tgrid)); %Convolution
B2_increased = interp1(0:resampling:resampling*(length(conv_j)-1),conv_j,tseg,'linear'); %Interpolate

Ct_k3increased = ((1-vB_fin)*K1/(sqrt((k2+k3*1.01*k4)^2-4*k2*k4)))*((k3*1.01+k4-(k2+k3*1.01+k4-sqrt((k2+k3*1.01*k4)^2-4*k2*k4))/2)*B1_increased'+...
    ((k2+k3*1.01+k4+sqrt((k2+k3*1.01*k4)^2-4*k2*k4))/2-k3*1.01-k4)*B2_increased')+vB_fin*Cb;
% Ct_k3increased = ((1-vB_fin)*K1/(a2_fin-a1_fin))*((k3*1.01+k4-a1_fin)*B1(i_min,:)'+ (a2_fin-k3*1.01-k4)*B2(j_min,:)')+vB_fin*Cb;

conv_i = resampling*conv(Cp_grid,exp(-((k2+k3+k4*1.01-sqrt((k2+k3*k4*1.01)^2-4*k2*k4*1.01))/2)*tgrid)); %Convolution
B1_increased = interp1(0:resampling:resampling*(length(conv_i)-1),conv_i,tseg,'linear'); %Interpolate
conv_j = resampling*conv(Cp_grid,exp(-((k2+k3+k4*1.01+sqrt((k2+k3*k4*1.01)^2-4*k2*k4*1.01))/2)*tgrid)); %Convolution
B2_increased = interp1(0:resampling:resampling*(length(conv_j)-1),conv_j,tseg,'linear'); %Interpolate

Ct_k4increased = ((1-vB_fin)*K1/(sqrt((k2+k3*k4*1.01)^2-4*k2*k4*1.01)))*((k3+k4*1.01-(k2+k3+k4*1.01-sqrt((k2+k3*k4*1.01)^2-4*k2*k4*1.01))/2)*B1_increased'+...
    ((k2+k3+k4*1.01+sqrt((k2+k3*k4*1.01)^2-4*k2*k4*1.01))/2-k3-k4*1.01)*B2_increased')+vB_fin*Cb;
% Ct_k4increased = ((1-vB_fin)*K1/(a2_fin-a1_fin))*((k3+k4*1.01-a1_fin)*B1(i_min,:)'+ (a2_fin-k3-k4*1.01)*B2(j_min,:)')+vB_fin*Cb;

% Sensitivity functions
Sens=100*[(Ct_vBincreased-Ct_estimate)./Ct_estimate,(Ct_K1increased-Ct_estimate)./Ct_estimate,...
   (Ct_k2increased-Ct_estimate)./Ct_estimate,(Ct_k3increased-Ct_estimate)./Ct_estimate,(Ct_k4increased-Ct_estimate)./Ct_estimate];
Sens(1,:)=0; %Remove NaN

for i=1:1:5
    for j=1:1:5
        Aux=cumsum(diff(times,1,2).*(Sens(:,i).*Sens(:,j)));
        SM(i,j)=Aux(end);
    end
end

CovarMatrix=SM\eye(5);% inv(SM);
Variance_Param=diag(CovarMatrix);
CorrM_Param=CovarMatrix./(sqrt(Variance_Param*Variance_Param'));


% CorrM_Param=[];
% Variance_Param=[];
K1 = K1*(K1>0) + 0*(K1<=0);
results = {vB_fin,K1*60,k2*60,k3*60,k4*60,mean(NSME),CorrM_C(1,2),CorrM_Param,Variance_Param};

%k3/k4
%K_i=60*(K1*k3/(k2+k3)) 

 B={B1 B2};
 alpha={a1 a2};
 QMlastPreprocess.TwoTCM.B=B;
 QMlastPreprocess.TwoTCM.alpha=alpha; 
 QMlastPreprocess.TwoTCM.M=M;
 plots=Ct_estimate;
 
end