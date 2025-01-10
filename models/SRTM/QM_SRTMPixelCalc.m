function images = QM_SRTMPixelCalc(Cpet,Cr,k2amin,k2amax,threshold,M,B,numBF,waitbarhandle,options)
% QM_SRTMPixelCalc function
% function images = QM_SRTMPixelCalc(Cpet,Cr,k2amin,k2amax,threshold,M,B,numBF,waitbarhandle,options)
%
% QM_SRTMPixelCalc is a function with calculate parametric images
% applying the SRTM model at voxel level. The input data come 
% from the previous preprocessing and the parameters introduced by the user
%
% Parameters:
% Cpet -> Matrix with concentrations of PET study
% Cr -> Vector containing the TAC of the Reference region
% k2amin -> k2a minimum defined by the user
% k2amax -> k2a maximum defined by the user
% threshold -> threshold to apply in the parametric images
% M -> matrix with data of the last preprocessing
% B -> data of the basis functions used in the last preprocessing
% numBF -> number of basis functions used
% waitbarhandle -> handle of the current waitbar
% options -> options defined by the user previously in the GUI
%---------------------------------------------------------------------------------
% QM_SRTMPixelCalc is part of QModeling.
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

global QMmaindata

num_files = size(Cpet,1);
num_vox = size(Cpet,2);

W=eye(num_files);
w = ones(1,num_files);

RSS_menor=Inf(1,num_vox); 
i_min=ones(num_vox,1);
theta_1_fin=Inf(num_vox,1);
theta_2_fin=Inf(num_vox,1);

% Logarithmic scale
theta_3 = exp(log(k2amin):(log(k2amax)-log(k2amin))/(numBF-1):log(k2amax));
theta_3 = theta_3/60;

% Calculate cut
energy = diff(QMmaindata.times,1,2)'*(Cpet.^2);
cutoff = max(energy)*(threshold/100);
indcut = find(energy < cutoff);

% Calculations
% mitad = floor(num_vox/2);
for i=1:1:length(theta_3)
    
    theta=M(2*i-1:2*i,:)*W*Cpet;

%     Ct_i(:,1:mitad)=Cr*theta(1,1:mitad)+B(i,:)'*theta(2,1:mitad);   
%     Ct_i(:,mitad+1:num_vox)=Cr*theta(1,mitad+1:end)+B(i,:)'*theta(2,mitad+1:end);     
    
    A = [Cr B(i,:)'];
    Ct_i = A*theta;
    
    RSS_aux=(w*((Cpet-Ct_i).^2));
    clear Ct_i;
    
    indices_men=find(RSS_aux<RSS_menor); 

    RSS_menor(indices_men)=RSS_aux(indices_men);
    i_min(indices_men)=i;
    
    theta_1_fin(indices_men)=theta(1,indices_men);
    theta_2_fin(indices_men)=theta(2,indices_men);
    
    waitbar(i/(length(theta_3)),waitbarhandle)
    
end
theta_3_fin=theta_3(i_min); 

% Get parametric images
R1=theta_1_fin;
k2=60*(theta_1_fin'.*theta_3_fin+theta_2_fin');
k2a=theta_3_fin*60;
BPnd=(k2-k2a)./k2a;

%Perform the cut
BPnd(indcut)=0;
k2(indcut)=0;
k2a(indcut)=0;
R1(indcut)=0;

% Apply retrictions
if options(1,1)==1
    positions = BPnd < options(1,2);
    BPnd(positions) = options(1,2);
    positions = BPnd > options(1,3);
    BPnd(positions) = options(1,3);
end
    
if options(2,1)==1
    positions = k2 < options(2,2);
    k2(positions) = options(2,2);
    positions = k2 > options(2,3);
    k2(positions) = options(2,3);
end
    
if options(3,1)==1
    positions = R1 < options(3,2);
    R1(positions) = options(3,2);
    positions = R1 > options(3,3);
    R1(positions) = options(3,3);
end
    
if options(4,1)==1
    positions = k2a < options(4,2);
    k2a(positions) = options(4,2);
    positions = k2a > options(4,3);
    k2a(positions) = options(4,3);
end
      
%Save   
images(1,:)=BPnd;    
images(2,:)=k2;
images(3,:)=R1;
images(4,:)=k2a;

end

