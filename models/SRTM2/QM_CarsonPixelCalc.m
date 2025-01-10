function images = QM_CarsonPixelCalc(Cpet,B,th3,threshold,k2_p,options,waitbarhandle)
% QM_CarsonPixelCalc function
% function images = QM_CarsonPixelCalc(Cpet,B,th3,threshold,k2_p,options,waitbarhandle)
%
% QM_CarsonPixelCalc is a function with calculate parametric images
% applying the SRTM2 model at voxel level. The input data come 
% from the previous preprocessing and the parameters introduced by the user
%
% Parameters:
% Cpet -> Matrix with concentrations of PET study
% B -> data of the basis functions used in the last preprocessing
% th3 -> theta 3 calculated in the last preprocessing
% threshold -> threshold to apply in the parametric images
% k2_p -> k2' of the last preprocessing
% options -> options defined by the user previously in the GUI
% waitbarhandle -> handle of the current waitbar
%---------------------------------------------------------------------------------
% QM_CarsonPixelCalc is part of QModeling.
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
s=size(Cpet);
n=size(B);
num_vox=s(2);
w = ones(1,n(2));

RSS_menor=Inf(1,num_vox); %Initialize the minimum square error for the basis function i(num_vox x 1) 
i_min=ones(num_vox,1);
RI_temp=Inf(1,num_vox);

% Calculate the cut (weighted by the length of the frame)
energy = diff(QMmaindata.times,1,2)'*(Cpet.^2);
cutoff = max(energy)*(threshold/100);
indcut = find(energy < cutoff);

for i=1:1:length(th3) %for each basis function
    RI = (B(i,:)*Cpet)/sum(B(i,:).^2);    
    L=(Cpet-B(i,:)'*RI).^2; %Matrix to the square error (num_frames x 1)
    RSS_aux = w*L;
    
    indices_men=find(RSS_aux<RSS_menor); 
    %Take the minimum values of the mean square error for each voxel and the basis function's number which generate it
    RSS_menor(indices_men)=RSS_aux(indices_men);
    i_min(indices_men)=i;
    
    RI_temp(indices_men)=RI(indices_men); 

    waitbar(i/length(th3),waitbarhandle);
end

theta_3_fin=th3(i_min); %for each voxel, load the value of theta_3 which produces the minimum square error 

R1=RI_temp;
k2a=theta_3_fin*60;
%k2a=theta_3_fin;
k2=k2_p*R1;
BPnd=((R1./k2a)*k2_p) -1;

%Apply the cut(all under this limit is change by zero)
BPnd(indcut)=0;
k2(indcut)=0;
k2a(indcut)=0;
R1(indcut)=0;

%Save the images
images(1,:)=BPnd;
images(2,:)=R1;
images(3,:)=k2a;
images(4,:)=k2;

%Apply restrictions
if options(1,1) == 1 %to BPnd
    images(1,images(1,:)<options(1,2)) = options(1,2);
    images(1,images(1,:)>options(1,3)) = options(1,3);
end
if options(2,1) == 1 %to R1
    images(2,images(2,:)<options(2,2)) = options(2,2);
    images(2,images(2,:)>options(2,3)) = options(2,3);
end
if options(3,1) == 1 %to k2a
    images(3,images(3,:)<options(3,2)) = options(3,2);
    images(3,images(3,:)>options(3,3)) = options(3,3);
end
if options(4,1) == 1 %to k2a 
    images(4,images(4,:)<options(4,2)) = options(4,2);
    images(4,images(4,:)>options(4,3)) = options(4,3);
end