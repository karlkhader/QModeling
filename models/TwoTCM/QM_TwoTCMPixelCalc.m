function images = QM_TwoTCMPixelCalc(options,waitbarhandle)
% QM_TwoTCMPixelCalc function
% function images = QM_TwoTCMPixelCalc(options,waitbarhandle)
%
% QM_TwoTCMPixelCalc is a function which calculates parametric images
% applying the 2TCM model at voxel level. The input data comes 
% from the previous preprocessing and the parameters introduced by the user
%
% Parameters:
% options -> options defined by the user previously in the GUI
% waitbarhandle -> handle of the current waitbar
%---------------------------------------------------------------------------------
% QM_TwoTCMPixelCalc is part of QModeling.
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

global QMmaindata
global QMlastPreprocess

k=1;
Cpet=QMmaindata.Cpet;
times=QMmaindata.times;
Cp=QMmaindata.Crois{2,QMlastPreprocess.TwoTCM.idCplasma};
threshold=QMlastPreprocess.TwoTCM.threshold;
M=QMlastPreprocess.TwoTCM.M;
B1=QMlastPreprocess.TwoTCM.B{1};
B2=QMlastPreprocess.TwoTCM.B{2};
a1=QMlastPreprocess.TwoTCM.alpha{1};
a2=QMlastPreprocess.TwoTCM.alpha{2};
vB=QMlastPreprocess.TwoTCM.vB;
k4_checked=QMlastPreprocess.TwoTCM.k4_checked;

if isempty(QMlastPreprocess.TwoTCM.Cblood)
    Cb=Cp;
else
    Cb=QMlastPreprocess.TwoTCM.Cblood;
end

s=size(Cpet);
num_vox=s(2);
W=eye(s(1));
w = ones(1,s(1));
t_minutes = ((times(1:s(1),1)+times(1:s(1),2))/2)/60;

RSS_menor=Inf(1,num_vox); %Initialize the minimum square error for the basis function i(num_vox x 1) 
theta_1_fin=Inf(num_vox,1);
theta_2_fin=Inf(num_vox,1);
vB_fin=Inf(num_vox,1);
i_min=ones(num_vox,1);
j_min=ones(num_vox,1);
k_min=ones(num_vox,1);


% Calculate the cut (weighted by the length of the frame)
energy = diff(QMmaindata.times,1,2)'*(Cpet.^2);
cutoff = max(energy)*(threshold/100);
indcut = find(energy < cutoff);

for i=1:1:length(a1)
    for j=1:1:length(a2)
        A = [B1(i,:)' B2(j,:)' Cb];         

        if vB(1) %vB will be fitted
            sol=M(3*k-2:3*k,:)*W*Cpet;
        else
            sol(1:2,:)=M(3*k-2:3*k-1,:)*W*Cpet;
            sol(3,:)=vB(2)*ones(1,s(2));
        end

        Ct_ij= A*sol;

        RSS_aux=(w*((Cpet-Ct_ij).^2));
        clear Ct_ij;

        indices_men=find(RSS_aux<RSS_menor); 

        RSS_menor(indices_men)=RSS_aux(indices_men);
        i_min(indices_men)=i;
        j_min(indices_men)=j;
        k_min(indices_men)=k;

        theta_1_fin(indices_men)=sol(1,indices_men);
        theta_2_fin(indices_men)=sol(2,indices_men);
        
        vB_fin(indices_men)=sol(3,indices_men);
                
        k=k+1;
        waitbar(k/(length(a1)*length(a2)),waitbarhandle)
    end
    
end


%Get the parametric images
% a1_fin=a1(i_min);
if length(a2)==1
    a2_fin=a2;
else
    a2_fin=a2(j_min);
end
K1 = (theta_1_fin + theta_2_fin)./(1-vB_fin);
% K1= K1*(K1>0) + 0*(K1<=0);

if k4_checked %k4 fitted
    a1_fin=a1(i_min);
    k2 = (theta_1_fin .* a1_fin' + theta_2_fin .* a2_fin')./(theta_1_fin + theta_2_fin);
    k4 = (a1_fin .* a2_fin)'./ k2;
    k3 = (a1_fin + a2_fin)' - k2 -k4;
    k3divk4 = k3./k4;    
else %alpha1=0=k4
    a1_fin=zeros(size(K1)); %this is for nothing
    k2 = (theta_2_fin .* a2_fin')./(theta_1_fin + theta_2_fin);
    k4 = zeros(size(K1));
    k3 = a2_fin' - k2;
    k3divk4 = zeros(size(K1));
end


Ki = (K1.*k3)./(k2+k3);
K1divk2 = K1./k2;
Vs = K1divk2.*k3divk4;
Vt = K1divk2.*(1+k3divk4);


%Apply the cut(all under this limit is change by zero)
vB_fin(indcut)=0;
K1(indcut)=0;
k2(indcut)=0;
k3(indcut)=0;
k4(indcut)=0;
Ki(indcut)=0;
K1divk2(indcut)=0;
k3divk4(indcut)=0;
Vs(indcut)=0;
Vt(indcut)=0;
a1_fin(indcut)=0;
a2_fin(indcut)=0;

%time_to_peak
[~,i_max]=max(Cpet,[],1);
TTP_image=t_minutes(i_max);



%Save the images
images(1,:)=vB_fin;
images(2,:)=K1.*60;
images(3,:)=k2.*60;
images(4,:)=k3.*60;
images(5,:)=k4.*60;
images(6,:)=Ki*60;
images(7,:)=Vs;
images(8,:)=Vt;
images(9,:)=K1divk2;
images(10,:)=k3divk4;
images(11,:)=a1_fin.*60;
images(12,:)=a2_fin.*60;

images(13,:)=TTP_image';



%Apply restrictions
if options(1,1) == 1 %to vB
    images(1,images(1,:)<options(1,2)) = options(1,2);
    images(1,images(1,:)>options(1,3)) = options(1,3);
end
if options(2,1) == 1 %to K1
    images(2,images(2,:)<options(2,2)) = options(2,2);
    images(2,images(2,:)>options(2,3)) = options(2,3);
end
if options(3,1) == 1 %to k2
    images(3,images(3,:)<options(3,2)) = options(3,2);
    images(3,images(3,:)>options(3,3)) = options(3,3);
end
if options(4,1) == 1 %to k3 
    images(4,images(4,:)<options(4,2)) = options(4,2);
    images(4,images(4,:)>options(4,3)) = options(4,3);
end
if options(5,1) == 1 %to k4 
    images(5,images(5,:)<options(5,2)) = options(5,2);
    images(5,images(5,:)>options(5,3)) = options(5,3);
end
if options(6,1) == 1 %to Ki 
    images(6,images(6,:)<options(6,2)) = options(6,2);
    images(6,images(6,:)>options(6,3)) = options(6,3);
end
if options(7,1) == 1 %to Vs 
    images(7,images(7,:)<options(7,2)) = options(7,2);
    images(7,images(7,:)>options(7,3)) = options(7,3);
end
if options(8,1) == 1 %to Vt
    images(8,images(8,:)<options(8,2)) = options(8,2);
    images(8,images(8,:)>options(8,3)) = options(8,3);
end
if options(9,1) == 1 %to K1divk2 
    images(9,images(9,:)<options(9,2)) = options(9,2);
    images(9,images(9,:)>options(9,3)) = options(9,3);
end
if options(10,1) == 1 %to k3divk5 
    images(10,images(10,:)<options(10,2)) = options(10,2);
    images(10,images(10,:)>options(10,3)) = options(10,3);
end
if options(11,1) == 1 %to alpha1 
    images(11,images(11,:)<options(11,2)) = options(11,2);
    images(11,images(11,:)>options(11,3)) = options(11,3);
end
if options(12,1) == 1 %to alpha2
    images(12,images(12,:)<options(12,2)) = options(12,2);
    images(12,images(12,:)>options(12,3)) = options(12,3);
end