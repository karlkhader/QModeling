function images = QM_LoganPixelCalc(X,Cpet,time,threshold,options,fr,waitbarhandle)
% QM_LoganPixelCalc function
% function images = QM_LoganPixelCalc(X,Cpet,time,threshold,options,fr,waitbarhandle)
%
% QM_LoganPixelCalc is a function with calculate parametric images
% applying the Logan plot model at voxel level. The input data come 
% from the previous preprocessing and the parameters introduced by the
% user.
%
% Parameters:
% X -> some calculation from the last preprocessing
% Cpet -> Matrix containing the TACs of the ROIs
% time -> timetable of the study
% threshold -> threshold to apply in the parametric images
% options -> options defined by the user previously in the GUI
% fr -> the frame index of the t*
% waitbarhandle -> handle of the current waitbar
%---------------------------------------------------------------------------------
% QM_LoganPixelCalc is part of QModeling.
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


% Weighted
energy = diff(time,1,2)'*(Cpet.^2);
cutoff = max(energy)*(threshold/100);


% Variable start
nVoxels = size(Cpet,2);
num_files = size(Cpet);
RCpet = Cpet(fr:end,:);

images = zeros(2,nVoxels);
P = zeros(6,nVoxels);
y = zeros(num_files);

%lhs numerator of equation
for i=1:nVoxels
    y(:,i)=cumsum(Cpet(:,i).*diff(time,1,2)); 
end

%take values from t* obtained in the preprocess.
y=y(fr:end,:);
waitbar(2/5,waitbarhandle);

%Parametric image
for i = 1:nVoxels
  %Obtain only values differents of NaN or infinite
  nonz=find(RCpet(:,i));
  den=RCpet(nonz,i);
  
  %At least two points
  if(length(nonz)>=2)    
  Yc=y(nonz,i)./den;
  Xc=X(nonz)./den;
  
  P(1,i)=mean(Xc);
  P(2,i)=mean(Yc);
  P(3,i)=sum(Yc.*Xc)*length(nonz);
  P(4,i)=sum(Xc);
  P(5,i)=sum(Yc);
  P(6,i)=sum(Xc.^2)*length(nonz);
  
  end
end

% manual regresion 
images(1,:)=((P(3,:)-P(4,:).*P(5,:))./(P(6,:)-P(4,:).^2))-1;
images(2,:)=P(2,:)-(images(1,:)+1).*P(1,:);

waitbar(3/5,waitbarhandle);
  
%images with values introduced discriminatory by user
images(1,energy < cutoff) = 0;
images(2,energy < cutoff) = 0;

if options(1,1) == 1
    images(1,images(1,:)<options(1,2)) = options(1,2);
    images(1,images(1,:)>options(1,3)) = options(1,3);
end
if options(2,1) == 1
    images(2,images(2,:)<options(2,2)) = options(2,2);
    images(2,images(2,:)>options(2,3)) = options(2,3);
end
