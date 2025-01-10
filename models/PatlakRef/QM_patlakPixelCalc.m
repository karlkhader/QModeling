function images = QM_patlakPixelCalc(Cpet,Cr,tfit_ind,X,times,threshold,options,waitbarhandle)
% QM_patlakPixelCalc function
% function images = QM_patlakPixelCalc(Cpet,Cr,tfit_ind,X,times,threshold,options,waitbarhandle)
%
% QM_patlakPixelCalc is a function with calculate parametric images
% applying the Patlak model at voxel level. The input data come 
% from the previous preprocessing and the parameters introduced by the
% user.
%
% Parameters:
% Cpet -> Matrix containing the TACs of the ROIs
% Cr -> Vector containing the TAC of the Reference region
% tfit_ind -> t* introduced by the user
% X -> plots of last preprocessing
% times -> timetable of the study
% threshold -> threshold to apply in the parametric images
% options -> options defined by the user previously in the GUI
% waitbarhandle -> handle of the current waitbar
%---------------------------------------------------------------------------------
% QM_patlakPixelCalc is part of QModeling.
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


% Compute cut (ponderando por la longitud del frame)
energy = diff(times,1,2)'*(Cpet.^2);
cutoff = max(energy)*(threshold/100);

% Define voxel-TAC matrix
Yall = zeros(size(Cpet));

% Normalize data
for i = 1:length(Cr)   
    Yall(i,:) = Cpet(i,:)/Cr(i);
end

Yall(isnan(Yall)|isinf(Yall)) = 0;

waitbar(1/5,waitbarhandle)

% Clear the unuse variable
clear Cpet

Xf = X(tfit_ind:end);

% Construct Vandermonde matrix.
V(:,2) = ones(length(Xf),1,class(Xf));
V(:,1) = Xf.*V(:,2);

% Solve least squares problem.
[Q,R] = qr(V,0);

ws = warning('off','all'); 
images = R\(Q'*Yall(tfit_ind:end,:));    % Same as p = V\y;
warning(ws);

waitbar(2/5,waitbarhandle)

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
