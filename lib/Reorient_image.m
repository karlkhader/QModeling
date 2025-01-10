% Reorient_image function
%
% Reorient the selected images to correct turns and displacements
% and have the identity matrix as the affine transformation matrix
% of the images. It is a function executed by QModeling to reorient
% images when it is necesary.
%
% Parameters:
% PP -> paths of the images to reorient
%---------------------------------------------------------------------------------
% Reorient_image is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% Reorient_image is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with Reorient_image.  If not, see <http://www.gnu.org/licenses/>.
%---------------------------------------------------------------------------------
% Copyright (C) 2015 Jos� Paredes Pacheco

function Reorient_image(PP)
try 
    if nargin<1, PP = spm_select(Inf,'image','Select images to reorient'); end;
    VV = spm_vol(PP);

    for V=VV'
        d = V.dim(1:3);
        c = [	1    1    1    1
            1    1    d(3) 1
            1    d(2) 1    1
            1    d(2) d(3) 1
            d(1) 1    1    1
            d(1) 1    d(3) 1
            d(1) d(2) 1    1
            d(1) d(2) d(3) 1]'; %Coordenadas en vóxeles de los vértices del cubo

        %A esas coordenadas se les aplica la matriz vóxel mundo de la imagen
        %original
        tc = V.mat(1:3,1:4)*c;
        %Se encuentra el vértice que tiene la x más grande, el que
        %tiene la y más grande y e que tiene la z más grande, y se guardan
        %esas x, y y z.
        mx = round(max(tc,[],2)');

        %Se encuentra el vértice que tiene la x más pequeña, el que
        %tiene la y más pequeña y e que tiene la z más pequeña, y se guardan
        %esas x, y y z.
        mn = round(min(tc,[],2)');

        %se genera una matriz vóxel - mundo sin giros ni cizallas, con el 
        %origen de coordenadas en la combinación de x, y y z más pequeñas 
        %encontradas en el apartado anterior. Nosotros añadimos que el tamaño 
        %de vóxel sea  igual que el de la imagen original. En la versión de
        %Ashburner era 1 mm x 1 mm x 1 mm.

        %Sacamos info sobre tamaño de voxel, voxel origen y giros
        P = spm_imatrix(V.mat);

        mat = spm_matrix(cat(2,[mn-1],[0 0 0 P(1,7) P(1,8) P(1,9)]));
        %dim = mx-mn; Solo válido si de partida nuestros vóxeles son de 1x1x1mm
        dim = (mat\[mx 1]')';

        VO               = V;
        [lpath,name,ext] = fileparts(V.fname);
        VO.fname         = fullfile(lpath,['r' name ext]);
        VO.dim(1:3)      = round(abs(dim(1:3)));
        VO.mat           = mat;

        %VO = spm_create_image(VO);
        %N.Roe (26/11/2013) to adapt it to SPM8 routines:
        VO = spm_create_vol(VO);

        for i=1:VO.dim(3)
            % Partiendo del slice original inclinada, esta M lo mueve al slice final i, 
            % con tamaño de vóxel y voxel origen el que queremos indicado por V0.mat 
            % y deshacemos el giro con inv(V.mat)  
            M = inv(spm_matrix([0 0 -i])*inv(VO.mat)*V.mat); %==inv(V.mat)*VO.mat*spm_matrix([0 0 i])
            img = spm_slice_vol(V,M,VO.dim(1:2),1);
            spm_write_plane(VO,img,i);
        end
    end
catch 
    errordlg(char({'An error ocurred during the reorient, revise the PET study'}),'Reorient error','modal');
end
return