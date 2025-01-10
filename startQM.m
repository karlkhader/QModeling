% startQM function
%
% startQM shows a welcome message when QModeling is starting

function startQM()
%---------------------------------------------------------------------------------
% startQM is part of QModeling.
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
% Copyright (C) 2021 Francisco Javier Lopez Gonzalez, Jose Paredes Pacheco, 
% Karl-Khader Thurnhofer Hemsi

clear all
clc

actualpath = mfilename('fullpath');
QMpath = actualpath(1:end-8);
image=imread(strcat(QMpath,filesep,'icons',filesep,'qmodeling.png'));
image_size=size(image);
height=image_size(1);
width=image_size(2);
screen_size=get(0,'screensize'); %we get the current screen size
pos=[(screen_size(3)-width)/2 (screen_size(4)-height)/2 width height];
fig=figure('Menubar','none','Name','Welcome to QModeling v 2.1 (UIM CIMES 2021)','NumberTitle','off','ToolBar','none','Visible','off','Position',pos);
axes('position',[0 0 1 1]);
xlabel('UIM CIMES');
imshow(image);
set(fig,'Visible','on');


disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('                  __________     _           _ ')
disp('                /           \   | \         / |')
disp('                |           |   |  \       /  |')
disp('                |           |   |   \     /   |')
disp('                |           |   |    \   /    |')
disp('                |          \|   |     \_/     |')
disp('                |           \   |             |')
disp('                \___________/\  |             |')
disp(' ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' ')
disp('Welcome to QModeling: Software for compartimental modeling')
disp(' ')
disp('--------------------------------------------------------------------')
disp('QModeling is free software: you can redistribute it and/or modify')
disp('it under the terms of the GNU General Public License as published by')
disp('the Free Software Foundation, either version 3 of the License, or')
disp('(at your option) any later version.')
disp(' ')
disp('QModeling is distributed in the hope that it will be useful,')
disp('but WITHOUT ANY WARRANTY; without even the implied warranty of')
disp('MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the')
disp('GNU General Public License for more details.')
disp(' ')
disp('You should have received a copy of the GNU General Public License')
disp('along with QModeling.  If not, see <http://www.gnu.org/licenses/>.')
disp('--------------------------------------------------------------------')
disp('Copyright (C) 2014, 2021 Francisco Javier Lopez Gonzalez, Jose Paredes Pacheco,') 
disp('Karl-Khader Thurnhofer Hemsi, Nuria Roe Vellve, Antonio L. Gutierrez Cardo,')
disp('Manuel Enciso Garcia Oliveros, Carlos Rossi Jimenez')
disp(' ')

pause(3)
close(fig)
clear all
