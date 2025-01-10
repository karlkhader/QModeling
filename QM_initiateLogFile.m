function QM_initiateLogFile()
%---------------------------------------------------------------------------------
% QM_initiateLogFile is part of QModeling.
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
% Copyright (C) 2016 Francisco Javier Lopez Gonzalez, Karl-Khader Thurnhofer Hemsi

    global QMf_logfile
    global QMpath
    
    dir_logfile=dir(strcat(QMpath,filesep,'logfile.txt'));

    if ~isempty(dir_logfile) && (dir_logfile.bytes > 4000000) %it exists and its size is more than 4 Mb
        d=datestr(now);
        movefile(strcat(QMpath,filesep,dir_logfile.name),strcat(QMpath,filesep,'logfile_',d(1:end-9),'.txt')); %rename the old log file with the current date
    end

    QMf_logfile=fopen(strcat(QMpath,filesep,'logfile.txt'),'at'); %create or open the log file

    if QMf_logfile == -1
        errordlg(char({'Error opening the log file.'}),'Initializing error','modal');
        return;
    end
end