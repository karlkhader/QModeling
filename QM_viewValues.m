function varargout = QM_viewValues(varargin)
% QM_VIEWVALUES MATLAB code for QM_viewValues.fig
%
%      User interface which allows you to see the values of the current
%      TACs and save them into a txt or .mat file. It is executed from QModeling GUI

% Last Modified by GUIDE v2.5 30-Mar-2014 20:10:14
%---------------------------------------------------------------------------------
% QM_viewValues is part of QModeling.
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
% Copyright (C) 2014, 2016 Francisco Javier Lopez Gonzalez, Jose Paredes Pacheco, 
% Karl-Khader Thurnhofer Hemsi

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @QM_viewValues_OpeningFcn, ...
                   'gui_OutputFcn',  @QM_viewValues_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before QM_viewValues is made visible.
function QM_viewValues_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to QM_viewValues (see VARARGIN)

% Choose default command line output for QM_viewValues
handles.output = hObject;

%Save mainGUI handle
mainGuiInput = find(strcmp(varargin, 'mainhandles'));
handles.mainGUI = varargin{mainGuiInput+1};

% Update handles structure
guidata(hObject, handles);

%Showing results in text area
fillTextValues(handles)

% UIWAIT makes QM_viewValues wait for user response (see UIRESUME)
uiwait(handles.QM_viewValues_figure);


% --- Outputs from this function are returned to the command line.
function varargout = QM_viewValues_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

delete(hObject);


% --- Executes on button press in save_Pushbutton.
function save_Pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to save_Pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% strvalues = get(handles.values_Text,'String');

% if ~isempty(strvalues)
    
    [filename, pathname, filterindex] = uiputfile( ...
           {'*.txt','Text file (*.txt)'; ...
           '*.mat','MAT-files (*.mat)'},...
            'Save values as');

    if ~(isequal(filename,0) || isequal(pathname,0))

        waitbarhandle = waitbar(0,'Saving, please wait...');
        ver=version('-release');
        if (str2num(ver(1:4))<2014) || isequal(ver,'2014a')
            set(findobj(waitbarhandle,'type','patch'),'edgecolor','b','facecolor','b');
        else
            wbc = allchild(waitbarhandle);     %you need to get at a hidden child
            wbc(1).JavaPeer.setForeground( wbc(1).JavaPeer.getBackground.BLUE )
            wbc(1).JavaPeer.setStringPainted(true)
        end
        
        names=get(handles.tac_table,'ColumnName')'; %1 x numberOfColumns
        values=get(handles.tac_table,'Data'); %numberOfFrames x numberOfColumns
        nC = length(names);
        nF = length(values);
        
        if filterindex == 1 
            f = fopen(strcat(pathname,filename), 'wt');
            for i=1:nC
                fprintf(f,'%s\t',names{1,i});
            end
            fprintf(f,'\n');
            for j=1:nF
                for i=1:nC                
                    fprintf(f,'%3.10f\t',values{j,i});
                end
                fprintf(f,'\n');
                waitbar(j/nF,waitbarhandle);
            end
            fclose(f);
        else
            myfile = fullfile(pathname,filename);
            for j=1:nF
                for i=1:nC
                    tac_values(j,i)=sscanf(num2str(values{j,i}),'%f');
                end
                waitbar(j/nF,waitbarhandle);
            end
            save(myfile,'tac_values');
        end
        
                
                
%                 if filterindex == 1       
%             f = fopen(strcat(pathname,filename), 'wt');
%             for i = 1:n
%                 fprintf(f,'%s\n',strvalues{i,1});
%                 waitbar(i/n,waitbarhandle);
%             end
%             fclose(f);
%         else 
%             % Create a MAT-file
%            myfile = fullfile(pathname,filename);
%            %matObj = matfile(myfile,'Writable',true);
%             
%            for i = 1:n-3
%                column = sscanf(char(strvalues{i+2,1}),'%f');
%                values(i,:) = column';
%            end
%            % Save into a variable in the file
%            save(myfile,'values');
%            %matObj.savedVar = values;
%         end
        
        delete(waitbarhandle);

    end
    
% end



% --- Executes when user attempts to close QM_viewValues_figure.
function QM_viewValues_figure_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to QM_viewValues_figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
uiresume(handles.QM_viewValues_figure);


% --- Executes during object creation, after setting all properties.
function values_Text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to values_Text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function fillTextValues(handles)

%Obtain handles using GUIDATA with the caller's handle 
mainhandles = guidata(handles.mainGUI);

plots=get(mainhandles.axes,'UserData');
if ~isempty(plots)
    switch plots{1,1}
        case 'TACs' 
            if size(plots{1,2},2)==3 %Only one TAC plotted
                set(handles.tac_table,'Data',plots{1,2});
                set(handles.tac_table,'ColumnName',{'Frame Nº','Time (minutes)','TAC 1 (Orginal Units)'});                
            else  %Two TACs plotted
                set(handles.tac_table,'Data',plots{1,2});
                set(handles.tac_table,'ColumnName',{'Frame Nº','Time (minutes)','TAC 1 (Orginal Units)','TAC 2 (Original Units)'});
            end
        case 'SRTM*'
            set(handles.tac_table,'Data',plots{1,2});
            set(handles.tac_table,'ColumnName',{'Frame Nº','Time (minutes)','TAC 1 (Orginal Units)','Fit to TAC 1','TAC 2 (Original Units)'});
        case 'Patlak'
            set(handles.tac_table,'Data',plots{1,2});
            set(handles.tac_table,'ColumnName',{'Frame Nº','int(Cr)/Cr','Ct/Cr','Patlak fit (linear regression)'});
        case 'LoganPlot'
            set(handles.tac_table,'Data',plots{1,2});
            set(handles.tac_table,'ColumnName',{'Frame Nº',strcat('[int(Ct)+Ct/k2',char(39),']/Ct'),'int(Ct)/Ct','Logan fit (linear regression)'});
        case 'TwoTCM*'            
            set(handles.tac_table,'Data',plots{1,2});
            set(handles.tac_table,'ColumnName',{'Frame Nº','Time (minutes)','TAC 1 (Orginal Units)','Fit to TAC 1','TAC 2 (Original Units)'});
        case 'TwoTCM'            
            set(handles.tac_table,'Data',plots{1,2});
            set(handles.tac_table,'ColumnName',{'Frame Nº','Time (minutes)','TAC 1 (Orginal Units)','Fit to TAC 1','TAC 2 (Original Units)','Whole blood activity'});
    end
    
end




% --------------------------------------------------------------------
function Zoom_Callback(hObject, eventdata, handles)
% hObject    handle to Zoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Increase_Callback(hObject, eventdata, handles)
% hObject    handle to Increase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.values_Text,'FontSize',get(handles.values_Text,'FontSize')+1.0);


% --------------------------------------------------------------------
function Decrease_Callback(hObject, eventdata, handles)
% hObject    handle to Decrease (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.values_Text,'FontSize',get(handles.values_Text,'FontSize')-1.0);
