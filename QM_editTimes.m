function varargout = QM_editTimes(varargin)
% QM_EDITTIMES MATLAB code for QM_editTimes.fig
%
%   User interface that allows you to load a time table 
%   for the current PET study or save the current time 
%   table into a .txt file. It is executed from QModeling GUI

% Last Modified by GUIDE v2.5 27-Dec-2013 20:27:45
%---------------------------------------------------------------------------------
% QM_editTimes is part of QModeling.
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
% Copyright (C) 2013, 2016 Francisco Javier Lopez Gonzalez, Jose Paredes Pacheco, 
% Karl-Khader Thurnhofer Hemsi

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @QM_editTimes_OpeningFcn, ...
                   'gui_OutputFcn',  @QM_editTimes_OutputFcn, ...
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


% --- Executes just before QM_editTimes is made visible.
function QM_editTimes_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to QM_editTimes (see VARARGIN)

% Choose default command line output for QM_editTimes
handles.output = hObject;

mainGuiInput = find(strcmp(varargin, 'mainhandles'));
handles.mainGUI = varargin{mainGuiInput+1};

% Update handles structure
guidata(hObject, handles);

% Update actual PET times
global QMmaindata
try
    for i=1:(QMmaindata.num_files)
        rows{1,i} = sprintf('Frame %d: from %d sec to %d sec', i, QMmaindata.times(i,1), QMmaindata.times(i,2));
    end
    set(handles.frame_times_Listbox,'String',rows)
catch er
    errordlg(char({'An error ocurred reading times from PET studio';'Please, check it.'}),'Reading error','modal');
end

%Obtain handles using GUIDATA with the caller's handle 
mainhandles = guidata(handles.mainGUI);

%Inicialize the consistency checkbox
check = get(mainhandles.ConsistencyChecks,'Checked');
if strcmp(check,'on')
    set(handles.checkConsistency_Checkbox,'Value',1);
else
    set(handles.checkConsistency_Checkbox,'Value',0);
end

uiwait(handles.edit_time_figure);


% --- Outputs from this function are returned to the command line.
function varargout = QM_editTimes_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

delete(handles.edit_time_figure);


% --- Executes on button press in applyTimes_Pushbutton.
function applyTimes_Pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to applyTimes_Pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global QMmaindata

if ~isempty(get(handles.frame_times_Listbox,'Userdata'))
    table = get(handles.frame_times_Listbox,'Userdata');
    QMmaindata.times = table;
end

uiresume(handles.edit_time_figure);


% --- Executes on selection change in saveTypeFile_Popupmenu.
function saveTypeFile_Popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to saveTypeFile_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns saveTypeFile_Popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from saveTypeFile_Popupmenu


% --- Executes during object creation, after setting all properties.
function saveTypeFile_Popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to saveTypeFile_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in saveTimes_Pushbutton.
function saveTimes_Pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to saveTimes_Pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global QMmaindata

selectformat = get(handles.saveTypeFile_Popupmenu,'Value');
name = get(handles.nameFileTimes_TextBox,'String');

if isempty(name)==1
    warndlg(char({'The save name is empty';'Please, write it.'}),'Reading error','modal');
    return;
end

set(handles.edit_time_figure,'WindowStyle','normal');
directoryname = uigetdir('', 'Select the output directory');
set(handles.edit_time_figure,'WindowStyle','modal');
if(directoryname~=0)
    try
        if isdir(directoryname) & selectformat == 1
            file = strcat(directoryname,filesep,name,'.txt');
        else if isdir(directoryname) & selectformat == 2
            file = strcat(directoryname,filesep,name,'.xls'); 
            else
                file = strcat(directoryname,filesep,name,'.mat');
            end
        end
    catch er
        warndlg(char({'It is not a correct directory or you canceled save';'Please, check it.'}),'Reading error','modal');
        return;
    end

    text = get(handles.frame_times_Listbox, 'String');

    waitbarhandle = waitbar(0,'Saving, please wait...');
    set(findobj(waitbarhandle,'type','patch'),'edgecolor','b','facecolor','b');

    if selectformat == 1 % .txt
        f = fopen(file, 'wt');
        for i=1:(QMmaindata.num_files)
            line = sscanf(char(text(i,1)),['Frame %d: from ' '%d' ' sec to ' '%d' ' sec']);

            waitbar(i/(QMmaindata.num_files),waitbarhandle);
            fprintf(f,'%d\t%d',line(2:3,1));
            fprintf(f,'\n');
        end
        fclose(f);
    else if selectformat == 2
            if exist(file) ~= 0
                delete(file);
            end

            for i=1:(QMmaindata.num_files)
                lines(i,:) = sscanf(char(text(i,1)),['Frame %d: from ' '%d' ' sec to ' '%d' ' sec']);

                waitbar(i/((QMmaindata.num_files)*3),waitbarhandle);
            end
            try
                title = {'FRAME' 'TIEMPO INICIO (s)' 'TIEMPO FIN (s)'};
                xlswrite(file,title, 1, 'A1');
                waitbar(2/3,waitbarhandle);
                xlswrite(file,lines, 1, 'A2');
                waitbar(1,waitbarhandle);
            catch er
                errordlg(char({'An error ocurred creating xls file';'Please, see the help.'}),'Writing error','modal');
            end

        else % .mat
            for i=1:(QMmaindata.num_files)
                line(i,:) = sscanf(char(text(i,1)),['Frame %d: from ' '%d' ' sec to ' '%d' ' sec']);
                waitbar(i/(QMmaindata.num_files),waitbarhandle);
            end

            ttimes = line(:,2:3);
            save(file,'ttimes');
        end
    end
    delete(waitbarhandle);
end

% --- Executes on button press in CancelTimes_Pushbutton.
function CancelTimes_Pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to CancelTimes_Pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.edit_time_figure);

function nameFileTimes_TextBox_Callback(hObject, eventdata, handles)
% hObject    handle to nameFileTimes_TextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nameFileTimes_TextBox as text
%        str2double(get(hObject,'String')) returns contents of nameFileTimes_TextBox as a double


% --- Executes during object creation, after setting all properties.
function nameFileTimes_TextBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nameFileTimes_TextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes when user attempts to close edit_time_figure.
function edit_time_figure_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to edit_time_figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
uiresume(handles.edit_time_figure);


% --- Executes on selection change in loadTypeFile_Popupmenu.
function loadTypeFile_Popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to loadTypeFile_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns loadTypeFile_Popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from loadTypeFile_Popupmenu


% --- Executes during object creation, after setting all properties.
function loadTypeFile_Popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to loadTypeFile_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in frame_times_Listbox.
function frame_times_Listbox_Callback(hObject, eventdata, handles)
% hObject    handle to frame_times_Listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns frame_times_Listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from frame_times_Listbox


% --- Executes during object creation, after setting all properties.
function frame_times_Listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frame_times_Listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in loadTimes_Pushbutton.
function loadTimes_Pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to loadTimes_Pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global QMmaindata

selectformat = get(handles.loadTypeFile_Popupmenu,'Value');

if selectformat == 1 %txt
    time_table_dir = spm_select(1,'txt','Select .txt time table document ');
    if(~isempty(time_table_dir))
        try
            time_table = load(time_table_dir);
            set(handles.frame_times_Listbox,'Userdata',time_table(1:QMmaindata.num_files,:)); %guardo tabla actual

            for i=1:(QMmaindata.num_files)
                rows{1,i} = sprintf('Frame %d: from %d sec to %d sec', i, time_table(i,1), time_table(i,2));
            end
            set(handles.frame_times_Listbox,'String',rows);
        catch er
            errordlg(char({'An error ocurred reading time table from txt document';'Please, check it or see the help.'}),'Reading error','modal');
        end
    end
    
else if selectformat == 2 %xls
    time_table_dir = spm_select(1,'xls','Select .xls time table document ');
    if(~isempty(time_table_dir))

        [time,txt] = xlsread(time_table_dir,-1);
        time_table = [ time(1:size(time(:,1)),1) , time(1:size(time(:,1)),2)];
        set(handles.frame_times_Listbox,'Userdata',time_table); %guardo tabla actual

        rows{1,1} = sprintf('Frame %d: from %d sec to %d sec', 1, time_table(1,1), time_table(1,2));
        for i=2:(QMmaindata.num_files)
            rows{1,i} = sprintf('Frame %d: from %d sec to %d sec', i, time_table(i,1), time_table(i,2));
        end
        
        set(handles.frame_times_Listbox,'String',rows);
    end
    
else %mat
        time_table_dir = spm_select(1,'^*\.mat$','Select .mat time table');
        if(~isempty(time_table_dir))
            try
                struct_table = load(time_table_dir);
                time_table=struct2cell(struct_table); %Porque no puedo acceder directamente al contenido sin conocer el nombre de la variable, tendr�a que pedirlo

                for i=1:(QMmaindata.num_files)
                    time_table_copy(i,:)=[time_table{1,1}(i,1), time_table{1,1}(i,2)];
                    rows{1,i} = sprintf('Frame %d: from %d sec to %d sec', i, time_table{1,1}(i,1), time_table{1,1}(i,2)); 
                end

                set(handles.frame_times_Listbox,'Userdata',time_table_copy); %guardo tabla actual
                set(handles.frame_times_Listbox,'String',rows);
            catch er
                errordlg(char({'An error ocurred reading time table from mat file';'Please, check it or see the help.'}),'Reading error','modal');
            end  
        end
    end
end


% --- Executes on button press in checkConsistency_Checkbox.
function checkConsistency_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to checkConsistency_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkConsistency_Checkbox

global QMmaindata

%Obtain handles using GUIDATA with the caller's handle 
mainhandles = guidata(handles.mainGUI);

check = get(hObject,'Value');
if check == 0
    QMmaindata.flags.consistency = 0;
    set(mainhandles.ConsistencyChecks,'Checked','off');
else
    QMmaindata.flags.consistency = 1;
    set(mainhandles.ConsistencyChecks,'Checked','on');
end
