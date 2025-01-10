function varargout = QM_mapParamLogan(varargin)
% %QM_MAPPARAMLOGAN MATLAB code file for QM_mapParamLogan.fig
%
%   QM_mapParamLogan generates a GUI to choose the parametric
%   images which you want to calculate from the previous preprocessing
%   of Logan Plot model. This code must be executed from QModeling GUI.
%
% Last Modified by GUIDE v2.5 27-Dec-2017 17:29:46
%---------------------------------------------------------------------------------
% QM_mapParamLogan is part of QModeling.
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

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @QM_mapParamLogan_OpeningFcn, ...
                   'gui_OutputFcn',  @QM_mapParamLogan_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before QM_mapParamLogan is made visible.
function QM_mapParamLogan_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for QM_mapParamLogan
handles.output = hObject;

mainGuiInput = find(strcmp(varargin, 'mainhandles'));
handles.mainGUI = varargin{mainGuiInput+1};

% Update handles structure
guidata(hObject, handles);

%UIWAIT makes QM_mapParamLogan wait for user response (see UIRESUME)
uiwait(handles.MapParLoganFigure);


% --- Outputs from this function are returned to the command line.
function varargout = QM_mapParamLogan_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
delete(handles.MapParLoganFigure);

% --- Executes on button press in BPndCheckbox.
function BPndCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to BPndCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of BPndCheckbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
    rest = get(handles.BPndRestrictCheckbox,'Value');
    sav = get(handles.BPndSaveCheckbox,'Value');
    if rest == 1
        set(handles.BPndLowerTextbox,'Enable','on');
        set(handles.BPndUpperTextbox,'Enable','on');
    end
    if sav == 1
        set(handles.BPndImage,'Enable','on');
        set(handles.BPndExtPopupmenu,'Enable','on');
    end
    set(handles.RunPixelWisePushbutton,'Enable','on')
else
    enable = 'off';
    set(handles.BPndLowerTextbox,'Enable','off');
    set(handles.BPndUpperTextbox,'Enable','off');
    set(handles.BPndImage,'Enable','off');
    set(handles.BPndExtPopupmenu,'Enable','off');
    if ~get(handles.InterceptCheckbox,'Value')
        set(handles.RunPixelWisePushbutton,'Enable','off')
    end
end
set(handles.BPndRestrictCheckbox,'Enable',enable);
set(handles.BPndSaveCheckbox,'Enable',enable);



function BPndLowerTextbox_Callback(hObject, eventdata, handles)
% hObject    handle to BPndLowerTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BPndLowerTextbox as text
%        str2double(get(hObject,'String')) returns contents of BPndLowerTextbox as a double
value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('Lower value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function BPndLowerTextbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BPndLowerTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function InterLowerTextbox_Callback(hObject, eventdata, handles)
% hObject    handle to InterLowerTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of InterLowerTextbox as text
%        str2double(get(hObject,'String')) returns contents of InterLowerTextbox as a double
value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('Lower value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function InterLowerTextbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to InterLowerTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function InterImage_Callback(hObject, eventdata, handles)
% hObject    handle to InterImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of InterImage as text
%        str2double(get(hObject,'String')) returns contents of InterImage as a double


% --- Executes during object creation, after setting all properties.
function InterImage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to InterImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function BPndUpperTextbox_Callback(hObject, eventdata, handles)
% hObject    handle to BPndUpperTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BPndUpperTextbox as text
%        str2double(get(hObject,'String')) returns contents of BPndUpperTextbox as a double
value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('Upper value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function BPndUpperTextbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BPndUpperTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function InterUpperTextbox_Callback(hObject, eventdata, handles)
% hObject    handle to InterUpperTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of InterUpperTextbox as text
%        str2double(get(hObject,'String')) returns contents of InterUpperTextbox as a double
value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('Upper value must be a positive number','Warning','modal');
    uicontrol(hObject);
end


% --- Executes during object creation, after setting all properties.
function InterUpperTextbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to InterUpperTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in InterceptCheckbox.
function InterceptCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to InterceptCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of InterceptCheckbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
    rest = get(handles.InterRestrictCheckbox,'Value');
    sav = get(handles.InterSaveCheckbox,'Value');
    if rest == 1
        set(handles.InterLowerTextbox,'Enable','on');
        set(handles.InterUpperTextbox,'Enable','on');
    end
    if sav == 1
        set(handles.InterImage,'Enable','on');
        set(handles.InterExtPopupmenu,'Enable','on');
    end
    
    set(handles.RunPixelWisePushbutton,'Enable','on');

else
    if ~get(handles.BPndCheckbox,'Value')
        set(handles.RunPixelWisePushbutton,'Enable','off')
    end
    enable ='off';
    set(handles.InterLowerTextbox,'Enable',enable);
    set(handles.InterUpperTextbox,'Enable',enable);
    set(handles.InterImage,'Enable',enable);
    set(handles.InterExtPopupmenu,'Enable',enable); 

    
end

set(handles.InterRestrictCheckbox,'Enable',enable);
set(handles.InterSaveCheckbox,'Enable',enable);




% --- Executes on button press in InterSaveCheckbox.
function InterSaveCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to InterSaveCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of InterSaveCheckbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
else
    enable = 'off';
end
set(handles.InterImage,'Enable',enable);
set(handles.InterExtPopupmenu,'Enable',enable);


% --- Executes on button press in BPndSaveCheckbox.
function BPndSaveCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to BPndSaveCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of BPndSaveCheckbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
else
    enable = 'off';
end
set(handles.BPndImage,'Enable',enable);
set(handles.BPndExtPopupmenu,'Enable',enable);


% --- Executes on button press in InterRestrictCheckbox.
function InterRestrictCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to InterRestrictCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of InterRestrictCheckbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
else
    enable = 'off';
end
set(handles.InterLowerTextbox,'Enable',enable);
set(handles.InterUpperTextbox,'Enable',enable);

% --- Executes on button press in BPndRestrictCheckbox.
function BPndRestrictCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to BPndRestrictCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of BPndRestrictCheckbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
else
    enable = 'off';
end
set(handles.BPndLowerTextbox,'Enable',enable);
set(handles.BPndUpperTextbox,'Enable',enable);

% --- Executes on selection change in BPndExtPopupmenu.
function BPndExtPopupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to BPndExtPopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns BPndExtPopupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from BPndExtPopupmenu


% --- Executes during object creation, after setting all properties.
function BPndExtPopupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BPndExtPopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in InterExtPopupmenu.
function InterExtPopupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to InterExtPopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns InterExtPopupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from InterExtPopupmenu


% --- Executes during object creation, after setting all properties.
function InterExtPopupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to InterExtPopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function BPndImage_Callback(hObject, eventdata, handles)
% hObject    handle to BPndImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BPndImage as text
%        str2double(get(hObject,'String')) returns contents of BPndImage as a double


% --- Executes during object creation, after setting all properties.
function BPndImage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BPndImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in RunPixelWisePushbutton.
function RunPixelWisePushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to RunPixelWisePushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global QMmaindata
global QMlastPreprocess
global QMpath
global QMf_logfile


BPndLow = 0;
BPndUpp = 0;
interlow = 0;
interupp = 0;

%Obtain handles using GUIDATA with the caller's handle 
mainhandles = guidata(handles.mainGUI);

%Delete the last parametric images created in temp file
mainhandles.deleteTempImages(mainhandles);

outputDir = get(handles.OutputDirTextbox,'String');
if ~isdir(outputDir)    %Check directory
    warndlg('Output directory does not exist','Pixel-wise error','modal');
    uicontrol(handles.OutputDirTextbox);
    return
end

if get(handles.BPndSaveCheckbox,'Value') == 1 || get(handles.InterSaveCheckbox,'Value') == 1
    [s dirAtrib ex] = fileattrib(outputDir);
    if s == 1 && dirAtrib.UserWrite == 0  %Check if directory is writable
        warndlg('Output directory does not writable','Pixel-wise error','modal');
        uicontrol(handles.outputDirTextbox);
        return
    end
end

BPndcheck = get(handles.BPndCheckbox,'Value');
intercheck = get(handles.InterceptCheckbox,'Value');
BPndrest = get(handles.BPndRestrictCheckbox,'Value');
interrest = get(handles.InterRestrictCheckbox,'Value');

%Checking errors
if BPndcheck == 1 && BPndrest == 1
    
    BPndLow = str2double(get(handles.BPndLowerTextbox,'String'));
    if isnan(BPndLow)
        warndlg('Lower restrict must be a positive number','Pixel-wise error','modal');
        uicontrol(handles.BPndLowerTextbox);
        return
    end
    
    BPndUpp = str2double(get(handles.BPndUpperTextbox,'String'));
    if isnan(BPndUpp)
        warndlg('Upper restrict must be a positive number','Pixel-wise error','modal');
        uicontrol(handles.BPndUpperTextbox);
        return
    end
    
    if BPndLow > BPndUpp
        warndlg('Upper restrict must be bigger than Lower restrict','Pixel-wise error','modal');
        uicontrol(handles.BPndUpperTextbox);
        return
    end
    
end

if intercheck == 1 && interrest == 1
    
    interlow = str2double(get(handles.InterLowerTextbox,'String'));
    if isnan(interlow)
        warndlg('Lower restrict must be a positive number','Pixel-wise error','modal');
        uicontrol(handles.IntercLowerTextbox);
        return
    end
    
    interupp = str2double(get(handles.InterUpperTextbox,'String'));
    if isnan(interupp)
        warndlg('Upper restrict must be a positive number','Pixel-wise error','modal');
        uicontrol(handles.IntertUpperTextbox);
        return
    end
    
    if interlow > interupp
        warndlg('Upper restrict must be bigger than Lower restrict','Pixel-wise error','modal');
        uicontrol(handles.InterUpperTextbox);
        return
    end
    
end

waitbarhandle = waitbar(0,'Preparing parametric images, please wait...');
ver=version('-release');
if (str2num(ver(1:4))<2014) || isequal(ver,'2014a')
    set(findobj(waitbarhandle,'type','patch'),'edgecolor','b','facecolor','b');
else
    wbc = allchild(waitbarhandle);     %you need to get at a hidden child
    wbc(1).JavaPeer.setForeground( wbc(1).JavaPeer.getBackground.BLUE)
    wbc(1).JavaPeer.setStringPainted(true)
end

%Update the logfile
aux_date=datestr(now);
fprintf(QMf_logfile,strcat('\n-Logan Plot-\n'));
fprintf(QMf_logfile,strcat(aux_date(end-8:end),' Generating parametric images (t*=',num2str(QMlastPreprocess.LoganPlot.t),', threshold=',num2str(QMlastPreprocess.LoganPlot.threshold),')\n'));

waitbar(1/5,waitbarhandle);

%Image calculation
images = QM_LoganPixelCalc(QMlastPreprocess.LoganPlot.numX(QMlastPreprocess.LoganPlot.l:end),QMmaindata.Cpet,...
            QMmaindata.times,QMlastPreprocess.LoganPlot.threshold,[BPndrest BPndLow BPndUpp ; interrest interlow interupp],...
            QMlastPreprocess.LoganPlot.l,...
            waitbarhandle);
waitbar(4/5,waitbarhandle);

%Update the logfile
aux_date=datestr(now);
fprintf(QMf_logfile,strcat(aux_date(end-8:end),' Parametric images generated \n'));        
 
studyPaths = get(mainhandles.study_Textbox,'UserData');
num_files = size(studyPaths,1);
if num_files == 1     
    multifile = studyPaths(1,1:max(strfind(studyPaths,','))-1);      % To raise ',1'
    aux=spm_vol(multifile);
    Vtemp = aux(1);
else
    file = studyPaths(1,:);
    Vtemp = spm_vol(file);
end

imgTempPaths = {};
imgTempFilenames = {};
indTemp=1;
if BPndcheck == 1
    BPndfile = get(handles.BPndImage,'String');
    
    if get(handles.BPndExtPopupmenu,'Value') == 1
        BPndFileExt = '.img';
    else
        BPndFileExt = '.nii';
    end

    VBPnd = struct('fname',strcat(BPndfile,BPndFileExt),'dim',Vtemp.dim,'mat',Vtemp.mat,...
                'dt',[spm_type('float64') spm_platform('bigend')],'n',[1 1]);

    % Make image
    val = zeros(VBPnd.dim);
    val(:) = images(1,:);
    

    if get(handles.BPndSaveCheckbox,'Value') == 1
        VBPnd.fname = strcat(outputDir,filesep,BPndfile,BPndFileExt);
        spm_write_vol(VBPnd,val);
        fprintf(QMf_logfile,strcat(['  -BPnd image saved in ',outputDir,'\n']));
    end
    
    VBPnd.fname = strcat(QMpath,filesep,'temp',filesep,BPndfile,BPndFileExt);
    spm_write_vol(VBPnd,val);
    
    imgTempPaths{1,indTemp} = VBPnd.fname;
    imgTempFilenames{1,indTemp} = 'BPnd (binding potential)'; 
    indTemp = indTemp + 1;
    
end        
  waitbar(5/5,waitbarhandle);


if intercheck == 1
    interfile = get(handles.InterImage,'String');

    if get(handles.InterExtPopupmenu,'Value') == 1
        interFileExt = '.img';
    else
        interFileExt = '.nii';
    end

    VI = struct('fname',strcat(interfile,interFileExt),'dim',Vtemp.dim,'mat',Vtemp.mat,...
                'dt',[spm_type('float64') spm_platform('bigend')],'n',[1 1]);
    
    % Make image
    val2 = zeros(VI.dim);
    val2(:) = images(2,:);
    
    if get(handles.InterSaveCheckbox,'Value') == 1
        VI.fname = strcat(outputDir,filesep,interfile,interFileExt);
        
        spm_write_vol(VI,val2);
        fprintf(QMf_logfile,strcat(['  -Intercept image saved in ',outputDir,'\n']));
    end
    
    VI.fname = strcat(QMpath,filesep,'temp',filesep,interfile,interFileExt);
    spm_write_vol(VI,val2);
    
    imgTempPaths{1,indTemp} = VI.fname;
    imgTempFilenames{1,indTemp} = 'Intercept (intercept of the linear regression)'; 
    

end


set(mainhandles.selectImage_Popupmenu,'UserData',imgTempPaths);
set(mainhandles.selectImage_Popupmenu,'String',imgTempFilenames);

set(mainhandles.setMapParam_Pushbutton,'UserData',1);

delete(waitbarhandle)

uiresume(handles.MapParLoganFigure);



% --- Executes on button press in SetDefaultTextbox.
function SetDefaultTextbox_Callback(hObject, eventdata, handles)
% hObject    handle to SetDefaultTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.BPndCheckbox,'Value',1);
set(handles.BPndRestrictCheckbox,'Value',1);
set(handles.BPndRestrictCheckbox,'Enable','on');
set(handles.BPndLowerTextbox,'String','0.0');
set(handles.BPndLowerTextbox,'Enable','on');
set(handles.BPndUpperTextbox,'String','60.0');
set(handles.BPndUpperTextbox,'Enable','on');
set(handles.BPndSaveCheckbox,'Value',0);
set(handles.BPndSaveCheckbox,'Enable','on');
set(handles.BPndImage,'String','BPnd_image');
set(handles.BPndImage,'Enable','off');
set(handles.BPndExtPopupmenu,'Value',1);
set(handles.BPndExtPopupmenu,'Enable','off');

set(handles.InterceptCheckbox,'Value',1);
set(handles.InterRestrictCheckbox,'Value',0);
set(handles.InterRestrictCheckbox,'Enable','on');
set(handles.InterLowerTextbox,'String','0.0');
set(handles.InterLowerTextbox,'Enable','off');
set(handles.InterUpperTextbox,'String','1.0');
set(handles.InterUpperTextbox,'Enable','off');
set(handles.InterSaveCheckbox,'Value',0);
set(handles.InterSaveCheckbox,'Enable','on');
set(handles.InterImage,'String','Intercept_image');
set(handles.InterImage,'Enable','off');
set(handles.InterExtPopupmenu,'Value',1);
set(handles.InterExtPopupmenu,'Enable','off');

set(handles.RunPixelWisePushbutton,'Enable','on')

dir = pwd;
set(handles.OutputDirTextbox,'String',dir)




function OutputDirTextbox_Callback(hObject, eventdata, handles)
% hObject    handle to OutputDirTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of OutputDirTextbox as text
%        str2double(get(hObject,'String')) returns contents of OutputDirTextbox as a double
dir = get(hObject,'String');
if ~isdir(dir)
    warndlg('This directory does not exist','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function OutputDirTextbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to OutputDirTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
dir = pwd;
set(hObject,'String',dir)


% --- Executes on button press in SelectDirTextbox.
function SelectDirTextbox_Callback(hObject, eventdata, handles)
% hObject    handle to SelectDirTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.MapParLoganFigure,'WindowStyle','normal');
directoryname = uigetdir('', 'Select the output directory');
set(handles.MapParLoganFigure,'WindowStyle','modal');

if directoryname ~= 0
     set(handles.OutputDirTextbox,'String',directoryname);
     
end


% --- Executes when user attempts to close MapParLoganFigure.
function MapParLoganFigure_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to MapParLoganFigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
% delete(hObject);
uiresume(handles.MapParLoganFigure);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over SetDefaultTextbox.
function SetDefaultTextbox_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to SetDefaultTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
