function varargout = QM_mapParamPatlak(varargin)
% QM_MAPPARAMPATLAK MATLAB code for QM_mapParamPatlak.fig
%      
%   QM_mapParamPatlak generates a GUI to choose the parametric
%   images which you want to calculate from the previous preprocessing
%   of Patlak model. This code must be executed from QModeling GUI.
%
% Last Modified by GUIDE v2.5 18-Oct-2013 17:53:44
%---------------------------------------------------------------------------------
% QM_mapParamPatlak is part of QModeling.
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

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @QM_mapParamPatlak_OpeningFcn, ...
                   'gui_OutputFcn',  @QM_mapParamPatlak_OutputFcn, ...
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


% --- Executes just before QM_mapParamPatlak is made visible.
function QM_mapParamPatlak_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to QM_mapParamPatlak (see VARARGIN)

% Choose default command line output for QM_mapParamPatlak
handles.output = hObject;

mainGuiInput = find(strcmp(varargin, 'mainhandles'));
handles.mainGUI = varargin{mainGuiInput+1};

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes QM_mapParamPatlak wait for user response (see UIRESUME)
uiwait(handles.MapParPatlak_figure);


% --- Outputs from this function are returned to the command line.
function varargout = QM_mapParamPatlak_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

delete(handles.MapParPatlak_figure);


% --- Executes on button press in runPixelCalc_Pushbutton.
function runPixelCalc_Pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to runPixelCalc_Pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global QMmaindata
global QMlastPreprocess
global QMpath
global QMf_logfile

Klow = 0;
Kupp = 0;
interlow = 0;
interupp = 0;

%Obtain handles using GUIDATA with the caller's handle 
mainhandles = guidata(handles.mainGUI);

%Delete the last parametric images created in temp file
mainhandles.deleteTempImages(mainhandles);

outputDir = get(handles.outputDir_Textbox,'String');
if ~isdir(outputDir)    %Check existing directory
    warndlg('Output directory does not exist','Pixel-wise error','modal');
    uicontrol(handles.outputDir_Textbox);
    return
end

if get(handles.saveK_Checkbox,'Value') == 1 || get(handles.saveIntercept_Checkbox,'Value') == 1
    [s dirAtrib ex] = fileattrib(outputDir);
    if s == 1 && dirAtrib.UserWrite == 0  %Check if directory is writable
        warndlg('Output directory does not writable','Pixel-wise error','modal');
        uicontrol(handles.outputDir_Textbox);
        return
    end
end


Kcheck = get(handles.K_Checkbox,'Value');
intercheck = get(handles.intercept_Checkbox,'Value');
Krest = get(handles.Krestrict_Checkbox,'Value');
interrest = get(handles.interceptRest_Checkbox,'Value');

%Checking errors
if Kcheck == 1 && Krest == 1
    
    Klow = str2double(get(handles.Klower_Textbox,'String'));
    if isnan(Klow)
        warndlg('Lower restriction must be a positive number','Pixel-wise error','modal');
        uicontrol(handles.Klower_Textbox);
        return
    end
    
    Kupp = str2double(get(handles.Kupper_Textbox,'String'));
    if isnan(Kupp)
        warndlg('Upper restriction must be a positive number','Pixel-wise error','modal');
        uicontrol(handles.Kupper_Textbox);
        return
    end
    
    if Klow > Kupp
        warndlg('Upper restriction must be bigger than Lower restriction','Pixel-wise error','modal');
        uicontrol(handles.Kupper_Textbox);
        return
    end
    
end

if intercheck == 1 && interrest == 1
    
    interlow = str2double(get(handles.interceptLower_Textbox,'String'));
    if isnan(interlow)
        warndlg('Lower restriction must be a positive number','Pixel-wise error','modal');
        uicontrol(handles.interceptLower_Textbox);
        return
    end
    
    interupp = str2double(get(handles.interceptUpper_Textbox,'String'));
    if isnan(interupp)
        warndlg('Upper restriction must be a positive number','Pixel-wise error','modal');
        uicontrol(handles.interceptUpper_Textbox);
        return
    end
    
    if interlow > interupp
        warndlg('Upper restriction must be bigger than Lower restriction','Pixel-wise error','modal');
        uicontrol(handles.interceptUpper_Textbox);
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
fprintf(QMf_logfile,strcat('\n-Patlak Reference-\n'));
fprintf(QMf_logfile,strcat(aux_date(end-8:end),' Generating parametric images (t*=',num2str(QMlastPreprocess.PatlakRef.t),', threshold=',num2str(QMlastPreprocess.PatlakRef.threshold),')\n'));

%Image calculation
images = QM_patlakPixelCalc(QMmaindata.Cpet,QMmaindata.Crois{2,QMlastPreprocess.PatlakRef.idCr},...
            find(QMmaindata.times(:,1)/60 == QMlastPreprocess.PatlakRef.t),QMlastPreprocess.PatlakRef.plots(:,1),...
            QMmaindata.times,QMlastPreprocess.PatlakRef.threshold,[Krest Klow Kupp ; interrest interlow interupp],...
            waitbarhandle);

%Update the logfile
aux_date=datestr(now);
fprintf(QMf_logfile,strcat(aux_date(end-8:end),' Parametric images generated \n'));        
        
waitbar(3/5,waitbarhandle)

%Reference header to make the new images
% studyPaths = get(mainhandles.mask_Textbox,'UserData');
% Vtemp = spm_vol(studyPaths{1});

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

if Kcheck == 1
    Kfile = get(handles.Kfilename_Textbox,'String');
    
    if get(handles.KfileExt_Popupmenu,'Value') == 1
        KfileExt = '.img';
    else
        KfileExt = '.nii';
    end

    VK = struct('fname',strcat(Kfile,KfileExt),'dim',Vtemp.dim,'mat',Vtemp.mat,...
                'dt',[spm_type('float64') spm_platform('bigend')],'n',[1 1]);

    % Make image
    val = zeros(VK.dim);
    val(:) = images(1,:)*60;
    

    if get(handles.saveK_Checkbox,'Value') == 1
        VK.fname = strcat(outputDir,filesep,Kfile,KfileExt);
        spm_write_vol(VK,val);
        fprintf(QMf_logfile,strcat(['  -K image saved in ',outputDir,'\n']));
    end
    
    VK.fname = strcat(QMpath,filesep,'temp',filesep,Kfile,KfileExt);
    spm_write_vol(VK,val);
    
    imgTempPaths{1,indTemp} = VK.fname;
    imgTempFilenames{1,indTemp} = 'K (slope of the linear regression)'; %strcat(Kfile,KfileExt);
    indTemp = indTemp + 1;
    
end        

waitbar(4/5,waitbarhandle)

if intercheck == 1
    interfile = get(handles.interceptFilename_Textbox,'String');

    if get(handles.interceptFileExt_Popupmenu,'Value') == 1
        interFileExt = '.img';
    else
        interFileExt = '.nii';
    end

    VI = struct('fname',strcat(interfile,interFileExt),'dim',Vtemp.dim,'mat',Vtemp.mat,...
                'dt',[spm_type('float64') spm_platform('bigend')],'n',[1 1]);
    
    % Make image
    val = zeros(VI.dim);
    val(:) = images(2,:);
    
    if get(handles.saveIntercept_Checkbox,'Value') == 1
        VI.fname = strcat(outputDir,filesep,interfile,interFileExt);
        spm_write_vol(VI,val);
        fprintf(QMf_logfile,strcat(['  -BPnd image saved in ',outputDir,'\n']));
    end
    
    VI.fname = strcat(QMpath,filesep,'temp',filesep,interfile,interFileExt);
    spm_write_vol(VI,val);
    
    imgTempPaths{1,indTemp} = VI.fname;
    imgTempFilenames{1,indTemp} = 'Intercept (intercept of the linear regression)'; %strcat(interfile,interFileExt);
    
end

waitbar(5/5,waitbarhandle)

set(mainhandles.selectImage_Popupmenu,'UserData',imgTempPaths);
set(mainhandles.selectImage_Popupmenu,'String',imgTempFilenames);

set(mainhandles.setMapParam_Pushbutton,'UserData',1);

delete(waitbarhandle)

uiresume(handles.MapParPatlak_figure);


% --- Executes on button press in setDefault_Pushbutton.
function setDefault_Pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to setDefault_Pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.K_Checkbox,'Value',1);
set(handles.Krestrict_Checkbox,'Value',0);
set(handles.Krestrict_Checkbox,'Enable','on');
set(handles.Klower_Textbox,'String','0.0');
set(handles.Klower_Textbox,'Enable','off');
set(handles.Kupper_Textbox,'String','1.0');
set(handles.Kupper_Textbox,'Enable','off');
set(handles.saveK_Checkbox,'Value',1);
set(handles.saveK_Checkbox,'Enable','on');
set(handles.Kfilename_Textbox,'String','K_image');
set(handles.KfileExt_Popupmenu,'Value',1);
set(handles.Kfilename_Textbox,'Enable','on');
set(handles.KfileExt_Popupmenu,'Enable','on');

set(handles.intercept_Checkbox,'Value',0);
set(handles.interceptRest_Checkbox,'Value',0);
set(handles.interceptRest_Checkbox,'Enable','off');
set(handles.interceptLower_Textbox,'String','0.0');
set(handles.interceptLower_Textbox,'Enable','off');
set(handles.interceptUpper_Textbox,'String','1.0');
set(handles.interceptUpper_Textbox,'Enable','off');
set(handles.saveIntercept_Checkbox,'Value',1);
set(handles.saveIntercept_Checkbox,'Enable','off');
set(handles.interceptFilename_Textbox,'String','Intercept_image');
set(handles.interceptFileExt_Popupmenu,'Value',1);
set(handles.interceptFilename_Textbox,'Enable','off');
set(handles.interceptFileExt_Popupmenu,'Enable','off');

dir = pwd;
set(handles.outputDir_Textbox,'String',dir)


% --- Executes on button press in K_Checkbox.
function K_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to K_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of K_Checkbox

checked = get(hObject,'Value');
if checked
    enable = 'on';
    rest = get(handles.Krestrict_Checkbox,'Value');
    sav = get(handles.saveK_Checkbox,'Value');
    if rest == 1
        set(handles.Klower_Textbox,'Enable','on');
        set(handles.Kupper_Textbox,'Enable','on');
    end
    if sav == 1
        set(handles.Kfilename_Textbox,'Enable','on');
        set(handles.KfileExt_Popupmenu,'Enable','on');
    end
    set(handles.runPixelCalc_Pushbutton,'Enable','on')
else
    enable = 'off';
    set(handles.Klower_Textbox,'Enable','off');
    set(handles.Kupper_Textbox,'Enable','off');
    set(handles.Kfilename_Textbox,'Enable','off');
    set(handles.KfileExt_Popupmenu,'Enable','off');
    if ~get(handles.intercept_Checkbox,'Value')
        set(handles.runPixelCalc_Pushbutton,'Enable','off')
    end
end
set(handles.Krestrict_Checkbox,'Enable',enable);
set(handles.saveK_Checkbox,'Enable',enable);




% --- Executes on button press in intercept_Checkbox.
function intercept_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to intercept_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of intercept_Checkbox

checked = get(hObject,'Value');
if checked
    enable = 'on';
    rest = get(handles.interceptRest_Checkbox,'Value');
    sav = get(handles.saveIntercept_Checkbox,'Value');
    if rest == 1
        set(handles.interceptLower_Textbox,'Enable','on');
        set(handles.interceptUpper_Textbox,'Enable','on');
    end
    if sav == 1
        set(handles.interceptFilename_Textbox,'Enable','on');
        set(handles.interceptFileExt_Popupmenu,'Enable','on');
    end
    set(handles.runPixelCalc_Pushbutton,'Enable','on')
else
    enable = 'off';
    set(handles.interceptLower_Textbox,'Enable','off');
    set(handles.interceptUpper_Textbox,'Enable','off');
    set(handles.interceptFilename_Textbox,'Enable','off');
    set(handles.interceptFileExt_Popupmenu,'Enable','off');
    if ~get(handles.K_Checkbox,'Value')
        set(handles.runPixelCalc_Pushbutton,'Enable','off')
    end
end
set(handles.interceptRest_Checkbox,'Enable',enable);
set(handles.saveIntercept_Checkbox,'Enable',enable);





% --- Executes on button press in Krestrict_Checkbox.
function Krestrict_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to Krestrict_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Krestrict_Checkbox

checked = get(hObject,'Value');
if checked
    enable = 'on';
else
    enable = 'off';
end
set(handles.Klower_Textbox,'Enable',enable);
set(handles.Kupper_Textbox,'Enable',enable);


function Klower_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to Klower_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Klower_Textbox as text
%        str2double(get(hObject,'String')) returns contents of Klower_Textbox as a double
value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('Lower value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function Klower_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Klower_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Kupper_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to Kupper_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Kupper_Textbox as text
%        str2double(get(hObject,'String')) returns contents of Kupper_Textbox as a double
value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('Upper value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function Kupper_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Kupper_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in interceptRest_Checkbox.
function interceptRest_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to interceptRest_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of interceptRest_Checkbox

checked = get(hObject,'Value');
if checked
    enable = 'on';
else
    enable = 'off';
end
set(handles.interceptLower_Textbox,'Enable',enable);
set(handles.interceptUpper_Textbox,'Enable',enable);


function interceptLower_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to interceptLower_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of interceptLower_Textbox as text
%        str2double(get(hObject,'String')) returns contents of interceptLower_Textbox as a double
value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('Lower value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function interceptLower_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to interceptLower_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function interceptUpper_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to interceptUpper_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of interceptUpper_Textbox as text
%        str2double(get(hObject,'String')) returns contents of interceptUpper_Textbox as a double
value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('Upper value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function interceptUpper_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to interceptUpper_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in saveK_Checkbox.
function saveK_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to saveK_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of saveK_Checkbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
else
    enable = 'off';
end
set(handles.Kfilename_Textbox,'Enable',enable);
set(handles.KfileExt_Popupmenu,'Enable',enable);


function Kfilename_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to Kfilename_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Kfilename_Textbox as text
%        str2double(get(hObject,'String')) returns contents of Kfilename_Textbox as a double


% --- Executes during object creation, after setting all properties.
function Kfilename_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Kfilename_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in saveIntercept_Checkbox.
function saveIntercept_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to saveIntercept_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of saveIntercept_Checkbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
else
    enable = 'off';
end
set(handles.interceptFilename_Textbox,'Enable',enable);
set(handles.interceptFileExt_Popupmenu,'Enable',enable);

function interceptFilename_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to interceptFilename_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of interceptFilename_Textbox as text
%        str2double(get(hObject,'String')) returns contents of interceptFilename_Textbox as a double


% --- Executes during object creation, after setting all properties.
function interceptFilename_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to interceptFilename_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function outputDir_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to outputDir_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of outputDir_Textbox as text
%        str2double(get(hObject,'String')) returns contents of outputDir_Textbox as a double
dir = get(hObject,'String');
if ~isdir(dir)
    warndlg('This directory does not exist','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function outputDir_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to outputDir_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

dir = pwd;
set(hObject,'String',dir)


% --- Executes on button press in selectDir_Pushbutton.
function selectDir_Pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to selectDir_Pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.MapParPatlak_figure,'WindowStyle','normal');
directoryname = uigetdir('', 'Select the output directory');
set(handles.MapParPatlak_figure,'WindowStyle','modal');

if directoryname ~= 0
     set(handles.outputDir_Textbox,'String',directoryname);
end


% --- Executes on selection change in KfileExt_Popupmenu.
function KfileExt_Popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to KfileExt_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns KfileExt_Popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from KfileExt_Popupmenu


% --- Executes during object creation, after setting all properties.
function KfileExt_Popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to KfileExt_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in interceptFileExt_Popupmenu.
function interceptFileExt_Popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to interceptFileExt_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns interceptFileExt_Popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from interceptFileExt_Popupmenu


% --- Executes during object creation, after setting all properties.
function interceptFileExt_Popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to interceptFileExt_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close MapParPatlak_figure.
function MapParPatlak_figure_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to MapParPatlak_figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
uiresume(handles.MapParPatlak_figure);
