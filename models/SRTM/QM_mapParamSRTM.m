function varargout = QM_mapParamSRTM(varargin)
% QM_MAPPARAMSRTM MATLAB code for QM_mapParamSRTM.fig
%    
%   QM_mapParamSRTM generates a GUI to choose the parametric
%   images which you want to calculate from the previous preprocessing
%   of SRTM model. This code must be executed from QModeling GUI.
%
% Last Modified by GUIDE v2.5 28-Oct-2013 19:03:03
%---------------------------------------------------------------------------------
% QM_mapParamSRTM is part of QModeling.
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
                   'gui_OpeningFcn', @QM_mapParamSRTM_OpeningFcn, ...
                   'gui_OutputFcn',  @QM_mapParamSRTM_OutputFcn, ...
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


% --- Executes just before QM_mapParamSRTM is made visible.
function QM_mapParamSRTM_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to QM_mapParamSRTM (see VARARGIN)

% Choose default command line output for QM_mapParamSRTM
handles.output = hObject;

mainGuiInput = find(strcmp(varargin, 'mainhandles'));
handles.mainGUI = varargin{mainGuiInput+1};

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes QM_mapParamSRTM wait for user response (see UIRESUME)
uiwait(handles.MapParSRTM_figure);


% --- Outputs from this function are returned to the command line.
function varargout = QM_mapParamSRTM_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

delete(handles.MapParSRTM_figure);


% --- Executes on button press in runPixelCalc_Pushbutton.
function runPixelCalc_Pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to runPixelCalc_Pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global QMmaindata
global QMlastPreprocess
global QMpath

global QMf_logfile

waitbarhandle = waitbar(0,'Creating parametric images, please wait...');
ver=version('-release');
if (str2num(ver(1:4))<2014) || isequal(ver,'2014a')
    set(findobj(waitbarhandle,'type','patch'),'edgecolor','b','facecolor','b');
else
    wbc = allchild(waitbarhandle);     %you need to get at a hidden child
    wbc(1).JavaPeer.setForeground( wbc(1).JavaPeer.getBackground.BLUE)
    wbc(1).JavaPeer.setStringPainted(true)
end

%Obtain handles using GUIDATA with the caller's handle 
mainhandles = guidata(handles.mainGUI);

%Delete the last parametric images created in temp file
mainhandles.deleteTempImages(mainhandles);

outputDir = get(handles.outputDir_Textbox,'String');
if ~isdir(outputDir)
    warndlg('Output directory does not exist','Pixel-wise error','modal');
    uicontrol(handles.outputDir_Textbox);
    return
end

%Get checked options
BPndcheck = get(handles.BPnd_Checkbox,'Value');
k2acheck = get(handles.k2a_Checkbox,'Value');
k2check = get(handles.k2_Checkbox,'Value');
R1check = get(handles.R1_Checkbox,'Value');

BPndrest = get(handles.BPndRestrict_Checkbox,'Value');
k2arest = get(handles.k2aRestrict_Checkbox,'Value');
k2rest = get(handles.k2Restrict_Checkbox,'Value');
R1rest = get(handles.R1Restrict_Checkbox,'Value');

BPndsave = get(handles.saveBPnd_Checkbox,'Value');
k2asave = get(handles.savek2a_Checkbox,'Value');
k2save = get(handles.savek2_Checkbox,'Value');
R1save = get(handles.saveR1_Checkbox,'Value');

BPndlow = 0;
BPndupp = 0;
k2alow = 0;
k2aupp = 0;
k2low = 0;
k2upp = 0;
R1low = 0;
R1upp = 0;

%Verify conditions
if BPndcheck == 1 && BPndrest == 1 
    BPndlow = str2double(get(handles.BPndLower_Textbox,'String'));
    if isnan(BPndlow) || (BPndlow < 0)
        warndlg('Lower restrict must be a positive number','Pixel-wise error','modal');
        uicontrol(handles.BPndLower_Textbox);
        return
    end
    
    BPndupp = str2double(get(handles.BPndUpper_Textbox,'String'));
    if isnan(BPndupp) || (BPndupp < 0)
        warndlg('Upper restrict must be a positive number','Pixel-wise error','modal');
        uicontrol(handles.BPndUpper_Textbox);
        return
    end
    
    if BPndupp < BPndlow
        warndlg('Upper restrict must be bigger than lower restrict','Pixel-wise error','modal');
        uicontrol(handles.BPndUpper_Textbox);
        return
    end
end

if k2acheck == 1 && k2arest == 1
    k2alow = str2double(get(handles.k2aLower_Textbox,'String'));
    if isnan(k2alow)
        warndlg('Lower restrict must be a positive number','Pixel-wise error','modal');
        uicontrol(handles.k2aLower_Textbox);
        return
    end
    
    k2aupp = str2double(get(handles.k2aUpper_Textbox,'String'));
    if isnan(k2aupp)
        warndlg('Upper restrict must be a positive number','Pixel-wise error','modal');
        uicontrol(handles.k2aUpper_Textbox);
        return
    end
    
    if k2aupp < k2alow
        warndlg('Upper restrict must be bigger than lower restrict','Pixel-wise error','modal');
        uicontrol(handles.k2aUpper_Textbox);
        return
    end
end

if k2check == 1 && k2rest == 1
    k2low = str2double(get(handles.k2Lower_Textbox,'String'));
    if isnan(k2low)
        warndlg('Lower restrict must be a positive number','Pixel-wise error','modal');
        uicontrol(handles.k2Lower_Textbox);
        return
    end
    
    k2upp = str2double(get(handles.k2Upper_Textbox,'String'));
    if isnan(k2upp)
        warndlg('Upper restrict must be a positive number','Pixel-wise error','modal');
        uicontrol(handles.k2Upper_Textbox);
        return
    end
    
    if k2upp < k2low
        warndlg('Upper restrict must be bigger than lower restrict','Pixel-wise error','modal');
        uicontrol(handles.k2Upper_Textbox);
        return
    end
end

if R1check == 1 && R1rest == 1
    R1low = str2double(get(handles.R1Lower_Textbox,'String'));
    if isnan(R1low)
        warndlg('Lower restrict must be a positive number','Pixel-wise error','modal');
        uicontrol(handles.R1Lower_Textbox);
        return
    end
    
    R1upp = str2double(get(handles.R1Upper_Textbox,'String'));
    if isnan(R1upp)
        warndlg('Upper restrict must be a positive number','Pixel-wise error','modal');
        uicontrol(handles.R1Upper_Textbox);
        return
    end
    
    if R1upp < R1low
        warndlg('Upper restrict must be bigger than lower restrict','Pixel-wise error','modal');
        uicontrol(handles.R1Upper_Textbox);
        return
    end
end

%Update the logfile
aux_date=datestr(now);
fprintf(QMf_logfile,strcat('\n-SRTM-\n'));
fprintf(QMf_logfile,strcat(aux_date(end-8:end),' Generating parametric images (BF= ',num2str(QMlastPreprocess.SRTM.numBF),', threshold= ',num2str(QMlastPreprocess.SRTM.threshold),')\n'));

%Image calculation
images = QM_SRTMPixelCalc(QMmaindata.Cpet,QMmaindata.Crois{2,QMlastPreprocess.SRTM.idCr},...
            QMlastPreprocess.SRTM.k2amin,QMlastPreprocess.SRTM.k2amax,QMlastPreprocess.SRTM.threshold,QMlastPreprocess.SRTM.M,...
            QMlastPreprocess.SRTM.B,QMlastPreprocess.SRTM.numBF,waitbarhandle,...
            [BPndrest BPndlow BPndupp BPndcheck; k2rest k2low k2upp k2check; R1rest R1low R1upp R1check; k2arest k2alow k2aupp k2acheck]);

%Update the logfile
aux_date=datestr(now);
fprintf(QMf_logfile,strcat(aux_date(end-8:end),' Parametric images generated \n'));
        
% %Reference header to make the new images
% studyPaths = get(mainhandles.mask_Textbox,'UserData');
% Vtemp = spm_vol(studyPaths{1});
% %Vtemp.mat(1:3,4) = 0;

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
    BPndfile = get(handles.BPndFilename_Textbox,'String');
    
    if get(handles.BPndFileExt_Popupmenu,'Value') == 1
        BPndfileExt = '.img';
    else
        BPndfileExt = '.nii';
    end
    
    VBPnd = struct('fname',strcat(BPndfile,BPndfileExt),'dim',Vtemp.dim,'mat',Vtemp.mat,...
            'dt',[spm_type('float64') spm_platform('bigend')],'n',[1 1]);

    %Create image
    aux = zeros(VBPnd.dim);
    aux(:) = images(1,:);
    
    if get(handles.saveBPnd_Checkbox,'Value') == 1
        VBPnd.fname = strcat(outputDir,filesep,BPndfile,BPndfileExt);
        spm_write_vol(VBPnd,aux);
        fprintf(QMf_logfile,strcat(['  -BPnd image saved in ',char(32),outputDir,'\n']));
    end
    VBPnd.fname = strcat(QMpath,filesep,'temp',filesep,BPndfile,BPndfileExt);
    spm_write_vol(VBPnd,aux); %Save image in hard disk
        
    imgTempPaths{1,indTemp} = VBPnd.fname;
    imgTempFilenames{1,indTemp} = 'BPnd (estimated binding potencial)';
    indTemp=indTemp+1;
end

if k2check == 1 
    k2file = get(handles.k2Filename_Textbox,'String');

    if get(handles.k2FileExt_Popupmenu,'Value') == 1
        k2fileExt = '.img';
    else
        k2fileExt = '.nii';
    end

    Vk2 = struct('fname',strcat(k2file,k2fileExt),'dim',Vtemp.dim,'mat',Vtemp.mat,...
            'dt',[spm_type('float64') spm_platform('bigend')],'n',[1 1]);
    
    %Create image
    aux = zeros(Vk2.dim);
    aux(:) = images(2,:);
    
    if get(handles.savek2_Checkbox,'Value') == 1
        Vk2.fname = strcat(outputDir,filesep,k2file,k2fileExt);
        spm_write_vol(Vk2,aux);
        fprintf(QMf_logfile,strcat(['  -k2 image saved in ',outputDir,'\n']));
    end        
    Vk2.fname = strcat(QMpath,filesep,'temp',filesep,k2file,k2fileExt);
    spm_write_vol(Vk2,aux);  %Save image in hard disk
        
    imgTempPaths{1,indTemp} = Vk2.fname;
    imgTempFilenames{1,indTemp} = 'k2 (estimated efflux rate constant)';
    indTemp=indTemp+1;
end

if R1check == 1 
    R1file = get(handles.R1Filename_Textbox,'String');

    if get(handles.R1FileExt_Popupmenu,'Value') == 1
        R1fileExt = '.img';
    else
        R1fileExt = '.nii';
    end

    VR1 = struct('fname',strcat(R1file,R1fileExt),'dim',Vtemp.dim,'mat',Vtemp.mat,...
            'dt',[spm_type('float64') spm_platform('bigend')],'n',[1 1]);
    
    %Create image
    aux = zeros(VR1.dim);
    aux(:) = images(3,:);
    
    if get(handles.saveR1_Checkbox,'Value') == 1
        VR1.fname = strcat(outputDir,filesep,R1file,R1fileExt);
        spm_write_vol(VR1,aux);
        fprintf(QMf_logfile,strcat(['  -R1 image saved in ',outputDir,'\n']));
    end
    VR1.fname = strcat(QMpath,filesep,'temp',filesep,R1file,R1fileExt);
    spm_write_vol(VR1,aux);  %Save image in hard disk
        
    imgTempPaths{1,indTemp} = VR1.fname;
    imgTempFilenames{1,indTemp} = 'R1 (ratio of tracer delivery)';
    indTemp=indTemp+1;
end

if k2acheck == 1 
    k2afile = get(handles.k2aFilename_Textbox,'String');
    
    if get(handles.k2aFileExt_Popupmenu,'Value') == 1
        k2afileExt = '.img';
    else
        k2afileExt = '.nii';
    end
    
    Vk2a = struct('fname',strcat(k2afile,k2afileExt),'dim',Vtemp.dim,'mat',Vtemp.mat,...
            'dt',[spm_type('float64') spm_platform('bigend')],'n',[1 1]);

    %Create image
    aux = zeros(Vk2a.dim);
    aux(:) = images(4,:);
    
    if get(handles.savek2a_Checkbox,'Value') == 1
        Vk2a.fname = strcat(outputDir,filesep,k2afile,k2afileExt);
        spm_write_vol(Vk2a,aux);
        fprintf(QMf_logfile,strcat(['  -k2a image saved in ',outputDir,'\n']));
    end
    Vk2a.fname = strcat(QMpath,filesep,'temp',filesep,k2afile,k2afileExt);
    spm_write_vol(Vk2a,aux);  %Save image in hard disk
        
    imgTempPaths{1,indTemp} = Vk2a.fname;
    imgTempFilenames{1,indTemp} = 'k2a (the best least squares fit)';
end

delete(waitbarhandle);  %Need to delete before ploting
set(mainhandles.selectImage_Popupmenu,'UserData',imgTempPaths);
set(mainhandles.selectImage_Popupmenu,'String',imgTempFilenames);

set(mainhandles.setMapParam_Pushbutton,'UserData',1);

uiresume(handles.MapParSRTM_figure);


% --- Executes on button press in setDefault_Pushbutton.
function setDefault_Pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to setDefault_Pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Update components
set(handles.BPndLower_Textbox,'String','0.0','Enable','on');
set(handles.BPndUpper_Textbox,'String','10.0','Enable','on');
set(handles.k2Lower_Textbox,'String','0.0','Enable','on');
set(handles.k2Upper_Textbox,'String','1.0','Enable','on');
set(handles.R1Lower_Textbox,'String','0.0','Enable','on');
set(handles.R1Upper_Textbox,'String','5.0','Enable','on');
set(handles.k2aLower_Textbox,'String','0.006','Enable','off');
set(handles.k2aUpper_Textbox,'String','0.6','Enable','off');

set(handles.BPndFilename_Textbox,'String','BPnd_image','Enable','off');
set(handles.k2Filename_Textbox,'String','k2_image','Enable','off');
set(handles.R1Filename_Textbox,'String','R1_image','Enable','off');
set(handles.k2aFilename_Textbox,'String','k2a_image','Enable','off');

set(handles.BPnd_Checkbox,'Value',1.0);
set(handles.k2_Checkbox,'Value',1.0);
set(handles.R1_Checkbox,'Value',1.0);
set(handles.k2a_Checkbox,'Value',1.0);

set(handles.BPndRestrict_Checkbox,'Enable','on');
set(handles.k2Restrict_Checkbox,'Enable','on');
set(handles.R1Restrict_Checkbox,'Enable','on');
set(handles.k2aRestrict_Checkbox,'Enable','on');

set(handles.BPndRestrict_Checkbox,'Value',1.0);
set(handles.k2Restrict_Checkbox,'Value',1.0);
set(handles.R1Restrict_Checkbox,'Value',1.0);
set(handles.k2aRestrict_Checkbox,'Value',0.0);

set(handles.saveBPnd_Checkbox,'Enable','on');
set(handles.savek2_Checkbox,'Enable','on');
set(handles.saveR1_Checkbox,'Enable','on');
set(handles.savek2a_Checkbox,'Enable','on');

set(handles.savek2a_Checkbox,'Value',0.0);
set(handles.saveBPnd_Checkbox,'Value',0.0);
set(handles.savek2_Checkbox,'Value',0.0);
set(handles.saveR1_Checkbox,'Value',0.0);

set(handles.BPndFileExt_Popupmenu,'Value',1,'Enable','off');
set(handles.k2FileExt_Popupmenu,'Value',1,'Enable','off');
set(handles.R1FileExt_Popupmenu,'Value',1,'Enable','off');
set(handles.k2aFileExt_Popupmenu,'Value',1,'Enable','off');

dir = pwd;
set(handles.outputDir_Textbox,'String',dir);
set(handles.runPixelCalc_Pushbutton,'Enable','on')

% --- Executes on button press in BPnd_Checkbox.
function BPnd_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to BPnd_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of BPnd_Checkbox

checked = get(hObject,'Value');
if checked
    enable = 'on';
    set(handles.runPixelCalc_Pushbutton,'Enable',enable)
else
    enable = 'off';
    if (~get(handles.R1_Checkbox,'Value')) && (~get(handles.k2a_Checkbox,'Value'))&&(~get(handles.k2a_Checkbox,'Value'))
        set(handles.runPixelCalc_Pushbutton,'Enable','off')
    end
end
set(handles.BPndRestrict_Checkbox,'Enable',enable);
set(handles.saveBPnd_Checkbox,'Enable',enable);
% set(handles.BPndfilename_Textbox,'Enable',enable);
% set(handles.BPndfileExt_Popupmenu,'Enable',enable);
set(handles.BPndFilename_Textbox,'Enable','off');
set(handles.BPndFileExt_Popupmenu,'Enable','off');
set(handles.BPndLower_Textbox,'Enable',enable);
set(handles.BPndUpper_Textbox,'Enable',enable);


% --- Executes on button press in k2a_Checkbox.
function k2a_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to k2a_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of k2a_Checkbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
    set(handles.runPixelCalc_Pushbutton,'Enable',enable)
else
    enable = 'off';
    if (~get(handles.BPnd_Checkbox,'Value')) && (~get(handles.k2_Checkbox,'Value'))&&(~get(handles.R1_Checkbox,'Value'))
        set(handles.runPixelCalc_Pushbutton,'Enable','off')
    end
end
set(handles.k2aRestrict_Checkbox,'Enable',enable);
set(handles.savek2a_Checkbox,'Enable',enable);
% set(handles.k2afilename_Textbox,'Enable',enable);
% set(handles.k2afileExt_Popupmenu,'Enable',enable);
set(handles.k2aFilename_Textbox,'Enable','off');
set(handles.k2aFileExt_Popupmenu,'Enable','off');
set(handles.k2aLower_Textbox,'Enable',enable);
set(handles.k2aUpper_Textbox,'Enable',enable);


% --- Executes on button press in BPndRestrict_Checkbox.
function BPndRestrict_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to BPndRestrict_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of BPndRestrict_Checkbox

checked = get(hObject,'Value');
if checked
    enable = 'on';
else
    enable = 'off';
end
set(handles.BPndLower_Textbox,'Enable',enable);
set(handles.BPndUpper_Textbox,'Enable',enable);


function BPndLower_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to BPndLower_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BPndLower_Textbox as text
%        str2double(get(hObject,'String')) returns contents of BPndLower_Textbox as a double
value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('Lower value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function BPndLower_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BPndLower_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function BPndUpper_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to BPndUpper_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BPndUpper_Textbox as text
%        str2double(get(hObject,'String')) returns contents of BPndUpper_Textbox as a double
value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('Upper value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function BPndUpper_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BPndUpper_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in k2aRestrict_Checkbox.
function k2aRestrict_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to k2aRestrict_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of k2aRestrict_Checkbox

checked = get(hObject,'Value');
if checked
    enable = 'on';
else
    enable = 'off';
end
set(handles.k2aLower_Textbox,'Enable',enable);
set(handles.k2aUpper_Textbox,'Enable',enable);


function k2aLower_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to k2aLower_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k2aLower_Textbox as text
%        str2double(get(hObject,'String')) returns contents of k2aLower_Textbox as a double
value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('Lower value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function k2aLower_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k2aLower_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function k2aUpper_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to k2aUpper_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k2aUpper_Textbox as text
%        str2double(get(hObject,'String')) returns contents of k2aUpper_Textbox as a double
value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('Upper value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function k2aUpper_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k2aUpper_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in saveBPnd_Checkbox.
function saveBPnd_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to saveBPnd_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of saveBPnd_Checkbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
else
    enable = 'off';
end
set(handles.BPndFilename_Textbox,'Enable',enable);
set(handles.BPndFileExt_Popupmenu,'Enable',enable);


function BPndFilename_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to BPndFilename_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BPndFilename_Textbox as text
%        str2double(get(hObject,'String')) returns contents of BPndFilename_Textbox as a double


% --- Executes during object creation, after setting all properties.
function BPndFilename_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BPndFilename_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in savek2a_Checkbox.
function savek2a_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to savek2a_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of savek2a_Checkbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
else
    enable = 'off';
end
set(handles.k2aFilename_Textbox,'Enable',enable);
set(handles.k2aFileExt_Popupmenu,'Enable',enable);


function k2aFilename_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to k2aFilename_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k2aFilename_Textbox as text
%        str2double(get(hObject,'String')) returns contents of k2aFilename_Textbox as a double


% --- Executes during object creation, after setting all properties.
function k2aFilename_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k2aFilename_Textbox (see GCBO)
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

set(handles.MapParSRTM_figure,'WindowStyle','normal');
directoryname = uigetdir('', 'Select the output directory');
set(handles.MapParSRTM_figure,'WindowStyle','modal');

if directoryname ~= 0
     set(handles.outputDir_Textbox,'String',directoryname);
end


% --- Executes on selection change in BPndFileExt_Popupmenu.
function BPndFileExt_Popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to BPndFileExt_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns BPndFileExt_Popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from BPndFileExt_Popupmenu


% --- Executes during object creation, after setting all properties.
function BPndFileExt_Popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BPndFileExt_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in k2aFileExt_Popupmenu.
function k2aFileExt_Popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to k2aFileExt_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns k2aFileExt_Popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from k2aFileExt_Popupmenu


% --- Executes during object creation, after setting all properties.
function k2aFileExt_Popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k2aFileExt_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in k2_Checkbox.
function k2_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to k2_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of k2_Checkbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
    set(handles.runPixelCalc_Pushbutton,'Enable',enable)
else
    enable = 'off';
    if (~get(handles.BPnd_Checkbox,'Value')) && (~get(handles.k2a_Checkbox,'Value'))&&(~get(handles.R1_Checkbox,'Value'))
        set(handles.runPixelCalc_Pushbutton,'Enable','off')
    end
end
set(handles.k2Restrict_Checkbox,'Enable',enable);
set(handles.savek2_Checkbox,'Enable',enable);
% set(handles.k2filename_Textbox,'Enable',enable);
% set(handles.k2fileExt_Popupmenu,'Enable',enable);
set(handles.k2Filename_Textbox,'Enable','off');
set(handles.k2FileExt_Popupmenu,'Enable','off');
set(handles.k2Lower_Textbox,'Enable',enable);
set(handles.k2Upper_Textbox,'Enable',enable);


% --- Executes on button press in k2Restrict_Checkbox.
function k2Restrict_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to k2Restrict_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of k2Restrict_Checkbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
else
    enable = 'off';
end
set(handles.k2Lower_Textbox,'Enable',enable);
set(handles.k2Upper_Textbox,'Enable',enable);


% --- Executes on button press in savek2_Checkbox.
function savek2_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to savek2_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of savek2_Checkbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
else
    enable = 'off';
end
set(handles.k2Filename_Textbox,'Enable',enable);
set(handles.k2FileExt_Popupmenu,'Enable',enable);



function k2Lower_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to k2Lower_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k2Lower_Textbox as text
%        str2double(get(hObject,'String')) returns contents of k2Lower_Textbox as a double
value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('Lower value must be a positive number','Warning','modal');
    uicontrol(hObject);
end


% --- Executes during object creation, after setting all properties.
function k2Lower_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k2Lower_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function k2Filename_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to k2Filename_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k2Filename_Textbox as text
%        str2double(get(hObject,'String')) returns contents of k2Filename_Textbox as a double


% --- Executes during object creation, after setting all properties.
function k2Filename_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k2Filename_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function k2Upper_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to k2Upper_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k2Upper_Textbox as text
%        str2double(get(hObject,'String')) returns contents of k2Upper_Textbox as a double
value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('Upper value must be a positive number','Warning','modal');
    uicontrol(hObject);
end


% --- Executes during object creation, after setting all properties.
function k2Upper_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k2Upper_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in k2FileExt_Popupmenu.
function k2FileExt_Popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to k2FileExt_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns k2FileExt_Popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from k2FileExt_Popupmenu


% --- Executes during object creation, after setting all properties.
function k2FileExt_Popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k2FileExt_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in R1_Checkbox.
function R1_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to R1_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of R1_Checkbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
    set(handles.runPixelCalc_Pushbutton,'Enable',enable)
else
    enable = 'off';
    if (~get(handles.BPnd_Checkbox,'Value')) && (~get(handles.k2a_Checkbox,'Value'))&&(~get(handles.k2a_Checkbox,'Value'))
        set(handles.runPixelCalc_Pushbutton,'Enable','off')
    end
end
set(handles.R1Restrict_Checkbox,'Enable',enable);
set(handles.saveR1_Checkbox,'Enable',enable);
% set(handles.R1filename_Textbox,'Enable',enable);
% set(handles.R1fileExt_Popupmenu,'Enable',enable);
set(handles.R1Filename_Textbox,'Enable','off');
set(handles.R1FileExt_Popupmenu,'Enable','off');
set(handles.R1Lower_Textbox,'Enable',enable);
set(handles.R1Upper_Textbox,'Enable',enable);


% --- Executes on button press in R1Restrict_Checkbox.
function R1Restrict_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to R1Restrict_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of R1Restrict_Checkbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
else
    enable = 'off';
end
set(handles.R1Lower_Textbox,'Enable',enable);
set(handles.R1Upper_Textbox,'Enable',enable);


% --- Executes on button press in saveR1_Checkbox.
function saveR1_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to saveR1_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of saveR1_Checkbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
else
    enable = 'off';
end
set(handles.R1Filename_Textbox,'Enable',enable);
set(handles.R1FileExt_Popupmenu,'Enable',enable);



function R1Lower_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to R1Lower_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of R1Lower_Textbox as text
%        str2double(get(hObject,'String')) returns contents of R1Lower_Textbox as a double
value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('Lower value must be a positive number','Warning','modal');
    uicontrol(hObject);
end


% --- Executes during object creation, after setting all properties.
function R1Lower_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to R1Lower_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function R1Filename_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to R1Filename_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of R1Filename_Textbox as text
%        str2double(get(hObject,'String')) returns contents of R1Filename_Textbox as a double


% --- Executes during object creation, after setting all properties.
function R1Filename_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to R1Filename_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function R1Upper_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to R1Upper_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of R1Upper_Textbox as text
%        str2double(get(hObject,'String')) returns contents of R1Upper_Textbox as a double
value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('Upper value must be a positive number','Warning','modal');
    uicontrol(hObject);
end


% --- Executes during object creation, after setting all properties.
function R1Upper_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to R1Upper_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in R1FileExt_Popupmenu.
function R1FileExt_Popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to R1FileExt_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns R1FileExt_Popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from R1FileExt_Popupmenu


% --- Executes during object creation, after setting all properties.
function R1FileExt_Popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to R1FileExt_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close MapParSRTM_figure.
function MapParSRTM_figure_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to MapParSRTM_figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
uiresume(handles.MapParSRTM_figure);
