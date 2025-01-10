function varargout = QM_mapParamCarson(varargin)
% QM_MAPPARAMCARSON MATLAB code for QM_mapParamCarson.fig
%      
%   QM_mapParamCarson generates a GUI to choose the parametric
%   images which you want to calculate from the previous preprocessing
%   of SRTM2 model. This code must be executed from QModeling GUI.
%
% Last Modified by GUIDE v2.5 29-Oct-2013 12:58:49
%---------------------------------------------------------------------------------
% QM_mapParamCarson is part of QModeling.
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
                   'gui_OpeningFcn', @QM_mapParamCarson_OpeningFcn, ...
                   'gui_OutputFcn',  @QM_mapParamCarson_OutputFcn, ...
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


% --- Executes just before QM_mapParamCarson is made visible.
function QM_mapParamCarson_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to QM_mapParamCarson (see VARARGIN)

% Choose default command line output for QM_mapParamCarson
handles.output = hObject;

mainGuiInput = find(strcmp(varargin, 'mainhandles'));
handles.mainGUI = varargin{mainGuiInput+1};

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes QM_mapParamCarson wait for user response (see UIRESUME)
uiwait(handles.MapParCarson_figure);


% --- Outputs from this function are returned to the command line.
function varargout = QM_mapParamCarson_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

delete(handles.MapParCarson_figure);


% --- Executes on button press in runPixelCalc_Pushbutton.
function runPixelCalc_Pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to runPixelCalc_Pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global QMmaindata
global QMlastPreprocess
global QMpath

global QMf_logfile

BPndlow = 0;
BPndupp = 0;
R1low = 0;
R1upp = 0;
k2alow=0;
k2aupp=0;
k2low=0;
k2upp=0;

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

BPndcheck = get(handles.BPnd_Checkbox,'Value');
R1check = get(handles.R1_Checkbox,'Value');
k2acheck = get(handles.k2a_Checkbox,'Value');
k2check = get(handles.k2_Checkbox,'Value');

BPndrest = get(handles.BPndrestrict_Checkbox,'Value');
R1rest = get(handles.R1restrict_Checkbox,'Value');
k2arest = get(handles.k2arestrict_Checkbox,'Value');
k2rest = get(handles.k2restrict_Checkbox,'Value');

%Checking errors
if BPndcheck == 1 && BPndrest == 1    
    BPndlow = str2double(get(handles.BPndlower_Textbox,'String'));
    if isnan(BPndlow)
        warndlg('Lower restrict must be a positive number','Pixel-wise error','modal');
        uicontrol(handles.BPndlower_Textbox);
        return
    end
    
    BPndupp = str2double(get(handles.BPndupper_Textbox,'String'));
    if isnan(BPndupp)
        warndlg('Upper restrict must be a positive number','Pixel-wise error','modal');
        uicontrol(handles.BPndupper_Textbox);
        return
    end  
    
    if BPndlow > BPndupp
        warndlg('Upper restrict for BPnd value must be bigger than Lower restrict','Pixel-wise error','modal');
        uicontrol(handles.BPndupper_Textbox);
        return
    end 
end

if R1check == 1 && R1rest == 1    
    R1low = str2double(get(handles.R1lower_Textbox,'String'));
    if isnan(R1low)
        warndlg('Lower restrict must be a positive number','Pixel-wise error','modal');
        uicontrol(handles.R1lower_Textbox);
        return
    end
    
    R1upp = str2double(get(handles.R1upper_Textbox,'String'));
    if isnan(R1upp)
        warndlg('Upper restrict must be a positive number','Pixel-wise error','modal');
        uicontrol(handles.R1upper_Textbox);
        return
    end
    
    if R1low > R1upp
        warndlg('Upper restrict for R1 value must be bigger than Lower restrict','Pixel-wise error','modal');
        uicontrol(handles.R1upper_Textbox);
        return
    end 
end

if k2acheck == 1 && k2arest == 1    
    k2alow = str2double(get(handles.k2alower_Textbox,'String'));
    if isnan(k2alow)
        warndlg('Lower restrict must be a positive number','Pixel-wise error','modal');
        uicontrol(handles.k2alower_Textbox);
        return
    end
    
    k2aupp = str2double(get(handles.k2aupper_Textbox,'String'));
    if isnan(k2aupp)
        warndlg('Upper restrict must be a positive number','Pixel-wise error','modal');
        uicontrol(handles.k2aupper_Textbox);
        return
    end  
    
    if k2alow > k2aupp
        warndlg('Upper restrict for k2a value must be bigger than Lower restrict','Pixel-wise error','modal');
        uicontrol(handles.k2aupper_Textbox);
        return
    end 
end

if k2check == 1 && k2rest == 1    
    k2low = str2double(get(handles.k2lower_Textbox,'String'));
    if isnan(k2low)
        warndlg('Lower restrict must be a positive number','Pixel-wise error','modal');
        uicontrol(handles.k2lower_Textbox);
        return
    end
    
    k2upp = str2double(get(handles.k2upper_Textbox,'String'));
    if isnan(k2upp)
        warndlg('Upper restrict must be a positive number','Pixel-wise error','modal');
        uicontrol(handles.k2upper_Textbox);
        return
    end  
    
    if k2low > k2upp
        warndlg('Upper restrict for k2 value must be bigger than Lower restrict','Pixel-wise error','modal');
        uicontrol(handles.k2upper_Textbox);
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
fprintf(QMf_logfile,strcat('\n-SRTM2-\n'));
fprintf(QMf_logfile,strcat(aux_date(end-8:end),' Generating parametric images (BF= ',num2str(length(QMlastPreprocess.SRTM2.th3)),', threshold= ',num2str(QMlastPreprocess.SRTM2.threshold),', k2',char(39),'=',...
    num2str(QMlastPreprocess.SRTM2.results{1}),')\n'));

%Image calculation
images = QM_CarsonPixelCalc(QMmaindata.Cpet,QMlastPreprocess.SRTM2.B,QMlastPreprocess.SRTM2.th3,QMlastPreprocess.SRTM2.threshold,QMlastPreprocess.SRTM2.results{1},...
    [BPndrest BPndlow BPndupp ; R1rest R1low R1upp ; k2arest k2alow k2aupp ; k2rest k2low k2upp], waitbarhandle);

%Update the logfile
aux_date=datestr(now);
fprintf(QMf_logfile,strcat(aux_date(end-8:end),' Parametric images generated \n'));

delete(waitbarhandle);  %Need to delete before ploting

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

if BPndcheck==1
    BPndfile=get(handles.BPndfilename_Textbox,'String');    
   
    if get(handles.BPndfileExt_Popupmenu,'Value') == 1
        BPndfileExt = '.img';
    else
        BPndfileExt = '.nii';
    end
    
    VBPnd = struct('fname',strcat(BPndfile,BPndfileExt),'dim',Vtemp.dim,'mat',Vtemp.mat,...
                'dt',[spm_type('float64') spm_platform('bigend')],'n',[1 1]); 
            
    aux = zeros(VBPnd.dim);
    aux(:) = images(1,:); 
    
    if get(handles.saveBPnd_Checkbox,'Value') == 1
        VBPnd.fname = strcat(outputDir,filesep,BPndfile,BPndfileExt);
        spm_write_vol(VBPnd,aux);
        fprintf(QMf_logfile,strcat(['  -BPnd image saved in ',outputDir,'\n']));
    end
        
    VBPnd.fname = strcat(QMpath,filesep,'temp',filesep,BPndfile,BPndfileExt);
    spm_write_vol(VBPnd,aux);
    
    imgTempPaths{1,indTemp} = VBPnd.fname;
    imgTempFilenames{1,indTemp} = 'BPnd (estimated binding potencial)';
    indTemp=indTemp+1;    
    
end

if R1check==1    
    R1file=get(handles.R1filename_Textbox,'String');    
    
    if get(handles.R1fileExt_Popupmenu,'Value') == 1
        R1fileExt = '.img';
    else
        R1fileExt = '.nii';
    end       
    
    VR1 = struct('fname',strcat(R1file,R1fileExt),'dim',Vtemp.dim,'mat',Vtemp.mat,...
                'dt',[spm_type('float64') spm_platform('bigend')],'n',[1 1]); 
            
    aux = zeros(VR1.dim);
    aux(:) = images(2,:); 
    
    if get(handles.saveR1_Checkbox,'Value') == 1
        VR1.fname = strcat(outputDir,filesep,R1file,R1fileExt); 
        spm_write_vol(VR1,aux);
        fprintf(QMf_logfile,strcat(['  -R1 image saved in ',outputDir,'\n']));
    end
    
    VR1.fname = strcat(QMpath,filesep,'temp',filesep,R1file,R1fileExt);
    spm_write_vol(VR1,aux);
    
    imgTempPaths{1,indTemp} = VR1.fname;
    imgTempFilenames{1,indTemp} = 'R1(ratio of tracer delivery)';
    indTemp=indTemp+1;      
end

if k2acheck==1
    k2afile=get(handles.k2afilename_Textbox,'String');    
        
    if get(handles.k2afileExt_Popupmenu,'Value') == 1
        k2afileExt = '.img';
    else
        k2afileExt = '.nii';
    end
    
    Vk2a = struct('fname',strcat(k2afile,k2afileExt),'dim',Vtemp.dim,'mat',Vtemp.mat,...
                'dt',[spm_type('float64') spm_platform('bigend')],'n',[1 1]); 
            
    aux = zeros(Vk2a.dim);
    aux(:) = images(3,:);
    
    if get(handles.savek2a_Checkbox,'Value') == 1
        Vk2a.fname = strcat(outputDir,filesep,k2afile,k2afileExt);
        spm_write_vol(Vk2a,aux); 
        fprintf(QMf_logfile,strcat(['  -k2a image saved in ',outputDir,'\n']));
    end
    
    Vk2a.fname =strcat(QMpath,filesep,'temp',filesep,k2afile,k2afileExt);       
    spm_write_vol(Vk2a,aux);
       
    imgTempPaths{1,indTemp} = Vk2a.fname;
    imgTempFilenames{1,indTemp} = 'K2a (the best least squares fit)';
    indTemp=indTemp+1; 
end

if k2check==1
    k2file=get(handles.k2filename_Textbox,'String');   
            
    if get(handles.k2fileExt_Popupmenu,'Value') == 1
        k2fileExt = '.img';
    else
        k2fileExt = '.nii';
    end
    
     Vk2 = struct('fname',strcat(k2file,k2fileExt),'dim',Vtemp.dim,'mat',Vtemp.mat,...
                'dt',[spm_type('float64') spm_platform('bigend')],'n',[1 1]); 
            
    aux = zeros(Vk2.dim);
    aux(:) = images(4,:);
    
    if get(handles.savek2_Checkbox,'Value') == 1
        Vk2.fname = strcat(outputDir,filesep,k2file,k2fileExt); 
        spm_write_vol(Vk2,aux);
        fprintf(QMf_logfile,strcat(['  -k2 image saved in ',outputDir,'\n']));
    end
    
    Vk2.fname = strcat(QMpath,filesep,'temp',filesep, k2file,k2fileExt);
    spm_write_vol(Vk2,aux);  
    
    imgTempPaths{1,indTemp} = Vk2.fname;
    imgTempFilenames{1,indTemp} = 'K2 (estimated efflux rate constant k2)';
    
end

set(mainhandles.selectImage_Popupmenu,'UserData',imgTempPaths);
set(mainhandles.selectImage_Popupmenu,'String',imgTempFilenames);

set(mainhandles.setMapParam_Pushbutton,'UserData',1);

uiresume(handles.MapParCarson_figure);


% --- Executes on button press in setDefault_Pushbutton.
function setDefault_Pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to setDefault_Pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.BPndlower_Textbox,'String','0.0','Enable','on');
set(handles.BPndupper_Textbox,'String','20.0','Enable','on');
set(handles.R1lower_Textbox,'String','0.0','Enable','on');
set(handles.R1upper_Textbox,'String','10.0','Enable','on');
set(handles.k2alower_Textbox,'String','0.0','Enable','on');
set(handles.k2aupper_Textbox,'String','1.0','Enable','on');
set(handles.k2lower_Textbox,'String','0.0','Enable','on');
set(handles.k2upper_Textbox,'String','1.0','Enable','on');

set(handles.BPndfilename_Textbox,'String','BPnd_image','Enable','off');
set(handles.R1filename_Textbox,'String','R1_image','Enable','off');
set(handles.k2afilename_Textbox,'String','k2a_image','Enable','off');
set(handles.k2filename_Textbox,'String','k2_image','Enable','off');

set(handles.BPnd_Checkbox,'Value',1);
set(handles.R1_Checkbox,'Value',1);
set(handles.k2a_Checkbox,'Value',1);
set(handles.k2_Checkbox,'Value',1);

set(handles.BPndrestrict_Checkbox,'Value',1,'Enable','on');
set(handles.R1restrict_Checkbox,'Value',1,'Enable','on');
set(handles.k2arestrict_Checkbox,'Value',1,'Enable','on');
set(handles.k2restrict_Checkbox,'Value',1,'Enable','on');
set(handles.saveBPnd_Checkbox,'Value',0,'Enable','on');
set(handles.saveR1_Checkbox,'Value',0,'Enable','on');
set(handles.savek2a_Checkbox,'Value',0,'Enable','on');
set(handles.savek2_Checkbox,'Value',0,'Enable','on');

set(handles.BPndfileExt_Popupmenu,'Value',1,'Enable','off');
set(handles.R1fileExt_Popupmenu,'Value',1,'Enable','off');
set(handles.k2afileExt_Popupmenu,'Value',1,'Enable','off');
set(handles.k2fileExt_Popupmenu,'Value',1,'Enable','off');

dir=pwd;
set(handles.outputDir_Textbox,'String', dir);
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
set(handles.BPndrestrict_Checkbox,'Enable',enable);
set(handles.saveBPnd_Checkbox,'Enable',enable);
% set(handles.BPndfilename_Textbox,'Enable',enable);
% set(handles.BPndfileExt_Popupmenu,'Enable',enable);
set(handles.BPndfilename_Textbox,'Enable','off');
set(handles.BPndfileExt_Popupmenu,'Enable','off');
set(handles.BPndlower_Textbox,'Enable',enable);
set(handles.BPndupper_Textbox,'Enable',enable);


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
set(handles.R1restrict_Checkbox,'Enable',enable);
set(handles.saveR1_Checkbox,'Enable',enable);
% set(handles.R1filename_Textbox,'Enable',enable);
% set(handles.R1fileExt_Popupmenu,'Enable',enable);
set(handles.R1filename_Textbox,'Enable','off');
set(handles.R1fileExt_Popupmenu,'Enable','off');
set(handles.R1lower_Textbox,'Enable',enable);
set(handles.R1upper_Textbox,'Enable',enable);

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
set(handles.k2restrict_Checkbox,'Enable',enable);
set(handles.savek2_Checkbox,'Enable',enable);
% set(handles.k2filename_Textbox,'Enable',enable);
% set(handles.k2fileExt_Popupmenu,'Enable',enable);
set(handles.k2filename_Textbox,'Enable','off');
set(handles.k2fileExt_Popupmenu,'Enable','off');
set(handles.k2lower_Textbox,'Enable',enable);
set(handles.k2upper_Textbox,'Enable',enable);


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
set(handles.k2arestrict_Checkbox,'Enable',enable);
set(handles.savek2a_Checkbox,'Enable',enable);
% set(handles.k2afilename_Textbox,'Enable',enable);
% set(handles.k2afileExt_Popupmenu,'Enable',enable);
set(handles.k2afilename_Textbox,'Enable','off');
set(handles.k2afileExt_Popupmenu,'Enable','off');
set(handles.k2alower_Textbox,'Enable',enable);
set(handles.k2aupper_Textbox,'Enable',enable);



% --- Executes on button press in BPndrestrict_Checkbox.
function BPndrestrict_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to BPndrestrict_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of BPndrestrict_Checkbox

checked = get(hObject,'Value');
if checked
    enable = 'on';
else
    enable = 'off';
end
set(handles.BPndlower_Textbox,'Enable',enable);
set(handles.BPndupper_Textbox,'Enable',enable);


function BPndlower_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to BPndlower_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BPndlower_Textbox as text
%        str2double(get(hObject,'String')) returns contents of BPndlower_Textbox as a double
value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('Lower value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function BPndlower_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BPndlower_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function BPndupper_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to BPndupper_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BPndupper_Textbox as text
%        str2double(get(hObject,'String')) returns contents of BPndupper_Textbox as a double
value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('Upper value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function BPndupper_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BPndupper_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in R1restrict_Checkbox.
function R1restrict_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to R1restrict_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of R1restrict_Checkbox

checked = get(hObject,'Value');
if checked
    enable = 'on';
else
    enable = 'off';
end
set(handles.R1lower_Textbox,'Enable',enable);
set(handles.R1upper_Textbox,'Enable',enable);


function R1lower_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to R1lower_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of R1lower_Textbox as text
%        str2double(get(hObject,'String')) returns contents of R1lower_Textbox as a double
value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('Lower value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function R1lower_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to R1lower_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function R1upper_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to R1upper_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of R1upper_Textbox as text
%        str2double(get(hObject,'String')) returns contents of R1upper_Textbox as a double
value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('Upper value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function R1upper_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to R1upper_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in k2arestrict_Checkbox.
function k2arestrict_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to k2arestrict_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of k2arestrict_Checkbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
else
    enable = 'off';
end
set(handles.k2alower_Textbox,'Enable',enable);
set(handles.k2aupper_Textbox,'Enable',enable);

function k2alower_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to k2alower_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k2alower_Textbox as text
%        str2double(get(hObject,'String')) returns contents of k2alower_Textbox as a double
value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('Lower value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function k2alower_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k2alower_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function k2aupper_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to k2aupper_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k2aupper_Textbox as text
%        str2double(get(hObject,'String')) returns contents of k2aupper_Textbox as a double
value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('Upper value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function k2aupper_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k2aupper_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in k2restrict_Checkbox.
function k2restrict_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to k2restrict_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of k2restrict_Checkbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
else
    enable = 'off';
end
set(handles.k2lower_Textbox,'Enable',enable);
set(handles.k2upper_Textbox,'Enable',enable);

function k2lower_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to k2lower_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k2lower_Textbox as text
%        str2double(get(hObject,'String')) returns contents of k2lower_Textbox as a double
value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('Lower value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function k2lower_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k2lower_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function k2upper_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to k2upper_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k2upper_Textbox as text
%        str2double(get(hObject,'String')) returns contents of k2upper_Textbox as a double
value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('Upper value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function k2upper_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k2upper_Textbox (see GCBO)
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
set(handles.BPndfilename_Textbox,'Enable',enable);
set(handles.BPndfileExt_Popupmenu,'Enable',enable);

function BPndfilename_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to BPndfilename_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BPndfilename_Textbox as text
%        str2double(get(hObject,'String')) returns contents of BPndfilename_Textbox as a double


% --- Executes during object creation, after setting all properties.
function BPndfilename_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BPndfilename_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


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
set(handles.R1filename_Textbox,'Enable',enable);
set(handles.R1fileExt_Popupmenu,'Enable',enable);


function R1filename_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to R1filename_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of R1filename_Textbox as text
%        str2double(get(hObject,'String')) returns contents of R1filename_Textbox as a double


% --- Executes during object creation, after setting all properties.
function R1filename_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to R1filename_Textbox (see GCBO)
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

set(handles.MapParCarson_figure,'WindowStyle','normal');
directoryname = uigetdir('', 'Select the output directory');
set(handles.MapParCarson_figure,'WindowStyle','modal');

if directoryname ~= 0
     set(handles.outputDir_Textbox,'String',directoryname);
end


% --- Executes on selection change in BPndfileExt_Popupmenu.
function BPndfileExt_Popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to BPndfileExt_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns BPndfileExt_Popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from BPndfileExt_Popupmenu


% --- Executes during object creation, after setting all properties.
function BPndfileExt_Popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BPndfileExt_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in R1fileExt_Popupmenu.
function R1fileExt_Popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to R1fileExt_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns R1fileExt_Popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from R1fileExt_Popupmenu


% --- Executes during object creation, after setting all properties.
function R1fileExt_Popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to R1fileExt_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
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
set(handles.k2afilename_Textbox,'Enable',enable);
set(handles.k2afileExt_Popupmenu,'Enable',enable);

function k2afilename_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to k2afilename_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k2afilename_Textbox as text
%        str2double(get(hObject,'String')) returns contents of k2afilename_Textbox as a double


% --- Executes during object creation, after setting all properties.
function k2afilename_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k2afilename_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


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
set(handles.k2filename_Textbox,'Enable',enable);
set(handles.k2fileExt_Popupmenu,'Enable',enable);

function k2filename_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to k2filename_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k2filename_Textbox as text
%        str2double(get(hObject,'String')) returns contents of k2filename_Textbox as a double


% --- Executes during object creation, after setting all properties.
function k2filename_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k2filename_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in k2afileExt_Popupmenu.
function k2afileExt_Popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to k2afileExt_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns k2afileExt_Popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from k2afileExt_Popupmenu


% --- Executes during object creation, after setting all properties.
function k2afileExt_Popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k2afileExt_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in k2fileExt_Popupmenu.
function k2fileExt_Popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to k2fileExt_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns k2fileExt_Popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from k2fileExt_Popupmenu


% --- Executes during object creation, after setting all properties.
function k2fileExt_Popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k2fileExt_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close MapParCarson_figure.
function MapParCarson_figure_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to MapParCarson_figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
uiresume(handles.MapParCarson_figure);
