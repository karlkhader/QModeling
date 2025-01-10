function varargout = QM_mapParamTwoTCM(varargin)
% QM_MAPPARAMTWOTCM MATLAB code for QM_mapParamTwoTCM.fig
%    
%   QM_mapParamTwoTCM generates a GUI to choose the parametric
%   images which you want to calculate from the previous preprocessing
%   of the Two Tissue Compartment model. This code must be executed from 
%   QModeling GUI.
%
% Last Modified by GUIDE v2.5 21-Jun-2018 18:32:42

%---------------------------------------------------------------------------------
% QM_mapParamTwoTCM is part of QModeling.
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
% Copyright (C) 2015, 2018 Francisco Javier Lopez Gonzalez, Jose Paredes Pacheco, 
% Karl-Khader Thurnhofer Hemsi

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @QM_mapParamTwoTCM_OpeningFcn, ...
                   'gui_OutputFcn',  @QM_mapParamTwoTCM_OutputFcn, ...
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


% --- Executes just before QM_mapParamTwoTCM is made visible.
function QM_mapParamTwoTCM_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to QM_mapParamTwoTCM (see VARARGIN)
global QMlastPreprocess
% Choose default command line output for QM_mapParamTwoTCM
handles.output = hObject;

mainGuiInput = find(strcmp(varargin, 'mainhandles'));
handles.mainGUI = varargin{mainGuiInput+1};

% Update handles structure
guidata(hObject, handles);

if ~QMlastPreprocess.TwoTCM.k4_checked %k4 fixed
    set(handles.k4_Checkbox,'Value',0.0);
    set(handles.k4_Checkbox,'Enable','off');
    set(handles.restrict_k4_Checkbox,'Enable','off');
    set(handles.save_k4_Checkbox,'Enable','off');
    set(handles.k4_Filename_Textbox,'Enable','off');
    set(handles.k4_FileExt_Popupmenu,'Enable','off');
    set(handles.k4_LowerTextbox,'Enable','off');
    set(handles.k4_UpperTextbox,'Enable','off');
    
    set(handles.k3divk4_Checkbox,'Value',0.0);
    set(handles.k3divk4_Checkbox,'Enable','off');
    set(handles.restrict_k3divk4_Checkbox,'Enable','off');
    set(handles.save_k3divk4_Checkbox,'Enable','off');
    set(handles.k3divk4_Filename_Textbox,'Enable','off');
    set(handles.k3divk4_FileExt_Popupmenu,'Enable','off');
    set(handles.k3divk4_LowerTextbox,'Enable','off');
    set(handles.k3divk4_UpperTextbox,'Enable','off');
    
    set(handles.alpha1_Checkbox,'Value',0.0);
    set(handles.alpha1_Checkbox,'Enable','off');
    set(handles.restrict_alpha1_Checkbox,'Enable','off');
    set(handles.save_alpha1_Checkbox,'Enable','off');
    set(handles.alpha1_Filename_Textbox,'Enable','off');
    set(handles.alpha1_FileExt_Popupmenu,'Enable','off');
    set(handles.alpha1_LowerTextbox,'Enable','off');
    set(handles.alpha1_UpperTextbox,'Enable','off');
    
end

% UIWAIT makes QM_mapParamTwoTCM wait for user response (see UIRESUME)
uiwait(handles.MapParTwoTCM_figure);




% --- Outputs from this function are returned to the command line.
function varargout = QM_mapParamTwoTCM_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

delete(handles.MapParTwoTCM_figure);

% --- Executes when user attempts to close MapParTwoTCM_figure.
function MapParTwoTCM_figure_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to MapParTwoTCM_figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
uiresume(handles.MapParTwoTCM_figure);


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
vBcheck = get(handles.vB_Checkbox,'Value');
K1check = get(handles.K1_Checkbox,'Value');
k2check = get(handles.k2_Checkbox,'Value');
k3check = get(handles.k3_Checkbox,'Value');
k4check = get(handles.k4_Checkbox,'Value');
Kicheck = get(handles.Ki_Checkbox,'Value');
Vscheck = get(handles.Vs_Checkbox,'Value');
Vtcheck = get(handles.Vt_Checkbox,'Value');
K1divk2check = get(handles.K1divk2_Checkbox,'Value');
k3divk4check = get(handles.k3divk4_Checkbox,'Value');
alpha1check = get(handles.alpha1_Checkbox,'Value');
alpha2check = get(handles.alpha2_Checkbox,'Value');
ttpcheck = get(handles.ttp_Checkbox,'Value');

vBrest = get(handles.restrict_vB_Checkbox,'Value');
K1rest = get(handles.restrict_K1_Checkbox,'Value');
k2rest = get(handles.restrict_k2_Checkbox,'Value');
k3rest = get(handles.restrict_k3_Checkbox,'Value');
k4rest = get(handles.restrict_k4_Checkbox,'Value');
Kirest = get(handles.restrict_Ki_Checkbox,'Value');
Vsrest = get(handles.restrict_Vs_Checkbox,'Value');
Vtrest = get(handles.restrict_Vt_Checkbox,'Value');
K1divk2rest = get(handles.restrict_K1divk2_Checkbox,'Value');
k3divk4rest = get(handles.restrict_k3divk4_Checkbox,'Value');
alpha1rest = get(handles.restrict_alpha1_Checkbox,'Value');
alpha2rest = get(handles.restrict_alpha2_Checkbox,'Value');
ttprest = get(handles.restrict_ttp_Checkbox,'Value');

% vBsave = get(handles.save_vB_Checkbox,'Value');
% K1save = get(handles.save_K1_Checkbox,'Value');
% k2save = get(handles.save_k2_Checkbox,'Value');
% k3save = get(handles.save_k3_Checkbox,'Value');
% k4save = get(handles.save_k4_Checkbox,'Value');
% Kisave = get(handles.save_Ki_Checkbox,'Value');
% Vssave = get(handles.save_Vs_Checkbox,'Value');
% Vtsave = get(handles.save_Vt_Checkbox,'Value');
% K1divk2save = get(handles.save_K1divk2_Checkbox,'Value');
% k3ddivk4save = get(handles.save_k3divk4_Checkbox,'Value');
% alpha1save = get(handles.save_alpha1_Checkbox,'Value');
% alpha2save = get(handles.save_alpha2_Checkbox,'Value');
% ttpsave = get(handles.save_ttp_Checkbox,'Value');

vBlow = 0;
vBupp = 0;
K1low = 0;
K1upp = 0;
k2low = 0;
k2upp = 0;
k3low = 0;
k3upp = 0;
k4low = 0;
k4upp = 0;
Kilow = 0;
Kiupp = 0;
Vslow = 0;
Vsupp = 0;
Vtlow = 0;
Vtupp = 0;
K1divk2low = 0;
K1divk2upp = 0;
k3divk4low = 0;
k3divk4upp = 0;
alpha1low = 0;
alpha1upp = 0;
alpha2low = 0;
alpha2upp = 0;
ttplow = 0;
ttpupp = 0;

%Checking errors and conditions
if vBcheck == 1 && vBrest == 1 
    vBlow = str2double(get(handles.vB_LowerTextbox,'String'));
    if isnan(vBlow) || (vBlow < 0)
        warndlg('Lower restrict must be a positive number','Pixel-wise error','modal');
        uicontrol(handles.vB_LowerTextbox);
        return
    end
    
    vBupp = str2double(get(handles.vB_UpperTextbox,'String'));
    if isnan(vBupp) || (vBupp < 0)
        warndlg('Upper restrict must be a positive number','Pixel-wise error','modal');
        uicontrol(handles.vB_UpperTextbox);
        return
    end
    
    if vBupp < vBlow
        warndlg('Upper restrict must be bigger than lower restrict','Pixel-wise error','modal');
        uicontrol(handles.vB_UpperTextbox);
        return
    end
end

if K1check == 1 && K1rest == 1 
    K1low = str2double(get(handles.K1_LowerTextbox,'String'));
    if isnan(K1low) || (K1low < 0)
        warndlg('Lower restrict must be a positive number','Pixel-wise error','modal');
        uicontrol(handles.K1_LowerTextbox);
        return
    end
    
    K1upp = str2double(get(handles.K1_UpperTextbox,'String'));
    if isnan(K1upp) || (K1upp < 0)
        warndlg('Upper restrict must be a positive number','Pixel-wise error','modal');
        uicontrol(handles.K1_UpperTextbox);
        return
    end
    
    if K1upp < K1low
        warndlg('Upper restrict must be bigger than lower restrict','Pixel-wise error','modal');
        uicontrol(handles.K1_UpperTextbox);
        return
    end
end

if k2check == 1 && k2rest == 1 
    k2low = str2double(get(handles.k2_LowerTextbox,'String'));
    if isnan(k2low) || (k2low < 0)
        warndlg('Lower restrict must be a positive number','Pixel-wise error','modal');
        uicontrol(handles.k2_LowerTextbox);
        return
    end
    
   k2upp = str2double(get(handles.k2_UpperTextbox,'String'));
    if isnan(k2upp) || (k2upp < 0)
        warndlg('Upper restrict must be a positive number','Pixel-wise error','modal');
        uicontrol(handles.k2_UpperTextbox);
        return
    end
    
    if k2upp < k2low
        warndlg('Upper restrict must be bigger than lower restrict','Pixel-wise error','modal');
        uicontrol(handles.k2_UpperTextbox);
        return
    end
end

if k3check == 1 && k3rest == 1 
    k3low = str2double(get(handles.k3_LowerTextbox,'String'));
    if isnan(k3low) || (k3low < 0)
        warndlg('Lower restrict must be a positive number','Pixel-wise error','modal');
        uicontrol(handles.k3_LowerTextbox);
        return
    end
    
   k3upp = str2double(get(handles.k3_UpperTextbox,'String'));
    if isnan(k3upp) || (k3upp < 0)
        warndlg('Upper restrict must be a positive number','Pixel-wise error','modal');
        uicontrol(handles.k3_UpperTextbox);
        return
    end
    
    if k3upp < k3low
        warndlg('Upper restrict must be bigger than lower restrict','Pixel-wise error','modal');
        uicontrol(handles.k3_UpperTextbox);
        return
    end
end

if k4check == 1 && k4rest == 1 
    k4low = str2double(get(handles.k4_LowerTextbox,'String'));
    if isnan(k4low) || (k4low < 0)
        warndlg('Lower restrict must be a positive number','Pixel-wise error','modal');
        uicontrol(handles.k4_LowerTextbox);
        return
    end
    
   k4upp = str2double(get(handles.k4_UpperTextbox,'String'));
    if isnan(k4upp) || (k4upp < 0)
        warndlg('Upper restrict must be a positive number','Pixel-wise error','modal');
        uicontrol(handles.k4_UpperTextbox);
        return
    end
    
    if k4upp < k4low
        warndlg('Upper restrict must be bigger than lower restrict','Pixel-wise error','modal');
        uicontrol(handles.k4_UpperTextbox);
        return
    end
end

if Kicheck == 1 && Kirest == 1 
    Kilow = str2double(get(handles.Ki_LowerTextbox,'String'));
    if isnan(Kilow) || (Kilow < 0)
        warndlg('Lower restrict must be a positive number','Pixel-wise error','modal');
        uicontrol(handles.Ki_LowerTextbox);
        return
    end
    
   Kiupp = str2double(get(handles.Ki_UpperTextbox,'String'));
    if isnan(Kiupp) || (Kiupp < 0)
        warndlg('Upper restrict must be a positive number','Pixel-wise error','modal');
        uicontrol(handles.Ki_UpperTextbox);
        return
    end
    
    if Kiupp < Kilow
        warndlg('Upper restrict must be bigger than lower restrict','Pixel-wise error','modal');
        uicontrol(handles.Ki_UpperTextbox);
        return
    end
end

if Vscheck == 1 && Vsrest == 1 
    Vslow = str2double(get(handles.Vs_LowerTextbox,'String'));
    if isnan(Vslow) || (Vslow < 0)
        warndlg('Lower restrict must be a positive number','Pixel-wise error','modal');
        uicontrol(handles.Vs_LowerTextbox);
        return
    end
    
   Vsupp = str2double(get(handles.Vs_UpperTextbox,'String'));
    if isnan(Vsupp) || (Vsupp < 0)
        warndlg('Upper restrict must be a positive number','Pixel-wise error','modal');
        uicontrol(handles.Vs_UpperTextbox);
        return
    end
    
    if Vsupp < Vslow
        warndlg('Upper restrict must be bigger than lower restrict','Pixel-wise error','modal');
        uicontrol(handles.Vs_UpperTextbox);
        return
    end
end

if Vtcheck == 1 && Vtrest == 1 
    Vtlow = str2double(get(handles.Vt_LowerTextbox,'String'));
    if isnan(Vtlow) || (Vtlow < 0)
        warndlg('Lower restrict must be a positive number','Pixel-wise error','modal');
        uicontrol(handles.Vt_LowerTextbox);
        return
    end
    
   Vtupp = str2double(get(handles.Vt_UpperTextbox,'String'));
    if isnan(Vtupp) || (Vtupp < 0)
        warndlg('Upper restrict must be a positive number','Pixel-wise error','modal');
        uicontrol(handles.Vt_UpperTextbox);
        return
    end
    
    if Vtupp < Vtlow
        warndlg('Upper restrict must be bigger than lower restrict','Pixel-wise error','modal');
        uicontrol(handles.Vt_UpperTextbox);
        return
    end
end

if K1divk2check == 1 && K1divk2rest == 1 
    K1divk2low = str2double(get(handles.K1divk2_LowerTextbox,'String'));
    if isnan(K1divk2low) || (K1divk2low < 0)
        warndlg('Lower restrict must be a positive number','Pixel-wise error','modal');
        uicontrol(handles.K1divk2_LowerTextbox);
        return
    end
    
   K1divk2upp = str2double(get(handles.K1divk2_UpperTextbox,'String'));
    if isnan(K1divk2upp) || (K1divk2upp < 0)
        warndlg('Upper restrict must be a positive number','Pixel-wise error','modal');
        uicontrol(handles.K1divk2_UpperTextbox);
        return
    end
    
    if K1divk2upp < K1divk2low
        warndlg('Upper restrict must be bigger than lower restrict','Pixel-wise error','modal');
        uicontrol(handles.K1divk2_UpperTextbox);
        return
    end
end

if k3divk4check == 1 && k3divk4rest == 1 
    k3divk4low = str2double(get(handles.k3divk4_LowerTextbox,'String'));
    if isnan(k3divk4low) || (k3divk4low < 0)
        warndlg('Lower restrict must be a positive number','Pixel-wise error','modal');
        uicontrol(handles.k3divk4_LowerTextbox);
        return
    end
    
   k3divk4upp = str2double(get(handles.k3divk4_UpperTextbox,'String'));
    if isnan(k3divk4upp) || (k3divk4upp < 0)
        warndlg('Upper restrict must be a positive number','Pixel-wise error','modal');
        uicontrol(handles.k3divk4_UpperTextbox);
        return
    end
    
    if k3divk4upp < k3divk4low
        warndlg('Upper restrict must be bigger than lower restrict','Pixel-wise error','modal');
        uicontrol(handles.k3divk4_UpperTextbox);
        return
    end
end

if alpha1check == 1 && alpha1rest == 1 
    alpha1low = str2double(get(handles.alpha1_LowerTextbox,'String'));
    if isnan(alpha1low) || (alpha1low < 0)
        warndlg('Lower restrict must be a positive number','Pixel-wise error','modal');
        uicontrol(handles.alpha1_LowerTextbox);
        return
    end
    
   alpha1upp = str2double(get(handles.alpha1_UpperTextbox,'String'));
    if isnan(alpha1upp) || (alpha1upp < 0)
        warndlg('Upper restrict must be a positive number','Pixel-wise error','modal');
        uicontrol(handles.alpha1_UpperTextbox);
        return
    end
    
    if alpha1upp < alpha1low
        warndlg('Upper restrict must be bigger than lower restrict','Pixel-wise error','modal');
        uicontrol(handles.alpha1_UpperTextbox);
        return
    end
end

if alpha2check == 1 && alpha2rest == 1 
    alpha2low = str2double(get(handles.alpha2_LowerTextbox,'String'));
    if isnan(alpha2low) || (alpha2low < 0)
        warndlg('Lower restrict must be a positive number','Pixel-wise error','modal');
        uicontrol(handles.alpha2_LowerTextbox);
        return
    end
    
   alpha2upp = str2double(get(handles.alpha2_UpperTextbox,'String'));
    if isnan(alpha2upp) || (alpha2upp < 0)
        warndlg('Upper restrict must be a positive number','Pixel-wise error','modal');
        uicontrol(handles.alpha2_UpperTextbox);
        return
    end
    
    if alpha2upp < alpha2low
        warndlg('Upper restrict must be bigger than lower restrict','Pixel-wise error','modal');
        uicontrol(handles.alpha2_UpperTextbox);
        return
    end
end

if ttpcheck == 1 && ttprest == 1 
    ttplow = str2double(get(handles.ttp_LowerTextbox,'String'));
    if isnan(ttplow) || (ttplow < 0)
        warndlg('Lower restrict must be a positive number','Pixel-wise error','modal');
        uicontrol(handles.ttp_LowerTextbox);
        return
    end
    
   ttpupp = str2double(get(handles.ttp_UpperTextbox,'String'));
    if isnan(ttpupp) || (ttpupp < 0)
        warndlg('Upper restrict must be a positive number','Pixel-wise error','modal');
        uicontrol(handles.ttp_UpperTextbox);
        return
    end
    
    if ttpupp < ttplow
        warndlg('Upper restrict must be bigger than lower restrict','Pixel-wise error','modal');
        uicontrol(handles.ttp_UpperTextbox);
        return
    end
end

%Update the logfile
aux_date=datestr(now);
fprintf(QMf_logfile,strcat('\n-2TCM-\n'));
fprintf(QMf_logfile,strcat(aux_date(end-8:end),' Generating parametric images (BF= ',num2str(length(QMlastPreprocess.TwoTCM.alpha(1))),', threshold= ',num2str(QMlastPreprocess.TwoTCM.threshold),')\n'));

%Image calculation
images = QM_TwoTCMPixelCalc([vBrest vBlow vBupp vBcheck; K1rest K1low K1upp K1check; k2rest k2low k2upp k2check; k3rest k3low k3upp k3check; k4rest k4low k4upp k4check;...
                             Kirest Kilow Kiupp Kicheck; Vsrest Vslow Vsupp Vscheck; Vtrest Vtlow Vtupp Vtcheck; K1divk2rest K1divk2low K1divk2upp K1divk2check;...
                             k3divk4rest k3divk4low k3divk4upp k3divk4check; alpha1rest alpha1low alpha1upp alpha1check; alpha2rest alpha2low alpha2upp alpha2check; ...
                             ttprest ttplow ttpupp ttpcheck], waitbarhandle);
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

if vBcheck == 1 
    vBfile = get(handles.vB_Filename_Textbox,'String');
    
    if get(handles.vB_FileExt_Popupmenu,'Value') == 1
        vBfileExt = '.img';
    else
        vBfileExt = '.nii';
    end
    
    VvB = struct('fname',strcat(vBfile,vBfileExt),'dim',Vtemp.dim,'mat',Vtemp.mat,...
            'dt',[spm_type('float64') spm_platform('bigend')],'n',[1 1]);

    %Create image
    aux = zeros(VvB.dim);
    aux(:) = images(1,:);
    
    if get(handles.save_vB_Checkbox,'Value') == 1
        VvB.fname = strcat(outputDir,filesep,vBfile,vBfileExt);
        spm_write_vol(VvB,aux);
        fprintf(QMf_logfile,strcat(['  -vB image saved in ',char(32),outputDir,'\n']));
    end
    VvB.fname = strcat(QMpath,filesep,'temp',filesep,vBfile,vBfileExt);
    spm_write_vol(VvB,aux); %Save image in hard disk
        
    imgTempPaths{1,indTemp} = VvB.fname;
    imgTempFilenames{1,indTemp} = 'vB (Blood volume fraction)';
    indTemp=indTemp+1;
end

if K1check == 1 
    K1file = get(handles.K1_Filename_Textbox,'String');
    
    if get(handles.K1_FileExt_Popupmenu,'Value') == 1
        K1fileExt = '.img';
    else
        K1fileExt = '.nii';
    end
    
    VK1 = struct('fname',strcat(K1file,K1fileExt),'dim',Vtemp.dim,'mat',Vtemp.mat,...
            'dt',[spm_type('float64') spm_platform('bigend')],'n',[1 1]);

    %Create image
    aux = zeros(VK1.dim);
    aux(:) = images(2,:);
    
    if get(handles.save_K1_Checkbox,'Value') == 1
        VK1.fname = strcat(outputDir,filesep,K1file,K1fileExt);
        spm_write_vol(VK1,aux);
        fprintf(QMf_logfile,strcat(['  -K1 image saved in ',char(32),outputDir,'\n']));
    end
    VK1.fname = strcat(QMpath,filesep,'temp',filesep,K1file,K1fileExt);
    spm_write_vol(VK1,aux); %Save image in hard disk
        
    imgTempPaths{1,indTemp} = VK1.fname;
    imgTempFilenames{1,indTemp} = 'K1 (Rate constant)';
    indTemp=indTemp+1;
end

if k2check == 1 
    k2file = get(handles.k2_Filename_Textbox,'String');
    
    if get(handles.k2_FileExt_Popupmenu,'Value') == 1
        k2fileExt = '.img';
    else
        k2fileExt = '.nii';
    end
    
    Vk2 = struct('fname',strcat(k2file,k2fileExt),'dim',Vtemp.dim,'mat',Vtemp.mat,...
            'dt',[spm_type('float64') spm_platform('bigend')],'n',[1 1]);

    %Create image
    aux = zeros(Vk2.dim);
    aux(:) = images(3,:);
    
    if get(handles.save_k2_Checkbox,'Value') == 1
        Vk2.fname = strcat(outputDir,filesep,k2file,k2fileExt);
        spm_write_vol(Vk2,aux);
        fprintf(QMf_logfile,strcat(['  -k2 image saved in ',char(32),outputDir,'\n']));
    end
    Vk2.fname = strcat(QMpath,filesep,'temp',filesep,k2file,k2fileExt);
    spm_write_vol(Vk2,aux); %Save image in hard disk
        
    imgTempPaths{1,indTemp} = Vk2.fname;
    imgTempFilenames{1,indTemp} = 'k2 (Rate constant)';
    indTemp=indTemp+1;
end

if k3check == 1 
    k3file = get(handles.k3_Filename_Textbox,'String');
    
    if get(handles.k3_FileExt_Popupmenu,'Value') == 1
        k3fileExt = '.img';
    else
        k3fileExt = '.nii';
    end
    
    Vk3 = struct('fname',strcat(k3file,k3fileExt),'dim',Vtemp.dim,'mat',Vtemp.mat,...
            'dt',[spm_type('float64') spm_platform('bigend')],'n',[1 1]);

    %Create image
    aux = zeros(Vk3.dim);
    aux(:) = images(4,:);
    
    if get(handles.save_k3_Checkbox,'Value') == 1
        Vk3.fname = strcat(outputDir,filesep,k3file,k3fileExt);
        spm_write_vol(Vk3,aux);
        fprintf(QMf_logfile,strcat(['  -k3 image saved in ',char(32),outputDir,'\n']));
    end
    Vk3.fname = strcat(QMpath,filesep,'temp',filesep,k3file,k3fileExt);
    spm_write_vol(Vk3,aux); %Save image in hard disk
        
    imgTempPaths{1,indTemp} = Vk3.fname;
    imgTempFilenames{1,indTemp} = 'k3 (Rate constant)';
    indTemp=indTemp+1;
end

if k4check == 1 
    k4file = get(handles.k4_Filename_Textbox,'String');
    
    if get(handles.k4_FileExt_Popupmenu,'Value') == 1
        k4fileExt = '.img';
    else
        k4fileExt = '.nii';
    end
    
    Vk4 = struct('fname',strcat(k4file,k4fileExt),'dim',Vtemp.dim,'mat',Vtemp.mat,...
            'dt',[spm_type('float64') spm_platform('bigend')],'n',[1 1]);

    %Create image
    aux = zeros(Vk4.dim);
    aux(:) = images(5,:);
    
    if get(handles.save_k4_Checkbox,'Value') == 1
        Vk4.fname = strcat(outputDir,filesep,k4file,k4fileExt);
        spm_write_vol(Vk4,aux);
        fprintf(QMf_logfile,strcat(['  -k4 image saved in ',char(32),outputDir,'\n']));
    end
    Vk4.fname = strcat(QMpath,filesep,'temp',filesep,k4file,k4fileExt);
    spm_write_vol(Vk4,aux); %Save image in hard disk
        
    imgTempPaths{1,indTemp} = Vk4.fname;
    imgTempFilenames{1,indTemp} = 'k4 (Rate constant)';
    indTemp=indTemp+1;
end

if Kicheck == 1 
    Kifile = get(handles.Ki_Filename_Textbox,'String');
    
    if get(handles.Ki_FileExt_Popupmenu,'Value') == 1
        KifileExt = '.img';
    else
        KifileExt = '.nii';
    end
    
    VKi = struct('fname',strcat(Kifile,KifileExt),'dim',Vtemp.dim,'mat',Vtemp.mat,...
            'dt',[spm_type('float64') spm_platform('bigend')],'n',[1 1]);

    %Create image
    aux = zeros(VKi.dim);
    aux(:) = images(6,:);
    
    if get(handles.save_Ki_Checkbox,'Value') == 1
        VKi.fname = strcat(outputDir,filesep,Kifile,KifileExt);
        spm_write_vol(VKi,aux);
        fprintf(QMf_logfile,strcat(['  -Ki or Flux image saved in ',char(32),outputDir,'\n']));
    end
    VKi.fname = strcat(QMpath,filesep,'temp',filesep,Kifile,KifileExt);
    spm_write_vol(VKi,aux); %Save image in hard disk
        
    imgTempPaths{1,indTemp} = VKi.fname;
    imgTempFilenames{1,indTemp} = 'Ki or Flux (Influx of the racer)';
    indTemp=indTemp+1;
end

if Vscheck == 1 
    Vsfile = get(handles.Vs_Filename_Textbox,'String');
    
    if get(handles.Vs_FileExt_Popupmenu,'Value') == 1
        VsfileExt = '.img';
    else
        VsfileExt = '.nii';
    end
    
    VVs = struct('fname',strcat(Vsfile,VsfileExt),'dim',Vtemp.dim,'mat',Vtemp.mat,...
            'dt',[spm_type('float64') spm_platform('bigend')],'n',[1 1]);

    %Create image
    aux = zeros(VVs.dim);
    aux(:) = images(7,:);
    
    if get(handles.save_Vs_Checkbox,'Value') == 1
        VVs.fname = strcat(outputDir,filesep,Vsfile,VsfileExt);
        spm_write_vol(VVs,aux);
        fprintf(QMf_logfile,strcat(['  -Vs image saved in ',char(32),outputDir,'\n']));
    end
    VVs.fname = strcat(QMpath,filesep,'temp',filesep,Vsfile,VsfileExt);
    spm_write_vol(VVs,aux); %Save image in hard disk
        
    imgTempPaths{1,indTemp} = VVs.fname;
    imgTempFilenames{1,indTemp} = 'Vs (Distribution volume of the second compartment)';
    indTemp=indTemp+1;
end

if Vtcheck == 1 
    Vtfile = get(handles.Vt_Filename_Textbox,'String');
    
    if get(handles.Vt_FileExt_Popupmenu,'Value') == 1
        VtfileExt = '.img';
    else
        VtfileExt = '.nii';
    end
    
    VVt = struct('fname',strcat(Vtfile,VtfileExt),'dim',Vtemp.dim,'mat',Vtemp.mat,...
            'dt',[spm_type('float64') spm_platform('bigend')],'n',[1 1]);

    %Create image
    aux = zeros(VVt.dim);
    aux(:) = images(8,:);
    
    if get(handles.save_Vt_Checkbox,'Value') == 1
        VVt.fname = strcat(outputDir,filesep,Vtfile,VtfileExt);
        spm_write_vol(VVt,aux);
        fprintf(QMf_logfile,strcat(['  -Vt image saved in ',char(32),outputDir,'\n']));
    end
    VVt.fname = strcat(QMpath,filesep,'temp',filesep,Vtfile,VtfileExt);
    spm_write_vol(VVt,aux); %Save image in hard disk
        
    imgTempPaths{1,indTemp} = VVt.fname;
    imgTempFilenames{1,indTemp} = 'Vt (Distribution volume)';
    indTemp=indTemp+1;
end

if K1divk2check == 1 
    K1divk2file = get(handles.K1divk2_Filename_Textbox,'String');
    
    if get(handles.K1divk2_FileExt_Popupmenu,'Value') == 1
        K1divk2fileExt = '.img';
    else
        K1divk2fileExt = '.nii';
    end
    
    VK1divk2 = struct('fname',strcat(K1divk2file,K1divk2fileExt),'dim',Vtemp.dim,'mat',Vtemp.mat,...
            'dt',[spm_type('float64') spm_platform('bigend')],'n',[1 1]);

    %Create image
    aux = zeros(VK1divk2.dim);
    aux(:) = images(9,:);
    
    if get(handles.save_K1divk2_Checkbox,'Value') == 1
        VK1divk2.fname = strcat(outputDir,filesep,K1divk2file,K1divk2fileExt);
        spm_write_vol(VK1divk2,aux);
        fprintf(QMf_logfile,strcat(['  -K1/k2 image saved in ',char(32),outputDir,'\n']));
    end
    VK1divk2.fname = strcat(QMpath,filesep,'temp',filesep,K1divk2file,K1divk2fileExt);
    spm_write_vol(VK1divk2,aux); %Save image in hard disk
        
    imgTempPaths{1,indTemp} = VK1divk2.fname;
    imgTempFilenames{1,indTemp} = 'K1/k2 (Distribution volume of the non-displaceable compartment)';
    indTemp=indTemp+1;
end

if k3divk4check == 1 
    k3divk4file = get(handles.k3divk4_Filename_Textbox,'String');
    
    if get(handles.k3divk4_FileExt_Popupmenu,'Value') == 1
        k3divk4fileExt = '.img';
    else
        k3divk4fileExt = '.nii';
    end
    
    Vk3divk4 = struct('fname',strcat(k3divk4file,k3divk4fileExt),'dim',Vtemp.dim,'mat',Vtemp.mat,...
            'dt',[spm_type('float64') spm_platform('bigend')],'n',[1 1]);

    %Create image
    aux = zeros(Vk3divk4.dim);
    aux(:) = images(10,:);
    
    if get(handles.save_k3divk4_Checkbox,'Value') == 1
        Vk3divk4.fname = strcat(outputDir,filesep,k3divk4file,k3divk4fileExt);
        spm_write_vol(Vk3divk4,aux);
        fprintf(QMf_logfile,strcat(['  -k3/k4 image saved in ',char(32),outputDir,'\n']));
    end
    Vk3divk4.fname = strcat(QMpath,filesep,'temp',filesep,k3divk4file,k3divk4fileExt);
    spm_write_vol(Vk3divk4,aux); %Save image in hard disk
        
    imgTempPaths{1,indTemp} = Vk3divk4.fname;
    imgTempFilenames{1,indTemp} = 'k3/k4 (Binding potential of receptor tracers)';
    indTemp=indTemp+1;
end

if alpha1check == 1 
    alpha1file = get(handles.alpha1_Filename_Textbox,'String');
    
    if get(handles.alpha1_FileExt_Popupmenu,'Value') == 1
        alpha1fileExt = '.img';
    else
        alpha1fileExt = '.nii';
    end
    
    Valpha1 = struct('fname',strcat(alpha1file,alpha1fileExt),'dim',Vtemp.dim,'mat',Vtemp.mat,...
            'dt',[spm_type('float64') spm_platform('bigend')],'n',[1 1]);

    %Create image
    aux = zeros(Valpha1.dim);
    aux(:) = images(11,:);
    
    if get(handles.save_alpha1_Checkbox,'Value') == 1
        Valpha1.fname = strcat(outputDir,filesep,alpha1file,alpha1fileExt);
        spm_write_vol(Valpha1,aux);
        fprintf(QMf_logfile,strcat(['  -alpha1 image saved in ',char(32),outputDir,'\n']));
    end
    Valpha1.fname = strcat(QMpath,filesep,'temp',filesep,alpha1file,alpha1fileExt);
    spm_write_vol(Valpha1,aux); %Save image in hard disk
        
    imgTempPaths{1,indTemp} = Valpha1.fname;
    imgTempFilenames{1,indTemp} = 'alpha1';
    indTemp=indTemp+1;
end

if alpha2check == 1 
    alpha2file = get(handles.alpha2_Filename_Textbox,'String');
    
    if get(handles.alpha2_FileExt_Popupmenu,'Value') == 1
        alpha2fileExt = '.img';
    else
        alpha2fileExt = '.nii';
    end
    
    Valpha2 = struct('fname',strcat(alpha2file,alpha2fileExt),'dim',Vtemp.dim,'mat',Vtemp.mat,...
            'dt',[spm_type('float64') spm_platform('bigend')],'n',[1 1]);

    %Create image
    aux = zeros(Valpha2.dim);
    aux(:) = images(12,:);
    
    if get(handles.save_alpha2_Checkbox,'Value') == 1
        Valpha2.fname = strcat(outputDir,filesep,alpha2file,alpha2fileExt);
        spm_write_vol(Valpha2,aux);
        fprintf(QMf_logfile,strcat(['  -alpha2 image saved in ',char(32),outputDir,'\n']));
    end
    Valpha2.fname = strcat(QMpath,filesep,'temp',filesep,alpha2file,alpha2fileExt);
    spm_write_vol(Valpha2,aux); %Save image in hard disk
        
    imgTempPaths{1,indTemp} = Valpha2.fname;
    imgTempFilenames{1,indTemp} = 'alpha2';  
    indTemp=indTemp+1;
end

if ttpcheck == 1 %image time to peak
    ttpfile = get(handles.ttp_Filename_Textbox,'String');    
    
    if get(handles.ttp_FileExt_Popupmenu,'Value') == 1
        ttpfileExt = '.img';
    else
        ttpfileExt = '.nii';
    end
    
    Vttp = struct('fname',strcat(ttpfile,ttpfileExt),'dim',Vtemp.dim,'mat',Vtemp.mat,...
            'dt',[spm_type('float64') spm_platform('bigend')],'n',[1 1]);

    %Create image
    aux = zeros(Vttp.dim);
    aux(:) = images(13,:);
    
    if get(handles.save_ttp_Checkbox,'Value') == 1
        Vttp.fname = strcat(outputDir,filesep,ttpfile,ttpfileExt);
        spm_write_vol(Vttp,aux);
        fprintf(QMf_logfile,strcat(['  -ttp image saved in ',char(32),outputDir,'\n']));
    end
    Vttp.fname = strcat(QMpath,filesep,'temp',filesep,ttpfile,ttpfileExt);
    spm_write_vol(Vttp,aux); %Save image in hard disk
        
    imgTempPaths{1,indTemp} = Vttp.fname;
    imgTempFilenames{1,indTemp} = 'Time to peak';
    
end

delete(waitbarhandle);  %Need to delete before ploting
set(mainhandles.selectImage_Popupmenu,'UserData',imgTempPaths);
set(mainhandles.selectImage_Popupmenu,'String',imgTempFilenames);

set(mainhandles.setMapParam_Pushbutton,'UserData',1);

uiresume(handles.MapParTwoTCM_figure);

% --- Executes on button press in setDefault_Pushbutton.
function setDefault_Pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to setDefault_Pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global QMlastPreprocess
% Update components
set(handles.vB_LowerTextbox,'String','0.0','Enable','on');
set(handles.vB_UpperTextbox,'String','1.0','Enable','on');
set(handles.K1_LowerTextbox,'String','0.0','Enable','on');
set(handles.K1_UpperTextbox,'String','2.0','Enable','on');
set(handles.k2_LowerTextbox,'String','0.0','Enable','on');
set(handles.k2_UpperTextbox,'String','2.0','Enable','on');
set(handles.k3_LowerTextbox,'String','0.0','Enable','on');
set(handles.k3_UpperTextbox,'String','1.0','Enable','on');
set(handles.k4_LowerTextbox,'String','0.0','Enable','on');
set(handles.k4_UpperTextbox,'String','1.0','Enable','on');
set(handles.Ki_LowerTextbox,'String','0.0','Enable','off');
set(handles.Ki_UpperTextbox,'String','30.0','Enable','off');
set(handles.Vs_LowerTextbox,'String','0.0','Enable','on');
set(handles.Vs_UpperTextbox,'String','8.0','Enable','on');
set(handles.Vt_LowerTextbox,'String','0.0','Enable','off');
set(handles.Vt_UpperTextbox,'String','10.0','Enable','off');
set(handles.K1divk2_LowerTextbox,'String','0.0','Enable','off');
set(handles.K1divk2_UpperTextbox,'String','10.0','Enable','off');
set(handles.k3divk4_LowerTextbox,'String','0.0','Enable','off');
set(handles.k3divk4_UpperTextbox,'String','10.0','Enable','off');
set(handles.alpha1_LowerTextbox,'String','0.0','Enable','off');
set(handles.alpha1_UpperTextbox,'String','0.1','Enable','off');
set(handles.alpha2_LowerTextbox,'String','0.0','Enable','off');
set(handles.alpha2_UpperTextbox,'String','1.0','Enable','off');
set(handles.ttp_LowerTextbox,'String','0.0','Enable','off');
set(handles.ttp_UpperTextbox,'String','60.0','Enable','off');

set(handles.vB_Filename_Textbox,'String','vB_Image','Enable','off');
set(handles.K1_Filename_Textbox,'String','K1_Image','Enable','off');
set(handles.k2_Filename_Textbox,'String','k2_Image','Enable','off');
set(handles.k3_Filename_Textbox,'String','k3_Image','Enable','off');
set(handles.k4_Filename_Textbox,'String','k4_Image','Enable','off');
set(handles.Ki_Filename_Textbox,'String','Ki_Image','Enable','off');
set(handles.Vs_Filename_Textbox,'String','Vs_Image','Enable','off');
set(handles.Vt_Filename_Textbox,'String','Vt_Image','Enable','off');
set(handles.K1divk2_Filename_Textbox,'String','K1overk2_Image','Enable','off');
set(handles.k3divk4_Filename_Textbox,'String','k3overk4_Image','Enable','off');
set(handles.alpha1_Filename_Textbox,'String','alpha1_Image','Enable','off');
set(handles.alpha2_Filename_Textbox,'String','alpha2_Image','Enable','off');
set(handles.ttp_Filename_Textbox,'String','time_to_peak_Image','Enable','off');

set(handles.vB_Checkbox,'Value',1.0);
set(handles.K1_Checkbox,'Value',1.0);
set(handles.k2_Checkbox,'Value',1.0);
set(handles.k3_Checkbox,'Value',1.0);
set(handles.k4_Checkbox,'Value',1.0);
set(handles.Vs_Checkbox,'Value',1.0);
set(handles.Ki_Checkbox,'Value',0.0);
set(handles.Vt_Checkbox,'Value',0.0);
set(handles.K1divk2_Checkbox,'Value',0.0);
set(handles.k3divk4_Checkbox,'Value',0.0);
set(handles.alpha1_Checkbox,'Value',0.0);
set(handles.alpha2_Checkbox,'Value',0.0);
set(handles.ttp_Checkbox,'Value',0.0);

set(handles.restrict_vB_Checkbox,'Enable','on','Value',1.0);
set(handles.restrict_K1_Checkbox,'Enable','on','Value',1.0);
set(handles.restrict_k2_Checkbox,'Enable','on','Value',1.0);
set(handles.restrict_k3_Checkbox,'Enable','on','Value',1.0);
set(handles.restrict_k4_Checkbox,'Enable','on','Value',1.0);
set(handles.restrict_Vs_Checkbox,'Enable','on','Value',1.0);
set(handles.restrict_Ki_Checkbox,'Enable','off','Value',1.0);
set(handles.restrict_Vt_Checkbox,'Enable','off','Value',1.0);
set(handles.restrict_K1divk2_Checkbox,'Enable','off','Value',1.0);
set(handles.restrict_k3divk4_Checkbox,'Enable','off','Value',1.0);
set(handles.restrict_alpha1_Checkbox,'Enable','off','Value',1.0);
set(handles.restrict_alpha2_Checkbox,'Enable','off','Value',1.0);
set(handles.restrict_ttp_Checkbox,'Enable','off','Value',1.0);

set(handles.save_vB_Checkbox,'Enable','on','Value',0.0);
set(handles.save_K1_Checkbox,'Enable','on','Value',0.0);
set(handles.save_k2_Checkbox,'Enable','on','Value',0.0);
set(handles.save_k3_Checkbox,'Enable','on','Value',0.0);
set(handles.save_k4_Checkbox,'Enable','on','Value',0.0);
set(handles.save_Ki_Checkbox,'Enable','off','Value',0.0);
set(handles.save_Vs_Checkbox,'Enable','on','Value',0.0);
set(handles.save_Vt_Checkbox,'Enable','off','Value',0.0);
set(handles.save_K1divk2_Checkbox,'Enable','off','Value',0.0);
set(handles.save_k3divk4_Checkbox,'Enable','off','Value',0.0);
set(handles.save_alpha1_Checkbox,'Enable','off','Value',0.0);
set(handles.save_alpha2_Checkbox,'Enable','off','Value',0.0);
set(handles.save_ttp_Checkbox,'Enable','off','Value',0.0);

set(handles.vB_FileExt_Popupmenu,'Value',1,'Enable','off');
set(handles.K1_FileExt_Popupmenu,'Value',1,'Enable','off');
set(handles.k2_FileExt_Popupmenu,'Value',1,'Enable','off');
set(handles.k3_FileExt_Popupmenu,'Value',1,'Enable','off');
set(handles.k4_FileExt_Popupmenu,'Value',1,'Enable','off');
set(handles.Ki_FileExt_Popupmenu,'Value',1,'Enable','off');
set(handles.Vs_FileExt_Popupmenu,'Value',1,'Enable','off');
set(handles.Vt_FileExt_Popupmenu,'Value',1,'Enable','off');
set(handles.K1divk2_FileExt_Popupmenu,'Value',1,'Enable','off');
set(handles.k3divk4_FileExt_Popupmenu,'Value',1,'Enable','off');
set(handles.alpha1_FileExt_Popupmenu,'Value',1,'Enable','off');
set(handles.alpha2_FileExt_Popupmenu,'Value',1,'Enable','off');
set(handles.ttp_FileExt_Popupmenu,'Value',1,'Enable','off');

dir = pwd;
set(handles.outputDir_Textbox,'String',dir);
set(handles.runPixelCalc_Pushbutton,'Enable','on')

if ~QMlastPreprocess.TwoTCM.k4_checked %k4 fixed
    set(handles.k4_Checkbox,'Value',0.0);
    set(handles.k4_Checkbox,'Enable','off');
    set(handles.restrict_k4_Checkbox,'Enable','off');
    set(handles.save_k4_Checkbox,'Enable','off');
    set(handles.k4_Filename_Textbox,'Enable','off');
    set(handles.k4_FileExt_Popupmenu,'Enable','off');
    set(handles.k4_LowerTextbox,'Enable','off');
    set(handles.k4_UpperTextbox,'Enable','off');
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
set(handles.MapParTwoTCM_figure,'WindowStyle','normal');
directoryname = uigetdir('', 'Select the output directory');
set(handles.MapParTwoTCM_figure,'WindowStyle','modal');

if directoryname ~= 0
     set(handles.outputDir_Textbox,'String',directoryname);
end


% --- Executes on button press in FAIL_vB_Checkbox.
function vB_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to FAIL_vB_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of FAIL_vB_Checkbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
    set(handles.runPixelCalc_Pushbutton,'Enable',enable)
else
    enable = 'off';
    if (~get(handles.vB_Checkbox,'Value') && ~get(handles.K1_Checkbox,'Value') && ~get(handles.k2_Checkbox,'Value')...
            && ~get(handles.k3_Checkbox,'Value') && ~get(handles.k4_Checkbox,'Value') && ~get(handles.Ki_Checkbox,'Value')...
            && ~get(handles.Vs_Checkbox,'Value') && ~get(handles.Vt_Checkbox,'Value') && ~get(handles.K1divk2_Checkbox,'Value')...
            && ~get(handles.k3divk4_Checkbox,'Value') && ~get(handles.alpha1_Checkbox,'Value') && ~get(handles.alpha2_Checkbox,'Value'))    
        set(handles.runPixelCalc_Pushbutton,'Enable','off')
    end
end
set(handles.restrict_vB_Checkbox,'Enable',enable);
set(handles.save_vB_Checkbox,'Enable',enable);
set(handles.vB_Filename_Textbox,'Enable','off');
set(handles.vB_FileExt_Popupmenu,'Enable','off');
set(handles.vB_LowerTextbox,'Enable',enable);
set(handles.vB_UpperTextbox,'Enable',enable);


% --- Executes on button press in restrict_vB_Checkbox.
function restrict_vB_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to restrict_vB_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of restrict_vB_Checkbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
else
    enable = 'off';
end
set(handles.vB_LowerTextbox,'Enable',enable);
set(handles.vB_UpperTextbox,'Enable',enable);


function vB_LowerTextbox_Callback(hObject, eventdata, handles)
% hObject    handle to vB_LowerTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vB_LowerTextbox as text
%        str2double(get(hObject,'String')) returns contents of vB_LowerTextbox as a double
value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('Lower value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function vB_LowerTextbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vB_LowerTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function vB_UpperTextbox_Callback(hObject, eventdata, handles)
% hObject    handle to vB_UpperTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vB_UpperTextbox as text
%        str2double(get(hObject,'String')) returns contents of vB_UpperTextbox as a double
value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('Upper value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function vB_UpperTextbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vB_UpperTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in save_vB_Checkbox.
function save_vB_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to save_vB_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of save_vB_Checkbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
else
    enable = 'off';
end
set(handles.vB_Filename_Textbox,'Enable',enable);
set(handles.vB_FileExt_Popupmenu,'Enable',enable);


function vB_Filename_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to vB_Filename_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vB_Filename_Textbox as text
%        str2double(get(hObject,'String')) returns contents of vB_Filename_Textbox as a double


% --- Executes during object creation, after setting all properties.
function vB_Filename_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vB_Filename_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in vB_FileExt_Popupmenu.
function vB_FileExt_Popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to vB_FileExt_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns vB_FileExt_Popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from vB_FileExt_Popupmenu


% --- Executes during object creation, after setting all properties.
function vB_FileExt_Popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vB_FileExt_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in K1_Checkbox.
function K1_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to K1_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of K1_Checkbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
    set(handles.runPixelCalc_Pushbutton,'Enable',enable)
else
    enable = 'off';
    if (~get(handles.vB_Checkbox,'Value') && ~get(handles.K1_Checkbox,'Value') && ~get(handles.k2_Checkbox,'Value')...
            && ~get(handles.k3_Checkbox,'Value') && ~get(handles.k4_Checkbox,'Value') && ~get(handles.Ki_Checkbox,'Value')...
            && ~get(handles.Vs_Checkbox,'Value') && ~get(handles.Vt_Checkbox,'Value') && ~get(handles.K1divk2_Checkbox,'Value')...
            && ~get(handles.k3divk4_Checkbox,'Value') && ~get(handles.alpha1_Checkbox,'Value') && ~get(handles.alpha2_Checkbox,'Value'))    
        set(handles.runPixelCalc_Pushbutton,'Enable','off')
    end
end
set(handles.restrict_K1_Checkbox,'Enable',enable);
set(handles.save_K1_Checkbox,'Enable',enable);
set(handles.K1_Filename_Textbox,'Enable','off');
set(handles.K1_FileExt_Popupmenu,'Enable','off');
set(handles.K1_LowerTextbox,'Enable',enable);
set(handles.K1_UpperTextbox,'Enable',enable);

% --- Executes on button press in restrict_K1_Checkbox.
function restrict_K1_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to restrict_K1_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of restrict_K1_Checkbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
else
    enable = 'off';
end
set(handles.K1_LowerTextbox,'Enable',enable);
set(handles.K1_UpperTextbox,'Enable',enable);


function K1_LowerTextbox_Callback(hObject, eventdata, handles)
% hObject    handle to K1_LowerTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of K1_LowerTextbox as text
%        str2double(get(hObject,'String')) returns contents of K1_LowerTextbox as a double
value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('Lower value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function K1_LowerTextbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to K1_LowerTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function K1_UpperTextbox_Callback(hObject, eventdata, handles)
% hObject    handle to K1_UpperTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of K1_UpperTextbox as text
%        str2double(get(hObject,'String')) returns contents of K1_UpperTextbox as a double
value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('Upper value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function K1_UpperTextbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to K1_UpperTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in save_K1_Checkbox.
function save_K1_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to save_K1_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of save_K1_Checkbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
else
    enable = 'off';
end
set(handles.K1_Filename_Textbox,'Enable',enable);
set(handles.K1_FileExt_Popupmenu,'Enable',enable);


function K1_Filename_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to K1_Filename_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of K1_Filename_Textbox as text
%        str2double(get(hObject,'String')) returns contents of K1_Filename_Textbox as a double


% --- Executes during object creation, after setting all properties.
function K1_Filename_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to K1_Filename_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in K1_FileExt_Popupmenu.
function K1_FileExt_Popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to K1_FileExt_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns K1_FileExt_Popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from K1_FileExt_Popupmenu


% --- Executes during object creation, after setting all properties.
function K1_FileExt_Popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to K1_FileExt_Popupmenu (see GCBO)
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
    if (~get(handles.vB_Checkbox,'Value') && ~get(handles.K1_Checkbox,'Value') && ~get(handles.k2_Checkbox,'Value')...
            && ~get(handles.k3_Checkbox,'Value') && ~get(handles.k4_Checkbox,'Value') && ~get(handles.Ki_Checkbox,'Value')...
            && ~get(handles.Vs_Checkbox,'Value') && ~get(handles.Vt_Checkbox,'Value') && ~get(handles.K1divk2_Checkbox,'Value')...
            && ~get(handles.k3divk4_Checkbox,'Value') && ~get(handles.alpha1_Checkbox,'Value') && ~get(handles.alpha2_Checkbox,'Value'))    
        set(handles.runPixelCalc_Pushbutton,'Enable','off')
    end
end
set(handles.restrict_k2_Checkbox,'Enable',enable);
set(handles.save_k2_Checkbox,'Enable',enable);
set(handles.k2_Filename_Textbox,'Enable','off');
set(handles.k2_FileExt_Popupmenu,'Enable','off');
set(handles.k2_LowerTextbox,'Enable',enable);
set(handles.k2_UpperTextbox,'Enable',enable);

% --- Executes on button press in restrict_k2_Checkbox.
function restrict_k2_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to restrict_k2_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of restrict_k2_Checkbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
else
    enable = 'off';
end
set(handles.k2_LowerTextbox,'Enable',enable);
set(handles.k2_UpperTextbox,'Enable',enable);


function k2_LowerTextbox_Callback(hObject, eventdata, handles)
% hObject    handle to k2_LowerTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k2_LowerTextbox as text
%        str2double(get(hObject,'String')) returns contents of k2_LowerTextbox as a double
value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('Lower value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function k2_LowerTextbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k2_LowerTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function k2_UpperTextbox_Callback(hObject, eventdata, handles)
% hObject    handle to k2_UpperTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k2_UpperTextbox as text
%        str2double(get(hObject,'String')) returns contents of k2_UpperTextbox as a double
value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('Upper value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function k2_UpperTextbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k2_UpperTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in save_k2_Checkbox.
function save_k2_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to save_k2_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of save_k2_Checkbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
else
    enable = 'off';
end
set(handles.k2_Filename_Textbox,'Enable',enable);
set(handles.k2_FileExt_Popupmenu,'Enable',enable);


function k2_Filename_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to k2_Filename_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k2_Filename_Textbox as text
%        str2double(get(hObject,'String')) returns contents of k2_Filename_Textbox as a double


% --- Executes during object creation, after setting all properties.
function k2_Filename_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k2_Filename_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in k2_FileExt_Popupmenu.
function k2_FileExt_Popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to k2_FileExt_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns k2_FileExt_Popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from k2_FileExt_Popupmenu


% --- Executes during object creation, after setting all properties.
function k2_FileExt_Popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k2_FileExt_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in k3_Checkbox.
function k3_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to k3_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of k3_Checkbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
    set(handles.runPixelCalc_Pushbutton,'Enable',enable)
else
    enable = 'off';
    if (~get(handles.vB_Checkbox,'Value') && ~get(handles.K1_Checkbox,'Value') && ~get(handles.k2_Checkbox,'Value')...
            && ~get(handles.k3_Checkbox,'Value') && ~get(handles.k4_Checkbox,'Value') && ~get(handles.Ki_Checkbox,'Value')...
            && ~get(handles.Vs_Checkbox,'Value') && ~get(handles.Vt_Checkbox,'Value') && ~get(handles.K1divk2_Checkbox,'Value')...
            && ~get(handles.k3divk4_Checkbox,'Value') && ~get(handles.alpha1_Checkbox,'Value') && ~get(handles.alpha2_Checkbox,'Value'))    
        set(handles.runPixelCalc_Pushbutton,'Enable','off')
    end
end
set(handles.restrict_k3_Checkbox,'Enable',enable);
set(handles.save_k3_Checkbox,'Enable',enable);
set(handles.k3_Filename_Textbox,'Enable','off');
set(handles.k3_FileExt_Popupmenu,'Enable','off');
set(handles.k3_LowerTextbox,'Enable',enable);
set(handles.k3_UpperTextbox,'Enable',enable);

% --- Executes on button press in restrict_k3_Checkbox.
function restrict_k3_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to restrict_k3_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of restrict_k3_Checkbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
else
    enable = 'off';
end
set(handles.k3_LowerTextbox,'Enable',enable);
set(handles.k3_UpperTextbox,'Enable',enable);


function k3_LowerTextbox_Callback(hObject, eventdata, handles)
% hObject    handle to k3_LowerTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k3_LowerTextbox as text
%        str2double(get(hObject,'String')) returns contents of k3_LowerTextbox as a double
value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('Lower value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function k3_LowerTextbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k3_LowerTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function k3_UpperTextbox_Callback(hObject, eventdata, handles)
% hObject    handle to k3_UpperTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k3_UpperTextbox as text
%        str2double(get(hObject,'String')) returns contents of k3_UpperTextbox as a double
value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('Upper value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function k3_UpperTextbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k3_UpperTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in save_k3_Checkbox.
function save_k3_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to save_k3_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of save_k3_Checkbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
else
    enable = 'off';
end
set(handles.k3_Filename_Textbox,'Enable',enable);
set(handles.k3_FileExt_Popupmenu,'Enable',enable);


function k3_Filename_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to k3_Filename_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k3_Filename_Textbox as text
%        str2double(get(hObject,'String')) returns contents of k3_Filename_Textbox as a double


% --- Executes during object creation, after setting all properties.
function k3_Filename_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k3_Filename_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in k3_FileExt_Popupmenu.
function k3_FileExt_Popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to k3_FileExt_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns k3_FileExt_Popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from k3_FileExt_Popupmenu


% --- Executes during object creation, after setting all properties.
function k3_FileExt_Popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k3_FileExt_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in k4_Checkbox.
function k4_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to k4_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of k4_Checkbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
    set(handles.runPixelCalc_Pushbutton,'Enable',enable)
else
    enable = 'off';
    if (~get(handles.vB_Checkbox,'Value') && ~get(handles.K1_Checkbox,'Value') && ~get(handles.k2_Checkbox,'Value')...
            && ~get(handles.k3_Checkbox,'Value') && ~get(handles.k4_Checkbox,'Value') && ~get(handles.Ki_Checkbox,'Value')...
            && ~get(handles.Vs_Checkbox,'Value') && ~get(handles.Vt_Checkbox,'Value') && ~get(handles.K1divk2_Checkbox,'Value')...
            && ~get(handles.k3divk4_Checkbox,'Value') && ~get(handles.alpha1_Checkbox,'Value') && ~get(handles.alpha2_Checkbox,'Value'))    
        set(handles.runPixelCalc_Pushbutton,'Enable','off')
    end
end
set(handles.restrict_k4_Checkbox,'Enable',enable);
set(handles.save_k4_Checkbox,'Enable',enable);
set(handles.k4_Filename_Textbox,'Enable','off');
set(handles.k4_FileExt_Popupmenu,'Enable','off');
set(handles.k4_LowerTextbox,'Enable',enable);
set(handles.k4_UpperTextbox,'Enable',enable);

% --- Executes on button press in restrict_k4_Checkbox.
function restrict_k4_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to restrict_k4_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of restrict_k4_Checkbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
else
    enable = 'off';
end
set(handles.k4_LowerTextbox,'Enable',enable);
set(handles.k4_UpperTextbox,'Enable',enable);


function k4_LowerTextbox_Callback(hObject, eventdata, handles)
% hObject    handle to k4_LowerTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k4_LowerTextbox as text
%        str2double(get(hObject,'String')) returns contents of k4_LowerTextbox as a double
value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('Lower value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function k4_LowerTextbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k4_LowerTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function k4_UpperTextbox_Callback(hObject, eventdata, handles)
% hObject    handle to k4_UpperTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k4_UpperTextbox as text
%        str2double(get(hObject,'String')) returns contents of k4_UpperTextbox as a double
value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('Upper value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function k4_UpperTextbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k4_UpperTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in save_k4_Checkbox.
function save_k4_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to save_k4_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of save_k4_Checkbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
else
    enable = 'off';
end
set(handles.k4_Filename_Textbox,'Enable',enable);
set(handles.k4_FileExt_Popupmenu,'Enable',enable);


function k4_Filename_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to k4_Filename_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k4_Filename_Textbox as text
%        str2double(get(hObject,'String')) returns contents of k4_Filename_Textbox as a double


% --- Executes during object creation, after setting all properties.
function k4_Filename_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k4_Filename_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in k4_FileExt_Popupmenu.
function k4_FileExt_Popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to k4_FileExt_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns k4_FileExt_Popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from k4_FileExt_Popupmenu


% --- Executes during object creation, after setting all properties.
function k4_FileExt_Popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k4_FileExt_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Ki_Checkbox.
function Ki_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to Ki_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Ki_Checkbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
    set(handles.runPixelCalc_Pushbutton,'Enable',enable)
else
    enable = 'off';
    if (~get(handles.vB_Checkbox,'Value') && ~get(handles.K1_Checkbox,'Value') && ~get(handles.k2_Checkbox,'Value')...
            && ~get(handles.k3_Checkbox,'Value') && ~get(handles.k4_Checkbox,'Value') && ~get(handles.Ki_Checkbox,'Value')...
            && ~get(handles.Vs_Checkbox,'Value') && ~get(handles.Vt_Checkbox,'Value') && ~get(handles.K1divk2_Checkbox,'Value')...
            && ~get(handles.k3divk4_Checkbox,'Value') && ~get(handles.alpha1_Checkbox,'Value') && ~get(handles.alpha2_Checkbox,'Value'))    
        set(handles.runPixelCalc_Pushbutton,'Enable','off')
    end
end
set(handles.restrict_Ki_Checkbox,'Enable',enable);
set(handles.save_Ki_Checkbox,'Enable',enable);
set(handles.Ki_Filename_Textbox,'Enable','off');
set(handles.Ki_FileExt_Popupmenu,'Enable','off');
set(handles.Ki_LowerTextbox,'Enable',enable);
set(handles.Ki_UpperTextbox,'Enable',enable);

% --- Executes on button press in restrict_Ki_Checkbox.
function restrict_Ki_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to restrict_Ki_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of restrict_Ki_Checkbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
else
    enable = 'off';
end
set(handles.Ki_LowerTextbox,'Enable',enable);
set(handles.Ki_UpperTextbox,'Enable',enable);


function Ki_LowerTextbox_Callback(hObject, eventdata, handles)
% hObject    handle to Ki_LowerTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ki_LowerTextbox as text
%        str2double(get(hObject,'String')) returns contents of Ki_LowerTextbox as a double
value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('Lower value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function Ki_LowerTextbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ki_LowerTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Ki_UpperTextbox_Callback(hObject, eventdata, handles)
% hObject    handle to Ki_UpperTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ki_UpperTextbox as text
%        str2double(get(hObject,'String')) returns contents of Ki_UpperTextbox as a double
value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('Upper value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function Ki_UpperTextbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ki_UpperTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in save_Ki_Checkbox.
function save_Ki_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to save_Ki_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of save_Ki_Checkbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
else
    enable = 'off';
end
set(handles.Ki_Filename_Textbox,'Enable',enable);
set(handles.Ki_FileExt_Popupmenu,'Enable',enable);


function Ki_Filename_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to Ki_Filename_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ki_Filename_Textbox as text
%        str2double(get(hObject,'String')) returns contents of Ki_Filename_Textbox as a double


% --- Executes during object creation, after setting all properties.
function Ki_Filename_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ki_Filename_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Ki_FileExt_Popupmenu.
function Ki_FileExt_Popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to Ki_FileExt_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Ki_FileExt_Popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Ki_FileExt_Popupmenu


% --- Executes during object creation, after setting all properties.
function Ki_FileExt_Popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ki_FileExt_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Vs_Checkbox.
function Vs_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to Vs_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Vs_Checkbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
    set(handles.runPixelCalc_Pushbutton,'Enable',enable)
else
    enable = 'off';
    if (~get(handles.vB_Checkbox,'Value') && ~get(handles.K1_Checkbox,'Value') && ~get(handles.k2_Checkbox,'Value')...
            && ~get(handles.k3_Checkbox,'Value') && ~get(handles.k4_Checkbox,'Value') && ~get(handles.Ki_Checkbox,'Value')...
            && ~get(handles.Vs_Checkbox,'Value') && ~get(handles.Vt_Checkbox,'Value') && ~get(handles.K1divk2_Checkbox,'Value')...
            && ~get(handles.k3divk4_Checkbox,'Value') && ~get(handles.alpha1_Checkbox,'Value') && ~get(handles.alpha2_Checkbox,'Value'))    
        set(handles.runPixelCalc_Pushbutton,'Enable','off')
    end
end
set(handles.restrict_Vs_Checkbox,'Enable',enable);
set(handles.save_Vs_Checkbox,'Enable',enable);
set(handles.Vs_Filename_Textbox,'Enable','off');
set(handles.Vs_FileExt_Popupmenu,'Enable','off');
set(handles.Vs_LowerTextbox,'Enable',enable);
set(handles.Vs_UpperTextbox,'Enable',enable);

% --- Executes on button press in restrict_Vs_Checkbox.
function restrict_Vs_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to restrict_Vs_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of restrict_Vs_Checkbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
else
    enable = 'off';
end
set(handles.Vs_LowerTextbox,'Enable',enable);
set(handles.Vs_UpperTextbox,'Enable',enable);


function Vs_LowerTextbox_Callback(hObject, eventdata, handles)
% hObject    handle to Vs_LowerTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Vs_LowerTextbox as text
%        str2double(get(hObject,'String')) returns contents of Vs_LowerTextbox as a double
value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('Lower value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function Vs_LowerTextbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Vs_LowerTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Vs_UpperTextbox_Callback(hObject, eventdata, handles)
% hObject    handle to Vs_UpperTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Vs_UpperTextbox as text
%        str2double(get(hObject,'String')) returns contents of Vs_UpperTextbox as a double
value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('Upper value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function Vs_UpperTextbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Vs_UpperTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in save_Vs_Checkbox.
function save_Vs_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to save_Vs_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of save_Vs_Checkbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
else
    enable = 'off';
end
set(handles.Vs_Filename_Textbox,'Enable',enable);
set(handles.Vs_FileExt_Popupmenu,'Enable',enable);


function Vs_Filename_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to Vs_Filename_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Vs_Filename_Textbox as text
%        str2double(get(hObject,'String')) returns contents of Vs_Filename_Textbox as a double


% --- Executes during object creation, after setting all properties.
function Vs_Filename_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Vs_Filename_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Vs_FileExt_Popupmenu.
function Vs_FileExt_Popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to Vs_FileExt_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Vs_FileExt_Popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Vs_FileExt_Popupmenu


% --- Executes during object creation, after setting all properties.
function Vs_FileExt_Popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Vs_FileExt_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Vt_Checkbox.
function Vt_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to Vt_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Vt_Checkbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
    set(handles.runPixelCalc_Pushbutton,'Enable',enable)
else
    enable = 'off';
    if (~get(handles.vB_Checkbox,'Value') && ~get(handles.K1_Checkbox,'Value') && ~get(handles.k2_Checkbox,'Value')...
            && ~get(handles.k3_Checkbox,'Value') && ~get(handles.k4_Checkbox,'Value') && ~get(handles.Ki_Checkbox,'Value')...
            && ~get(handles.Vs_Checkbox,'Value') && ~get(handles.Vt_Checkbox,'Value') && ~get(handles.K1divk2_Checkbox,'Value')...
            && ~get(handles.k3divk4_Checkbox,'Value') && ~get(handles.alpha1_Checkbox,'Value') && ~get(handles.alpha2_Checkbox,'Value'))    
        set(handles.runPixelCalc_Pushbutton,'Enable','off')
    end
end
set(handles.restrict_Vt_Checkbox,'Enable',enable);
set(handles.save_Vt_Checkbox,'Enable',enable);
set(handles.Vt_Filename_Textbox,'Enable','off');
set(handles.Vt_FileExt_Popupmenu,'Enable','off');
set(handles.Vt_LowerTextbox,'Enable',enable);
set(handles.Vt_UpperTextbox,'Enable',enable);

% --- Executes on button press in restrict_Vt_Checkbox.
function restrict_Vt_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to restrict_Vt_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of restrict_Vt_Checkbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
else
    enable = 'off';
end
set(handles.Vt_LowerTextbox,'Enable',enable);
set(handles.Vt_UpperTextbox,'Enable',enable);


function Vt_LowerTextbox_Callback(hObject, eventdata, handles)
% hObject    handle to Vt_LowerTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Vt_LowerTextbox as text
%        str2double(get(hObject,'String')) returns contents of Vt_LowerTextbox as a double
value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('Lower value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function Vt_LowerTextbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Vt_LowerTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Vt_UpperTextbox_Callback(hObject, eventdata, handles)
% hObject    handle to Vt_UpperTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Vt_UpperTextbox as text
%        str2double(get(hObject,'String')) returns contents of Vt_UpperTextbox as a double
value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('Upper value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function Vt_UpperTextbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Vt_UpperTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in save_Vt_Checkbox.
function save_Vt_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to save_Vt_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of save_Vt_Checkbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
else
    enable = 'off';
end
set(handles.Vt_Filename_Textbox,'Enable',enable);
set(handles.Vt_FileExt_Popupmenu,'Enable',enable);


function Vt_Filename_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to Vt_Filename_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Vt_Filename_Textbox as text
%        str2double(get(hObject,'String')) returns contents of Vt_Filename_Textbox as a double


% --- Executes during object creation, after setting all properties.
function Vt_Filename_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Vt_Filename_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Vt_FileExt_Popupmenu.
function Vt_FileExt_Popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to Vt_FileExt_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Vt_FileExt_Popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Vt_FileExt_Popupmenu


% --- Executes during object creation, after setting all properties.
function Vt_FileExt_Popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Vt_FileExt_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in K1divk2_Checkbox.
function K1divk2_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to K1divk2_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of K1divk2_Checkbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
    set(handles.runPixelCalc_Pushbutton,'Enable',enable)
else
    enable = 'off';
    if (~get(handles.vB_Checkbox,'Value') && ~get(handles.K1_Checkbox,'Value') && ~get(handles.k2_Checkbox,'Value')...
            && ~get(handles.k3_Checkbox,'Value') && ~get(handles.k4_Checkbox,'Value') && ~get(handles.Ki_Checkbox,'Value')...
            && ~get(handles.Vs_Checkbox,'Value') && ~get(handles.Vt_Checkbox,'Value') && ~get(handles.K1divk2_Checkbox,'Value')...
            && ~get(handles.k3divk4_Checkbox,'Value') && ~get(handles.alpha1_Checkbox,'Value') && ~get(handles.alpha2_Checkbox,'Value'))    
        set(handles.runPixelCalc_Pushbutton,'Enable','off')
    end
end
set(handles.restrict_K1divk2_Checkbox,'Enable',enable);
set(handles.save_K1divk2_Checkbox,'Enable',enable);
set(handles.K1divk2_Filename_Textbox,'Enable','off');
set(handles.K1divk2_FileExt_Popupmenu,'Enable','off');
set(handles.K1divk2_LowerTextbox,'Enable',enable);
set(handles.K1divk2_UpperTextbox,'Enable',enable);

% --- Executes on button press in restrict_K1divk2_Checkbox.
function restrict_K1divk2_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to restrict_K1divk2_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of restrict_K1divk2_Checkbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
else
    enable = 'off';
end
set(handles.K1divk2_LowerTextbox,'Enable',enable);
set(handles.K1divk2_UpperTextbox,'Enable',enable);


function K1divk2_LowerTextbox_Callback(hObject, eventdata, handles)
% hObject    handle to K1divk2_LowerTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of K1divk2_LowerTextbox as text
%        str2double(get(hObject,'String')) returns contents of K1divk2_LowerTextbox as a double
value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('Lower value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function K1divk2_LowerTextbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to K1divk2_LowerTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function K1divk2_UpperTextbox_Callback(hObject, eventdata, handles)
% hObject    handle to K1divk2_UpperTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of K1divk2_UpperTextbox as text
%        str2double(get(hObject,'String')) returns contents of K1divk2_UpperTextbox as a double
value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('Upper value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function K1divk2_UpperTextbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to K1divk2_UpperTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in save_K1divk2_Checkbox.
function save_K1divk2_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to save_K1divk2_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of save_K1divk2_Checkbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
else
    enable = 'off';
end
set(handles.K1divk2_Filename_Textbox,'Enable',enable);
set(handles.K1divk2_FileExt_Popupmenu,'Enable',enable);


function K1divk2_Filename_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to K1divk2_Filename_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of K1divk2_Filename_Textbox as text
%        str2double(get(hObject,'String')) returns contents of K1divk2_Filename_Textbox as a double


% --- Executes during object creation, after setting all properties.
function K1divk2_Filename_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to K1divk2_Filename_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in K1divk2_FileExt_Popupmenu.
function K1divk2_FileExt_Popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to K1divk2_FileExt_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns K1divk2_FileExt_Popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from K1divk2_FileExt_Popupmenu


% --- Executes during object creation, after setting all properties.
function K1divk2_FileExt_Popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to K1divk2_FileExt_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in k3divk4_Checkbox.
function k3divk4_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to k3divk4_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of k3divk4_Checkbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
    set(handles.runPixelCalc_Pushbutton,'Enable',enable)
else
    enable = 'off';
    if (~get(handles.vB_Checkbox,'Value') && ~get(handles.K1_Checkbox,'Value') && ~get(handles.k2_Checkbox,'Value')...
            && ~get(handles.k3_Checkbox,'Value') && ~get(handles.k4_Checkbox,'Value') && ~get(handles.Ki_Checkbox,'Value')...
            && ~get(handles.Vs_Checkbox,'Value') && ~get(handles.Vt_Checkbox,'Value') && ~get(handles.K1divk2_Checkbox,'Value')...
            && ~get(handles.k3divk4_Checkbox,'Value') && ~get(handles.alpha1_Checkbox,'Value') && ~get(handles.alpha2_Checkbox,'Value'))    
        set(handles.runPixelCalc_Pushbutton,'Enable','off')
    end
end
set(handles.restrict_k3divk4_Checkbox,'Enable',enable);
set(handles.save_k3divk4_Checkbox,'Enable',enable);
set(handles.k3divk4_Filename_Textbox,'Enable','off');
set(handles.k3divk4_FileExt_Popupmenu,'Enable','off');
set(handles.k3divk4_LowerTextbox,'Enable',enable);
set(handles.k3divk4_UpperTextbox,'Enable',enable);

% --- Executes on button press in restrict_k3divk4_Checkbox.
function restrict_k3divk4_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to restrict_k3divk4_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of restrict_k3divk4_Checkbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
else
    enable = 'off';
end
set(handles.k3divk4_LowerTextbox,'Enable',enable);
set(handles.k3divk4_UpperTextbox,'Enable',enable);


function k3divk4_LowerTextbox_Callback(hObject, eventdata, handles)
% hObject    handle to k3divk4_LowerTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k3divk4_LowerTextbox as text
%        str2double(get(hObject,'String')) returns contents of k3divk4_LowerTextbox as a double
value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('Lower value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function k3divk4_LowerTextbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k3divk4_LowerTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function k3divk4_UpperTextbox_Callback(hObject, eventdata, handles)
% hObject    handle to k3divk4_UpperTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k3divk4_UpperTextbox as text
%        str2double(get(hObject,'String')) returns contents of k3divk4_UpperTextbox as a double
value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('Upper value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function k3divk4_UpperTextbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k3divk4_UpperTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in save_k3divk4_Checkbox.
function save_k3divk4_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to save_k3divk4_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of save_k3divk4_Checkbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
else
    enable = 'off';
end
set(handles.k3divk4_Filename_Textbox,'Enable',enable);
set(handles.k3divk4_FileExt_Popupmenu,'Enable',enable);


function k3divk4_Filename_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to k3divk4_Filename_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k3divk4_Filename_Textbox as text
%        str2double(get(hObject,'String')) returns contents of k3divk4_Filename_Textbox as a double


% --- Executes during object creation, after setting all properties.
function k3divk4_Filename_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k3divk4_Filename_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in k3divk4_FileExt_Popupmenu.
function k3divk4_FileExt_Popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to k3divk4_FileExt_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns k3divk4_FileExt_Popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from k3divk4_FileExt_Popupmenu


% --- Executes during object creation, after setting all properties.
function k3divk4_FileExt_Popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k3divk4_FileExt_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in alpha1_Checkbox.
function alpha1_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to alpha1_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of alpha1_Checkbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
    set(handles.runPixelCalc_Pushbutton,'Enable',enable)
else
    enable = 'off';
    if (~get(handles.vB_Checkbox,'Value') && ~get(handles.K1_Checkbox,'Value') && ~get(handles.k2_Checkbox,'Value')...
            && ~get(handles.k3_Checkbox,'Value') && ~get(handles.k4_Checkbox,'Value') && ~get(handles.Ki_Checkbox,'Value')...
            && ~get(handles.Vs_Checkbox,'Value') && ~get(handles.Vt_Checkbox,'Value') && ~get(handles.K1divk2_Checkbox,'Value')...
            && ~get(handles.k3divk4_Checkbox,'Value') && ~get(handles.alpha1_Checkbox,'Value') && ~get(handles.alpha2_Checkbox,'Value'))    
        set(handles.runPixelCalc_Pushbutton,'Enable','off')
    end
end
set(handles.restrict_alpha1_Checkbox,'Enable',enable);
set(handles.save_alpha1_Checkbox,'Enable',enable);
set(handles.alpha1_Filename_Textbox,'Enable','off');
set(handles.alpha1_FileExt_Popupmenu,'Enable','off');
set(handles.alpha1_LowerTextbox,'Enable',enable);
set(handles.alpha1_UpperTextbox,'Enable',enable);

% --- Executes on button press in restrict_alpha1_Checkbox.
function restrict_alpha1_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to restrict_alpha1_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of restrict_alpha1_Checkbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
else
    enable = 'off';
end
set(handles.alpha1_LowerTextbox,'Enable',enable);
set(handles.alpha1_UpperTextbox,'Enable',enable);


function alpha1_LowerTextbox_Callback(hObject, eventdata, handles)
% hObject    handle to alpha1_LowerTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alpha1_LowerTextbox as text
%        str2double(get(hObject,'String')) returns contents of alpha1_LowerTextbox as a double
value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('Lower value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function alpha1_LowerTextbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alpha1_LowerTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function alpha1_UpperTextbox_Callback(hObject, eventdata, handles)
% hObject    handle to alpha1_UpperTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alpha1_UpperTextbox as text
%        str2double(get(hObject,'String')) returns contents of alpha1_UpperTextbox as a double
value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('Upper value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function alpha1_UpperTextbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alpha1_UpperTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in save_alpha1_Checkbox.
function save_alpha1_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to save_alpha1_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of save_alpha1_Checkbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
else
    enable = 'off';
end
set(handles.alpha1_Filename_Textbox,'Enable',enable);
set(handles.alpha1_FileExt_Popupmenu,'Enable',enable);


function alpha1_Filename_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to alpha1_Filename_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alpha1_Filename_Textbox as text
%        str2double(get(hObject,'String')) returns contents of alpha1_Filename_Textbox as a double


% --- Executes during object creation, after setting all properties.
function alpha1_Filename_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alpha1_Filename_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in alpha1_FileExt_Popupmenu.
function alpha1_FileExt_Popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to alpha1_FileExt_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns alpha1_FileExt_Popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from alpha1_FileExt_Popupmenu


% --- Executes during object creation, after setting all properties.
function alpha1_FileExt_Popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alpha1_FileExt_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in alpha2_Checkbox.
function alpha2_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to alpha2_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of alpha2_Checkbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
    set(handles.runPixelCalc_Pushbutton,'Enable',enable)
else
    enable = 'off';
    if (~get(handles.vB_Checkbox,'Value') && ~get(handles.K1_Checkbox,'Value') && ~get(handles.k2_Checkbox,'Value')...
            && ~get(handles.k3_Checkbox,'Value') && ~get(handles.k4_Checkbox,'Value') && ~get(handles.Ki_Checkbox,'Value')...
            && ~get(handles.Vs_Checkbox,'Value') && ~get(handles.Vt_Checkbox,'Value') && ~get(handles.K1divk2_Checkbox,'Value')...
            && ~get(handles.k3divk4_Checkbox,'Value') && ~get(handles.alpha1_Checkbox,'Value') && ~get(handles.alpha2_Checkbox,'Value'))    
        set(handles.runPixelCalc_Pushbutton,'Enable','off')
    end
end
set(handles.restrict_alpha2_Checkbox,'Enable',enable);
set(handles.save_alpha2_Checkbox,'Enable',enable);
set(handles.alpha2_Filename_Textbox,'Enable','off');
set(handles.alpha2_FileExt_Popupmenu,'Enable','off');
set(handles.alpha2_LowerTextbox,'Enable',enable);
set(handles.alpha2_UpperTextbox,'Enable',enable);

% --- Executes on button press in restrict_alpha2_Checkbox.
function restrict_alpha2_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to restrict_alpha2_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of restrict_alpha2_Checkbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
else
    enable = 'off';
end
set(handles.alpha2_LowerTextbox,'Enable',enable);
set(handles.alpha2_UpperTextbox,'Enable',enable);


function alpha2_LowerTextbox_Callback(hObject, eventdata, handles)
% hObject    handle to alpha2_LowerTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alpha2_LowerTextbox as text
%        str2double(get(hObject,'String')) returns contents of alpha2_LowerTextbox as a double
value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('Lower value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function alpha2_LowerTextbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alpha2_LowerTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function alpha2_UpperTextbox_Callback(hObject, eventdata, handles)
% hObject    handle to alpha2_UpperTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alpha2_UpperTextbox as text
%        str2double(get(hObject,'String')) returns contents of alpha2_UpperTextbox as a double
value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('Upper value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function alpha2_UpperTextbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alpha2_UpperTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in save_alpha2_Checkbox.
function save_alpha2_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to save_alpha2_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of save_alpha2_Checkbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
else
    enable = 'off';
end
set(handles.alpha2_Filename_Textbox,'Enable',enable);
set(handles.alpha2_FileExt_Popupmenu,'Enable',enable);


function alpha2_Filename_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to alpha2_Filename_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alpha2_Filename_Textbox as text
%        str2double(get(hObject,'String')) returns contents of alpha2_Filename_Textbox as a double


% --- Executes during object creation, after setting all properties.
function alpha2_Filename_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alpha2_Filename_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in alpha2_FileExt_Popupmenu.
function alpha2_FileExt_Popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to alpha2_FileExt_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns alpha2_FileExt_Popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from alpha2_FileExt_Popupmenu


% --- Executes during object creation, after setting all properties.
function alpha2_FileExt_Popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alpha2_FileExt_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ttp_Checkbox.
function ttp_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to ttp_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ttp_Checkbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
    set(handles.runPixelCalc_Pushbutton,'Enable',enable)
else
    enable = 'off';
    if (~get(handles.vB_Checkbox,'Value') && ~get(handles.K1_Checkbox,'Value') && ~get(handles.k2_Checkbox,'Value')...
            && ~get(handles.k3_Checkbox,'Value') && ~get(handles.k4_Checkbox,'Value') && ~get(handles.Ki_Checkbox,'Value')...
            && ~get(handles.Vs_Checkbox,'Value') && ~get(handles.Vt_Checkbox,'Value') && ~get(handles.K1divk2_Checkbox,'Value')...
            && ~get(handles.k3divk4_Checkbox,'Value') && ~get(handles.alpha1_Checkbox,'Value') && ~get(handles.alpha2_Checkbox,'Value'))    
        set(handles.runPixelCalc_Pushbutton,'Enable','off')
    end
end
set(handles.restrict_ttp_Checkbox,'Enable',enable);
set(handles.save_ttp_Checkbox,'Enable',enable);
set(handles.ttp_Filename_Textbox,'Enable','off');
set(handles.ttp_FileExt_Popupmenu,'Enable','off');
set(handles.ttp_LowerTextbox,'Enable',enable);
set(handles.ttp_UpperTextbox,'Enable',enable);

% --- Executes on button press in restrict_ttp_Checkbox.
function restrict_ttp_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to restrict_ttp_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of restrict_ttp_Checkbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
else
    enable = 'off';
end
set(handles.ttp_LowerTextbox,'Enable',enable);
set(handles.ttp_UpperTextbox,'Enable',enable);


function ttp_LowerTextbox_Callback(hObject, eventdata, handles)
% hObject    handle to ttp_LowerTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ttp_LowerTextbox as text
%        str2double(get(hObject,'String')) returns contents of ttp_LowerTextbox as a double
value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('Lower value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function ttp_LowerTextbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ttp_LowerTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ttp_UpperTextbox_Callback(hObject, eventdata, handles)
% hObject    handle to ttp_UpperTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ttp_UpperTextbox as text
%        str2double(get(hObject,'String')) returns contents of ttp_UpperTextbox as a double
value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('Upper value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function ttp_UpperTextbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ttp_UpperTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in save_ttp_Checkbox.
function save_ttp_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to save_ttp_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of save_ttp_Checkbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
else
    enable = 'off';
end
set(handles.ttp_Filename_Textbox,'Enable',enable);
set(handles.ttp_FileExt_Popupmenu,'Enable',enable);


function ttp_Filename_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to ttp_Filename_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ttp_Filename_Textbox as text
%        str2double(get(hObject,'String')) returns contents of ttp_Filename_Textbox as a double


% --- Executes during object creation, after setting all properties.
function ttp_Filename_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ttp_Filename_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ttp_FileExt_Popupmenu.
function ttp_FileExt_Popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to ttp_FileExt_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ttp_FileExt_Popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ttp_FileExt_Popupmenu


% --- Executes during object creation, after setting all properties.
function ttp_FileExt_Popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ttp_FileExt_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
