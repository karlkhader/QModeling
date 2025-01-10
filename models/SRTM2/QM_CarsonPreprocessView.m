function varargout = QM_CarsonPreprocessView(varargin)
% QM_CARSONPREPROCESSVIEW MATLAB code for QM_CarsonPreprocessView.fig
%  
%   QM_CarsonPreprocessView generates a GUI to set the parameters
%   used in the preprocessing of SRTM2 model. This code must be 
%   executed from QModeling GUI.  
%
% Last Modified by GUIDE v2.5 29-Nov-2013 13:45:16
%---------------------------------------------------------------------------------
% QM_CarsonPreprocessView is part of QModeling.
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
                   'gui_OpeningFcn', @QM_CarsonPreprocessView_OpeningFcn, ...
                   'gui_OutputFcn',  @QM_CarsonPreprocessView_OutputFcn, ...
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


% --- Executes just before QM_CarsonPreprocessView is made visible.
function QM_CarsonPreprocessView_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to QM_CarsonPreprocessView (see VARARGIN)
global QMmaindata

% Choose default command line output for QM_CarsonPreprocessView
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

%Save mainGUI handle
mainGuiInput = find(strcmp(varargin, 'mainhandles'));
handles.mainGUI = varargin{mainGuiInput+1};

%Set the popupmenus
num_rois = size(QMmaindata.Crois,2);
strmenu = cell(1,num_rois);
for i = 1:num_rois
    %strmenu{1,i}  = strcat('TAC ',num2str(i));
    strmenu{1,i} = QMmaindata.Crois{1,i};
end

set(handles.CroiTac_Popupmenu,'String',strmenu);
set(handles.CroiTac_Popupmenu,'Enable','on');
set(handles.CroiTac_Popupmenu,'Value',QMmaindata.idCt);

set(handles.CrefTac_Popupmenu,'String',strmenu);
set(handles.CrefTac_Popupmenu,'Value',QMmaindata.idCr);
set(handles.CrefTac_Popupmenu,'Enable','on');

%Enable preprocess button
set(handles.preprocessCarson_Pushbutton,'Enable','on');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes QM_CarsonPreprocessView wait for user response (see UIRESUME)
uiwait(handles.CarsonPreprocess_figure);


% --- Outputs from this function are returned to the command line.
function varargout = QM_CarsonPreprocessView_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
delete(hObject);


% --- Executes on selection change in CroiTac_Popupmenu.
function CroiTac_Popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to CroiTac_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns CroiTac_Popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from CroiTac_Popupmenu

idCt = get(hObject,'Value');
idCr = get(handles.CrefTac_Popupmenu,'Value');

if idCt == idCr
    warndlg(char({'Selecting the same TAC for both ROI and REF tac may be incorrect';'Please, change one TAC.'}),'Warning','modal');
    uicontrol(hObject);
end


% --- Executes during object creation, after setting all properties.
function CroiTac_Popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CroiTac_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in CrefTac_Popupmenu.
function CrefTac_Popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to CrefTac_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns CrefTac_Popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from CrefTac_Popupmenu

idCr = get(hObject,'Value');
idCt = get(handles.CroiTac_Popupmenu,'Value');

if idCt == idCr
    warndlg(char({'Selecting the same TAC for both ROI and REF tac may be incorrect';'Please, change one TAC.'}),'Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function CrefTac_Popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CrefTac_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in preprocessCarson_Pushbutton.
function preprocessCarson_Pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to preprocessCarson_Pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global QMmaindata
global QMlastPreprocess
global QMf_logfile

idCt = get(handles.CroiTac_Popupmenu,'Value');
idCr = get(handles.CrefTac_Popupmenu,'Value');

if idCt == idCr
    warndlg(char({'For a good preprocessing, the region of interest (TAC1) must be different than reference region (TAC2)';...
                    'Please, change one.'}),'Preprocessing error','modal');
    uicontrol(handles.CrefTac_Popupmenu);
    return
end

k2_pChecked = get(handles.k2_pfit_Checkbox,'Value');
k2_pFit = str2double(get(handles.k2_pfit_Textbox,'String'));
if isnan(k2_pFit)
    warndlg(strcat('k2',char(39),' must be a positive number'),'Preprocessing error','modal');
    uicontrol(handles.k2_pfit_Textbox);
    return
end

resampling = str2double(get(handles.resampling_Textbox,'String'));
if isnan(resampling)
    warndlg('Resampling must be a positive number','Preprocessing error','modal');
    uicontrol(handles.resampling_Textbox);
    return
end

threshold = str2double(get(handles.threshold_Textbox,'String'));
if isnan(threshold)|| (threshold < 0) || (threshold > 100)
    warndlg('Threshold must be a positive number and less than 100','Preprocessing error','modal');
    uicontrol(handles.threshold_Textbox);
    return
end

basisF = str2double(get(handles.basis_Textbox,'String'));
if isnan(basisF)
    warndlg('Basis must be a positive number','Preprocessing error','modal');
    uicontrol(handles.basis_Textbox);
    return
end

k2a_max = str2double(get(handles.k2amax_Textbox,'String'));
if isnan(k2a_max)
    warndlg('k2a max must be a positive number','Preprocessing error','modal');
    uicontrol(handles.k2amax_Textbox);
    return
end

k2a_min = str2double(get(handles.k2amin_Textbox,'String'));
if isnan(k2a_min)
    warndlg('k2a min must be a positive number','Preprocessing error','modal');
    uicontrol(handles.k2amin_Textbox);
    return
end

waitbarhandle = waitbar(0,'Preparing TACs, please wait...');
ver=version('-release');
if (str2num(ver(1:4))<2014) || isequal(ver,'2014a')
    set(findobj(waitbarhandle,'type','patch'),'edgecolor','b','facecolor','b');
else
    wbc = allchild(waitbarhandle);     %you need to get at a hidden child
    wbc(1).JavaPeer.setForeground( wbc(1).JavaPeer.getBackground.BLUE)
    wbc(1).JavaPeer.setStringPainted(true)
end

error = 0;

%Update the logfile
aux_date=datestr(now);
fprintf(QMf_logfile,strcat('\n-SRTM2-\n'));
if k2_pChecked
    fprintf(QMf_logfile,strcat(aux_date(end-8:end),' Preprocessing model (k2amin= ',num2str(k2a_min),...
    ', k2amax= ',num2str(k2a_max),', BF= ',num2str(basisF),', resampling= ',num2str(resampling),')\n'));
else
    fprintf(QMf_logfile,strcat(aux_date(end-8:end),' Preprocessing model (k2amin= ',num2str(k2a_min),...
    ', k2amax= ',num2str(k2a_max),', BF= ',num2str(basisF),', resampling= ',num2str(resampling),', k2',char(39),' =',num2str(k2_pFit),')\n'));
end
    

try
    [results, plots, B, th3] = QM_preprocessCarson(QMmaindata.Crois{2,idCt},QMmaindata.Crois{2,idCr},QMmaindata.times,...
                                            [k2_pChecked k2_pFit],resampling,basisF,k2a_max,k2a_min,waitbarhandle);
catch ME
    errordlg(char({'An unexpected error ocurred in the preprocessing';'Please, check logfile.txt'}),'Preprocessing error','modal');
    error = 1;
    fprintf(QMf_logfile,strcat('  ERROR: Preprocessing did not finish succesfully','\n','   Error message: ',ME.message,'\n',...
        '   File: ',ME.stack(1).file,'\n','   Line: ',num2str(ME.stack(1).line),'\n'));
end

delete(waitbarhandle);  %Need to delete before ploting
    

if error == 0
    %Update the logfile
    aux_date=datestr(now);
    fprintf(QMf_logfile,strcat(aux_date(end-8:end),' Preprocessing finished \n'));
    
    %Obtain handles using GUIDATA with the caller's handle 
    mainhandles = guidata(handles.mainGUI);
    
    %Plotting the preprocess   
    tavg = ((QMmaindata.times(:,1)+QMmaindata.times(:,2))/2);
    tmin = tavg/60;
    
    plot(mainhandles.axes,tmin,QMmaindata.Crois{2,idCt},'b.-',tmin,QMmaindata.Crois{2,idCr},'r*-',tmin,plots,'g+-');
    xlabel(mainhandles.axes,'Time (minutes)');
    ylabel(mainhandles.axes,'Original units');
    title(mainhandles.axes,'SRTM2 fit')
    legend(mainhandles.axes,'TAC1 from receptor-rich region', 'TAC2 from reference region', 'Fit to TAC1');
    
    %Save values (it is used, for example, to show the values of the TACs when the "view Values" option is selected
%     set(mainhandles.axes,'UserData',{'SRTM*',[tmin QMmaindata.Crois{2,idCt} plots' QMmaindata.Crois{2,idCr}]});
    t_aux=table([1:1:length(plots)]',tmin,QMmaindata.Crois{2,idCt},plots',QMmaindata.Crois{2,idCr});
    set(mainhandles.axes,'UserData',{'SRTM*',table2cell(t_aux)});
    
    %Showing results in text area
    strresults = cell(1,6);
    strresults{1,1} = 'Carson preprocessing:';
    strresults{1,2} = '-----------------------------------';
    if ~k2_pChecked
        strresults{1,3}  = [strcat('k2',char(39),':'),'     ',num2str(results{1},'%.5f')];
    else
        strresults{1,3}  = [strcat('k2',char(39),':'),'     ',num2str(results{1},'%.5f'), '     (SE: ',num2str(sqrt(results{8}(1)),'%1.2e'),')'];
    end
    strresults{1,4}  = ['R1:      ',num2str(results{2},'%.5f'), '     (SE: ',num2str(sqrt(results{8}(2)),'%1.2e'),')'];
    strresults{1,5}  = ['k2a:    ',num2str(results{3},'%.5f'), '     (SE: ',num2str(sqrt(results{8}(3)),'%1.2e'),')'];
    strresults{1,6}  = ['BPnd: ',num2str(results{4},'%.5f')];
    strresults2 = cell(1,4);
    strresults2{1,1} = 'Goodness of fit:';
    strresults2{1,2} = '-----------------------------------';
    strresults2{1,3} = ['NMSE:           ',num2str(results{5},'%.5f')];
    strresults2{1,4} = ['Corr. Coef.:  ',num2str(results{6},'%.5f')];
    
    set(mainhandles.results_Text,'String',strresults);
    set(mainhandles.results_Text2,'String',strresults2);
    
    %Save last preprocess
    QMlastPreprocess.SRTM2.idCt = idCt;
    QMlastPreprocess.SRTM2.idCr = idCr;    
    QMlastPreprocess.SRTM2.B=B;
    QMlastPreprocess.SRTM2.th3=th3;
    QMlastPreprocess.SRTM2.threshold=threshold;
    QMlastPreprocess.SRTM2.plots = plots;
    QMlastPreprocess.SRTM2.results=results;
    QMlastPreprocess.SRTM2.k2_pChecked=k2_pChecked;
       
    %Restart same GUI components
    set(mainhandles.plottedCt_Popupmenu,'Value',idCt);
    set(mainhandles.plottedCr_Popupmenu,'Value',idCr);
    set(mainhandles.plottedCt_Popupmenu,'Enable','off');
    set(mainhandles.plottedCr_Popupmenu,'Enable','off');
    set(mainhandles.showTacs_Radiobutton,'Value',0)
    set(mainhandles.showPreprocess_Radiobutton,'Value',1)
    set(mainhandles.showTacs_Radiobutton,'Visible','on')
    set(mainhandles.showPreprocess_Radiobutton,'Visible','on')
        
    set(mainhandles.selectImage_Popupmenu,'String','Parametric Images');
    set(mainhandles.selectImage_Popupmenu,'Enable','off');
    set(mainhandles.viewImage_Pushbutton,'Enable','off');    
    set(mainhandles.setMapParam_Pushbutton,'Enable','off');   

    %Detele the last parametric images created in temp file
    mainhandles.deleteTempImages(mainhandles);
    
    %Restart flag
    set(mainhandles.setMapParam_Pushbutton,'UserData',0);
    
    %Set panel focus
    mainhandles.setFocus(mainhandles,'PlotResultsPanels');
    
    if QMmaindata.flags.minWindow == 1
        minfig(handles.CarsonPreprocess_figure,1);
    end
        
end

% --- Executes on button press in k2_pfit_Checkbox.
function k2_pfit_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to k2_pfit_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of k2_pfit_Checkbox

checked = get(hObject,'Value');
if checked
    enable = 'off';
else
    enable = 'on';
end
set(handles.k2_pfit_Textbox,'Enable',enable);


    
function k2_pfit_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to k2_pfit_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k2_pfit_Textbox as text
%        str2double(get(hObject,'String')) returns contents of k2_pfit_Textbox as a double

value = str2double(get(hObject,'String'));
if isnan(value)   
    warndlg(strcat(strcat('k2',char(39)),' value must be a positive number'),'Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function k2_pfit_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k2_pfit_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function resampling_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to resampling_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of resampling_Textbox as text
%        str2double(get(hObject,'String')) returns contents of resampling_Textbox as a double
value = str2double(get(hObject,'String'));
if isnan(value)| (value < 0)
    warndlg('Resampling value must be a positive number','Warning','modal');
    uicontrol(hObject);
end


% --- Executes during object creation, after setting all properties.
function resampling_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to resampling_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in setDefaults_Pushbutton.
function setDefaults_Pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to setDefaults_Pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.CroiTac_Popupmenu,'Value',1);
set(handles.CrefTac_Popupmenu,'Value',2);
set(handles.k2amin_Textbox,'String','0.006');
set(handles.k2amax_Textbox,'String', '0.6');
set(handles.basis_Textbox,'String','100');
set(handles.resampling_Textbox,'String','5.0');
set(handles.threshold_Textbox,'String','3.0');
set(handles.k2_pfit_Textbox,'String','0.128205');
set(handles.k2_pfit_Textbox, 'Enable','off');
set(handles.k2_pfit_Checkbox,'Value',1)

% --- Executes on button press in help_Pushbutton.
function help_Pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to help_Pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global QMpath

% open(strcat(QMpath,filesep,'help',filesep,'PatlakRef_help.pdf'))
open(strcat(QMpath,filesep,'help',filesep,'html',filesep,'help_Carson.html'))

% --- Executes during object creation, after setting all properties.
function help_Pushbutton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to help_Pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

global QMpath
[a,map]=imread(strcat(QMpath,filesep,'icons',filesep,'csh_icon.png')); 

% Determina el tama�o de la imagen en las tres dimensiones
[r,c,d]=size(a); 
% redondea el tama�o de la matriz r / 38
x=ceil(r/38);
% redondea el tama�o de la matriz c / 51
y=ceil(c/51);
% crea una matriz g, con elementos de la matriz "a", pero considerando intervalos x e y
g=a(1:x:end,1:y:end,:); 
g(g==0)=240; %Para poner el borde de la imagen al mismo color del fondo que tengo puesto
set(hObject,'CData',g); 


% --- Executes on button press in OK_Pushbutton.
function OK_Pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to OK_Pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiresume(handles.CarsonPreprocess_figure);


% --- Executes when user attempts to close CarsonPreprocess_figure.
function CarsonPreprocess_figure_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to CarsonPreprocess_figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
uiresume(handles.CarsonPreprocess_figure);


function k2amax_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to k2amax_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k2amax_Textbox as text
value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('k2a max value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function k2amax_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k2amax_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function k2amin_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to k2amin_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k2amin_Textbox as text
value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('k2a min value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function k2amin_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k2amin_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function basis_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to basis_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of basis_Textbox as text
value = str2double(get(hObject,'String'));
if isnan(value)| (value < 1)
    warndlg('Basis value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function basis_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to basis_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function threshold_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to threshold_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of threshold_Textbox as text
value = str2double(get(hObject,'String'));
if isnan(value)| (value < 0)
    warndlg('Threshold value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function threshold_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to threshold_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
