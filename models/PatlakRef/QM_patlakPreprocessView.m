function varargout = QM_patlakPreprocessView(varargin)
% QM_PATLAKPREPROCESSVIEW MATLAB code for QM_patlakPreprocessView.fig
%  
%   QM_patlakPreprocessView generates a GUI to set the parameters
%   used in the preprocessing of Patlak model. This code must be 
%   executed from QModeling GUI.  
%
% Last Modified by GUIDE v2.5 23-Dec-2013 20:21:10
%---------------------------------------------------------------------------------
% QM_patlakPreprocessView is part of QModeling.
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
                   'gui_OpeningFcn', @QM_patlakPreprocessView_OpeningFcn, ...
                   'gui_OutputFcn',  @QM_patlakPreprocessView_OutputFcn, ...
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


% --- Executes just before QM_patlakPreprocessView is made visible.
function QM_patlakPreprocessView_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to QM_patlakPreprocessView (see VARARGIN)
global QMmaindata

% Choose default command line output for QM_patlakPreprocessView
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
    strmenu{1,i}  = QMmaindata.Crois{1,i};
end
set(handles.CroiTac_Popupmenu,'String',strmenu);
set(handles.CroiTac_Popupmenu,'Enable','on');
set(handles.CroiTac_Popupmenu,'Value',QMmaindata.idCt);

set(handles.CrefTac_Popupmenu,'String',strmenu);
set(handles.CrefTac_Popupmenu,'Value',QMmaindata.idCr);
set(handles.CrefTac_Popupmenu,'Enable','on');

% Enable preprocess button
set(handles.preprocessPatlak_Pushbutton,'Enable','on');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes QM_patlakPreprocessView wait for user response (see UIRESUME)
uiwait(handles.PatlakPreprocess_figure);


% --- Outputs from this function are returned to the command line.
function varargout = QM_patlakPreprocessView_OutputFcn(hObject, eventdata, handles) 
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


% --- Executes on button press in preprocessPatlak_Pushbutton.
function preprocessPatlak_Pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to preprocessPatlak_Pushbutton (see GCBO)
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

tchecked = get(handles.tfit_Checkbox,'Value');
tfit = str2double(get(handles.tfit_Textbox,'String'));
if isnan(tfit)
    warndlg('t* must be a positive number','Preprocessing error','modal');
    uicontrol(handles.tfit_Textbox);
    return
end
% if (tfit > QMmaindata.times(:,2)/60)
%     warndlg('t* must be lower than PET duration','Preprocessing error','modal');
%     uicontrol(handles.tfit_Textbox);
%     return
% end

minR = 0;
maxR = 0;
restrict = get(handles.restrict_Checkbox,'Value');
if restrict == 1
    minR = str2double(get(handles.lowerRest_Textbox,'String'));
    if isnan(minR)
        warndlg('Lower restriction must be a positive number','Preprocessing error','modal');
        uicontrol(handles.lowerRest_Textbox);
        return
    end
    
    maxR = str2double(get(handles.upperRest_Textbox,'String'));
    if isnan(maxR)
        warndlg('Upper restriction must be a positive number','Preprocessing error','modal');
        uicontrol(handles.upperRest_Textbox);
        return
    end
    
    if minR > maxR
        warndlg('Upper restriction must be bigger than Lower restriction','Preprocessing error','modal');
        uicontrol(handles.upperRest_Textbox);
        return
    end
end

Max_err = str2double(get(handles.maxErr_Textbox,'String'));
if isnan(Max_err)
    warndlg('MaxErr must be a positive number','Preprocessing error','modal');
    uicontrol(handles.maxErr_Textbox);
    return
end

threshold = str2double(get(handles.threshold_Textbox,'String'));
if isnan(threshold) || (threshold < 0) || (threshold > 100)
    warndlg('threshold must be a positive number and less than 100','Preprocessing error','modal');
    uicontrol(handles.threshold_Textbox);
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
fprintf(QMf_logfile,strcat('\n-Patlak Reference-\n'));
if tchecked && ~restrict
    fprintf(QMf_logfile,strcat(aux_date(end-8:end),' Preprocessing model (t* checked= ',num2str(tfit),...
    ' without restrictions ',', Max_err= ',num2str(Max_err),')\n'));
elseif tchecked
    fprintf(QMf_logfile,strcat(aux_date(end-8:end),' Preprocessing model (t* checked= ',num2str(tfit),...
    ' with restrictions -',' MinRestrict= ',num2str(minR),', MaxRestrict= ',num2str(minR),'-, Max_err= ',num2str(Max_err),')\n'));
else
    fprintf(QMf_logfile,strcat(aux_date(end-8:end),' Preprocessing model (t* not checked ',', Max_err= ',num2str(Max_err),')\n'));
end


try
    [results plots] = QM_preprocessPatlak(QMmaindata.Crois{2,idCt},QMmaindata.Crois{2,idCr},QMmaindata.times,...
                                            [tchecked tfit],[restrict minR maxR],Max_err,waitbarhandle);
catch ME
    if strcmp(ME.identifier,'PreprocessPatlak:MaxErrCondition')
        errordlg(char({'Not full fill MaxErr condition';'Please, change MaxErr value.'}),'Preprocessing error','modal');
        uicontrol(handles.maxErr_Textbox);
        fprintf(QMf_logfile,strcat('  ERROR: Not full fill MaxErr condition','\n'));
    elseif strcmp(ME.identifier,'PreprocessPatlak:NumberDataPointsFitCondition')
        errordlg(char({'The t* value only takes the last frame';'Please, change t* value.'}),'Preprocessing error','modal');
        uicontrol(handles.maxErr_Textbox);
        fprintf(QMf_logfile,strcat('  The t* value only takes the last frame','\n'));
    else
        errordlg(char({'An unexpected error ocurred in the preprocessing';'Please, check logfile.txt'}),'Preprocessing error','modal');
        fprintf(QMf_logfile,strcat('  ERROR: Preprocessing did not finish succesfully','\n','   Error message: ',ME.message,'\n',...
        '   File: ',ME.stack(1).file,'\n','   Line: ',num2str(ME.stack(1).line),'\n'));
    end
    error = 1;
end

delete(waitbarhandle);  % Need to delete before ploting

if error == 0
    %Update the logfile
    aux_date=datestr(now);
    fprintf(QMf_logfile,strcat(aux_date(end-8:end),' Preprocessing finished \n'));
    
    %Obtain handles using GUIDATA with the caller's handle 
    mainhandles = guidata(handles.mainGUI);
    
    %Plotting the preprocess
    plot(mainhandles.axes,plots(:,1),plots(:,2),'ro',plots(:,1),plots(:,3),'b');
    xlabel(mainhandles.axes,'int(Cr)/Cr')
    ylabel(mainhandles.axes,'Ct/Cr')
    title(mainhandles.axes,'Patlak fit')
    
     %Save values (it is used, for example, to show the values of the TACs when the "view Values" option is selected
%     set(mainhandles.axes,'UserData',{'Patlak',[plots(:,1) plots(:,2) plots(:,3)]});
    t_aux=table([1:1:length(plots(:,1))]',plots(:,1),plots(:,2),plots(:,3));
    set(mainhandles.axes,'UserData',{'Patlak',table2cell(t_aux)});
    
    %Showing results in text area
    strresults = cell(1,7);
    strresults{1,1} = 'Patlak preprocessing:';
    strresults{1,2} = '-----------------------------------';
    strresults{1,3}  = ['t*:                ',num2str(results{1})];
    strresults{1,4}  = ['K:                 ',num2str(results{2},'%.5f'), '     (SE: ',num2str(sqrt(results{8}(1)),'%1.2e'),')'];
    strresults{1,5}  = ['intercept:   ',num2str(results{3},'%.5f'),'     (SE: ',num2str(sqrt(results{8}(2)),'%1.2e'),')'];
%     strresults{1,4}  = ['K:                 ',num2str(results{2},'%.5f'), '     (Var: ',num2str(results{8}(1),'%1.2e'),')'];
%     strresults{1,5}  = ['intercept:   ',num2str(results{3},'%.5f'),'     (Var: ',num2str(results{8}(2),'%1.2e'),')'];
    strresults{1,6}  = ['Start:          ',num2str(results{4},'%.5f')];
    strresults{1,7}  = ['Max. Error: ',num2str(results{5},'%.5f')];
    strresults2 = cell(1,4);
    strresults2{1,1} = 'Goodness of fit:';
    strresults2{1,2} = '-----------------------------------';
    strresults2{1,3} = ['NMSE:           ',num2str(results{6},'%.5f')];
    strresults2{1,4} = ['Corr. Coef.:  ',num2str(results{7},'%.5f')];
    
    set(mainhandles.results_Text,'String',strresults);
    set(mainhandles.results_Text2,'String',strresults2);
    
    %NUEVO
    if (tchecked == 0) && (results{1} ~= tfit)
        warndlg(char({'t* has been selected discarding some inital frames because they provoke 0 or NaN'}),'Preprocessing warning','modal');
        fprintf(QMf_logfile,strcat('  WARNING: t* has been selected discarding some inital frames because they provoke 0 or NaN','\n'));
    end
    
    %Save last preprocess
    QMlastPreprocess.PatlakRef.idCt = idCt;
    QMlastPreprocess.PatlakRef.idCr = idCr;
    QMlastPreprocess.PatlakRef.t = results(1);
    QMlastPreprocess.PatlakRef.threshold = threshold;
    QMlastPreprocess.PatlakRef.plots = plots;
    QMlastPreprocess.PatlakRef.results = results;
    
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
    
    %Minimize the current GUI
    if QMmaindata.flags.minWindow == 1
        minfig(handles.PatlakPreprocess_figure,1);
    end
%     jFig = get(handles.PatlakPreprocess_figure, 'JavaFrame');
%     jFig.setMinimized(true);
%     set(handles.PatlakPreprocess_figure,'Position',[0 0 0 0]);
    
end

% --- Executes on button press in tfit_Checkbox.
function tfit_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to tfit_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tfit_Checkbox

checked = get(hObject,'Value');
if checked
    enable = 'on';
    set(handles.tfit_Textbox,'Enable','off');
else
    enable = 'off';
    set(handles.tfit_Textbox,'Enable','on');
    set(handles.restrict_Checkbox,'Value',0);
    set(handles.lowerRest_Textbox,'Enable',enable);
    set(handles.upperRest_Textbox,'Enable',enable);
end
set(handles.restrict_Checkbox,'Enable',enable);


function tfit_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to tfit_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tfit_Textbox as text
%        str2double(get(hObject,'String')) returns contents of tfit_Textbox as a double

value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('t* value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function tfit_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tfit_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in maxErr_Checkbox.
function maxErr_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to maxErr_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of maxErr_Checkbox



function maxErr_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to maxErr_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxErr_Textbox as text
%        str2double(get(hObject,'String')) returns contents of maxErr_Textbox as a double
value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('Max Error value must be a positive number','Warning','modal');
    uicontrol(hObject);
end


% --- Executes during object creation, after setting all properties.
function maxErr_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxErr_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in restrict_Checkbox.
function restrict_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to restrict_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of restrict_Checkbox

checked = get(hObject,'Value');
if checked
    enable = 'on';
else
    enable = 'off';
end
set(handles.lowerRest_Textbox,'Enable',enable);
set(handles.upperRest_Textbox,'Enable',enable);


function lowerRest_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to lowerRest_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lowerRest_Textbox as text
%        str2double(get(hObject,'String')) returns contents of lowerRest_Textbox as a double

value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('Lower value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function lowerRest_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lowerRest_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function upperRest_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to upperRest_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of upperRest_Textbox as text
%        str2double(get(hObject,'String')) returns contents of upperRest_Textbox as a double

value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('Upper value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function upperRest_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to upperRest_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in setDefault_Pushbutton.
function setDefault_Pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to setDefault_Pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.CroiTac_Popupmenu,'Value',1);
set(handles.CrefTac_Popupmenu,'Value',2);

set(handles.tfit_Checkbox,'Value',1);
set(handles.tfit_Textbox,'String','40.0');
set(handles.tfit_Textbox,'Enable','off');

set(handles.restrict_Checkbox,'Value',0);
set(handles.restrict_Checkbox,'Enable','on');

set(handles.lowerRest_Textbox,'String','0.0');
set(handles.lowerRest_Textbox,'Enable','off');
set(handles.upperRest_Textbox,'String','60.0');
set(handles.upperRest_Textbox,'Enable','off');

set(handles.maxErr_Textbox,'String','10.0');


% --- Executes on button press in help_Pushbutton.
function help_Pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to help_Pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global QMpath

open(strcat(QMpath,filesep,'help',filesep,'html',filesep,'help_patlak.html'))


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
% crea una matriz g, con elementos de la matriz "a", pero considerando intervalos x y y
g=a(1:x:end,1:y:end,:); 
g(g==0)=240; %Para poner el borde de la imagen al mismo color del fondo que tengo puesto
set(hObject,'CData',g); 


% --- Executes on button press in OK_Pushbutton.
function OK_Pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to OK_Pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiresume(handles.PatlakPreprocess_figure);


% --- Executes when user attempts to close PatlakPreprocess_figure.
function PatlakPreprocess_figure_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to PatlakPreprocess_figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure

uiresume(handles.PatlakPreprocess_figure);


% --- Executes on button press in threshold_Checkbox.
function threshold_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to threshold_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of threshold_Checkbox



function threshold_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to threshold_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of threshold_Textbox as text
%        str2double(get(hObject,'String')) returns contents of threshold_Textbox as a double

threshold = str2double(get(hObject,'String'));
if isnan(threshold) || (threshold < 0) || (threshold > 100)
    warndlg('Threshold must be a positive number and less than 100','Warning','modal');
    uicontrol(hObject);
    return
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
