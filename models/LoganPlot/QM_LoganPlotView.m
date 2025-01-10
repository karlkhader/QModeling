function varargout = QM_LoganPlotView(varargin)
% QM_LOGANPLOTVIEW MATLAB code for QM_LoganPlotView.fig
%
%   QM_LoganPlotView generates a GUI to set the parameters
%   used in the preprocessing of Logan Plot model. This code must be 
%   executed from QModeling GUI.  
%
% Last Modified by GUIDE v2.5 29-Aug-2017 12:13:15

% QM_LoganPlotView is part of QModeling.
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
                   'gui_OpeningFcn', @QM_LoganPlotView_OpeningFcn, ...
                   'gui_OutputFcn',  @QM_LoganPlotView_OutputFcn, ...
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


% --- Executes just before QM_LoganPlotView is made visible.
function QM_LoganPlotView_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to QM_LoganPlotView (see VARARGIN)
global QMmaindata
% Choose default command line output for QM_LoganPlotView
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

%Save mainGUI handle
mainGuiInput = find(strcmp(varargin, 'mainhandles'));
handles.mainGUI = varargin{mainGuiInput+1};


%Set the popupmenus
num_rois = size(QMmaindata.Crois,2);
strmenu = cell(1,num_rois);
%copia los nombres de las ROIs en strmenu{}
for i = 1:num_rois
    strmenu{1,i}  = QMmaindata.Crois{1,i};
end

set(handles.PopupMenuROI,'String',strmenu);
set(handles.PopupMenuROI,'Enable','on');
set(handles.PopupMenuROI,'Value',QMmaindata.idCt);

set(handles.PopupMenuREF,'String',strmenu);
set(handles.PopupMenuREF,'Value',QMmaindata.idCr);
set(handles.PopupMenuREF,'Enable','on');

% Enable preprocess button
set(handles.LoganPreprocess,'Enable','on');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes QM_LoganPlotView for user response (see UIRESUME)
uiwait(handles.LoganPrepFigure);

% --- Outputs from this function are returned to the command line.
function varargout = QM_LoganPlotView_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

delete(hObject);


% --- Executes on selection change in listbox1.
function PopupMenuROI_Callback(hObject, eventdata, handles)
% hObject    handle to PopupMenuROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1

idCt = get(hObject,'Value');
idCr = get(handles.PopupMenuREF,'Value');
%Si se selecciona un idCt=idCr lanza un error porque no pueden ser iguales.
if idCt == idCr
    warndlg(char({'Selecting the same TAC for both ROI and REF tac may be incorrect';'Please, change one TAC.'}),'Warning','modal');
    uicontrol(hObject);
end
% --- Executes during object creation, after setting all properties.
function PopupMenuROI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in PopupMenuREF.
function PopupMenuREF_Callback(hObject, eventdata, handles)
% hObject    handle to PopupMenuREF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns PopupMenuREF contents as cell array
%        contents{get(hObject,'Value')} returns selected item from PopupMenuREF
idCt = get(hObject,'Value');
idCr = get(handles.PopupMenuROI,'Value');
%Si se selecciona un idCt=idCr lanza un error porque no pueden ser iguales.
if idCt == idCr
    warndlg(char({'Selecting the same TAC for both ROI and REF tac may be incorrect';'Please, change one TAC.'}),'Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function PopupMenuREF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PopupMenuREF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in LoganPreprocess.
function LoganPreprocess_Callback(hObject, eventdata, handles)
% hObject    handle to LoganPreprocess (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


global QMmaindata
global QMlastPreprocess
global QMf_logfile

idCt = get(handles.PopupMenuROI,'Value');
idCr = get(handles.PopupMenuREF,'Value');

%check del t*
tCh=get(handles.checkboxT,'Value');

k2=str2double(get(handles.k2,'String'));

if idCt == idCr
    warndlg(char({'For a good preprocessing, the region of interest (TAC1) must be different than reference region (TAC2)';...
                    'Please, change one.'}),'Preprocessing error','modal');
    uicontrol(handles.CrefTac_Popupmenu);
    return
end

%check of t upper and lower
tRC=get(handles.checkboxR,'Value');

%check of max error
mTCheck=get(handles.maxErrorCheckbox,'value');
if mTCheck
    gError = str2double(get(handles.MaxError,'String'));
else
    gError=100;
end



tRest=str2double(get(handles.tRestrict_Textbox,'String'));
tLwr=str2double(get(handles.tLwr_Textbox,'String'));
tUpr=str2double(get(handles.tUpr_Textbox,'String'));


if isnan(gError)
    warndlg('Max error must be a positive number','Preprocessing error','modal');
    uicontrol(handles.MaxError);
    return
end
if isnan(tRest)
    warndlg('t* must be a positive number','Preprocessing error','modal');
    uicontrol(handles.tRestrict_Textbox);
    return
end
if tRC == 1
    if isnan(tLwr)
        warndlg('Lower restriction must be a positive number','Preprocessing error','modal');
        uicontrol(handles.tLwr_Textbox);
        return
    end
    

    if isnan(tUpr)
        warndlg('Upper restriction must be a positive number','Preprocessing error','modal');
        uicontrol(handles.tUpr_Textbox);
        return
    end
    
    if tLwr > tUpr
        warndlg('The upper restriction must be bigger than the lower restriction','Preprocessing error','modal');
        uicontrol(handles.tUpr_Textbox);
        return
    end
end


if isnan(gError)
    warndlg('Max Error must be a positive number','Preprocessing error','modal');
    uicontrol(handles.MaxError);
    return
end

threshold = str2double(get(handles.thresholdTextbox,'String'));
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

times=QMmaindata.times;
Cref=QMmaindata.Crois{2,idCr};
Ct=QMmaindata.Crois{2,idCt};


trYerror=0;
    %Update the logfile
    aux_date=datestr(now);
    fprintf(QMf_logfile,strcat('\n-Logan Plot-\n'));
    if tCh && ~tRC
        fprintf(QMf_logfile,strcat(aux_date(end-8:end),' Preprocessing model (t* checked= ',num2str(tRest),...
        ' without restrictions ',', Max error= ',num2str(gError),')\n'));
    elseif tCh
        fprintf(QMf_logfile,strcat(aux_date(end-8:end),' Preprocessing model (t* checked= ',num2str(tRest),...
        ' with restrictions -',' MinRestrict= ',num2str(tLwr),', MaxRestrict= ',num2str(tLwr),'-, Max_err= ',num2str(gError),')\n'));
    else
        fprintf(QMf_logfile,strcat(aux_date(end-8:end),' Preprocessing model (t* not checked ',', Max_err= ',num2str(gError),')\n'));
    end
    

try
    [results, plots] = QM_preprocessLoganPlot(Ct, Cref, times, gError, k2,[tRest tLwr tUpr],[tRC tCh]);
catch ME

    if strcmp(ME.identifier, 'PreprocessLoganPlot:MaxErrCondition')
        errordlg(char({'Not full fill Max Error condition';'Please, change Max error value.'}),'Preprocessing error','modal');
        fprintf(QMf_logfile,strcat('  ERROR: Not full fill Max error condition','\n'));
        
    elseif strcmp(ME.identifier,'PreprocessLoganPlot:NumberDataPointsFitCondition')
        errordlg(char({'The t* value only takes the last frame';'Please, change t* value.'}),'Preprocessing error','modal');
        fprintf(QMf_logfile,strcat('  The t* value only takes the last frame','\n'));
        
    elseif strcmp(ME.identifier,'PreprocessLoganPlot:LowerTimeIncorrect')
        errordlg(char({'Lower time restriction only take the last frame';'Please, change lower time restriction value.'}),'Preprocessing error','modal');
        fprintf(QMf_logfile,strcat('Lower time restriction only take the last frame','\n'));

     else
       errordlg(char({'An unexpected error ocurred in the preprocessing';'Please, check logfile.txt'}),'Preprocessing error','modal');
       fprintf(QMf_logfile,strcat('  ERROR: Preprocessing did not finish succesfully','\n','   Error message: ',ME.message,'\n',...
        '   File: ',ME.stack(1).file,'\n','   Line: ',num2str(ME.stack(1).line),'\n'));
    end
    trYerror = 1;
end

delete(waitbarhandle);

  if trYerror == 0
    %Update the logfile
    aux_date=datestr(now);
    fprintf(QMf_logfile,strcat(aux_date(end-8:end),' Preprocessing finished \n'));
   
    
    %OBTAIN RESULTS OF PREPROCESS
    BPnd = results{1};
    intercept = results{2};
    tr=results{3};
%     start=results(4);
%     CCoef=results(5);
%     NMSE=results(6);
    CCoef=results{4};
    NMSE=results{5};
    
    %Obtain handles using GUIDATA with the caller's handle 
    mainhandles = guidata(handles.mainGUI);
    

    plot(mainhandles.axes,plots(:,1),plots(:,2),'ro',plots(:,1),plots(:,3),'b');
    set(mainhandles.axes,'UserData',{'LoganPlot',[plots(:,1) plots(:,2)]});
    xlabel(mainhandles.axes,strcat('[int(Ct)+Ct/k2',char(39),']/Ct'))
    ylabel(mainhandles.axes,'int(Ct)/Ct')
    title(mainhandles.axes,'Logan Fit')
    
    %Save values (it is used, for example, to show the values of the TACs when the "view Values" option is selected
%     set(mainhandles.axes,'UserData',{'LoganPlot',[plots(:,1) plots(:,2) plots(:,3)]});
    t_aux=table([1:1:length(plots(:,1))]',plots(:,1),plots(:,2),plots(:,3));
    set(mainhandles.axes,'UserData',{'LoganPlot',table2cell(t_aux)});
    
    

    
    %SHOW RESULTS
    strresults = cell(1,6);
    strresults{1,1} = 'Logan Plot preprocessing';
    strresults{1,2} = '-----------------------------------';
    strresults{1,3} = ['BPnd:        ',num2str(BPnd,'%.5f'), '             (SE: ',num2str(sqrt(results{6}(1)),'%1.2e'),')'];
    strresults{1,4} = ['intercept: ',num2str(intercept,'%.5f'), '     (SE: ',num2str(sqrt(results{6}(2)),'%1.2e'),')'];
%     strresults{1,3} = ['BPnd:        ',num2str(BPnd,'%.5f'), '             (Var: ',num2str(results{6}(1),'%1.2e'),')'];
%     strresults{1,4} = ['intercept: ',num2str(intercept,'%.5f'), '     (Var: ',num2str(results{6}(2),'%1.2e'),')'];
    strresults{1,5} = ['k2',char(39),':            ' , num2str(k2)]; 
    strresults{1,6} = ['t*:              ',num2str(tr)];
    %strresults{1,7} = ['start : ',num2str(start),' [min]'];
    strresults2 = cell(1,4);
    strresults2{1,1} = 'Goodness of Fit';  
    strresults2{1,2} = '-----------------------------------';
    strresults2{1,3} = ['Corr. Coef.:  ', num2str(CCoef,'%.5f')];
    strresults2{1,4} = ['NMSE:           ', num2str(NMSE,'%.5f')];

    set(mainhandles.results_Text2,'String',strresults2);
    set(mainhandles.results_Text,'String',strresults);
    
    %SAVE THE LAST PREPROCESS
    QMlastPreprocess.LoganPlot.k2 = k2;
    QMlastPreprocess.LoganPlot.idCt = idCt;
    QMlastPreprocess.LoganPlot.idCr = idCr;
%     QMlastPreprocess.LoganPlot.t = results(4);
    QMlastPreprocess.LoganPlot.t = results{3};
    QMlastPreprocess.LoganPlot.threshold = threshold;
    QMlastPreprocess.LoganPlot.plots = plots;
    QMlastPreprocess.LoganPlot.results = results;

    
    %Restart same GUI components
    set(mainhandles.plottedCt_Popupmenu,'Value',idCt);
    set(mainhandles.plottedCr_Popupmenu,'Value',idCr);
    set(mainhandles.plottedCt_Popupmenu,'Enable','off');
    set(mainhandles.plottedCr_Popupmenu,'Enable','off');
    set(mainhandles.showTacs_Radiobutton,'Value',0);
    set(mainhandles.showPreprocess_Radiobutton,'Value',1);
    set(mainhandles.showTacs_Radiobutton,'Visible','on');
    set(mainhandles.showPreprocess_Radiobutton,'Visible','on');
    
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
        minfig(handles.LoganPrepFigure,1);
    end
   end

function MaxError_Callback(hObject, eventdata, handles)
% hObject    handle to MaxError (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MaxError as text
%        str2double(get(hObject,'String')) returns contents of MaxError as a double
maxE = str2double(get(hObject,'String'));
if isnan(maxE)
    warndlg('Max Error value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function MaxError_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaxError (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function k2_Callback(hObject, eventdata, handles)
% hObject    handle to k2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k2 as text
%        str2double(get(hObject,'String')) returns contents of k2 as a double


% --- Executes during object creation, after setting all properties.
function k2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tRestrict_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to tRestrict_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tRestrict_Textbox as text
%        str2double(get(hObject,'String')) returns contents of tRestrict_Textbox as a double


% --- Executes during object creation, after setting all properties.
function tRestrict_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tRestrict_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton.
function radiobutton_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton


% --- Executes on button press in checkboxT.
function checkboxT_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxT
ch = get(hObject,'Value');
if ch
    set(handles.tLwr_Textbox,'Enable','off');
    set(handles.tUpr_Textbox,'Enable','off');
    set(handles.tRestrict_Textbox,'Enable','off');
    set(handles.checkboxR,'Enable','on');    
else
    set(handles.tLwr_Textbox,'Enable','off');
    set(handles.tUpr_Textbox,'Enable','off');
    set(handles.checkboxR,'Value',0);
    set(handles.tRestrict_Textbox,'Enable','on');
    set(handles.checkboxR,'Enable','off');
    
end

function tLwr_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to tLwr_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tLwr_Textbox as text
%        str2double(get(hObject,'String')) returns contents of tLwr_Textbox as a double


% --- Executes during object creation, after setting all properties.
function tLwr_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tLwr_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tUpr_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to tUpr_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tUpr_Textbox as text
%        str2double(get(hObject,'String')) returns contents of tUpr_Textbox as a double
value = str2double(get(hObject,'String'));
if isnan(value)
    warndlg('Upper value must be a positive number','Warning','modal');
    uicontrol(hObject);
end



% --- Executes during object creation, after setting all properties.
function tUpr_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tUpr_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkboxR.
function checkboxR_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxR
ch = get(hObject,'Value');
if ch
    enable = 'on';
    
else
    enable = 'off';
    
end
set(handles.tLwr_Textbox,'Enable',enable);
set(handles.tUpr_Textbox,'Enable',enable);


function thresholdTextbox_Callback(hObject, eventdata, handles)
% hObject    handle to thresholdTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thresholdTextbox as text
%        str2double(get(hObject,'String')) returns contents of thresholdTextbox as a double


% --- Executes during object creation, after setting all properties.
function thresholdTextbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thresholdTextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in thresholdCheckbox.
function thresholdCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to thresholdCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of thresholdCheckbox



% --- Executes on button press in maxErrorCheckbox.
function maxErrorCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to maxErrorCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of maxErrorCheckbox



% --- Executes on button press in setDefaultPushbutton.
function setDefaultPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to setDefaultPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.PopupMenuROI,'value',1);
set(handles.PopupMenuREF,'value',2);

set(handles.MaxError,'enable','on');
set(handles.MaxError,'String','10.0');
set(handles.maxErrorCheckbox,'value',1);

set(handles.k2,'String','0.15');

set(handles.checkboxT,'value',1);
set(handles.tRestrict_Textbox,'String','40');
set(handles.tRestrict_Textbox,'Enable','off');

set(handles.checkboxR,'value',0);
set(handles.checkboxR,'enable','on');

set(handles.tLwr_Textbox,'String','0.0');
set(handles.tLwr_Textbox,'Enable','off');
set(handles.tUpr_Textbox,'String','60.0');
set(handles.tUpr_Textbox,'Enable','off');
set(handles.thresholdCheckbox,'value',1);
set(handles.thresholdTextbox,'Enable','on');
set(handles.thresholdTextbox,'String','3.0');


% --- Executes on button press in OKpushbutton.
function OKpushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to OKpushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.LoganPrepFigure);


% --- Executes when user attempts to close LoganPrepFigure.
function LoganPrepFigure_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to LoganPrepFigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
uiresume(handles.LoganPrepFigure);


% --- Executes on button press in HelpPushbutton.
function HelpPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to HelpPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global QMpath

open(strcat(QMpath,filesep,'help',filesep,'html',filesep,'help_LoganPlot.html'))
