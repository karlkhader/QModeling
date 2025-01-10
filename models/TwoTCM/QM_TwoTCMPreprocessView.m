function varargout = QM_TwoTCMPreprocessView(varargin)
% QM_TWOTCMPREPROCESSVIEW MATLAB code for QM_TwoTCMPreprocessView.fig
%
%   QM_TwoTCMPreprocessView generates a GUI to set the parameters
%   used in the preprocessing of 2-tissue compartment model. 
%   This code must be executed from QModeling GUI.
%
% Last Modified by GUIDE v2.5 11-Jun-2018 14:49:04
%---------------------------------------------------------------------------------
% QM_TwoTCMPreprocessView is part of QModeling.
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
                   'gui_OpeningFcn', @QM_TwoTCMPreprocessView_OpeningFcn, ...
                   'gui_OutputFcn',  @QM_TwoTCMPreprocessView_OutputFcn, ...
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


% --- Executes just before QM_TwoTCMPreprocessView is made visible.
function QM_TwoTCMPreprocessView_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to QM_TwoTCMPreprocessView (see VARARGIN)
global QMmaindata

% Choose default command line output for QM_TwoTCMPreprocessView
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
    strmenu{1,i} = QMmaindata.Crois{1,i};
end

set(handles.CroiTac_Popupmenu,'String',strmenu);
set(handles.CroiTac_Popupmenu,'Enable','on');
set(handles.CroiTac_Popupmenu,'Value',QMmaindata.idCt);

set(handles.PlasmaTac_Popupmenu,'String',strmenu);
set(handles.PlasmaTac_Popupmenu,'Value',QMmaindata.idCr);
set(handles.PlasmaTac_Popupmenu,'Enable','on');

%Enable preprocess button
set(handles.preprocessTwoTCM_Pushbutton,'Enable','on');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes QM_TwoTCMView wait for user response (see UIRESUME)
uiwait(handles.TwoTCMPreprocess_figure);


% --- Outputs from this function are returned to the command line.
function varargout = QM_TwoTCMPreprocessView_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

delete(hObject);

% --- Executes on button press in help_Pushbutton.
function help_Pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to help_Pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global QMpath
% CUANDO ESTÉ HECHO EL HTML DE LA AYUDA, HABILITAR ESTA LINEA
% open(strcat(QMpath,filesep,'help',filesep,'html',filesep,'help_TwoTCM.html'))

% --- Executes during object creation, after setting all properties.
function help_Pushbutton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to help_Pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
global QMpath
%Create a button
[a,map]=imread(strcat(QMpath,filesep,'icons',filesep,'csh_icon.png'));
[r,c,d]=size(a); 
x=ceil(r/38);
y=ceil(c/51);
g=a(1:x:end,1:y:end,:); 
g(g==0)=240; 
set(hObject,'CData',g);



% --- Executes on button press in preprocessTwoTCM_Pushbutton.
function preprocessTwoTCM_Pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to preprocessTwoTCM_Pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global QMmaindata
global QMlastPreprocess
global QMf_logfile

idCt = get(handles.CroiTac_Popupmenu,'Value');
idCplasma = get(handles.PlasmaTac_Popupmenu,'Value');

if idCt == idCplasma
    warndlg(char({'For a good preprocessing, the region of interest (TAC1) must be different than Plasma region (TAC2)';...
                    'Please, change one.'}),'Preprocessing error','modal');
    uicontrol(handles.CrefTac_Popupmenu);
    return
end

alpha1_rest_checkbox=get(handles.alpha1_restrict_Checkbox,'Value');
if alpha1_rest_checkbox
    alpha1_lower = str2double(get(handles.alpha1_lowerRest_Textbox,'String'));
    if isnan(alpha1_lower) || (alpha1_lower < 0)
        warndlg('The lower value of alpha 1 must be a positive number','Preprocessing error','modal');
        uicontrol(handles.alpha1_lowerRest_Textbox);
        return
    end

    alpha1_upper = str2double(get(handles.alpha1_upperRest_Textbox,'String'));
    if isnan(alpha1_upper) || (alpha1_upper < 0)
        warndlg('The upper value of alpha 1 must be a positive number','Preprocessing error','modal');
        uicontrol(handles.alpha1_upperRest_Textbox);
        return
    end
    
    if alpha1_upper < alpha1_lower
        warndlg('The upper value of alpha 1 must be bigger than the lower value of alpha 1','Preprocessing error','modal');
        uicontrol(handles.alpha1_upperRest_Textbox);
        return
    end
else
    alpha1 = str2double(get(handles.alpha1_Textbox,'String'));    
    if isnan(alpha1) || (alpha1 < 0)
        warndlg('alpha 1 value must be a positive number','Preprocessing error','modal');
        uicontrol(handles.alpha1_Textbox);
        return
    else
        alpha1_lower = alpha1;
        alpha1_upper = alpha1;
    end
end
    

alpha2_rest_checkbox=get(handles.alpha2_restrict_Checkbox,'Value');
if alpha2_rest_checkbox
    alpha2_lower = str2double(get(handles.alpha2_lowerRest_Textbox,'String'));
    if isnan(alpha2_lower) || (alpha2_lower < 0)
        warndlg('The lower value of alpha 2 must be a positive number','Preprocessing error','modal');
        uicontrol(handles.alpha2_lowerRest_Textbox);
        return
    end

    alpha2_upper = str2double(get(handles.alpha2_upperRest_Textbox,'String'));
    if isnan(alpha2_upper) || (alpha2_upper < 0)
        warndlg('The upper value of alpha 2 must be a positive number','Preprocessing error','modal');
        uicontrol(handles.alpha2_upperRest_Textbox);
        return
    end
    
    if alpha2_upper < alpha2_lower
        warndlg('The upper value of alpha 2 must be bigger than the lower value of alpha 2','Preprocessing error','modal');
        uicontrol(handles.alpha2_upperRest_Textbox);
        return
    end
else
    alpha2 = str2double(get(handles.alpha2_Textbox,'String'));    
    if isnan(alpha2) || (alpha2 < 0)
        warndlg('alpha 2 value must be a positive number','Preprocessing error','modal');
        uicontrol(handles.alpha2_Textbox);
        return
    else
        alpha2_lower = alpha2;
        alpha2_upper = alpha2;
    end
end

vB_checked=get(handles.vB_Checkbox,'Value');
if ~vB_checked %vB is fixed
    vB = str2num(get(handles.vB_Textbox,'String'));
    if isnan(vB) || (vB < 0)
        warndlg('vB must be a positive integer','Preprocessing error','modal');
        uicontrol(handles.basis_Textbox);
        return
    end
else %vB will be fitted
    vB=Inf;
end

k4_checked=get(handles.k4_Checkbox,'Value');

basisF = str2num(get(handles.basis_Textbox,'String'));
if isnan(basisF) || (basisF < 1)
    warndlg('The number of basis functions must be a positive integer','Preprocessing error','modal');
    uicontrol(handles.basis_Textbox);
    return
end

resampling = str2double(get(handles.resampling_Textbox,'String'));
if isnan(resampling) || (resampling < 0)
    warndlg('Resampling must be a positive number','Preprocessing error','modal');
    uicontrol(handles.resampling_Textbox);
    return
end

threshold = str2double(get(handles.threshold_Textbox,'String'));
if isnan(threshold) || (threshold < 0) || (threshold > 100)
    warndlg('Threshold must be a positive number less than 100','Preprocessing error','modal');
    uicontrol(handles.threshold_Textbox);
    return
end

waitbarhandle = waitbar(0,'Preprocessing, please wait...');
ver=version('-release');
if (str2num(ver(1:4))<2014) || isequal(ver,'2014a')
    set(findobj(waitbarhandle,'type','patch'),'edgecolor','b','facecolor','b');
else
    wbc = allchild(waitbarhandle);     %you need to get at a hidden child
    wbc(1).JavaPeer.setForeground( wbc(1).JavaPeer.getBackground.BLUE )
    wbc(1).JavaPeer.setStringPainted(true)
end

error = 0;

%Update the logfile
aux_date=datestr(now);
fprintf(QMf_logfile,strcat('\n-2TCM-\n'));
if ~vB_checked
    vB_value=strcat('= ',num2str(vB));
else
    vB_value=' will be fitted';
end
if ~k4_checked
    k4_value='0';
else
    k4_value=' will be fitted';
end
fprintf(QMf_logfile,strcat(aux_date(end-8:end),' Preprocessing model (alpha1_min= ',num2str(alpha1_lower),...
    ', alpha1_max= ',num2str(alpha1_upper),', alpha2_min= ',num2str(alpha2_lower),...
    ', alpha2_max= ',num2str(alpha2_upper),', vB ',vB_value,', k4 ',k4_value,', BF= ',num2str(basisF),', resampling= ',num2str(resampling),')\n'));

try
    [results, plots] = QM_preprocessTwoTCM(QMmaindata.Crois{2,idCt},QMmaindata.Crois{2,idCplasma},QMmaindata.times,...
                                            [vB_checked vB],k4_checked,[alpha1_rest_checkbox,alpha1_lower,alpha1_upper],...
                                            [alpha2_rest_checkbox,alpha2_lower,alpha2_upper],resampling,basisF,waitbarhandle);
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
    
    if isempty(QMlastPreprocess.TwoTCM.Cblood)
        plot(mainhandles.axes,tmin,QMmaindata.Crois{2,idCt},'b.-',tmin,QMmaindata.Crois{2,idCplasma},'r*-',tmin,plots,'g+-');
        xlabel(mainhandles.axes,'Time (minutes)');
        ylabel(mainhandles.axes,'Original units');
        title(mainhandles.axes,'2-Tissue Compartment Model fit')
        legend(mainhandles.axes,'TAC1 from receptor-rich region', 'TAC2 from plasma region', 'Fit to TAC1');
        %Save values (it is used, for example, to show the values of the
        %TACs when the "view Values" option is selected
%         set(mainhandles.axes,'UserData',{'TwoTCM',[tmin QMmaindata.Crois{2,idCt} plots QMmaindata.Crois{2,idCplasma}]});
        t_aux=table([1:1:length(plots)]',tmin,QMmaindata.Crois{2,idCt},plots,QMmaindata.Crois{2,idCplasma});
        set(mainhandles.axes,'UserData',{'TwoTCM*',table2cell(t_aux)});
    else %optionally, a whole blood activity TAC has been included
        plot(mainhandles.axes,tmin,QMmaindata.Crois{2,idCt},'b.-',tmin,QMmaindata.Crois{2,idCplasma},'r*-',tmin,plots,'g+-',tmin,QMlastPreprocess.TwoTCM.Cblood,'y--');
        xlabel(mainhandles.axes,'Time (minutes)');
        ylabel(mainhandles.axes,'Original units');
        title(mainhandles.axes,'2-Tissue Compartment Model fit')
        legend(mainhandles.axes,'TAC1 from receptor-rich region', 'TAC2 from plasma region', 'Fit to TAC1','Whole blood activity');
        %Save values (it is used, for example, to show the values of the
        %TACs when the "view Values" option is selected
%         set(mainhandles.axes,'UserData',{'TwoTCM',[tmin QMmaindata.Crois{2,idCt} plots QMmaindata.Crois{2,idCplasma} QMlastPreprocess.TwoTCM.Cblood]});
        t_aux=table([1:1:length(plots(:,1))]',tmin,QMmaindata.Crois{2,idCt},plots,QMmaindata.Crois{2,idCplasma},QMlastPreprocess.TwoTCM.Cblood);
        set(mainhandles.axes,'UserData',{'TwoTCM',table2cell(t_aux)});
    end
        
     
    %Showing results in text area
    strresults = cell(1,6);
    strresults{1,1} = '2TCM preprocessing:';
    strresults{1,2} = '-----------------------------------';
    if ~vB_checked
        strresults{1,3}  = ['vB: ',num2str(results{1},'%.5f'), '     (fixed)'];
    else
        strresults{1,3}  = ['vB: ',num2str(results{1},'%.5f'), '     (SE: ',num2str(sqrt(results{9}(1)),'%1.2e'),')'];
    end
    strresults{1,4}  = ['K1: ',num2str(results{2},'%.5f'), '     (SE: ',num2str(sqrt(results{9}(2)),'%1.2e'),')'];
    strresults{1,5}  = ['k2: ',num2str(results{3},'%.5f'), '     (SE: ',num2str(sqrt(results{9}(3)),'%1.2e'),')'];
    strresults{1,6}  = ['k3: ',num2str(results{4},'%.5f'), '     (SE: ',num2str(sqrt(results{9}(4)),'%1.2e'),')'];
    if ~k4_checked
        strresults{1,7}  = ['k4: ',num2str(results{5},'%.5f'), '     (fixed)'];
    else
        strresults{1,7}  = ['k4: ',num2str(results{5},'%.5f'), '     (SE: ',num2str(sqrt(results{9}(5)),'%1.2e'),')'];
    end
%     strresults{1,7}  = ['Ki: ',num2str(results{6},'%.5f')];
    
    strresults2 = cell(1,4);
    strresults2{1,1} = 'Goodness of fit:';
    strresults2{1,2} = '-----------------------------------';
    strresults2{1,3} = ['NMSE:           ',num2str(results{6},'%.5f')];
    strresults2{1,4} = ['Corr. Coef.:  ',num2str(results{7},'%.5f')];
    
    set(mainhandles.results_Text,'String',strresults);
    set(mainhandles.results_Text2,'String',strresults2);
    
    %Save last preprocess
    QMlastPreprocess.TwoTCM.idCt = idCt;
    QMlastPreprocess.TwoTCM.idCplasma = idCplasma; 
    QMlastPreprocess.TwoTCM.vB=[vB_checked vB];
    QMlastPreprocess.TwoTCM.k4_checked=k4_checked;
%     QMlastPreprocess.TwoTCM.B=B;
%     QMlastPreprocess.TwoTCM.alpha=alpha;
    QMlastPreprocess.TwoTCM.threshold=threshold;
    QMlastPreprocess.TwoTCM.plots = plots;
    QMlastPreprocess.TwoTCM.results=results;
    
       
    %Restart same GUI components
    set(mainhandles.plottedCt_Popupmenu,'Value',idCt);
    set(mainhandles.plottedCr_Popupmenu,'Value',idCplasma);
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
        minfig(handles.TwoTCMPreprocess_figure,1);
    end
        
end

% --- Executes on button press in OK_Pushbutton.
function OK_Pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to OK_Pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDAT
uiresume(handles.TwoTCMPreprocess_figure);


% --- Executes when user attempts to close TwoTCMPreprocess_figure.
function TwoTCMPreprocess_figure_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to TwoTCMPreprocess_figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
uiresume(handles.TwoTCMPreprocess_figure);


% --- Executes on selection change in CroiTac_Popupmenu.
function CroiTac_Popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to CroiTac_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns CroiTac_Popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from CroiTac_Popupmenu

idCt = get(hObject,'Value');
idCplasma = get(handles.PlasmaTac_Popupmenu,'Value');

if idCt == idCplasma
    warndlg(char({'Selecting the same TAC for both ROI and Plasma tac may be incorrect';'Please, change one TAC.'}),'Warning','modal');
    uicontrol(hObject);
end

% --- Executes on selection change in PlasmaTac_Popupmenu.
function PlasmaTac_Popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to PlasmaTac_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns PlasmaTac_Popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from PlasmaTac_Popupmenu

idCplasma = get(hObject,'Value');
idCt = get(handles.CroiTac_Popupmenu,'Value');

if idCt == idCplasma
    warndlg(char({'Selecting the same TAC for both ROI and Plasma tac may be incorrect';'Please, change one TAC.'}),'Warning','modal');
    uicontrol(hObject);
end


% --- Executes on button press in setDefault_Pushbutton.
function setDefault_Pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to setDefault_Pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.CroiTac_Popupmenu,'Value',1);
set(handles.PlasmaTac_Popupmenu,'Value',2);
set(handles.optionalBloodTAC_Textbox,'String','');

set(handles.alpha1_Textbox,'String','0.005');
set(handles.alpha1_Textbox,'Enable','on');
set(handles.alpha2_Textbox,'String', '0.1');
set(handles.alpha2_Textbox,'Enable','on');

set(handles.alpha1_restrict_Checkbox,'Value',0);
set(handles.alpha1_lowerRest_Textbox,'String','0.0005');
set(handles.alpha1_upperRest_Textbox,'String','0.0015');
set(handles.alpha1_upperRest_Textbox,'Enable','off');
set(handles.alpha1_lowerRest_Textbox,'Enable','off');
set(handles.alpha1_lowerText,'Enable','off');
set(handles.alpha1_upperText,'Enable','off');
set(handles.alpha1_text,'Enable','off');

set(handles.alpha2_restrict_Checkbox,'Value',0);
set(handles.alpha2_lowerRest_Textbox,'String', '0.06');
set(handles.alpha2_upperRest_Textbox,'String', '0.6');
set(handles.alpha2_upperRest_Textbox,'Enable','off');
set(handles.alpha2_lowerRest_Textbox,'Enable','off');
set(handles.alpha2_lowerText,'Enable','off');
set(handles.alpha2_upperText,'Enable','off');
set(handles.alpha2_text,'Enable','off');

set(handles.vB_Checkbox,'Value',0)
set(handles.vB_Textbox,'String','0.05');
set(handles.vB_Textbox,'Visible','on');
set(handles.vBfitted_warnText,'Visible','off');

set(handles.k4_Checkbox,'Value',0)
set(handles.k4_Textbox,'String','0.0');
set(handles.k4_Textbox,'Visible','on');
set(handles.text_k4,'Visible','on');
set(handles.k4fitted_warnText,'Visible','off');

set(handles.basis_Textbox,'String','100');
set(handles.basis_Textbox,'Enable','off');
set(handles.resampling_Textbox,'String','1.0');
set(handles.threshold_Textbox,'String','3.0');



function optionalBloodTAC_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to optionalBloodTAC_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of optionalBloodTAC_Textbox as text
%        str2double(get(hObject,'String')) returns contents of optionalBloodTAC_Textbox as a double


% --- Executes during object creation, after setting all properties.
function optionalBloodTAC_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to optionalBloodTAC_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in loadOptionalMask_Pushbutton.
function loadOptionalMask_Pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to loadOptionalMask_Pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global QMmaindata
global QMlastPreprocess


file = spm_select([1,1],'any','Select TAC (.txt, .xls or .mat)');

if ~isempty(file)
    [~, name, ext] = fileparts(file);
    set(handles.optionalBloodTAC_Textbox,'String',strcat(name,ext));
    
    if strcmp(ext,'.mat')  
        try
            Cb_TAC=load(file);
            if (length(Cb_TAC.TACs)~=QMmaindata.num_files)
                warndlg(char({'Dimensions of the loaded PET study and'; 'the optional blood TAC dismatch'}),'Warning','modal');
                uicontrol(hObject);
            else
                QMlastPreprocess.TwoTCM.Cblood{1,1}=Cb_TAC.TACs_names{1};
                QMlastPreprocess.TwoTCM.Cblood{2,1}=Cb_TAC.TACs(:,1);
            end
        catch
            warndlg(char({'Something wrong in data text'}),'Warning','modal');
            uicontrol(hObject);            
        end            
    elseif strcmp(ext,'.txt')
        Cb_TAC=importdata(file);
        if (length(Cb_TAC.data)~=QMmaindata.num_files)
            warndlg(char({'Dimensions of the loaded PET study and'; 'the optional blood TAC dismatch'}),'Warning','modal');
            uicontrol(hObject);
        else           
            QMlastPreprocess.TwoTCM.Cblood{1,1}=Cb_TAC.textdata{1};
            QMlastPreprocess.TwoTCM.Cblood{2,1}=Cb_TAC.data(:,1);
        end         
        
    elseif (strcmp(ext,'.xls') || strcmp(ext,'.xlsx'))
        try
            [Cb_TAC, Cb_TAC_name, ~]=xlsread(file);
            if (length(Cb_TAC)~=QMmaindata.num_files)
                warndlg(char({'Dimensions of the loaded PET study and'; 'the optional blood TAC dismatch'}),'Warning','modal');
                uicontrol(hObject);
            else                
                QMlastPreprocess.TwoTCM.Cblood{1,1}=Cb_TAC_name;
                QMlastPreprocess.TwoTCM.Cblood{2,1}=Cb_TAC(:,1);
            end
        catch
            warndlg(char({'Something wrong in data text'}),'Warning','modal');
            uicontrol(hObject);            
        end            
                
    else
        warndlg(char({ 'TAC file format not accepted'}),'Warning','modal');
        uicontrol(hObject);
        set(handles.optionalBloodTAC_Textbox,'String','');
       
        
    end

else
     QMlastPreprocess.TwoTCM.Cblood=[];
end
        

function alpha1_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to alpha1_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alpha1_Textbox as text
value = str2double(get(hObject,'String'));
if isnan(value) | (value < 0)
    warndlg('alpha 1 value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function alpha1_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alpha1_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function alpha2_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to alpha2_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alpha2_Textbox as text
value = str2double(get(hObject,'String'));
if isnan(value) | (value < 0)
    warndlg('alpha 2 lower value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function alpha2_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alpha2_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in alpha1_restrict_Checkbox.
function alpha1_restrict_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to alpha1_restrict_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of alpha1_restrict_Checkbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
    set(handles.alpha1_Textbox,'Enable','off');
else
    enable = 'off';
    set(handles.alpha1_Textbox,'Enable','on');
end
set(handles.alpha1_lowerText,'Enable',enable);
set(handles.alpha1_upperText,'Enable',enable);
set(handles.alpha1_text,'Enable',enable);
set(handles.alpha1_lowerRest_Textbox,'Enable',enable);
set(handles.alpha1_upperRest_Textbox,'Enable',enable);

set(handles.basis_Textbox,'Enable',enable);




% --- Executes on button press in alpha2_restrict_Checkbox.
function alpha2_restrict_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to alpha2_restrict_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of alpha2_restrict_Checkbox
checked = get(hObject,'Value');
if checked
    enable = 'on';
    set(handles.alpha2_Textbox,'Enable','off');
else
    enable = 'off';
    set(handles.alpha2_Textbox,'Enable','on');
end
set(handles.alpha2_lowerText,'Enable',enable);
set(handles.alpha2_upperText,'Enable',enable);
set(handles.alpha2_text,'Enable',enable);
set(handles.alpha2_lowerRest_Textbox,'Enable',enable);
set(handles.alpha2_upperRest_Textbox,'Enable',enable);

set(handles.basis_Textbox,'Enable',enable);

function alpha1_lowerRest_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to alpha1_lowerRest_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alpha1_lowerRest_Textbox as text
value = str2double(get(hObject,'String'));
if isnan(value) | (value < 0)
    warndlg('alpha 1 lower value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function alpha1_lowerRest_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alpha1_lowerRest_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function alpha2_lowerRest_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to alpha2_lowerRest_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alpha2_lowerRest_Textbox as text
value = str2double(get(hObject,'String'));
if isnan(value) | (value < 0)
    warndlg('alpha 2 lower value must be a positive number','Warning','modal');
    uicontrol(hObject);
end


% --- Executes during object creation, after setting all properties.
function alpha2_lowerRest_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alpha2_lowerRest_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function alpha1_upperRest_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to alpha1_upperRest_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alpha1_upperRest_Textbox as text
value = str2double(get(hObject,'String'));
if isnan(value) | (value < 0)
    warndlg('alpha 1 upper value must be a positive number','Warning','modal');
    uicontrol(hObject);
end

% --- Executes during object creation, after setting all properties.
function alpha1_upperRest_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alpha1_upperRest_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function alpha2_upperRest_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to alpha2_upperRest_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alpha2_upperRest_Textbox as text
value = str2double(get(hObject,'String'));
if isnan(value) | (value < 0)
    warndlg('alpha 2 lower value must be a positive number','Warning','modal');
    uicontrol(hObject);
end


% --- Executes during object creation, after setting all properties.
function alpha2_upperRest_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alpha2_upperRest_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function vB_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to vB_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vB_Textbox as text
%        str2double(get(hObject,'String')) returns contents of vB_Textbox as a double


% --- Executes during object creation, after setting all properties.
function vB_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vB_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
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
    set(handles.k4_Textbox,'Visible','off');
    set(handles.text_k4,'Visible','off');
    set(handles.k4fitted_warnText,'Visible','on');    
else     
    set(handles.k4_Textbox,'Visible','on');
    set(handles.text_k4,'Visible','on');
%     set(handles.k4_Textbox,'Enable','off');
    set(handles.k4fitted_warnText,'Visible','off');
end



function k4_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to k4_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k4_Textbox as text
%        str2double(get(hObject,'String')) returns contents of k4_Textbox as a double


% --- Executes during object creation, after setting all properties.
function k4_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k4_Textbox (see GCBO)
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



function resampling_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to resampling_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of resampling_Textbox as text
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



function threshold_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to threshold_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of threshold_Textbox as text
value = str2double(get(hObject,'String'));
if isnan(value)|| (value < 0)  || (value > 100)
    warndlg('Threshold value must be a positive number less than 100','Warning','modal');
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


% --- Executes on button press in vB_Checkbox.
function vB_Checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to vB_Checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of vB_Checkbox
checked = get(hObject,'Value');
if checked
    set(handles.vB_Textbox,'Visible','off');
    set(handles.vBfitted_warnText,'Visible','on');    
else     
    set(handles.vB_Textbox,'Visible','on');
    %set(handles.vB_Textbox,'Enable','on');
    set(handles.vBfitted_warnText,'Visible','off');
end
