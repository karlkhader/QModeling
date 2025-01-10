function varargout = QModeling(varargin)
% QMODELING MATLAB code for QModeling.fig
%
% QModeling v 1.8 is an open source program, developed as a toolbox for SPM 
% (Stadistical Parametric Imaging), to make kinetic analysis for PET studies 
% based on compartmental models. 
%---------------------------------------------------------------------------------
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
% Copyright (C) 2014, 2018 Francisco Javier López González, Jose Paredes Pacheco, 
% Karl-Khader Thurnhofer Hemsi, Núria Roé Vellvé, Antonio L. Gutierrez
% Cardo, Manuel Enciso Garcia Oliveros, Carlos Rossi Jimenez

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @QModeling_OpeningFcn, ...
                   'gui_OutputFcn',  @QModeling_OutputFcn, ...
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


% --- Executes just before QModeling is made visible.
function QModeling_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to QModeling (see VARARGIN)
startQM()

global QMmaindata
global QMlastPreprocess
global QMpath

global QMf_logfile

%We save the temporary folder path
actualpath = mfilename('fullpath');
QMpath = actualpath(1:end-10);
if exist(strcat(QMpath,filesep,'temp'),'dir') ~= 7
    mkdir(QMpath,'temp');
end

%Create/open and initialize the log file
QM_initiateLogFile();

fprintf(QMf_logfile,strcat('\n','----------------------------------------------','\n'));
fprintf(QMf_logfile,date);
fprintf(QMf_logfile,strcat('\n','----------------------------------------------','\n'));
aux_date=datestr(now);
fprintf(QMf_logfile,strcat(aux_date(end-7:end),' Starting QModeling \n'));


%Add folders to path
addpath(QMpath,strcat(QMpath,filesep,'lib'),...
    strcat(QMpath,filesep,'lib',filesep,'nifti'),...
    strcat(QMpath,filesep,'lib',filesep,'wintools'),...
    strcat(QMpath,filesep,'lib',filesep,'CopyPaste'),...
    strcat(QMpath,filesep,'models'),...
    strcat(QMpath,filesep,'models',filesep,'SRTM'),...
    strcat(QMpath,filesep,'models',filesep,'SRTM2'),...
    strcat(QMpath,filesep,'models',filesep,'PatlakRef'),...
    strcat(QMpath,filesep,'models',filesep,'LoganPlot'),...
    strcat(QMpath,filesep,'models',filesep,'TwoTCM'),...
    strcat(QMpath,filesep,'icons'),...
    strcat(QMpath,filesep,'help'),...
    strcat(QMpath,filesep,'help',filesep,'html'),...
    strcat(QMpath,filesep,'command_line'));

%Publish focus function 
handles.setFocus = @setPanelFocus;
handles.deleteTempImages=@deleteTempImages;

% Choose default command line output for QModeling
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

format long
%Define global variable which saves main data for the toolbox working
QMmaindata = struct('num_files',0,'times',0,'Cpet',0,'Crois',0,'idCt',1,'idCr',2,'flags',struct('consistency',0,'minWindow',1,'studies',0));
QMlastPreprocess = struct('SRTM',struct('idCt',1,'idCr',2,'k2amin',0,'k2amax',0,'threshold',0,'M',[],'B',[],'numBF',0,'plots',[],'results',[]),...
    'SRTM2',struct('idCt',1,'idCr',2,'k2_pChecked',false,'B',0,'th3',0,'threshold',0,'plots',[],'results',[]),...
    'PatlakRef',struct('idCt',1,'idCr',2,'t',[],'plots',[],'results',[]),...
    'LoganPlot',struct('idCt',1,'idCr',2,'t',[],'plots',[],'results',[]),...
    'TwoTCM',struct('idCt',1,'idCplasma',2,'Cblood',[],'vB',[],'k4_checked',false,'M',[],'B',[],'alpha',[],'threshold',0,'plots',[],'results',[]));


% UIWAIT makes QModeling wait for user response (see UIRESUME)
% uiwait(handles.QModeling_figure);

% --- Outputs from this function are returned to the command line.
function varargout = QModeling_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

deleteTempImages(handles);

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in selectModel_Popupmenu.
function selectModel_Popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to selectModel_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns selectModel_Popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from selectModel_Popupmenu

global QMlastPreprocess
global QMmaindata
%global QMmaindata

selected = get(handles.selectModel_Popupmenu,'Value');
if selected == 1, model = 'SRTM'; end
if selected == 2, model = 'SRTM2'; end
if selected == 3, model = 'PatlakRef'; end
if selected == 4, model = 'LoganPlot'; end
if selected == 5, model = 'TwoTCM'; end

data = getfield(QMlastPreprocess,model);

if isempty(data.plots)
    
    Refresh_Tacs(handles,data); 
    
    set(handles.setMapParam_Pushbutton,'Enable','off')
    set(handles.selectImage_Popupmenu,'Enable','off');
    set(handles.viewImage_Pushbutton,'Enable','off');
    set(handles.showTacs_Radiobutton,'Visible','off')
    set(handles.showPreprocess_Radiobutton,'Visible','off')
    set(handles.results_Text,'String','');
    set(handles.results_Text2,'String','');
    set(handles.saveResults_Pushbutton,'Enable','off')
    set(handles.corrMatrix_Pushbutton,'Enable','off')
    
    %Set panel focus
    setPanelFocus(handles,'ModelingPanel');
    
else
    
    Refresh_Plots_Results(handles, selected);    
    set(handles.setMapParam_Pushbutton,'Enable','on')    
    
    if QMmaindata.flags.studies==selected
        set(handles.selectImage_Popupmenu,'Enable','on');
        set(handles.viewImage_Pushbutton,'Enable','on');
        
        %Set panel focus
        setPanelFocus(handles,'ViewParamImgPanel');
    else
        set(handles.selectImage_Popupmenu,'Enable','off');
        set(handles.viewImage_Pushbutton,'Enable','off');
        
        %Set panel focus
        setPanelFocus(handles,'PixelCalcPanel');
    end
        
    set(handles.showTacs_Radiobutton,'Visible','on')
    set(handles.showPreprocess_Radiobutton,'Visible','on')
    set(handles.showTacs_Radiobutton,'Value',0)
    set(handles.showPreprocess_Radiobutton,'Value',1)
    
end


% --- Executes during object creation, after setting all properties.
function selectModel_Popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to selectModel_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in preprocessModel_Pushbutton.
function preprocessModel_Pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to preprocessModel_Pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global QMmaindata
global QMf_logfile

error = 0;

%Set panel focus
setPanelFocus(handles,'ModelingPanel');

if size(QMmaindata.Crois,2) == 1
%     errordlg(char({'To use reference models you must have al least 2 diferent TACs';'View help.'}),'Preprocessing error','modal');
%     fprintf(QMf_logfile,strcat('  ERROR: Insufficient number of TACs to use reference models','\n'));
    errordlg(char({'To use reference models you must have at least 2 diferent TACs';...
        'In other case, the TAC for the region of interest and the blood/plasma radioactivity concentration are required';...
        'View help.'}),'Preprocessing error','modal');
    fprintf(QMf_logfile,strcat('  ERROR: Insufficient number of TACs','\n'));
    return;
end

set(handles.setMapParam_Pushbutton,'Enable','off');

set(hObject,'Enable','off');

selected = get(handles.selectModel_Popupmenu,'Value');
switch selected
    case 1,
        try
            QM_SRTMPreprocessView('mainhandles',handles.QModeling_figure);
            set(handles.corrMatrix_Pushbutton,'Enable','on');
        catch ME
            errordlg(char({'An error ocurred during the preprocess:';'Please, check logfile.txt'}),'Preprocessing error','modal');
            error = 1;
            fprintf(QMf_logfile,strcat('  ERROR: During the preprocess','\n','   Error message: ',ME.message,'\n',...
            '   File: ',ME.stack(1).file,'\n','   Line: ',num2str(ME.stack(1).line),'\n'));
        end
    case 2,
        try
            QM_CarsonPreprocessView('mainhandles',handles.QModeling_figure);
            set(handles.corrMatrix_Pushbutton,'Enable','on');
        catch ME
            errordlg(char({'An error ocurred during the preprocess:';'Please, check logfile.txt'}),'Preprocessing error','modal');
            error = 1;
            fprintf(QMf_logfile,strcat('  ERROR: During the preprocess','\n','   Error message: ',ME.message,'\n',...
            '   File: ',ME.stack(1).file,'\n','   Line: ',num2str(ME.stack(1).line),'\n'));
        end
    case 3,
        try
            QM_patlakPreprocessView('mainhandles',handles.QModeling_figure);
            set(handles.corrMatrix_Pushbutton,'Enable','off');
        catch ME
            errordlg(char({'An error ocurred during the preprocess:';'Please, check logfile.txt'}),'Preprocessing error','modal');
            error = 1;
            fprintf(QMf_logfile,strcat('  ERROR: During the preprocess','\n','   Error message: ',ME.message,'\n',...
            '   File: ',ME.stack(1).file,'\n','   Line: ',num2str(ME.stack(1).line),'\n'));
        end

    case 4,
         try
            QM_LoganPlotView('mainhandles',handles.QModeling_figure);
            set(handles.corrMatrix_Pushbutton,'Enable','off');
        catch ME
            errordlg(char({'An error ocurred during the preprocess:';'Please, check logfile.txt'}),'Preprocessing error','modal');
            error = 1;
            fprintf(QMf_logfile,strcat('  ERROR: During the preprocess','\n','   Error message: ',ME.message,'\n',...
            '   File: ',ME.stack(1).file,'\n','   Line: ',num2str(ME.stack(1).line),'\n'));
         end
         
    case 5,
         try
            QM_TwoTCMPreprocessView('mainhandles',handles.QModeling_figure);
            set(handles.corrMatrix_Pushbutton,'Enable','on');
        catch ME
            errordlg(char({'An error ocurred during the preprocess:';'Please, check logfile.txt'}),'Preprocessing error','modal');
            error = 1;
            fprintf(QMf_logfile,strcat('  ERROR: During the preprocess','\n','   Error message: ',ME.message,'\n',...
            '   File: ',ME.stack(1).file,'\n','   Line: ',num2str(ME.stack(1).line),'\n'));
         end
end

if error ~= 1
    
    %Set visible plots options
    set(handles.plottedCt_Popupmenu,'Enable','off');
    set(handles.plottedCr_Popupmenu,'Enable','off');
    
    if ~isempty(get(handles.results_Text,'String'))
        
        % If showTacs_Radiobutton selected, enable selection tacs
        if get(handles.showTacs_Radiobutton,'Value') == 1
            set(handles.plottedCt_Popupmenu,'Enable','on');
            set(handles.plottedCr_Popupmenu,'Enable','on');
        end
        
        %Set panel focus
        setPanelFocus(handles,'PixelCalcPanel');   
        %Set enable pixel-wise calculation button
        set(handles.setMapParam_Pushbutton,'Enable','on')
        set(handles.saveResults_Pushbutton,'Enable','on')
        
    else  % If nothing preprocessed, enable selection tacs
        set(handles.plottedCt_Popupmenu,'Enable','on');
        set(handles.plottedCr_Popupmenu,'Enable','on');
    end

end

set(hObject,'Enable','on');


% --- Executes during object creation, after setting all properties.
function LoadData_uipanel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LoadData_uipanel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

set(hObject,'BorderType','line','HighlightColor',[0.043;0.518;0.78],'BorderWidth',3);


% --- Executes on button press in loadStudy_Pushbutton.
function loadStudy_Pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to loadStudy_Pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global QMmaindata
% global reorient
global QMf_logfile
error = 0;
% reorient = false;

set(handles.setMapParam_Pushbutton,'UserData',0);
% set(handles.plottedCt_Popupmenu,'Enable','off');
% set(handles.plottedCr_Popupmenu,'Enable','off');
set(handles.editTime_Pushbutton,'Enable','off');
set(handles.loadStudy_Pushbutton,'Enable','off');
set(handles.selectMask_Pushbutton,'Enable','off');

%Set panel focus
setPanelFocus(handles,'LoadPanel');


%Load study
paths = spm_select([1,600], 'image', 'Select .img or .nii images from PET studio ','','','.*',[1 600]);

if ~isempty(paths)
    
    %Reset main GUI
    set(handles.prepareTacs_Pushbutton,'Enable','off');
    cla(handles.axes,'reset');
    set(handles.showTacs_Radiobutton,'Visible','off')
    set(handles.showPreprocess_Radiobutton,'Visible','off')
    set(handles.preprocessModel_Pushbutton,'Enable','off');
    set(handles.results_Text,'String','');
    set(handles.results_Text2,'String','');
    set(handles.setMapParam_Pushbutton,'Enable','off');
    set(handles.selectImage_Popupmenu,'Enable','off');
    set(handles.viewImage_Pushbutton,'Enable','off');
    set (handles.selectModel_Popupmenu,'Enable','off');

    num_files = size(paths,1);
    
    %Save in UserData from study_Textbox the PET paths
    set(handles.study_Textbox,'UserData',paths);
     
    if num_files ~= 0
        
        %Show PET name
        auxstr = '';
        auxpath = paths(1,:);
        auxpath = auxpath(find(auxpath == filesep,1,'last')+1:end-2);
        if length(auxpath)+11 > 34
            auxpath = strcat(auxpath(1:20),'...');
        end
        auxstr = strcat(auxstr,auxpath,' (',num2str(num_files),' files)');
        
    end

    set(handles.study_Textbox,'String',auxstr);

    %Update the logfile
    aux_date=datestr(now);
    fprintf(QMf_logfile,strcat(aux_date(end-8:end),' Loading PET study -',auxstr,'-\n'));
    
    waitbarhandle = waitbar(0,'Loading PET study, please wait...');
    ver=version('-release');
    if (str2num(ver(1:4))<2014) || isequal(ver,'2014a')
        set(findobj(waitbarhandle,'type','patch'),'edgecolor','b','facecolor','b');
    else
        wbc = allchild(waitbarhandle);     %you need to get at a hidden child
        wbc(1).JavaPeer.setForeground( wbc(1).JavaPeer.getBackground.BLUE )
        wbc(1).JavaPeer.setStringPainted(true)
    end
    
    if num_files ~= 0
        if num_files == 1 %Case multiframe

            %Open the multiframe PET study
            multifile = paths(1,1:max(strfind(paths,','))-1);      % To raise ',1'
            multifile_hdr = spm_vol(multifile);
            num_files = length(multifile_hdr);
            QMmaindata.num_files = num_files;

            %Initialize the PET times variable
            times = zeros(num_files,2);
            try
                for frame = 1:num_files
                    current_hdr = multifile_hdr(frame);
                    times(frame,1) = current_hdr.private.timing.toffset;
                    times(frame,2) = times(frame,1) + current_hdr.private.timing.tspace;
                    waitbar(frame/num_files,waitbarhandle)
                end
                QMmaindata.times = times;
            catch er
                %Checking if is an internal error
                if (isfield(current_hdr,'private') && isfield(current_hdr.private,'timing')...
                        && isfield(current_hdr.private.timing,'toffset') && isfield(current_hdr.private.timing,'tspace'))
                    
                    errordlg(char({'An error ocurred reading the PET times.';'Please, edit times.'}),'Reading error','modal');
                    uicontrol(handles.editTime_Pushbutton);
                    error = 1;
                    fprintf(QMf_logfile,strcat('  ERROR Reading the PET times','\n'));
                else        %No times found in study, we insert default timetable
                    actualpath = mfilename('fullpath');
                    actualpath = actualpath(1:end-10);
                    aux = load(strcat(actualpath,filesep,'TimeTable.txt'));
                    if length(aux)<num_files %length of Timetable is not sufficient. There are more frames than start and end times at the Timetable.txt
                        errordlg(char({'There are more frames at the study than start and end times at the default TimeTable file. Please, edit the default TimeTable file.'}),'Reading error','modal');
                        uicontrol(handles.editTime_Pushbutton);
                        error = 1;
                        fprintf(QMf_logfile,strcat('  ERROR There are more frames at the study than start and end times at the default TimeTable file.','\n'));
                    else
                        times = aux(1:num_files,:);
                        QMmaindata.times = times;
                        warndlg(char({'No times found in your study. Default times has been loaded';...
                                        'If you like to change it, select Change times'}),'Reading warning','modal');
                        fprintf(QMf_logfile,strcat('  WARNING: No times found in the study. Default times has been loaded','\n'));
                    end
                end
            end
            delete(waitbarhandle)
            
        else

            QMmaindata.num_files = num_files;
            %Initialize the PET times variable
            times = zeros(num_files,2);

            try
                for frame = 1:num_files
                    %Open the current frame for the PET study
                    file = paths(frame,:);
                    file_hdr = spm_vol(file);
%                     %Check if PET study is reoriented
%                     if(reorient == false && (~isequal(file_hdr.mat(1,2:3),[0,0]) ||...
%                        file_hdr.mat(2,1)~=0 || file_hdr.mat(2,3)~=0 || ~isequal(file_hdr.mat(3,1:2),[0,0]) ||...
%                        ~isequal(file_hdr.mat(4,1:3),[0,0,0])))
%                        warndlg(char({'The study is not oriented';'The parametric images will be show reoriented'}),'Reading warning','modal');
%                        reorient = true;
%                        fprintf(QMf_logfile,strcat('  WARNING: The study is not oriented. The parametric images will be show reoriented','\n'));
%                     end
                    
                    times(frame,1) = file_hdr.private.timing.toffset;
                    times(frame,2) = times(frame,1) + file_hdr.private.timing.tspace;
                    waitbar(frame/num_files,waitbarhandle)
                end
                QMmaindata.times = times;
            catch er
                %Checking if is an internal error
                if (isfield(file_hdr,'private') && isfield(file_hdr.private,'timing')...
                        && isfield(file_hdr.private.timing,'toffset') && isfield(file_hdr.private.timing,'tspace'))
                    
                    errordlg(char({'An error ocurred reading the PET times.';'Please, edit times.'}),'Reading error','modal');
                    uicontrol(handles.editTime_Pushbutton);
                    error = 1;
                    fprintf(QMf_logfile,strcat('  ERROR Reading the PET times','\n'));
                else        %No times in study, we insert default timetable
                    actualpath = mfilename('fullpath');
                    actualpath = actualpath(1:end-10);
                    aux = load(strcat(actualpath,filesep,'TimeTable.txt'));
                    times = aux(1:num_files,:);
                    QMmaindata.times = times;
                    warndlg(char({'No times found in your study. Default times has been loaded';...
                                    'If you like to change it, select Change times'}),'Reading warning','modal');
                    fprintf(QMf_logfile,strcat('  WARNING: No times found in the study. Default times has been loaded','\n'));
                end
            end
            delete(waitbarhandle)
        end
    end

    if error == 0 && ~isempty(get(handles.mask_Textbox,'String'))
        set(handles.prepareTacs_Pushbutton,'Enable','on');
    end
    set(handles.editTime_Pushbutton,'Enable','on');
    
    %Checking if the PET times are consistent
    checkTimesConsistency();
    
    if error == 0
        %Detele the last parametric images created in temp file
        deleteTempImages(handles)
        aux_date=datestr(now);
        fprintf(QMf_logfile,strcat(aux_date(end-8:end),' Study loaded \n'));
    else
        set(handles.selectImage_Popupmenu,'Enable','on');
    end
end

set(handles.loadStudy_Pushbutton,'Enable','on');
set(handles.selectMask_Pushbutton,'Enable','on');


function study_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to study_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of study_Textbox as text
%        str2double(get(hObject,'String')) returns contents of study_Textbox as a double

if isempty(get(handles.study_Textbox,'String'))
        set(handles.prepareTacs_Pushbutton,'Enable','off');
end


% --- Executes during object creation, after setting all properties.
function study_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to study_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in editTime_Pushbutton.
function editTime_Pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to editTime_Pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

QM_editTimes('mainhandles',handles.QModeling_figure);
checkTimesConsistency();


% --- Executes on button press in Load Mask/TACs_Pushbutton.
function selectMask_Pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to selectMask_Pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Set panel focus
setPanelFocus(handles,'LoadPanel');

file = spm_select([1,1],'any','Select group mask (.img or .nii) or TACs (.txt, .xls or .mat)');

if ~isempty(file)
    [path, name, ext]=fileparts(file);
    
%     if size(file,1) ~= 0
%         auxstr = file(find(file == filesep,1,'last')+1:end-2);
%     end
% 
%     set(handles.mask_Textbox,'String',auxstr);
    set(handles.mask_Textbox,'String',strcat(name,ext));
    
    if ~isempty(get(handles.study_Textbox,'String'))
        set(handles.prepareTacs_Pushbutton,'Enable','on');
    end
    
%     [path name ext] = fileparts(file);
    if (strcmp(ext,'.nii') || strcmp(ext,'.img')) %look for the file with ROI names
 
        codes = strcat(path,filesep,name,'.txt');
        if exist(codes,'file') == 0     %If file doesn't exist

            uiwait(msgbox({'The .txt with the names of the ROIs has not been found.';... 
                'After you close this window, a new window dialog will be open for select it.';...
                'The file format is:';...
                '';...
                   'region_A    id_A';...
                   'region_B    id_B';...
                   'region_C    id_C';...
                   'region_D    id_D';...
                   '...';...
                   '';...
                   'If the .txt file is not selected or an identifier doesnt exist,';...
                   'default names will be loaded for the ROIs.'},'Information','help','modal'));


            codes = spm_select([1,1],'txt','Select a .txt file with the ROIs names. Example per line: occipital   1');
        end
    elseif (strcmp(ext,'.mat')||strcmp(ext,'.txt')||strcmp(ext,'.xls')||strcmp(ext,'.xlsx'))
        codes=ext;
    else
        codes='ERROR';
    end
    
    set(handles.mask_Textbox,'UserData',{file;codes});
% file = spm_select([1,1],'image','Select group mask (.img or .nii)');
% 
% if ~isempty(file)
%     
%     if size(file,1) ~= 0
%         auxstr = file(find(file == filesep,1,'last')+1:end-2);
%     end
% 
%     set(handles.mask_Textbox,'String',auxstr);
%     
%     if ~isempty(get(handles.study_Textbox,'String'))
%         set(handles.prepareTacs_Pushbutton,'Enable','on');
%     end
%     
%     [path name ext] = fileparts(file);
%     codes = strcat(path,filesep,name,'.txt');
%     if exist(codes,'file') == 0     %If file doesn't exist
%         
%         uiwait(msgbox({'The .txt with the names of the ROIs has not been found.';... 
%             'After you close this window, a new window dialog will be open for select it.';...
%             'The file format is:';...
%             '';...
%                'region_A    id_A';...
%                'region_B    id_B';...
%                'region_C    id_C';...
%                'region_D    id_D';...
%                '...';...
%                '';...
%                'If the .txt file is not selected or an identifier doesnt exist,';...
%                'default names will be loaded for the ROIs.'},'Information','help','modal'));
%             
%                
%         codes = spm_select([1,1],'txt','Select a .txt file with the ROIs names. Example per line: occipital   1');
%     end
%     
%     set(handles.mask_Textbox,'UserData',{file;codes});
    
end


function mask_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to mask_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mask_Textbox as text
%        str2double(get(hObject,'String')) returns contents of mask_Textbox as a double

if isempty(get(handles.mask_Textbox,'String'))
        set(handles.prepareTacs_Pushbutton,'Enable','off');
end


% --- Executes during object creation, after setting all properties.
function mask_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mask_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in prepareTACs_Pushbutton.
function prepareTacs_Pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to prepareTacs_Pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global QMmaindata
global QMlastPreprocess
global QMf_logfile
error = 0;

QMlastPreprocess = struct('SRTM',struct('idCt',1,'idCr',2,'k2amin',0,'k2amax',0,'threshold',0,'M',[],'B',[],'numBF',0,'plots',[],'results',[]),...
    'SRTM2',struct('idCt',1,'idCr',2,'B',0,'th3',0,'threshold',0,'plots',[],'results',[]),...
    'PatlakRef',struct('idCt',1,'idCr',2,'t',[],'plots',[],'results',[]),...
    'LoganPlot',struct('idCt',1,'idCr',2,'t',[],'plots',[],'results',[]),...
    'TwoTCM',struct('idCt',1,'idCplasma',2,'Cblood',[],'vB',[],'k4_checked',false,'M',[],'B',[],'alpha',[],'threshold',0,'plots',[],'results',[]));

set(handles.setMapParam_Pushbutton,'UserData',0);

%Set panel focus
setPanelFocus(handles,'LoadPanel');

%Disable the follow steps
set(handles.preprocessModel_Pushbutton,'Enable','off');
set(handles.setMapParam_Pushbutton,'Enable','off');
set(handles.selectImage_Popupmenu,'Enable','off');
set(handles.viewImage_Pushbutton,'Enable','off');
set(handles.saveResults_Pushbutton,'Enable','off');
set(handles.corrMatrix_Pushbutton,'Enable','off');

%Update the logfile
aux_date=datestr(now);
fprintf(QMf_logfile,strcat(aux_date(end-8:end),' Preparing TACs \n'));
    
if QMmaindata.times == 0
    
    errordlg(char({'Cannot prepare TACs because have not times.';'Please, edit times.'}),'Reading error','modal');
    uicontrol(handles.editTime_Pushbutton);
    fprintf(QMf_logfile,strcat('  ERROR: No times loaded','\n'));
    
elseif size(QMmaindata.times,1) ~= QMmaindata.num_files
    
    errordlg(char({'Cannot prepare TACs because the number of frames is diferent than the number of times.';'Please, edit times.'}),'Reading error','modal');
    uicontrol(handles.editTime_Pushbutton);
    fprintf(QMf_logfile,strcat('  ERROR: Diferent number of frames and times','\n'));
    
else
    
    waitbarhandle = waitbar(0,'Preparing TACs, please wait...');
    ver=version('-release');
    if (str2num(ver(1:4))<2014) || isequal(ver,'2014a')
        set(findobj(waitbarhandle,'type','patch'),'edgecolor','b','facecolor','b');
    else
        wbc = allchild(waitbarhandle);     %you need to get at a hidden child
        wbc(1).JavaPeer.setForeground( wbc(1).JavaPeer.getBackground.BLUE )
        wbc(1).JavaPeer.setStringPainted(true)
    end
    
    try
        [QMmaindata.Cpet QMmaindata.Crois num_rois] = QM_prepareTACs(handles,waitbarhandle,QMmaindata.num_files);
    catch ME
        if strcmp(ME.identifier,'PrepareTACs:dimensionsDismatch')
            errordlg(char({'Cannot prepare TACs because the dimensions of the PET and the mask dismatch';'Please, select new mask.'}),'Reading error','modal');
            uicontrol(handles.selectMask_Pushbutton);
            fprintf(QMf_logfile,strcat('  ERROR: The dimensions of the PET and the mask dismatch','\n'));
        elseif strcmp(ME.identifier,'PrepareTACs:fileFormatIncorrect')
            errordlg(char({'Cannot prepare TACs because the file format is not accepted';'Please, select a new correct format TACs file.'}),'Reading error','modal');
            uicontrol(handles.selectMask_Pushbutton);
            fprintf(QMf_logfile,strcat('  ERROR: TACs file format not accepted','\n'));
        elseif strcmp(ME.identifier, 'PrepareTACs:badFormattedData')
            errordlg(char({'Cannot prepare TACs because data format is not correct';'Please, read the manual to ensure your data is well formatted.'}),'Reading error','modal');
            uicontrol(handles.selectMask_Pushbutton);
            fprintf(QMf_logfile,strcat('  ERROR: Data format for TACs data is not correct','\n'));
        else
            errordlg(char({'An unexpected error ocurred reading the PET';'Please, open help.'}),'Reading error','modal');
            fprintf(QMf_logfile,strcat('  ERROR: Reading the PET','\n'));
        end
        error = 1;
    end

    %Plot the TACs

    delete(waitbarhandle);  %Need to delete before ploting
    
    if error == 0
        %Update the logfile
        aux_date=datestr(now);
        fprintf(QMf_logfile,strcat(aux_date(end-8:end),' TACs prepared \n'));

        %Clean results panel
        set(handles.results_Text,'String','');
        set(handles.results_Text2,'String','');
        
        %Average time between frames and translated to minutes
        tavg = ((QMmaindata.times(:,1)+QMmaindata.times(:,2))/2);
        tmin = tavg/60;

        %Reset selected TACs
        QMmaindata.idCt = 1;
        if num_rois == 1
            QMmaindata.idCr = 1;
        else
            QMmaindata.idCr = 2;
        end

        %Plot TACs
%         plot(handles.axes,tmin,QMmaindata.Crois{2,QMmaindata.idCt},'b.-');
%         hold on
%         plot(handles.axes,tmin,QMmaindata.Crois{2,QMmaindata.idCr},'r.-');
%         hold off
%         plot(handles.axes,tmin,QMmaindata.Crois{2,QMmaindata.idCr},'r.-');
%         hold on
%         plot(handles.axes,tmin,QMmaindata.Crois{2,QMmaindata.idCt},'b.-'); 
%         hold off
        plot(handles.axes,tmin,QMmaindata.Crois{2,QMmaindata.idCt},'b.-',tmin,QMmaindata.Crois{2,QMmaindata.idCr},'r.-');
%         hold on
%         plot(handles.axes,tmin,QMmaindata.Crois{2,QMmaindata.idCt},'b.-'); 
%         hold off
        xlabel('Time (minutes)')
        ylabel('Original units')
%         aux = get(handles.mask_Textbox,'UserData');
%         [~, mask_name, mask_ext]=fileparts(aux{1});
%         mask_name = aux{1}(find(aux{1} == filesep,1,'last')+1:end-2);
%         title(['Time-activity curves for mask: ' strcat(mask_name,mask_ext)])
        title('Time-activity curves')
%         set(handles.axes,'UserData',{'TACs',[tmin QMmaindata.Crois{2,QMmaindata.idCt:QMmaindata.idCr}]});
        t_aux=table([1:1:length(tmin)]',tmin,QMmaindata.Crois{2,QMmaindata.idCt:QMmaindata.idCr});
        set(handles.axes,'UserData',{'TACs',table2cell(t_aux)});
        
        %Enable selection TACs popupmenus
        strmenu = cell(1,num_rois);
        for i = 1:num_rois
            strmenu{1,i}  = QMmaindata.Crois{1,i};  %strcat('TAC ',num2str(i));
        end

        set(handles.plottedCt_Popupmenu,'String',strmenu);
        set(handles.plottedCt_Popupmenu,'Enable','on');
        set(handles.plottedCt_Popupmenu,'Value',1);
        
        set(handles.plottedCr_Popupmenu,'String',strmenu);
        if num_rois == 1
            set(handles.plottedCr_Popupmenu,'Value',1);
        else
            set(handles.plottedCr_Popupmenu,'Value',2);
        end
        set(handles.plottedCr_Popupmenu,'Enable','on');
        
        set(handles.showTacs_Radiobutton,'Visible','off')
        set(handles.showPreprocess_Radiobutton,'Visible','off')

        %Set panel focus
        setPanelFocus(handles,'ModelingPanel');
        
        %Enable preprocess buttons
        set(handles.selectModel_Popupmenu,'Enable','on'); 
        set(handles.preprocessModel_Pushbutton,'Enable','on');
        
        %Detele the last parametric images created in temp file
        deleteTempImages(handles)
       
    end
    
end

% --------------------------------------------------------------------
function QModeling_MenuBar_Callback(hObject, eventdata, handles)
% hObject    handle to QModeling_MenuBar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Help_MenuBar_Callback(hObject, eventdata, handles)
% hObject    handle to Help_MenuBar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global QMpath
open(strcat(QMpath,filesep,'help',filesep,'html',filesep,'manual.html'))


% --- Executes on button press in setMapParam_Pushbutton.
function setMapParam_Pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to setMapParam_Pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global QMmaindata
global QMlastPreprocess
global QMf_logfile
error = 0;

%Set panel focus
setPanelFocus(handles,'PixelCalcPanel');

selected = get(handles.selectModel_Popupmenu,'Value');
if selected == 1, model = 'SRTM'; end
if selected == 2, model = 'SRTM2'; end
if selected == 3, model = 'PatlakRef'; end
if selected == 4, model = 'LoganPlot'; end
if selected == 5, model = 'TwoTCM'; end

data = getfield(QMlastPreprocess,model);

if isempty(data.plots)
   errordlg(char({'An error ocurred during the pixel-wise calculation:';'You do not have a preprocess with the model selected.';'View help.'}),'Pixel calculation','modal');
   fprintf(QMf_logfile,strcat('  ERROR: A preprocess with the model selected does not exist','\n'));
else
    switch selected
        case 1,
            %Call SRTM
            try
                QM_mapParamSRTM('mainhandles',handles.QModeling_figure);
            catch ME
                errordlg(char({'An error ocurred during the pixel-wise calculation:';'Please, check logfile.txt'}),'Pixel calculation','modal');
                error = 1;
                fprintf(QMf_logfile,strcat('  ERROR: During the pixel-wise calculation','\n','   Error message: ',ME.message,'\n',...
                '   File: ',ME.stack(1).file,'\n','   Line: ',num2str(ME.stack(1).line),'\n'));
            end
        case 2,
            %Call SRTM2
            try
                QM_mapParamCarson('mainhandles',handles.QModeling_figure);
            catch ME
               errordlg(char({'An error ocurred during the pixel-wise calculation:';'Please, check logfile.txt'}),'Pixel calculation','modal');
               error=1;  
               fprintf(QMf_logfile,strcat('  ERROR: During the pixel-wise calculation','\n','   Error message: ',ME.message,'\n',...
                '   File: ',ME.stack(1).file,'\n','   Line: ',num2str(ME.stack(1).line),'\n'));
            end
        case 3,
            %Call PatlakRef
            try
                QM_mapParamPatlak('mainhandles',handles.QModeling_figure);
            catch ME
                errordlg(char({'An error ocurred during the pixel-wise calculation:';'Please, check logfile.txt'}),'Pixel calculation','modal');
                error = 1;
                fprintf(QMf_logfile,strcat('  ERROR: During the pixel-wise calculation','\n','   Error message: ',ME.message,'\n',...
                '   File: ',ME.stack(1).file,'\n','   Line: ',num2str(ME.stack(1).line),'\n'));
            end
        case 4,
            %Call LoganPlot
            try
                QM_mapParamLogan('mainhandles',handles.QModeling_figure);
            catch ME
                errordlg(char({'An error ocurred during the pixel-wise calculation:';'Please, check logfile.txt'}),'Pixel calculation','modal');
                error = 1;
                fprintf(QMf_logfile,strcat('  ERROR: During the pixel-wise calculation','\n','   Error message: ',ME.message,'\n',...
                '   File: ',ME.stack(1).file,'\n','   Line: ',num2str(ME.stack(1).line),'\n'));
            end
            
        case 5,
            %Call TwoTCM
            try
                QM_mapParamTwoTCM('mainhandles',handles.QModeling_figure);
            catch ME
                errordlg(char({'An error ocurred during the pixel-wise calculation:';'Please, check logfile.txt'}),'Pixel calculation','modal');
                error = 1;
                fprintf(QMf_logfile,strcat('  ERROR: During the pixel-wise calculation','\n','   Error message: ',ME.message,'\n',...
                '   File: ',ME.stack(1).file,'\n','   Line: ',num2str(ME.stack(1).line),'\n'));
            end
    end

    executed=get(hObject, 'UserData');
    if error == 0 && executed == 1
        set(handles.selectImage_Popupmenu,'Enable','on');
        set(handles.viewImage_Pushbutton,'Enable','on');

        QMmaindata.flags.studies=selected;
        
        %Set panel focus
        setPanelFocus(handles,'ViewParamImgPanel');   
    else
        set(handles.selectImage_Popupmenu,'String','Parametric Images');
        set(handles.selectImage_Popupmenu,'Enable','off');
        set(handles.viewImage_Pushbutton,'Enable','off');
    end
end

% --- Executes on selection change in plottedCt_Popupmenu.
function plottedCt_Popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to plottedCt_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns plottedCt_Popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plottedCt_Popupmenu

global QMmaindata

QMmaindata.idCt = get(hObject,'Value');
QMmaindata.idCr = get(handles.plottedCr_Popupmenu,'Value');

tavg = ((QMmaindata.times(:,1)+QMmaindata.times(:,2))/2);
tmin = tavg/60;

plot(handles.axes,tmin,QMmaindata.Crois{2,QMmaindata.idCt},'b.-');
hold on
plot(handles.axes,tmin,QMmaindata.Crois{2,QMmaindata.idCr},'r.-');
hold off
xlabel('Time (minutes)')
ylabel('Original units')
% aux = get(handles.mask_Textbox,'UserData');
% mask_name = aux{1}(find(aux{1} == filesep,1,'last')+1:end-2);
% title(['Time-activity curves for mask: ' mask_name])
title('Time-activity curves');
% set(handles.axes,'UserData',{'TACs',[tmin QMmaindata.Crois{2,QMmaindata.idCt} QMmaindata.Crois{2,QMmaindata.idCr}]});
if QMmaindata.idCt==QMmaindata.idCr
    t_aux=table([1:1:length(tmin)]',tmin,QMmaindata.Crois{2,QMmaindata.idCt});
else
    t_aux=table([1:1:length(tmin)]',tmin,QMmaindata.Crois{2,QMmaindata.idCt}, QMmaindata.Crois{2,QMmaindata.idCr});
end

set(handles.axes,'UserData',{'TACs',table2cell(t_aux)});


% --- Executes during object creation, after setting all properties.
function plottedCt_Popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plottedCt_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in plottedCr_Popupmenu.
function plottedCr_Popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to plottedCr_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns plottedCr_Popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plottedCr_Popupmenu

global QMmaindata

QMmaindata.idCt = get(handles.plottedCt_Popupmenu,'Value');
QMmaindata.idCr = get(hObject,'Value');

tavg = ((QMmaindata.times(:,1)+QMmaindata.times(:,2))/2);
tmin = tavg/60;

plot(handles.axes,tmin,QMmaindata.Crois{2,QMmaindata.idCt},'b.-');
hold on
plot(handles.axes,tmin,QMmaindata.Crois{2,QMmaindata.idCr},'r.-');
hold off
xlabel('Time (minutes)')
ylabel('Original units')
% aux = get(handles.mask_Textbox,'UserData');
% mask_name = aux{1}(find(aux{1} == filesep,1,'last')+1:end-2);
% title(['Time-activity curves for mask: ' mask_name])
title('Time-activity curves')
% set(handles.axes,'UserData',{'TACs',[tmin QMmaindata.Crois{2,QMmaindata.idCt} QMmaindata.Crois{2,QMmaindata.idCr}]});
if QMmaindata.idCt==QMmaindata.idCr
    t_aux=table([1:1:length(tmin)]',tmin,QMmaindata.Crois{2,QMmaindata.idCt});
else
    t_aux=table([1:1:length(tmin)]',tmin,QMmaindata.Crois{2,QMmaindata.idCt}, QMmaindata.Crois{2,QMmaindata.idCr});
end
set(handles.axes,'UserData',{'TACs',table2cell(t_aux)});



% --- Executes during object creation, after setting all properties.
function plottedCr_Popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plottedCr_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in showTacs_Radiobutton.
function showTacs_Radiobutton_Callback(hObject, eventdata, handles)
% hObject    handle to showTacs_Radiobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showTacs_Radiobutton


% --- Executes on button press in showPreprocess_Radiobutton.
function showPreprocess_Radiobutton_Callback(hObject, eventdata, handles)
% hObject    handle to showPreprocess_Radiobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showPreprocess_Radiobutton


% --- Executes when selected object is changed in showingTacs_uibuttongroup.
function showingTacs_uibuttongroup_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in showingTacs_uibuttongroup 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
global QMmaindata
global QMlastPreprocess

modelID = get(handles.selectModel_Popupmenu,'Value');
if modelID == 1, model = 'SRTM'; end
if modelID == 2, model = 'SRTM2'; end
if modelID == 3, model = 'PatlakRef'; end
if modelID == 4, model = 'LoganPlot'; end
if modelID == 5, model = 'TwoTCM'; end

data = getfield(QMlastPreprocess,model);

tavg = ((QMmaindata.times(:,1)+QMmaindata.times(:,2))/2);
tmin = tavg/60;

switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'showTacs_Radiobutton'
        Refresh_Tacs(handles,data);
        
    case 'showPreprocess_Radiobutton'
        Refresh_Plots_Results(handles, modelID);       
end


% --- Executes on selection change in selectImage_Popupmenu.
function selectImage_Popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to selectImage_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns selectImage_Popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from selectImage_Popupmenu


% --- Executes during object creation, after setting all properties.
function selectImage_Popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to selectImage_Popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in viewImage_Pushbutton.
function viewImage_Pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to viewImage_Pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% global reorient 
global QMf_logfile
images = get(handles.selectImage_Popupmenu,'UserData');
imageSel = get(handles.selectImage_Popupmenu,'Value');

try
%     %Reorient the temporal parametric images before show them
%     if reorient == true
%         Reorient_image(images{1,imageSel});
%         [pth nam ext] = fileparts(images{1,imageSel});
%         rimg = [pth, filesep, 'r', nam, ext];
%         nii = load_nii(rimg);
%     else
        nii = load_nii(images{1,imageSel});
%     end
    
    opt.setviewpoint = size(nii.img)/2;
    opt.setcolorindex = 4;
    opt.useinterp = 1;
    view_nii(nii,opt);                    
catch ME
    errordlg(char({'An error ocurred: ';'Please, check logfile.txt'}),'Viewing error','modal');
    fprintf(QMf_logfile,strcat('  ERROR showing parametric image','\n','   Error message: ',ME.message,'\n',...
                '   File: ',ME.stack(1).file,'\n','   Line: ',num2str(ME.stack(1).line),'\n'));
end   


function deleteTempImages(handles)
global QMmaindata

QMmaindata.flags.studies=0;

images = get(handles.selectImage_Popupmenu,'UserData');

if ~isempty(images)
%     delreo=false; %Indicate if a temporal reoriented parametric images .img has been deleted 
    for i = 1:size(images,2)
        [~,~,ext] = fileparts(images{1,i});
        delete(images{1,i});
%         %delete temporal reoriented parametric image .img
%         if (exist([p,filesep,'r',n,ext]))
%             delete([p,filesep,'r',n,ext]); 
%             delreo=true;
%         end
        if strcmp(ext,'.img')
            delete(strcat(images{1,i}(1:end-3),'hdr'));
            %delete temporal reoriented parametric image .hdr
%             if (delreo)
%                 delete([p,filesep,'r',n,'.hdr']);
%                 delreo=false;
%             end
        end
    end
    set(handles.selectImage_Popupmenu,'UserData',{});
    set(handles.selectImage_Popupmenu,'Value',1);
    set(handles.selectImage_Popupmenu,'String','Parametric Images');
end




function setPanelFocus(handles,uipanel)

switch uipanel
    
    case 'LoadPanel'
        
        set(handles.LoadData_uipanel,'BorderType','line','HighlightColor',[0.043;0.518;0.78],'BorderWidth',3);
        set(handles.Modeling_uipanel,'BorderType','etchedin','HighlightColor',[1.0;1.0;1.0],'BorderWidth',1);
        set(handles.Plots_uipanel,'BorderType','etchedin','HighlightColor',[1.0;1.0;1.0],'BorderWidth',1);
        set(handles.Results_uipanel,'BorderType','etchedin','HighlightColor',[1.0;1.0;1.0],'BorderWidth',1);
        set(handles.PixelCalc_uipanel,'BorderType','etchedin','HighlightColor',[1.0;1.0;1.0],'BorderWidth',1);
        set(handles.ViewParamImg_uipanel,'BorderType','etchedin','HighlightColor',[1.0;1.0;1.0],'BorderWidth',1);
        
    case 'ModelingPanel'
        
        set(handles.LoadData_uipanel,'BorderType','etchedin','HighlightColor',[1.0;1.0;1.0],'BorderWidth',1);
        set(handles.Modeling_uipanel,'BorderType','line','HighlightColor',[0.043;0.518;0.78],'BorderWidth',3);
        set(handles.Plots_uipanel,'BorderType','etchedin','HighlightColor',[1.0;1.0;1.0],'BorderWidth',1);
        set(handles.Results_uipanel,'BorderType','etchedin','HighlightColor',[1.0;1.0;1.0],'BorderWidth',1);
        set(handles.PixelCalc_uipanel,'BorderType','etchedin','HighlightColor',[1.0;1.0;1.0],'BorderWidth',1);
        set(handles.ViewParamImg_uipanel,'BorderType','etchedin','HighlightColor',[1.0;1.0;1.0],'BorderWidth',1);
        
    case 'PixelCalcPanel'
        
        set(handles.LoadData_uipanel,'BorderType','etchedin','HighlightColor',[1.0;1.0;1.0],'BorderWidth',1);
        set(handles.Modeling_uipanel,'BorderType','etchedin','HighlightColor',[1.0;1.0;1.0],'BorderWidth',1);
        set(handles.Plots_uipanel,'BorderType','etchedin','HighlightColor',[1.0;1.0;1.0],'BorderWidth',1);
        set(handles.Results_uipanel,'BorderType','etchedin','HighlightColor',[1.0;1.0;1.0],'BorderWidth',1);
        set(handles.PixelCalc_uipanel,'BorderType','line','HighlightColor',[0.043;0.518;0.78],'BorderWidth',3);
        set(handles.ViewParamImg_uipanel,'BorderType','etchedin','HighlightColor',[1.0;1.0;1.0],'BorderWidth',1);
        
    case 'PlotResultsPanels'
        
        set(handles.LoadData_uipanel,'BorderType','etchedin','HighlightColor',[1.0;1.0;1.0],'BorderWidth',1);
        set(handles.Modeling_uipanel,'BorderType','etchedin','HighlightColor',[1.0;1.0;1.0],'BorderWidth',1);
        set(handles.Plots_uipanel,'BorderType','line','HighlightColor',[0.043;0.518;0.78],'BorderWidth',3);
        set(handles.Results_uipanel,'BorderType','line','HighlightColor',[0.043;0.518;0.78],'BorderWidth',3);
        set(handles.PixelCalc_uipanel,'BorderType','etchedin','HighlightColor',[1.0;1.0;1.0],'BorderWidth',1);
        set(handles.ViewParamImg_uipanel,'BorderType','etchedin','HighlightColor',[1.0;1.0;1.0],'BorderWidth',1);
        
    case 'ViewParamImgPanel'
        
        set(handles.LoadData_uipanel,'BorderType','etchedin','HighlightColor',[1.0;1.0;1.0],'BorderWidth',1);
        set(handles.Modeling_uipanel,'BorderType','etchedin','HighlightColor',[1.0;1.0;1.0],'BorderWidth',1);
        set(handles.Plots_uipanel,'BorderType','etchedin','HighlightColor',[1.0;1.0;1.0],'BorderWidth',1);
        set(handles.Results_uipanel,'BorderType','etchedin','HighlightColor',[1.0;1.0;1.0],'BorderWidth',1);
        set(handles.PixelCalc_uipanel,'BorderType','etchedin','HighlightColor',[1.0;1.0;1.0],'BorderWidth',1);
        set(handles.ViewParamImg_uipanel,'BorderType','line','HighlightColor',[0.043;0.518;0.78],'BorderWidth',3);
        
end


% --- Executes when user attempts to close QModeling_figure.
function QModeling_figure_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to QModeling_figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global QMf_logfile

deleteTempImages(handles);

aux_date=datestr(now);
fprintf(QMf_logfile,strcat(aux_date(end-7:end),' Finishing QModeling'));
fprintf(QMf_logfile,strcat('\n','----------------------------------------------','\n'));
fclose(QMf_logfile);

% Hint: delete(hObject) closes the figure
delete(hObject);

clear GLOBAL QMf_logfile QMmaindata QMlastPreprocess;


% --------------------------------------------------------------------
function ViewImage_MenuBar_Callback(hObject, eventdata, handles)
% hObject    handle to ViewImage_MenuBar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global QMf_logfile

file = spm_select(1,'image','Select an image (.img or .nii)');
file = file(1:end-2);

if ~isempty(file)
    try
    nii = load_nii(file,'','','','','',1);
    opt.setviewpoint = size(nii.img)/2;
    opt.setcolorindex = 4;
    opt.useinterp = 1;
    view_nii(nii,opt);                    
    catch ME
        errordlg(char({'An error ocurred: ';'Please, check logfile.txt'}),'Viewing error','modal');
        fprintf(QMf_logfile,strcat('   Error message: ',ME.message,'\n',...
                '   File: ',ME.stack(1).file,'\n','   Line: ',num2str(ME.stack(1).line),'\n'));
    end   
end


% --------------------------------------------------------------------
function SavePlot_AxesMenu_Callback(hObject, eventdata, handles)
% hObject    handle to SavePlot_AxesMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname, filterindex] = uiputfile( ...
       {'*.jpg','JPEG image (*.jpg)'; ...
       '*.fig','Figures (*.fig)'; ...
        '*.png','Portable Networks Graphics file (*.png)'},...
        'Save figure as');

if ~(isequal(filename,0) || isequal(pathname,0))

    % create a new figure for saving and printing
    fig = figure('Visible','off');
    % copy axes into the new figure
    newax = copyobj(handles.axes,fig);
    set(newax, 'units', 'normalized', 'position', [0.13 0.11 0.775 0.815]);

    saveas(fig,strcat(pathname,filename)) % save it
    close(fig) % clean up by closing it
    
end


% --------------------------------------------------------------------
function ViewValues_AxesMenu_Callback(hObject, eventdata, handles)
% hObject    handle to ViewValues_AxesMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

QM_viewValues('mainhandles',handles.QModeling_figure);




% --------------------------------------------------------------------
function AxesMenu_Callback(hObject, eventdata, handles)
% hObject    handle to AxesMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function WebResources_Callback(hObject, eventdata, handles)
% hObject    handle to WebResources (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function AboutQM_Callback(hObject, eventdata, handles)
% hObject    handle to AboutQM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
QM_about();

% --------------------------------------------------------------------
function WebCimes_Callback(hObject, eventdata, handles)
% hObject    handle to WebCimes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
web('http://www.cimes.es', '-browser');


% --------------------------------------------------------------------
function ContactUs_Callback(hObject, eventdata, handles)
% hObject    handle to ContactUs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

web('mailto:nroe@fguma.es')


% --------------------------------------------------------------------
function Configuration_Callback(hObject, eventdata, handles)
% hObject    handle to Configuration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function ConsistencyChecks_Callback(hObject, eventdata, handles)
% hObject    handle to ConsistencyChecks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global QMmaindata
global QMf_logfile

check = get(hObject,'Checked');
if strcmp(check,'on')
    QMmaindata.flags.consistency = 0;
    set(hObject,'Checked','off');
    aux_date=datestr(now);
    fprintf(QMf_logfile,strcat(aux_date(end-8:end),' CONFIGURATION: Check times consistency set OFF\n'));
else
    QMmaindata.flags.consistency = 1;
    set(hObject,'Checked','on');
    aux_date=datestr(now);
    fprintf(QMf_logfile,strcat(aux_date(end-8:end),' CONFIGURATION: Check times consistency set ON\n'));
end


% --- Executes on button press in CreateMask.
function CreateMask_Callback(hObject, eventdata, handles)
% hObject    handle to CreateMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

file = get(hObject,'UserData');

if ~isempty(file)
    run(file);
end


% --- Executes during object creation, after setting all properties.
function CreateMask_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CreateMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

set(hObject,'Enable','off');
set(hObject,'UserData','');

actualpath = mfilename('fullpath');
actualpath = actualpath(1:end-19);

strcat(actualpath,'wfu_pickatlas',filesep,'wfu_pickatlas.m');
if (isdir(strcat(actualpath,'wfu_pickatlas')) == 1 && exist(strcat(actualpath,'wfu_pickatlas',filesep,'wfu_pickatlas.m'),'file') == 2)
    addpath(strcat(actualpath,'wfu_pickatlas'));
    set(hObject,'Enable','on');
    set(hObject,'UserData',strcat(actualpath,'wfu_pickatlas',filesep,'wfu_pickatlas.m'));
end



function Refresh_Tacs(handles,data)

global QMmaindata

tavg = ((QMmaindata.times(:,1)+QMmaindata.times(:,2))/2);
tmin = tavg/60;

QMmaindata.idCt = data.idCt;
if isfield(data,'idCr')
    QMmaindata.idCr = data.idCr;
else %TwoTCM
    QMmaindata.idCr = data.idCplasma;
end

plot(handles.axes,tmin,QMmaindata.Crois{2,QMmaindata.idCt},'b.-');
hold on
plot(handles.axes,tmin,QMmaindata.Crois{2,QMmaindata.idCr},'r.-');
hold off
xlabel('Time (minutes)')
ylabel('Original units')
% aux = get(handles.mask_Textbox,'UserData');
% mask_name = aux{1}(find(aux{1} == filesep,1,'last')+1:end-2);
% title(['Time-activity curves for mask: ' mask_name])
title('Time-activity curves')

if QMmaindata.idCt==QMmaindata.idCr
    t_aux=table([1:1:length(tmin)]',tmin,QMmaindata.Crois{2,QMmaindata.idCt});
else
    t_aux=table([1:1:length(tmin)]',tmin,QMmaindata.Crois{2,QMmaindata.idCt}, QMmaindata.Crois{2,QMmaindata.idCr});
end
set(handles.axes,'UserData',{'TACs',table2cell(t_aux)});
% set(handles.axes,'UserData',{'TACs',[tmin QMmaindata.Crois{2,QMmaindata.idCt:QMmaindata.idCr}]});

set(handles.plottedCt_Popupmenu,'Enable','on');
set(handles.plottedCt_Popupmenu,'Value',QMmaindata.idCt);
set(handles.plottedCr_Popupmenu,'Enable','on');
set(handles.plottedCr_Popupmenu,'Value',QMmaindata.idCr);



function Refresh_Plots_Results(handles, modelID)

global QMmaindata
global QMlastPreprocess
tavg = ((QMmaindata.times(:,1)+QMmaindata.times(:,2))/2);
tmin = tavg/60;

%Set the popupmenus
num_rois = size(QMmaindata.Crois,2);
strmenu = cell(1,num_rois);

for i = 1:num_rois
    strmenu{1,i} = QMmaindata.Crois{1,i};
end
set(handles.plottedCt_Popupmenu,'String',strmenu);
set(handles.plottedCr_Popupmenu,'String',strmenu);

if modelID == 1
        %SRTM
        plot(handles.axes,tmin,QMlastPreprocess.SRTM.plots(:,1),'b.-',tmin,QMlastPreprocess.SRTM.plots(:,3),'r*-',tmin,QMlastPreprocess.SRTM.plots(:,2),'g+-');
        xlabel(handles.axes,'Time (minutes)')
        ylabel(handles.axes,'Original units')
        title(handles.axes,'SRTM fit')
        legend('TAC1 from receptor-rich region', 'TAC2 from reference region', 'Fit to TAC1');
%         set(handles.axes,'UserData',{'SRTM*',[tmin QMmaindata.Crois{2,QMmaindata.idCt} QMlastPreprocess.SRTM.plots(:,2) QMmaindata.Crois{2,QMmaindata.idCr}]});        
        
        t_aux=table([1:1:length(QMmaindata.Crois{2,QMmaindata.idCt})]',tmin,QMmaindata.Crois{2,QMmaindata.idCt},QMlastPreprocess.SRTM.plots(:,2),QMmaindata.Crois{2,QMmaindata.idCr});
        set(handles.axes,'UserData',{'SRTM*',table2cell(t_aux)});
        
        %Showing results in text area
        strresults1 = cell(1,6);
        strresults1{1,1} = 'SRTM preprocessing:';
        strresults1{1,2} = '-----------------------------------';
        strresults1{1,3}  = ['R1:      ',num2str(QMlastPreprocess.SRTM.results{1},'%.5f'), '     (SE: ',num2str(sqrt(QMlastPreprocess.SRTM.results{8}(1)),'%1.2e'),')'];
        strresults1{1,4}  = ['k2:      ',num2str(QMlastPreprocess.SRTM.results{2},'%.5f'), '     (SE: ',num2str(sqrt(QMlastPreprocess.SRTM.results{8}(2)),'%1.2e'),')'];
        strresults1{1,5}  = ['BPnd: ',num2str(QMlastPreprocess.SRTM.results{3},'%.5f'), '     (SE: ',num2str(sqrt(QMlastPreprocess.SRTM.results{8}(3)),'%1.2e'),')'];
        strresults1{1,6}  = ['k2a:    ',num2str(QMlastPreprocess.SRTM.results{4},'%.5f')];
        strresults2 = cell(1,4);
        strresults2{1,1} = 'Goodness of fit:';
        strresults2{1,2} = '-----------------------------------';
        strresults2{1,3} = ['NMSE:           ',num2str(QMlastPreprocess.SRTM.results{5},'%.5f')];
        strresults2{1,4} = ['Corr. Coef.:  ',num2str(QMlastPreprocess.SRTM.results{6},'%.5f')];    

        set(handles.results_Text,'String',strresults1);
        set(handles.results_Text2,'String',strresults2);
        
        set(handles.plottedCt_Popupmenu,'Value',QMlastPreprocess.SRTM.idCt);
        set(handles.plottedCr_Popupmenu,'Value',QMlastPreprocess.SRTM.idCr);

        set(handles.corrMatrix_Pushbutton,'Enable','on')
        
elseif modelID == 2
        %SRTM2
        plot(handles.axes,tmin,QMmaindata.Crois{2,QMlastPreprocess.SRTM2.idCt},'b.-',tmin,QMmaindata.Crois{2,QMlastPreprocess.SRTM2.idCr},'r*-',tmin,QMlastPreprocess.SRTM2.plots,'g+-');
        xlabel(handles.axes,'Time (minutes)')
        ylabel(handles.axes,'Original units')
        title(handles.axes,'SRTM2 fit')
        legend('TAC1 from receptor-rich region', 'TAC2 from reference region', 'Fit to TAC1');
%         set(handles.axes,'UserData',{'SRTM*',[tmin QMmaindata.Crois{2,QMmaindata.idCt} QMlastPreprocess.SRTM2.plots' QMmaindata.Crois{2,QMmaindata.idCr}]});
        
        t_aux=table([1:1:length(QMmaindata.Crois{2,QMmaindata.idCt})]',tmin,QMmaindata.Crois{2,QMmaindata.idCt},QMlastPreprocess.SRTM2.plots(:,2),QMmaindata.Crois{2,QMmaindata.idCr});
        set(handles.axes,'UserData',{'SRTM*',table2cell(t_aux)});
        
        %Showing results in text area
        strresults = cell(1,6);
        strresults{1,1} = 'Carson preprocessing:';
        strresults{1,2} = '-----------------------------------';
        if ~QMlastPreprocess.SRTM2.k2_pChecked
            strresults{1,3}  = [strcat('k2',char(39),':'),'     ',num2str(QMlastPreprocess.SRTM2.results{1},'%.5f')];
        else
            strresults{1,3}  = [strcat('k2',char(39),':'),'     ',num2str(QMlastPreprocess.SRTM2.results{1},'%.5f'), '     (SE: ',num2str(sqrt(QMlastPreprocess.SRTM2.results{8}(1)),'%1.2e'),')'];
        end
        strresults{1,4}  = ['R1:      ',num2str(QMlastPreprocess.SRTM2.results{2},'%.5f'), '     (SE: ',num2str(sqrt(QMlastPreprocess.SRTM2.results{8}(2)),'%1.2e'),')'];
        strresults{1,5}  = ['k2a:    ',num2str(QMlastPreprocess.SRTM2.results{3},'%.5f'), '     (SE: ',num2str(sqrt(QMlastPreprocess.SRTM2.results{8}(3)),'%1.2e'),')'];
        strresults{1,6}  = ['BPnd: ',num2str(QMlastPreprocess.SRTM2.results{4},'%.5f')];
        strresults2 = cell(1,4);
        strresults2{1,1} = 'Goodness of fit:';
        strresults2{1,2} = '-----------------------------------';
        strresults2{1,3} = ['NMSE:           ',num2str(QMlastPreprocess.SRTM2.results{5},'%.5f')];
        strresults2{1,4} = ['Corr. Coef.:  ',num2str(QMlastPreprocess.SRTM2.results{6},'%.5f')];

        set(handles.results_Text,'String',strresults);
        set(handles.results_Text2,'String',strresults2);
        
        set(handles.plottedCt_Popupmenu,'Value',QMlastPreprocess.SRTM2.idCt);
        set(handles.plottedCr_Popupmenu,'Value',QMlastPreprocess.SRTM2.idCr);

        set(handles.corrMatrix_Pushbutton,'Enable','on')
elseif modelID == 3
        %PatlakRef
        plot(handles.axes,QMlastPreprocess.PatlakRef.plots(:,1),QMlastPreprocess.PatlakRef.plots(:,2),'ro',QMlastPreprocess.PatlakRef.plots(:,1),QMlastPreprocess.PatlakRef.plots(:,3),'b');
        xlabel(handles.axes,'int(Cr)/Cr')
        ylabel(handles.axes,'Ct/Cr')
        title(handles.axes,'Patlak fit')
%         set(handles.axes,'UserData',{'Patlak',[QMlastPreprocess.PatlakRef.plots(:,1) QMlastPreprocess.PatlakRef.plots(:,2) QMlastPreprocess.PatlakRef.plots(:,3)]});
        
        t_aux=table([1:1:length(QMlastPreprocess.PatlakRef.plots(:,1))]',QMlastPreprocess.PatlakRef.plots(:,1),QMlastPreprocess.PatlakRef.plots(:,2),QMlastPreprocess.PatlakRef.plots(:,3));
        set(handles.axes,'UserData',{'Patlak',table2cell(t_aux)});
        
        %Showing results in text area
        strresults = cell(1,11);
        strresults{1,1} = 'Patlak preprocessing:';
        strresults{1,2} = '-----------------------------------';
        strresults{1,3}  = ['t*:                ',num2str(QMlastPreprocess.PatlakRef.results{1})];
        strresults{1,4}  = ['K:                 ',num2str(QMlastPreprocess.PatlakRef.results{2},'%.5f'), '     (Var: ',num2str(QMlastPreprocess.PatlakRef.results{8}(1),'%1.2e'),')'];
        strresults{1,5}  = ['intercept:   ',num2str(QMlastPreprocess.PatlakRef.results{3},'%.5f'), '     (Var: ',num2str(QMlastPreprocess.PatlakRef.results{8}(2),'%1.2e'),')'];
        strresults{1,6}  = ['Start:          ',num2str(QMlastPreprocess.PatlakRef.results{4},'%.5f')];
        strresults{1,7}  = ['Max. Error: ',num2str(QMlastPreprocess.PatlakRef.results{5},'%.5f')];
        strresults2 = cell(1,4);
        strresults2{1,1} = 'Goodness of fit:';
        strresults2{1,2} = '-----------------------------------';
        strresults2{1,3} = ['NMSE:           ',num2str(QMlastPreprocess.PatlakRef.results{6},'%.5f')];
        strresults2{1,4} = ['Corr. Coef.:  ',num2str(QMlastPreprocess.PatlakRef.results{7},'%.5f')];

        set(handles.results_Text,'String',strresults);
        set(handles.results_Text2,'String',strresults2);
        
        set(handles.plottedCt_Popupmenu,'Value',QMlastPreprocess.PatlakRef.idCt);
        set(handles.plottedCr_Popupmenu,'Value',QMlastPreprocess.PatlakRef.idCr);
        set(handles.corrMatrix_Pushbutton,'Enable','off');
elseif modelID == 4
	%LoganPlot
        plot(handles.axes,QMlastPreprocess.LoganPlot.plots(:,1),QMlastPreprocess.LoganPlot.plots(:,2),'ro',QMlastPreprocess.LoganPlot.plots(:,1),QMlastPreprocess.LoganPlot.plots(:,3),'b');
        %set(handles.axes,'UserData',{'LoganPlot',[QMlastPreprocess.LoganPlot.plots(:,1) QMlastPreprocess.LoganPlot.plots(:,2)]});
        xlabel(handles.axes,strcat('[int(Ct)+Ct/k2',char(39),']/Ct'))
        ylabel(handles.axes,'int(Ct)/Ct')
        title(handles.axes,'Logan Fit')
%         set(handles.axes,'UserData',{'LoganPlot',[QMlastPreprocess.LoganPlot.plots(:,1) QMlastPreprocess.LoganPlot.plots(:,2) QMlastPreprocess.LoganPlot.plots(:,3)]});

        t_aux=table([1:1:length(QMlastPreprocess.LoganPlot.plots(:,1))]',QMlastPreprocess.LoganPlot.plots(:,1),QMlastPreprocess.LoganPlot.plots(:,2),QMlastPreprocess.LoganPlot.plots(:,3));
        set(handles.axes,'UserData',{'LoganPlot',table2cell(t_aux)});
        
        %Showing results in text area
        strresults = cell(1,6);
        strresults{1,1} = 'Logan Plot preprocessing:';
        strresults{1,2} = '-----------------------------------';
        strresults{1,3} = ['BPnd:        ',num2str(QMlastPreprocess.LoganPlot.results{1},'%.5f'), '             (Var: ',num2str(QMlastPreprocess.LoganPlot.results{6}(1),'%1.2e'),')'];
        strresults{1,4} = ['intercept: ',num2str(QMlastPreprocess.LoganPlot.results{2},'%.5f'), '     (Var: ',num2str(QMlastPreprocess.LoganPlot.results{6}(2),'%1.2e'),')'];
        strresults{1,5} = ['k2',char(39),':            ' , num2str(QMlastPreprocess.LoganPlot.k2)];
        strresults{1,6} = ['t*:              ',num2str(QMlastPreprocess.LoganPlot.results{3})];
        strresults2 = cell(1,4);
        strresults2{1,1} = 'Goodness of Fit';
        strresults2{1,2} = '-----------------------------------';
        strresults2{1,3} = ['Corr. Coef.:  ', num2str(QMlastPreprocess.LoganPlot.results{4},'%.5f')];
        strresults2{1,4} = ['NMSE:           ', num2str(QMlastPreprocess.LoganPlot.results{5},'%.5f')];

        set(handles.results_Text,'String',strresults);
        set(handles.results_Text2,'String',strresults2);
            
        set(handles.plottedCt_Popupmenu,'Value',QMlastPreprocess.LoganPlot.idCt);
        set(handles.plottedCr_Popupmenu,'Value',QMlastPreprocess.LoganPlot.idCr);
        set(handles.corrMatrix_Pushbutton,'Enable','off');
elseif modelID == 5
  %TwoTCM
        if isempty(QMlastPreprocess.TwoTCM.Cblood)
            plot(handles.axes,tmin,QMmaindata.Crois{2,QMlastPreprocess.TwoTCM.idCt},'b.-',tmin,QMmaindata.Crois{2,QMlastPreprocess.TwoTCM.idCplasma},'r*-',tmin,QMlastPreprocess.TwoTCM.plots,'g+-');
            xlabel(handles.axes,'Time (minutes)')
            ylabel(handles.axes,'Original units')
            title(handles.axes,'2-Tissue Compartment Model fit')
            legend('TAC1 from receptor-rich region', 'TAC2 from plasma region', 'Fit to TAC1');
            %Save values
%             set(handles.axes,'UserData',{'TwoTCM*',[tmin QMmaindata.Crois{2,QMmaindata.idCt} QMlastPreprocess.TwoTCM.plots QMmaindata.Crois{2,QMmaindata.idCr}]});
            t_aux=table([1:1:length(QMlastPreprocess.TwoTCM.plots)]',tmin,QMmaindata.Crois{2,QMlastPreprocess.TwoTCM.idCt},QMlastPreprocess.TwoTCM.plots,QMmaindata.Crois{2,QMlastPreprocess.TwoTCM.idCplasma});
            set(handles.axes,'UserData',{'TwoTCM*',table2cell(t_aux)});
        else %optionally, a whole blood activity TAC has been included
            plot(handles.axes,tmin,QMmaindata.Crois{2,QMlastPreprocess.TwoTCM.idCt},'b.-',tmin,QMmaindata.Crois{2,QMlastPreprocess.TwoTCM.idCplasma},'r*-',tmin,QMlastPreprocess.TwoTCM.plots,'g+-',tmin,QMlastPreprocess.TwoTCM.Cblood,'y--');
            xlabel(handles.axes,'Time (minutes)')
            ylabel(handles.axes,'Original units')
            title(handles.axes,'2-Tissue Compartment Model fit')
            legend(mainhandles.axes,'TAC1 from receptor-rich region', 'TAC2 from plasma region', 'Fit to TAC1','Whole blood activity');
            %Save values
%             set(handles.axes,'UserData',{'TwoTCM',[tmin QMmaindata.Crois{2,QMmaindata.idCt} QMlastPreprocess.TwoTCM.plots QMmaindata.Crois{2,QMmaindata.idCr} QMlastPreprocess.TwoTCM.Cblood]});  
            t_aux=table([1:1:length(QMlastPreprocess.TwoTCM.plots)]',tmin,QMmaindata.Crois{2,QMlastPreprocess.TwoTCM.idCt},plots,QMmaindata.Crois{2,QMlastPreprocess.TwoTCM.idCplasma},QMlastPreprocess.TwoTCM.Cblood);
            set(handles.axes,'UserData',{'TwoTCM',table2cell(t_aux)});
        end
        
        
        %Showing results in text area   
        strresults = cell(1,6);
        strresults{1,1} = '2TCM preprocessing:';
        strresults{1,2} = '-----------------------------------';
        if ~(QMlastPreprocess.TwoTCM.vB(1)) %~vB_checked QMlastPreprocess.TwoTCM.vB=[vB_checked vB];
            strresults{1,3}  = ['vB: ',num2str(QMlastPreprocess.TwoTCM.results{1},'%.5f'), '     (fixed)'];
        else
            strresults{1,3}  = ['vB: ',num2str(QMlastPreprocess.TwoTCM.results{1},'%.5f'), '     (SE: ',num2str(sqrt(QMlastPreprocess.TwoTCM.results{9}(1)),'%1.2e'),')'];
        end        
        strresults{1,4}  = ['K1: ',num2str(QMlastPreprocess.TwoTCM.results{2},'%.5f'), '     (SE: ',num2str(sqrt(QMlastPreprocess.TwoTCM.results{9}(2)),'%1.2e'),')'];
        strresults{1,5}  = ['k2: ',num2str(QMlastPreprocess.TwoTCM.results{3},'%.5f'), '     (SE: ',num2str(sqrt(QMlastPreprocess.TwoTCM.results{9}(3)),'%1.2e'),')'];
        strresults{1,6}  = ['k3: ',num2str(QMlastPreprocess.TwoTCM.results{4},'%.5f'), '     (SE: ',num2str(sqrt(QMlastPreprocess.TwoTCM.results{9}(4)),'%1.2e'),')'];
        if ~QMlastPreprocess.TwoTCM.k4_checked
            strresults{1,7}  = ['k4: ',num2str(QMlastPreprocess.TwoTCM.results{5},'%.5f'), '     (fixed)'];
        else
            strresults{1,7}  = ['k4: ',num2str(QMlastPreprocess.TwoTCM.results{5},'%.5f'), '     (SE: ',num2str(sqrt(QMlastPreprocess.TwoTCM.results{9}(5)),'%1.2e'),')'];
        end
%         strresults{1,7}  = ['Ki: ',num2str(QMlastPreprocess.TwoTCM.results{5},'%.5f')];
        
        strresults2 = cell(1,4);
        strresults2{1,1} = 'Goodness of fit:';
        strresults2{1,2} = '-----------------------------------';
        strresults2{1,3} = ['NMSE:           ',num2str(QMlastPreprocess.TwoTCM.results{6},'%.5f')];
        strresults2{1,4} = ['Corr. Coef.:  ',num2str(QMlastPreprocess.TwoTCM.results{6},'%.5f')];

        set(handles.results_Text,'String',strresults);
        set(handles.results_Text2,'String',strresults2);
        
        set(handles.plottedCt_Popupmenu,'Value',QMlastPreprocess.TwoTCM.idCt);
        set(handles.plottedCr_Popupmenu,'Value',QMlastPreprocess.TwoTCM.idCplasma);

        set(handles.corrMatrix_Pushbutton,'Enable','on')
end

set(handles.saveResults_Pushbutton,'Enable','on');
set(handles.plottedCt_Popupmenu,'Enable','off');
set(handles.plottedCr_Popupmenu,'Enable','off');


   
function checkTimesConsistency()

global QMmaindata

if QMmaindata.flags.consistency == 1
    i = 1;
    error = 0;
    while ((i < QMmaindata.num_files) && (error == 0))    
        if QMmaindata.times(i,2) ~= QMmaindata.times(i+1,1)
            error = 1;
        end
        i = i+1;
    end
    if error == 1
        uiwait(warndlg({'A consistency check for the PET times has been executed and errors';...
                       'were found with respect to the end of a frame and the start of the next.';...
                       'If you dont want to get inconsistent results in the preprocessing and';...
                       'the parametric images, please modify the PET times.'},'Warning','modal'));     
    end
end


% --- Executes on button press in saveFigure_IconButton.
function saveFigure_IconButton_Callback(hObject, eventdata, handles)
% hObject    handle to saveFigure_IconButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname, filterindex] = uiputfile( ...
       {'*.jpg','JPEG image (*.jpg)'; ...
       '*.fig','Figures (*.fig)'; ...
        '*.png','Portable Networks Graphics file (*.png)'},...
        'Save figure as');

if ~(isequal(filename,0) || isequal(pathname,0))

    % create a new figure for saving and printing
    fig = figure('Visible','off');
    % copy axes into the new figure
    newax = copyobj(handles.axes,fig);
    set(newax, 'units', 'normalized', 'position', [0.13 0.11 0.775 0.815]);

    saveas(fig,strcat(pathname,filename)) % save it
    close(fig) % clean up by closing it
    
end


% --- Executes during object creation, after setting all properties.
function saveFigure_IconButton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to saveFigure_IconButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% global QMpath
% actualpath = mfilename('fullpath');
% QMpath = actualpath(1:end-10);
% 
% [a,map]=imread(strcat(QMpath,filesep,'icons',filesep,'file_save.png'));
% 
% % Determina el tama�o de la imagen en las tres dimensiones
% [r,c,d]=size(a); 
% % redondea el tama�o de la matriz r / 38
% x=ceil(r/38);
% % redondea el tama�o de la matriz c / 51
% y=ceil(c/51);
% % crea una matriz g, con elementos de la matriz "a", pero considerando intervalos x y y
% g=im2double(a(1:x:end,1:y:end,:));
% %g(g==0)=0.973; %Para poner el borde de la imagen al mismo color del fondo que tengo puesto
% set(hObject,'CData',g);


% --- Executes on button press in viewValues_IconButton.
function viewValues_IconButton_Callback(hObject, eventdata, handles)
% hObject    handle to viewValues_IconButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

QM_viewValues('mainhandles',handles.QModeling_figure);


% --- Executes during object creation, after setting all properties.
function viewValues_IconButton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to viewValues_IconButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% global QMpath
% 
% actualpath = mfilename('fullpath');
% QMpath = actualpath(1:end-10);
% 
% [a,map]=imread(strcat(QMpath,filesep,'icons',filesep,'notesicon.gif'));
% g = ind2rgb(a,map);
% 
% %[a,map]=imread(strcat(QMpath,filesep,'icons',filesep,'tool_data_cursor.gif'));
% % % Determina el tama�o de la imagen en las tres dimensiones
% % [r,c,d]=size(a); 
% % % redondea el tama�o de la matriz r / 38
% % x=ceil(r/38);
% % % redondea el tama�o de la matriz c / 51
% % y=ceil(c/51);
% % % crea una matriz g, con elementos de la matriz "a", pero considerando intervalos x y y
% % g=im2double(a(1:x:end,1:y:end,:));
% %g=a(1:x:end,1:y:end,:);
% %g(g==0.062)=0.973 %Para poner el borde de la imagen al mismo color del fondo que tengo puesto
% set(hObject,'CData',g);


% --------------------------------------------------------------------
function CopyResults_Callback(hObject, eventdata, handles)
% hObject    handle to CopyResults (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

copy(get(handles.results_Text,'String'));


% --------------------------------------------------------------------
function ResultsMenu_Callback(hObject, eventdata, handles)
% hObject    handle to ResultsMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in dataCursor_ToggleButton.
function dataCursor_ToggleButton_Callback(hObject, eventdata, handles)
% hObject    handle to dataCursor_ToggleButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dataCursor_ToggleButton

if get(hObject,'Value') == 1
    datacursormode on
else
    datacursormode off
end

% --- Executes during object creation, after setting all properties.
function dataCursor_ToggleButton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dataCursor_ToggleButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% global QMpath
% actualpath = mfilename('fullpath');
% QMpath = actualpath(1:end-10);
% 
% [a,map]=imread(strcat(QMpath,filesep,'icons',filesep,'tool_data_cursor.gif'));
% g = ind2rgb(a,map);
% set(hObject,'CData',g);


% --------------------------------------------------------------------
function Minimize_preprocessing_window_Callback(hObject, eventdata, handles)
% hObject    handle to Minimize_preprocessing_window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global QMmaindata
global QMf_logfile

check = get(hObject,'Checked');
if strcmp(check,'on')
    QMmaindata.flags.minWindow = 0;
    set(hObject,'Checked','off');
    aux_date=datestr(now);
    fprintf(QMf_logfile,strcat(aux_date(end-8:end),' CONFIGURATION: Minimize preprocessing window set OFF\n'));
else
    QMmaindata.flags.minWindow = 1;
    set(hObject,'Checked','on');
    aux_date=datestr(now);
    fprintf(QMf_logfile,strcat(aux_date(end-8:end),' CONFIGURATION: Minimize preprocessing window set ON\n'));
end

% --------------------------------------------------------------------
function WebUIMCIMES_Callback(hObject, eventdata, handles)
% hObject    handle to WebUIMCIMES (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
web('http://www.uimcimes.es', '-browser');


% --- Executes on button press in corrMatrix_Pushbutton.
function corrMatrix_Pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to corrMatrix_Pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global QMlastPreprocess
global QMf_logfile


%Set panel focus
setPanelFocus(handles,'PlotResultsPanels');

selected = get(handles.selectModel_Popupmenu,'Value');
if selected == 1, model = 'SRTM'; end
if selected == 2, model = 'SRTM2'; end
if selected == 3, model = 'PatlakRef'; end
if selected == 4, model = 'LoganPlot'; end
if selected == 5, model = 'TwoTCM';end

data = getfield(QMlastPreprocess,model);

if isempty(data.plots)
   errordlg(char({'An error ocurred during the pixel-wise calculation:';'You have not a preprocess with the model selected.';'View help.'}),'Pixel calculation','modal');
   fprintf(QMf_logfile,strcat('  ERROR: A preprocess with the model selected does not exist','\n'));
else
    
    switch selected
        case 1
            c_rName={'R1','k2','BPnd'};
            id_data=7;
        case 2
            c_rName={'R1',strcat('k2',char(39)),'k2a'};
            id_data=7;
        case 3 %the button must not be enable for this model            
            errordlg(char({'Correlation matrix is obtained only for compartmental models'}),'Correlation matrix','modal');
            fprintf(QMf_logfile,strcat('  ERROR: Correlation matrix is obtained only for compartmental models','\n'));
            
        case 4 % the button must not be enable for this model           
            errordlg(char({'Correlation matrix is obtained only for compartmental models'}),'Correlation matrix','modal');
            fprintf(QMf_logfile,strcat('  ERROR: Correlation matrix is obtained only for compartmental models','\n'));
            
        case 5
            c_rName={'vB','K1','k2','k3','k4'};
            id_data=8;
    end

    fig=figure('Menubar','none','Name','Correlation Matrix','NumberTitle','off','ToolBar','figure','Visible','off');
  
    t=uitable(fig,'Data',data.results{id_data},'RowName',c_rName,'ColumnName',c_rName);
    
    axis off
    title(model)
    t.Position(3:4)=t.Extent(3:4); %Set Width and Height to Accommodate Data Size
    fig.Position(3)=1.2*t.Extent(3); % figure width
    fig.Position(4)=1.5*t.Extent(4); % figure height
        
    set(fig,'Visible','on');
end

% --- Executes on button press in saveResults_Pushbutton.
function saveResults_Pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to saveResults_Pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% global QMlastPreprocess

[filename, pathname] = uiputfile( ...
           {'*.txt','Text file (*.txt)';},...
            'Save values as');
if ~(isequal(filename,0) || isequal(pathname,0))
    f = fopen(strcat(pathname,filename), 'wt');

    restext=get(handles.results_Text,'String');
    for i=1:size(restext,1)
        fprintf(f,'%s\n',restext{i,1});
    end
    fprintf(f,'\n');
    restext=get(handles.results_Text2,'String');
    for i=1:size(restext,1)
        fprintf(f,'%s\n',restext{i,1});
    end
    fclose(f);
end

