%% Model selection & Preprocessing
%% How do I obtain the kinetic parameters and the fitted TAC for a region of interest? 
% A PET dynamic study must be loaded (see the Load Data help section), 
% together with a mask or a set of TACs, to begin the kinetic analysis process. 
% After that, the Model selection and Preprocessing frame changes into blue, 
% and the kinetic model can be chosen. Then, you have to press the
% _Preprocess Model_ button, and a new configuration window appears,
% where the preprocessing parameters for the selected model can be specified.
% In this window there is also a new _Preprocess Model_ button. The calculation 
% starts when you click on it. The fitted and input TACs are shown in the _Plots_ frame 
% for comparison, and the kinetic parameters are presented in the _Results_ frame.
%
%% Which input parameters do I need to specify to perform the fitting?
% The input parameters are selected in the configuration window, which 
% opens when the user clicks on the _Preprocess Model_ button. Each kinetic model requires 
% different input parameters to be specified. You can see the information 
% for each model by clicking on the help icon or in the following links:
% <matlab:open('help_SRTM.html') SRTM model>, <matlab:open('help_Carson.html') Carson model>,
% <matlab:open('help_patlak.html') PatlakRef model> or <matlab:open('help_LoganPlot.html') LoganPlot model>.  
%
%% Can I perform more than one fit, to try different input parameters?
% Yes. After the fitting curve and the estimated parameters are shown,
% you can restore the minimized configuration window to choose another
% parameter combination. When the fit is considered correct, the configuration
% window can be closed py pressing the _OK_ button. The parameters used / obtained 
% for that particular configuration will be used in the subsequent processing steps. 
% If a new preprocessing fit is later needed, it can be performed by 
% clicking on _Preprocess Model_ in the main frame again.
%
%% Can I see and/or save the values of TACs?
% Yes, you just need to right-click on any point of the _Plots_ frame. 
% Two options are available: saving the plots _(Save plots)_ or showing 
% the numerical values of the curves _(View values)_. If _Save plots_ is chosen,
% the software will ask for a file name, a format and a save directory. 
% If _View values_ is selected, a new window will open, showing the values. 
% The values in this window can be saved by clicking on the _Save values_ button.
%
%% Can I see the measured TAC of any region again after performing a new fit?
% Yes. Two radio buttons will be enabled in the Plots frame after any fit. They allow
% switching between the display of the fitted curve with its associated TACs, 
% and the display of any pair of measured TACs _(Show TAC's)_,
% specified as region of interest and reference region.
%
%% Can I restore a preprocessing perfomed for a given model at any time?
% Yes, but only for the last preprocessing. If you have done a
% preprocessing for a given model, you just need to select that model in the _Model
% selection & Preprocessing_ frame and the results of last preprocessing will
% be loaded.