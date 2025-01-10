%% User manual
% 
% This manual explains the capabilities of QModeling, and the way to use
% it. It is divided in four main sections: _Graphical User Interface (GUI)_, 
% _Implemented models_, _Command line functions_ and _Include your own model_.
% 
% The _Graphical User Interface (GUI)_ section, which is organized in a questions
% and answers format, has four different subsections according to the four main 
% steps to be followed when using the toolbox. Each of them is carried out 
% by a specific submodule: Data loading and TAC generation; 
% model selection and preprocessing; generation of parametric
% images and display of parametric images. Thus, QModeling allows 
% generating parametric images for a particular dynamic PET study, by 
% fitting one of the implemented models (SRTM, SRTM2, Patlak or Logan) to the data.
% Brief information and references about the four compartmental models
% are included at the _Implemented models_ section.
% 
% You can also fit each of the four compartmental models to multiple studies 
% using the _Command line functions_. Three command line functions (to load studies, 
% preprocess them and generate whole brain images) are implemented for 
% that purpose. 
%
% At the _Include your own model_ section, it is described how to include
% in the QModeling software the implementation of a new kinetic model.
%% Graphical User Interface (GUI)
% The main window is divided into six different parts or subframes. The four 
% of them placed at the right and the bottom of the window represent the 
% four submodules of the toolbox. The other two show the graphical  
% and numerical results, at the right  of the main window. 
% When the program is started, all the subframes are inactive except for
% the first one. After introducing the relevant input data, the following
% frames will be activated.
%
% * <matlab:open('load_data.html') *Load Data*> 
%
% * <matlab:open('modeling.html') *Model selection & Preprocessing*> 
%
% * <matlab:open('pixel_calculation.html') *Pixel-wise Calculation*>
%
% * <matlab:open('view_images.html') *View Parametric Images*>
%
% In this <http://www.mathworks.com _flowchart_> you can see some pictures of the different steps.
%% Implemented models
% * <matlab:open('help_SRTM.html') *SRTM*>
% * <matlab:open('help_Carson.html') *SRTM2*>
% * <matlab:open('help_patlak.html') *Patlak Reference*>
% * <matlab:open('help_LoganPlot.html') *Logan Plot*>
%% Command line functions
% In order to fit one of the four kinetic models included in _QModeling_
% to multiple studies, three main command line functions have been
% implemented:
%
% *_QModeling_generateTACs.m_*: To load multiple PET studies (each one 
% with a reference and a target mask) and to obtain the TACs from all of
% them.
%
% *_QModeling_preprocessTACs.m_*: From the reference and target TACs of
% multiple studies, this function allows the user to fit one of the four kinetic
% models. 
%
% *_QModeling_parametric_images.m_*: To obtain the parametric images from 
% the fitting of one of the kinetic models implemented to multiple PET 
% studies (pixel-wise calculation).
% 
% To obtain more information about one of these functions, type _help_ and 
% the name of the function at the Matlab command window.
% Furthermore, a _demo script_ is included at the command line folder of the toolbox
% to show the user how the command line could be run. It is called _demoQModeling.m_
% and typing _help demoQModeling.m_ at the command window, more information
% will be displayed. Also, there are many comments on the script to help
% the user to understand all the steps.
%% Include your own model
% The aim of this section is to provide users with the necessary
% information to modify the QModeling software, in order to include in it
% additional kinetic models. These instructions only refer to the
% implementation of this new function via command line. Users are welcome
% to contact the authors for any further advice, of if they are interested
% in including a model in the GUI, after developing it
% for the command line. This information has not been included in this
% manual for the moment, as adding a new model to the GUI requires
% modifying some more routines in the package.
%   
% A new Matlab function has to be developed, under the name
% _QModeling_newModel_Name_. Here, _newModel_Name_ is the name that you
% want to give to the model. The three already implemented functions
% (_QModeling_SRTM.m_, _QModeling_SRTM2.m_ and _QModeling_PatlakRef.m_) can be
% taken as examples. The new function should be run from
% _QModeling_preprocessTACs.m_ and from _QModeling_parametric_images.m_. With
% the first one, the fit will be done at a ROI level. The second one will
% yield parametric images, after a voxel level fit.
%
% * <matlab:open('QM_preprocTACs.html') *QModeling_preprocessTACs.m*>
%
% * <matlab:open('QM_parametImages.html') *QModeling_parametric_images.m*>


