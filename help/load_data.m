%% Load Data
% 
%% How can I load a study? In which format? 
% To load a dynamic PET study, press the _Load Study_ button. A new auxiliary window
% will open, to allow searching and loading the study from your folder structure. 
% The supported formats are _Analyze (hdr/img)_ or _NIfTI_(nii). Single-file
% (multiframe) or multiple-file (one file for each frame of the PET
% study) studies can be loaded. Once the file/s are chosen, press 
% the _Done_ button and the PET study will be loaded.
%
%% What happens if my study does not have start and end frame times?
% Nothing. A warning window will appear, and default start and end times will
% be used, unless you specify values for your scan. See next question on
% how to specify these values.
%
%% Can I change or save the frame times of the loaded study? How can I do this?
% Yes. You can do this by clicking on the _Change Times_ button. A new auxiliary window will open. Then, 
% you can change or save current times, load new times from a file or check
% the consistency of the time values.
%
%% What is the format of the times file?
% You can load new frame times from a file. This file must have the following format:
% 
% * One row for each _frame_ of the study.
% * One column for the starting times and one for the ending times
% * Units must be seconds
%
% In the first line, the starting and ending times of the first _frame_ of the
% study, in the second line, the starting and ending times of the second _frame_
% of the study, and so on.
%
%% Can I check the times consistency? How can I do this?
% Yes, there are two ways:
%
% * By the _Change Times_ button on the main window.
% * By the toolbar of QModeling: Click _QModeling->Configuration_ and check
% the _Consistency check_ option. With this, the times consistency will be
% automatically checked. A warning will appear if they are not consistent.
%
%% How can I load a mask? What are the accepted formats?
% To load a mask with regions of interest (ROIs) 
% press the _Load Mask/TACs_ button.  It must include only one ROI for the
% reference region, and an indefinite number of ROIs representing
% specific-uptake areas for the preprocessing step. Only one of these shall
% be used for the preprocessing.
% When the button is pressed, a new auxiliary window opens, 
% which allows searching and loading the mask from the 
% folder structure. The allowed formats are _Analyze (hdr/img)_ or _NIfTI
% (nii)_. All the ROIs must be included in one single file. A numeric index
% (an integer) specifies each ROI. Clicking the _Done_ button, the mask
% will be loaded. 
% The name of the ROI associated with each numeric index can be specified
% in a plain text file. If such a file exists at the same path and with the 
% same name as the mask, it will be automatically loaded*. If not, 
% a new auxiliary window will appear asking for that text file*. 
% If no text file is chosen or detected, the TACs associated to each mask 
% will receive  generic names (_TAC1_, _TAC2_, ...) with numeric indexes following a correlative order with the indexes in the mask.
%
%
% *Specific format is required. Please, go to the "What is the format of the mask names file?" question.
%% What is the format of the mask names file?
% It must be a plain text file including the different regions on different
% rows. The name of each region must be at the first column and the corresponding integer
% at a second column. In order to be automatically loaded by the
% application, this text file and the mask file must be at the same path and have the same name.
% If they are in separate folders, the program will automatically open a window asking for the mask names
% file.
%
%% Can I directly load the time-activity curves?
% Yes. Instead of using a mask file, you can directly load the time-activity
% curves of the region(s) of interest and the reference region. The accepted
% formats are .txt, .xls, .xlsx and .mat. If a .mat file is used, it must
% contain a variable called _TACs_names_ with the names of the regions and
% another one called _TACs_ with the values of each region as columns. 
% If the file is in .txt, .xls, or .xlsx format, there must be a first row with the names of
% the regions.  The following rows will contain the values of the TACs for each
% frame. Each column will represent the TAC of a specific region.
% In the .txt case, you have to split the columns by the tabulator key.
% At /QModeling/help/template_TACs you have a template of TACs file for
% each format option.

%% Can I generate a mask for a loaded PET study? How can I do this?
% Yes. The _WFU_PickAtlas toolbox_ for SPM must be downloaded
% (http://www.nitrc.org/projects/wfu_pickatlas/ or
% http://fmri.wfubmc.edu/software/PickAtlas) and installed to enable _Create
% mask with pickatlas_ button in the main window.
% 
%% How can I see the TACs?
% Once you have succesfully loaded the PET study and the mask file, you can
% obtain the TACs by clicking on the _Prepare TACs_ button. 
% If the time-activity curves have been directly loaded, the same button has to be pressed. 
% If everything  is correct, the TACs will be displayed 
% in the _Plots_ panel of the main window; otherwise, an error message will
% be shown.
%
%% Can I see and save the numerical values of the TACs? Can I save the plotted TACs?
% Yes. Right-button click on the _Plots_ panel.
%
%% Can I change the reference and interest region?
% Yes, you can choose whatever you want on the appropriate drop-down list.
