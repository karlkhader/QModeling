%% _QModeling_parametric_images.m_
% 
% The common input parameters of this function are:
%
% * _studies_path_: Path to a directory containing one image folder for
% each study. The image file formats allowed are NIfTI or Analyze.
%
% * _reference_path_: Path to a directory which contains the reference
% masks or the TACs (with the same name as the study folder)
%
% * _interest_path_: Path to a directory where the masks or the TACs of the
% region of interest (with the same name as the study folder).
%
% *IMPORTANT:* The number of folders at the _studies_path_, the number of 
% files at the _reference_path_ and the number of files at the 
% _interest_path_ must be the same. The allowed formats for the masks are 
% NIfTI and Analyze.
%
% * _save_path_: Path to a directory where the results will be saved. If the
% paths is _d_, the current directory will be used.
%
% * _model_: A string indicating the name of the selected model. Throughout
% this manual, the new model to be included will be called _newModel_Name.
%
% * _parameters_: A vector containing the minimum information required to
% preprocess the selected model. Its length and values will depende on the
% model.
%
% * _timesfile_path_ (optional): Path to a .txt file containing the start and
% end times for each frame. Default times will be loaded if this optional
% file is not added and the studies do not include times series.
%
% We have written a few code lines and suggestions on the way to implement
% the new model step by step:
%% STEP 1: Include your model
% Change the _newModel_Name_ text by the name you will give to your model.
% Also, you have to change the _else_ condition.
% [Figure 1]
% And add a new case to the _if_ command:
% [Figure 2]
% After the first step, for each study including at the _studies_path_, the
% main function _QModeling_parametric_images.m_ loads this study in a
% variable called _Cpet_ and obtain the _target_ and _reference_ TACs.
% These TACs are saved to the _Ct_ and _Cr_ variables.
%% STEP 2: Check input parameters
% As an option, you could include a check for the number of input
% parameters, to make sure that it coincides with the number of parameters
% your model needs. It would also be useful to verify that the values are
% correct (non-negative values, correct numbers, between the allowed 
% maximum and minimum,...). In that case, the program should save all the
% errors detected to the log-file.
% [Figure 3]
%% STEP 3: Call your own model function
% [Figure 4 ]
% The implementation of the new model will be included in a new Matlab
% function, _QModeling_newModel_Name_. Each time a study is processed, the
% main program will call this function with the following input parameters:
%
% * the target and reference TACs: _Ct_ and _Cr_. Both of them are column
% vectors.
% * *_times_*: it contains the start and end times of each frame of the
% studies, in a _n_ by 2 matrix (_n_ being the number of
% frames of the PET studies). It is loaded from the header of the input
% studies; from the _timesfile_path_ input parameter if it is included or
% from a default times table if none of the previous options are succeeded.
% 
% * the input parameters the model needs from the user
%
% * *_Cpet_*: it contains the current study loaded. Each frame at one row
% and each voxel value at one column.
%
% * The time of execution will be controlled by the _tic-toc_ commands of
% Matlab.
%
%% STEP 4: Apply a "cut" to the generated images matrices (optional)
% It is possible to ask the user for a threshold input parameter to "cut"
% the parametric images. All the values of the parametric images under a
% certain value obtain from the threshold will be set to zero.
% [Figure 5]
%% STEP 5: Create and save the parametric images
% [Figure 6]
% Once you have run the main script, the images will be saved at the 
% _save_path_ directory.
%
%

