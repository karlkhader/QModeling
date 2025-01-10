%% _QModeling_preprocessTACs.m_
% 
% The common input parameters of this function are:
%
% * _TACsfile_: Path to .txt file containing the TACs to be analyzed. This
% file must be organized as follows. The first row must contain the
% mid-times of each frame. The following rows are organized in pairs, each
% pair corresponding to a different PET study. The first row of each pair
% contains the TACs of the reference region. The second row in the pair
% contains the TAC of the region of interest of the same PET study.
%
% * _TimeTable_: Path to a .txt file containing the start and end times for
% each frame.
%
% * _save_path_: Path to a directory where the results will be saved. If
% the path is "d", the current directory will be used.
%
% * _model_: A string indicating the name of the selected model. Throughout
% this manual, the new model bo be included will be called "newModel_Name".
%
% * _parameters_: A vector containing the minimum information required to
% preprocess the selected model. Its length and values will depend on the
% model.
%
% Near the end of the script we have written a few code lines and
% suggestions on the way to implement the new model step by step:
%
%% STEP 1: Include your model
% [Figure 1]
% Change the _newModel_Name_ text by the name you will give to your model.
%
%% STEP 2: Check input parameters
% [Figure 2]
% As an option, you could include a check for the number of input
% parameters, to make sure that it coincides with the number of parameters
% your model needs. It would also be useful to verify that the values are
% correct (non-negative values, between the allowed maximum and
% minimum,...). In that case, the program should display the values of the
% input parameters like it is done with the implemented models.
%% STEP 3: Call your own model function
% [Figure 3]
% At the start of the script, some variables are created. In this step, you
% will need the following ones: 
%
% * *_data_*: it contains the _TACs_ loaded from the _TACsfile_ input file.
%
% * *_num_sim_*: it is the number of studies to preprocess.
%
% * *_times_*: it contains the start and end times of each frame of the
% studies to preprocess, in a _n_ by 2 matrix (_n_ being the number of
% frames of the PET studies). It is loaded from the "TimeTable" input
% parameter.
%
% Furthermore you have to create a new variable to save the results of the
% preprocess: _results_model_matrix_. This variable will contain one row
% for each study to be preprocessed and one column for each parameter the
% model returns. Two additional columns are required to include the mean
% square error and the correlation coefficient obtained from the fit.
%
% The implementation of the new model will be included, as previously
% mentioned, in a new Matlab function _QModeling_newModel_Name_. Each time
% a study is preprocessed, the main program will call this function with
% the following input parameters:
%
% * the target and reference TACs: _Ct_, _Cr_. Both of them are column
% vectors.
%
% * the times: The same _times_ matrix described in the input parameters
% list for _QModeling_preprocessTACs.m_
%
% * the input parameters the model needs from the user.
%
% The time of execution will be controlled by the _tic-toc_ commands of
% Matlab.
%% STEP 4: Save the results
% [Figure 4]
% Once you have run the main script, the results should be saved in a text
% file called _PreprocessingModelName_Results_ at the _save_path_
% directory.
%
%

