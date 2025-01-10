%% Logan Reference Plot with fixed $k_2'$

%% Introduction
%
% The Logan Reference Tissue method is used to estimate the  
% distribution volume ratio (DVR) for PET studies with reversible  
% radioligands. $BP$ is then derived from $DVR$ as $BP = DVR-1$. 
% 
%
% The method does not depend on a specific model structure of the reference tissue. 
% Assuming the presence of reference region TAC $C_{T'}(t)$ with 
% an average tissue-to-plasma clearance $k_2'$, the target tissue TAC $C_T(t)$ 
% is transformed and plotted as a function of the 
% transformed reference TAC, which acts as a sort of "normalized time".
%
% The operational equation for the method is:
%
%
% $$ \frac{\int_{0}^{t}C_{T}(\tau)d\tau}{C_{T}(t)}=DVR \,
% \frac{\int_{0}^{t}C_{T'}(\tau)d\tau + C_{T'}(t)/k'_{2} }{C_{T}(t)} + b $$
%
% 
% This resembles a linear equation, where $DVR$ is the slope and $b$ is an error term 
% which decreases over time. This dependence also becomes negligible after
% some time t*. Then, if only the values obtained after t* are taken into account, 
% a linear regression allows estimating $DVR$ and $b$.
%
% As mentioned, $k_2'$ represents the average tissue to plasma efflux.  
% It has to be specified as an input for the model. A representative value  
% for the population based on previous studies can be entered, or,  
% alternatively, $k_2'$ can be estimated for each subject by applying,  
% for instance, the SRTM method, before starting the Logan analysis.

 

%% Input parameters
%
% * *TAC 1*: TAC from a region with specific uptake (typically a receptor
% rich area).
%
% * *TAC 2*: TAC from a reference region with no specific uptake
% (typically, a region devoid of target receptors).
%
% * *t**: The linear regression estimation should be restricted to a 
% range after an equilibration time. t* marks the beginning of 
% the range used in the linear regression analysis. 
% The graphical representation of the data provided by the program
% can allow the user to select appropriate values for t* based 
% on a visual analysis. This parameter can also be fitted automatically
% using the _Max. Err._ criterion. In this case, the lowest t* yielding 
% an error lower than _Max. Err._ (see below) will be chosen.
% Note that the t* is specified in acquisition time units. 
% The program later translates this time into the units of the x axis
% (the "normalized time") and shows it as the _Start_ parameter.
%
% * *Max. Err.*: Maximum relative error allowed between the linear 
% regression and the Logan Plot-transformed measurements in the 
% segment starting from t*.
%
% * *Threshold*: Discrimination threshold for background masking.
% 
% * $\bf{k_2'}$: an average value of the efflux rate constant from 
%   regions without receptors, which has been previously determined.
%  
%% Output parameters and goodness of fit
%
% * *BP*: Binding potential obtained as $BP = DVR-1$ from the 
% slope of the linear regression. 
%
% * *Intercept* :Intercept of the linear regression.
%
% * *Start*: The value of the x axis corresponding to the time point t*.
%   
% - *_Goodness of fit_*: 
%
% To show the goodness of fit at the preprocessing step, two parameters are given, together with the estimated parameters:
%
% * *Normalized Mean Squared Error (NMSE)*: 
% 
% $$ NMSE=\frac{||C_t - C_{t_{estimate}}||^2}{||C_t - mean(C_t)||^2} $$
%
% It measures the quality of the fit for the TACs. Values vary between $-\infty$ (bad fit) to 1 (perfect fit).
%
%
% * *Correlation coefficient (Corr. Coef.)*: The correlation coefficient between the values 
%                                            of the TAC of interest and the values of the TAC 
%                                            estimated by the model. Values closer to 1 are better.
%
%% Main references
% 
% [1] Logan J, Fowler JS, Volkow ND, Wang GJ, Ding YS, Alexoff DL (1996)
% Distribution volume ratios without blood sampling from graphical analysis 
% of PET data. _Journal of Cerebral Blood Flow and Metabolism_, 16(5):834-840. 
