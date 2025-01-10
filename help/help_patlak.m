%% Patlak Reference Plot 

%% Introduction
% 
% Tracers undergoing irreversible trapping can be analyzed with the _Patlak
% Reference Model_. The _Patlak plot_ [1] can be used as a reference model
% provided there is also some tissue where the tracer is not irreversibly
% bound. This is frequently applied in the analysis of [18F]-FDG studies,
% where a 2-tissue compartment model is used, with $k_4=0$
% 
% [Figure patlak_model.png]
% 
% However, the _Patlak Reference Model_ is also applicable to other
% compartment model situations, as long as there is one tissue compartment
% with irreversible trapping and a suitable _reference tissue_. Among other
% hypothesis, it is assumed that the distribution volume is the same for
% the _tissue of interest_ (tissue including irreversible trapping) and the
% _reference tissue_: $K_1/k_2=K_1'/k_2'$.
%
% The operational equation for the model is:
%
% $$ \frac{C_{tissue}(t)}{C_{ref}(t)}=K\,
% \frac{\int_{0}^{t}C_{ref}(u)du}{C_{ref}(t)} + V \textbf{   (1)} $$
%
% Where $C_{tissue}(t)$ is the TAC for the region including irreversible
% trapping and $C_{ref}$ is the TAC for the _reference region_.
%
% In this approach, the measured PET activity $C_{tissue}(t)$ is divided by
% the activity in the _reference tissue_. This value is plotted against the
% quotient at the other side of the equation: the integral of the TAC at
% the _reference region_ divided by the activity in this same region. The
% latter acts as a sort of "normalized time". If the system has an
% irreversible compartment, this relation becomes linear after some time t*
% after injection, which allows estimating the _slope_ and the _intercept_ with
% a linear fit. The interpretation of these parameters depends on the
% particular configuration of the system.
%
% When the compartment model in the figure is followed, the _slope_ is:
%
%
% $$ slope=K=\frac{k_2k_3}{k_2+k_3} \textbf{   (2)} $$
%
%
% The _intercept_ is also a function of the transfer constants of the
% system, but it is not usually taken into account for the analyses.


%% Input parameters
% 
% * *TAC 1*: TAC from a region of interest (region including irreversible trapping) 
%
% * *TAC 2*: TAC from a reference region
%
% * *Max. Err.*:  Maximum relative error allowed between the linear 
% regression and the Patlak-transformed measurements
% $(C_{tissue}(t)/C_{ref}(t) vs. \int_{0}^{t}C_{ref}(u)du/C_{ref}(t))$.
%
% * *t**:    The linear regression estimation should be restricted to a 
% range after an equilibration time. Parameter t* marks the beginning of 
% the range used in the linear regression analysis. It can
% be fitted based on the _Max. Err._ criterion. If so, the model is fit with
% an initial t* starting at the beginning of the scan, and the estimation
% error is compared to the specified _Max. Err._ value. If the obtained
% error is larger, then the next possible t* is tried. The process
% continues until the obtained error is lower than the specified _Max.
% Err._ Alternatively, t* can be directly specified by the user. Parameter
% t* is an actual acquisition time value, even though a "normalized time"
% is used at the rigth-hand side of the equation.
%
% * *Threshold*:    Discrimination threshold for background masking.
%  
%% Output parameters and goodness of fit
%
% * *t** : It marks the beginning of the range of times used in the linear
% regression analysis. If the user does not directly specify t*, it is
% estimated by the program on the basis of a maximum allowed error. This
% happens at the preprocessing stage, when the model is fit to the mean TAC
% of a region with irreversible trapping. The value obtained in this way is
% then used as an input for the calculation of the parametric images.
%
% * *K* :	Slope of the linear regression. 
%
% * *Intercept* :	Intercept of the linear regression.
%
% * *Start*: "Normalized time" that corresponds to t*.
%   
% * *Max. Error*: Maximum relative error obtained between the linear 
%           regression and the Patlak-transformed measurements in the 
%           segment starting from t*.
%
% - *_Goodness of fit_*: 
%
% To show the goodness of fit at the preprocessing step, two parameters are given, together with the estimated parameters:
%
% * *Normalized Mean Squared Error (NMSE)*:
%
% $$  NMSE=\frac{||C_t - C_{t_{estimate}}||^2}{||C_t - mean(C_t)||^2} $$
%
% It measures the quality of the fit for the TACs. Values vary between 
% $-\infty$ (bad fit) to 1 (perfect fit).
%
%
% * *Correlation coefficient (Corr. Coef.)*: The correlation coefficient
% between the values of the TAC of interest and the values of the TAC 
% estimated by the model. Values closer to 1 are better.
%
%% Main references
%
% [1] Patlak, C. S., Blasberg, R. G., & Fenstermacher, J. D. (1983). 
% Graphical Evaluation of Blood-to-Brain Transfer Constants from Multiple-Time 
% Uptake Data. _Journal of Cerebral Blood Flow and Metabolism_, 3(1), 1–7.
%
% [2] Patlak, C. S., & Blasberg, R. G. (1985). Graphical evaluation of 
% blood-to-brain transfer constants from multiple-time uptake data. 
% Generalizations. _Journal of Cerebral Blood Flow and Metabolism_, 
% 5(4), 584–590.
%
%%

