%% SRTM: Simplified Reference Tissue Model 

%% Introduction
% The _SRTM_ method is mostly used for receptor studies using reversibly 
% binding tracers. It was developed by Lammertsma [1], on the basis of the 
% _Full Reference Tissue Method_ (or _4 Parameter Reference Tissue Method_) [2].
%
% The formulation involves a _reference region_ devoid of specific binding,
% modeled with a one-tissue compartment, and a region with specific binding
% (_region of interest_), which is represented with a two-tissue
% compartment model. The rate constants are $K_1$ and $k_2$, representing
% the exchange of tracer between plasma and a free ligand compartment in
% the _region of interest_; $K_1'$, $k_2'$, which are the equivalent for the
% _reference region_; and $k_3$ and $k_4$, which represent the exchange of
% tracer between the free compartment and a specifically bound ligand
% compartment in the _region of interest_.
% 
% [Figure srtm_model.png]
%
% The model relies on the following assumptions:
%
% * The distribution volume is the same for the tissue of interest (tissue 
% including specific binding) and the reference tissue: $$ K_1/k_2=K_1'/k_2' $$ .
%
% * The kinetics in the tissue with specific binding is such that it is
% difficult to distinguish between the specific and the free/non-specific
% compartments. This happens when the exchange between the free and the
% specifically bound compartments is sufficiently rapid. In these cases the
% tissue region of interest may be approximated by a single compartment,
% with efflux constant $k_{2a}=k_2/(1+BP)$, where
% $BP=k_3/k_4$ is the binding potential.
%
% Defining $$ R_1=K_1/K_1' $$ as the ratio of tracer delivery, the 
% following operational equation can be derived for the measured TAC
% in the receptor-rich region:
%
% $$C_T(t)=R_1C_R(t)+R_1(k_2'-k_{2a})C_R(t)\otimes{}e^{-k_{2a}t}  \textbf{   (1)}$$ 
%
% The three unknowns, $R_1$, $k_2'$ and $k_{2a}$, in this equation can be fitted 
% using nonlinear regression techniques.
%
% *_Implementation_*
%
% The implementation of the _SRTM_ model developed by Gunn [3], using basis
% functions, was better suited for a pixel-wise application than the
% original approach. Equation *(1)* can be rewritten as:
%
% $$C_T(t)=\theta_1 C_R(t)+ \theta_2 C_R(t) \otimes e^{-\theta_3t}  \textbf{   (2)}$$
%
% Where $\theta_1=R_1$, $\theta_2=R_1(k_2'-k_{2a})$ and $\theta_3=k_{2a}$.
% Since this equation is linear on $\theta_1$ and $\theta_2$, these values
% can be estimated using standard linear least squares, when the value of
% $\theta_3$ is fixed. To obatin a solution for the nonlinear term, a
% discrete set of parameter values for $\theta_3$ can be chosen, to form
% the following basis functions:
%
% $$B_i(t)=C_R(t) \otimes e^{-\theta_{3,i}t}  \textbf{   (3)}$$
%
% Equation *(2)* can then be transformed into a linear equation for each
% basis function:
%
% $$ C_T(t)=\theta_1 C_R(t) + \theta_2 B_i(t) \textbf{   (4)}$$
%
% In this approach, which we have implemented in _QModeling_, equation *(4)* is
% solved using linear least squares for each $B_i$. After the index $i$,
% which minimizes the deviation between the TAC and the model curve is
% determined, the values of $\theta_1$, $\theta_2$ and $\theta_3$ are
% obtained. Values for $BP$, $R_1$ and $k_2$ are then easily deduced.
%
% The user can select a logarithmic range of values of $\theta_3$, to take
% into account all plausible values for this parameter.
%
%% Preprocessing algorithm   
%%
% # Calculation of the basis functions: convolution of the reference
%    TAC with decaying exponentials in the range $[k_{2a}min, k_{2a} max]$ 
%  and at the resolution (_Resampling_) selected by the user.
% # Least squares fit for each of the basis functions. This results in a
% set of optimal parameters and an estimated model curve for each basis
% function.
% # The fit with minimal deviation between receptor-rich TAC and model 
% curve is regarded as the result. The parameters of interest
% can be calculated from the fitted values.
    
%% Input parameters
%
% * *TAC 1*: TAC from a region with specific uptake (typically a receptor
% rich area).
%
% * *TAC 2*: TAC from a reference region with no specific uptake
% (typically, a region devoid of target receptors).
%
% * $\bf{k_{2a}min}$: Minimal value of $k_{2a}$ (slowest decay of exponential).
%  
% * $\bf{k_{2a}max}$: Maximal value of $k_{2a}$ (fastest decay of exponential).
%
% * *#Basis*: Number of basis functions. Each basis funtion is defined by
% its $k_{2a}$ value. The values of $k_{2a}$ for the basis functions will
% be in the range between $k_{2a} min$ and $k_{2a} max$, and there will be
% a total of #Basis different values of $k_{2a}$. Increments will be taken
% at logarithmic steps. This number is directly proportional to processing
% time: the bigger #Basis, the longer the processing time. If the number of
% basis functions is too low, the estimation may lack precision. However,
% increasing its value does not indefinitely improve the estimation.
%
% * *Resampling*: It specifies the interval at which the TACs will be 
% resampled, before convolving them with the exponentials to form the basis
% functions. This interval should be equal or smaller than the shortest 
% frame duration. Notice that the bigger this value, the shorter the 
% processing time. However, smaller intervals lead to more accurate 
% estimations. 
%  
% * *Threshold*: Discrimination threshold for background masking. All pixels 
%              with energy below the specified percentage of the maximal 
%              energy will be masked to zero.

%% Output parameters and goodness of fit
%
% * $\bf{BP}$: Binding Potential relative to non-displaceable uptake  $$ (BP = k_3/k_4) $$
%
% * $\bf{k_2}$: Efflux constant for the free/non-specific binding compartment in
% the region of interest
%
% * $\bf{R_1}$: Relative tracer delivery  $$ (R_1= K_1/K_1') $$ 
%
% * $\bf{k_{2a}}$ : Apparent efflux constant for the tissue of interest when
% considered as a single compartment. $$ (k_{2a}=k_2/(1+BP)) $$
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
% * *Correlation coefficient (Corr. Coef.)*: The correlation coefficient between the values 
%                                            of the TAC of interest and the values of the TAC 
%                                            estimated by the model. Values closer to 1 are better.
%
%% Image generation algorithm
%%
% # The basis functions are already calculated (see preprocessing
% algorithm).
% # Voxel-wise least squares fit for each of the basis functions. This
% results in a set of optimal parameters for each voxel.
% # The image for each of the selected parameters is written.
%
%% Main references
% 
% [1] Lammertsma, A. A., & Hume, S. P. (1996). Simplified Reference Tissue Model for PET
% Receptor Studies. _NeuroImage_ , 4(3), 153–158
% 
% [2] Lammertsma AA, Bech CJ, Hume SP, Osman S, Gunn K, Brooks DJ,
% Frackowiak RS (1996). Comparison of methods for analysis of clinical
% [11C]raclopride studies. _Journal of Cerebral Blood Flow and Metabolism_ , 
% 16(1):42-52
%
% [3] Gunn, R. N., Lammertsma, A. A., Hume, S. P., & Cunningham, V. J. (1997). Parametric
% Imaging of Ligand-Receptor Binding in PET Using a Simplified Reference Region Model.
% _NeuroImage_ , 6(4), 279–287.
%%

