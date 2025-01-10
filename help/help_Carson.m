%% Simplified Reference Tissue Model 2: Carson Model (Simplified Reference Tissue Model with fixed $k_2'$) 

%% Introduction
    %%
    % The basis for _SRTM2_ is largely the same as that for <http://www.mathworks.com _SRTM_> . 
    % Both are based on the _Full Reference Tissue Model_ (or _4
    % Parameter Reference Tissue Model_), with the hypothesis that:
    % 
    % * The distribution volume is the same for the _tissue of interest_
    % (tissue including specific binding) and the _reference tissue_ (tissue
    % with no specific binding): $K_1/k_2=K_1'/k_2'$
    % 
    % * The kinetics in the tissue with specific binding is such that it is
    % difficult to distinguish between the specific and the
    % free/non-specific compartments. This happens when the exchange
    % between the free and the specifically bound compartments is
    % sufficiently rapid. In these cases the tissue region of interest may
    % be approximated by a single compartment, with efflux constant
    % $k_{2a}=k_2/(1+BP)$, where $BP=k_3/k_4$ is the binding potential.
    % 
    % [Figure srtm_model.png]
    %
    % An additional restriction is added in the approach for _SRTM2_, in
    % order to reduce noise and accelerate the generation of parametric
    % images. The _SRTM_ operational equation is:
    %
    % $$C_T(t)=R_1C_R(t)+R_1(k_2'-k_{2a})C_R(t)\otimes{}e^{-k_{2a}t} \textbf{   (1)}$$ 
    %
    % *_Implementation_*
    %
    % Since $k_2'=k_2/R_1$ is the clearance rate constant from the
    % _reference region_, it has only one true value, and it should not be
    % necessary to estimate it for every pixel, as _SRTM_ does. The _SRTM2_
    % implementation [1] uses this property to reduce noise in the
    % calculations, by considering a fixed $k_2'$ parameter. In our
    % implementation, the value for $k_2'$ can be established in two ways:
    % either it is directly specified by the user, or it is estimated just
    % once, by applying _SRTM_ to the mean TAC of a region with high specific
    % uptake, using the appropriate TAC of the reference region as well.
    % From here on, the obtained $k_2'$ is considered fixed, and, thus, the
    % parametric images are generated with a fixed $k_2'$. The other
    % parameters, $k_{2a}$ and $R_1$, are estimated using the basis
    % functions approach.
    %
    % The basis functions can be defined as:
    %
    % $$ B_i=C_R(t)+(k_2'-k_{2a,i})C_R(t)\otimes{}e^{-k_{2a,i}t} \textbf{   (2)} $$
    % 
    % Thus, there is one basis function for each value of $k_{2a}$.
    % Equation *(1)* can now be written as: 
    %
    % $$ C_T(t)=R_1 B_i(t) \textbf{   (3)} $$
    %
    % Which shows that a least-squares fit allows calculating $R_1$ for
    % each basis function. After the index $i$, which minimizes the
    % deviation between the TAC and the model curve is determined, the
    % values of $k_{2a}$ and $R_1$ are obtained. Values for $BP$ and $k_2$,
    % are then easily deduced.
    %
    % The user can select a logarithmic range of values for $k_{2a}$, to
    % take into account all plausible values for this parameter.
%% Preprocessing algorithm   
    %%
    % # Calculation of the basis functions: convolution of the reference
    % TAC with decaying exponentials in the range $[k_{2a}min, k_{2a} max]$ 
    % and at the resolution (_Resampling_) selected by the user.
    % # Least squares fit of $R_1$ for each of the basis functions. This 
    % results in a set of optimal parameters and an estimated model curve 
    % for each basis function.
    % # The fit with minimal deviation between receptor-rich TAC and model 
    % curve is regarded as the result. The parameters of interest
    % can be calculated from the fitted values.
   
%% Input parameters
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
    % with energy below the specified percentage of the maximal energy will be masked to zero.
    %
    % * $\bf{k_2'}$: the clearance rate constant from the reference region. 
    
    
%% Output parameters and goodness of fit
    %%
    % * $\bf{BP}$: Binding Potential relative to non-displaceable uptake $(BP=k_3/k_4)$   
    %
    % * $\bf{k_2}$: Efflux constant for the free/non-specific binding
    % compartment in the region of interest.
    %
    % * $\bf{R_1}$: Relative tracer delivery $(R_1=K_1/K_1')$
    %
    % * $\bf{k_{2a}}$: Apparent efflux constant for the tissue of interest when
    % considered as a single compartment $(k_{2a}=k_2/(1+BP)$
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
% * *Correlation coefficient (Corr. Coef.)*: The correlation coefficient
% between the values of the TAC of interest and the values of the TAC 
% estimated by the model. Values closer to 1 are better.
%
%% Image generation algorithm
%%
% # The basis functions are already calculated (see preprocessing
% algorithm).
% # Voxel-wise least squares fit of $R_1$ for each of the basis functions. 
% This results in a set of optimal parameters for each voxel.
% # The image for each of the selected parameters is written.
%
%% Main references
% 
% [1] Wu, Y., & Carson, R. E. (2002). Noise reduction in the simplified 
% reference tissue model for neuroreceptor functional imaging. 
% _Journal of Cerebral Blood Flow and Metabolism_, 22(12), 1440–1452
%%

