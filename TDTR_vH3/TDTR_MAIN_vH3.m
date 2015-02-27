%TDTR_MAIN_vH3 - Given system & thermal parameters, do thermal modeling.
% This script will do thermal modeling and print/save a figure representing
% the fit. It can handle manual, automatic, and manual in-phase fitting.
% It will also make sensitivity plots and error bar calculations upon 
% request.
%
% Syntax:  run this script through an analyze_yymmdd script, not
% independently, unless you define the expected variables first. I'd make
% it a function, but it expects a lot of variables.
%
% Variables expected to exist in the workspace:
%    ii              - index from the analyze script's for loop.
%                      Used to pick out stack(ii,:) for individual samples.
%    r_pump, r_probe - pump and probe spot sizes, in meters
%    tau_rep         - Ti:sapphire pulse repetition period
%    f               - pump/EOM modulation frequency
%    TCR             - thermoreflectance coefficient
%    fname           - string indicating data file name.
%    datain          - string indicating path to data file.
%    datadir         - string indicating path to data folder.
%    stack           - cell string matrix from the analyze script.
%    LCTE            - 4xM thermal parameter matrix
%    Xguess          - used in FIT and MANUALFIT functions
%    XguessIJ        - used in FIT and MANUALFIT functions
%    Xij             - used in FIT and MANUALFIT functions
%    Zdelay          - optional, starting time delay in ps from which to 
%                      plot data and evaluate goodness of fit.
%    tdelay_min      - minimum of delay time to model
%    tdelay_max      - maximum of delay time to model
%    P,T             - pressure and temperature. Set to -1 to assume ambient.
%    BI              - TRUE if bidirectional heat flow
%    n_toplayer      - thermal modeling parameter for bidirectional heat flow.
%                      indicates number of layers above the heat deposition
%                      plane in the sample.
%    Voutlinfit      - boolean, TRUE if V(out) is replaced with a linear
%                      fit to V(out), made in the process script.
%    m1,m2,fitOK     - Used if Voutlinfit = TRUE; linear fit parameters.
%                      Vout = m1*tdelay + m2, fitOK judges quality of fit.
%                      m1 should be in units of uV/ps, m2 in units of uV.
%    intscheme       - integration scheme: 0 = Legendre-Gauss, 1 = Rombint,
%                      2 = Simpson integration.
%    nnodes          - number of nodes for Legendre-Gauss integration;
%                      affects numerical accuracy. Don't go below 35 nodes
%                      unless you know what you're doing.
%    options         - tolerances for auto-fitting
%    senseplot       - Generate Sensitivity Plot? This option is available
%                      dynamically in MANUALFIT.
%    ebar            - TRUE if calculating Errorbars (takes longer).
%    importdata      - TRUE if fitting data. FALSE if just running 
%                      sensitivity plots or error bar calculations.
%    manualfit       - TRUE if fitting manually, FALSE if auto-fitting.
%    sigfit          - 1 if fitting V(in), 2 if fitting V(out), any other
%                      value will fit the ratio.
%    doughnut        - TRUE if fitting beam offset data.
%
%
% Important products of this script:
%    Xsol - thermal parameter fit
%    fit result figure in .fig and .eps formats, saved to datadir.
%
% Other m-files required: TDTR_REFL_vH3.m, TDTR_FIT_vH3.m,
%                         TDTR_MANUALFIT_vH3.m, TDTR_TDEP_vH3.m,
%                         SS_Heating_vH3.m, TDTR_REFL_vH3.m,
%                         TDTR_TEMP_vH3.m, INITIALIZE_CELLPARAMS_vH3.m,
%                         senseplot_vH3.m, errorbars_vH3.m,
%                         errorbar_conditions_vH3.m, genLCTEtext.m,
%                         extract_interior.m, rombint_VV3.m, SimpsonInt.m,
%                         lgwt_V4.m, mtimesx package,
%                         TDTR_TEMP_DOUGHNUT_vH3.m, GetHankel_OffsetBeam.m
% Subfunctions: none
% MAT-files required: none
%
% See also: The TDTR_TTM_vH3 package

% Author: Gregory Hohensee
% Acknowledgement: code based on Joseph P. Feser's TDTR_MAIN_V4.m script.
% University of Illinois Urbana-Champaign
% email: hohense2@illinois.edu
% Website: n/a
% Revision history: 3/25/2014 - TDTR_MAIN_vH1.m
%                   4/8/2014  - comments, harmonized with TTM version
%                   6/13/2014 - hacked to support dR/dT calculation,
%                               see lines ~240-250.
%                   7/14/2014 - vH2.
%                   11/13/2014 - added minimal beam offset functionality.
%                   1/27/2015 - Cleaned up MAIN by splitting it into
%                               sub-scripts (MAIN_MODEL, _IMPORT, _FIT,
%                               _SAVE). Cleaned up beam offset
%                               functionality (still only Vout).
%                               Sensitivity plots should work, also
%                               errorbar.
%                   2/17/2015 - vH3.
%------------- BEGIN CODE --------------
%%
% I removed keepdir functionality for this version of MAIN, because this
% MAIN assumes the data folder and fnames are all generated by the
% analyze script.

%______PROGRAM OPTIONS______________________
% These are specified in the analyze script.
%
% If you want no laser heating *and* T not at room T, just give
% writeLCTE your different T, and still set T0 = T = -1.
%___________________________________________
%% Initialize

matparams = {LCTE aniso BI n_toplayer TCR doughnut};
sysparams = {tau_rep f r_pump r_probe};
Tparams = {T0, T_LCTE, A_pump, A_probe, absC, perpulse, jabs, jtrans};
% calparams assignment needs to wait until Zind is defined.
% datparams assignment needs to wait until tdelay and data are defined.
% If psc == TRUE (varying spot size), sysparams will be updated later,
%   after tdelay is defined.

% initialize variables relating to errorbar calculation (used in TDTR_FIT)
consider_error = 0;           % MAIN does not do error bars, errorbar_vH3 does.
LCTE_err = zeros(size(LCTE)); % MAIN does not do error bars, errorbar_vH3 does.

% define tolerances for automatic fitting in FIT and sensitivity calculations
if ~exist('options','var'); options = optimset('TolFun',1e-2,'TolX',1e-2); end

%% Create standard time delay and beam offset arrays for sensitivity and errorbar plots.
% DESIGN CHOICE: I calculate sensitivities and error bars using the
% thermal model generated best fit to the data, as specified by LCTE,
% NOT with respect to the data itself. This ensures that these 
% calculations do not second-guess my judgement of the fit, are not 
% influenced by random noise, and represent sensitivities and error 
% bars with respect ONLY to systematic or thermal parameters.

if senseplot || ebar
    % MAIN_MODEL generates models for tdelay, offset, deltaR, and ratio.
    % Also vectorizes r_pump based on psc variable. 
    % Rebuilds sysparams, builds calparams, datparams.
    TDTR_MAIN_MODEL_vH3;
end

%-------------Generate Sensitivity Plot--------------
if senseplot
    tic
    fprintf('Calculating Sensitivities...Please Wait.\n')  
    % senseplot calls REFL, requres cellparams
    senseplot_offset_vH3(datparams,sysparams, calparams, matparams, Tparams, offset); 
    toc
    
    fprintf('Sensitivities calculated.\n')
end

%--------------Compute Errorbars---------------------
if ebar
    tic
    options = optimset('TolFun',1e-1,'TolX',1e-1); % default tolerances for errorbars
    fprintf('Calculating Errorbar...Please Wait.\n')    
    % errorbars_vH3 calls FIT, requires all cellparams
    [kErr_perc, kErr_abs, ...
        ErrSummary_perc, ErrSummary_abs] = ...
        errorbars_vH3(Xijc,datparams,sysparams,...
                          calparams,matparams,Tparams,options); 
    
    % save Xsol and errorbar result to file
    dlmwrite(strcat(datadir,'Xsol_ErrSummary_perc_',filename,'.txt'),...
                    vertcat(Xsol,ErrSummary_perc));
    toc
    
    fprintf('Errorbar calculated.\n')
end


%% --------------Import Data---------------------------
% reads and extracts data matrix, generates Zind, checks that ratio is 
% positive, linearizes V(out) if desired.

if importdata
    TDTR_MAIN_IMPORT_vH3; % Input: tdelay, datadir, datain, mid, sigfit, 
                          %        psc, doughnut, r_pump scalar, Zdelay,
                          %        calparams, Vout_offset*, detV_raw*, 
                          %        Voutlinfit*, fitOK*, (*:optional).
                          % Output: signal_data and offset arrays,
                          %         r_pump array, datparams, calparams.
    
    TDTR_MAIN_FIT_vH3; % Uses signal_data and offset arrays to execute fit,
                       % creates Xsol, signal_model arrays; updates LCTE.
                       
else % not importing data, just doing errorbars, sensitivities, or making figures.
    Xsol = Xguess;
    Z = 0;
    T_adj = 0;
    
    % Assign tdelay, deltaR, ratio models to the model reference "data"
    % from the sensitivity and error bar calculations.
    TDTR_MAIN_MODEL_vH3;
end

%% Hack to include dR/dT estimation: outputs are dT0 and Vin0, vs. tdelay_data %
if get_dRdT == 1    
    % dT here removes conversion to dR (dR = TCR * dT), and
    % scales dT to be dT per unit pump power absorbed, which I can
    % then re-apply using A0 = 2(1-R)*A_pump/pi.
    dT0 = real(deltaR_model) / (TCR * A_pump);
    Vin0 = (Vin_data*1e-6) ./ (detV_data * 1e-3); % SI units, from uV/mV
    % Vin_data ./ (detV * dT0) is part of the dR/dT calculation.
    % You can index by tdelay_data for the appropriate delay time to apply
    % the dR/dT calculation.
    dlnRdT = Vin0 ./ dT0;
end
% % end hack % %

%% --------------Generate save data --------------------
%% saves last fit (tdelay_data or offset, ratio_model, deltaR_model)
if doughnut, xvals = offset*1e-6; % save in meters
else xvals = tdelay*1e12; end % save in picoseconds
dlmwrite('last_fit.txt',[xvals,ratio_model,deltaR_model],'delimiter','\t')


%% Create and save MATLAB figure and print eps of fit.
% uses saveres
% uses fname for results
% uses Vin_data, ratio_data
% saves a sensitivity plot, if you asked for it in the analyze script.

saveres=input('Want to save results? (0=no, 1=yes): ');
%saveres = 1;  % TRUE to save fit result figure and/or sensitivity figure.
if saveres
    TDTR_MAIN_SAVE_vH3;
end

%% ----------------------------------------------------
fprintf('Program Completed\n')
beep
pause(0.1);
beep
pause(0.1);
beep
%---------------- END OF CODE -----------------------