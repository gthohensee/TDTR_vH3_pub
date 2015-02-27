%% Type and conditions of modeling
P0 = 1;         % GPa initial. P = -1 assumes atmospheric pressure.

BI = 1;         % boolean: TRUE if using bidirectional heat flow model
n_toplayer = 3; % # of model layers above the point where heat is deposited.

Voutlinfit = 0; % Did you linearize the V(out) for any of the data?
if Voutlinfit
    Voutfitparams = dlmread(strcat(directory, datafolder, yymmdd, ...
                            '_Voutfitparams_', tagline, '.txt'));
end

%% Common sysparams - TDTR system settings
f=9.8e6; % pump laser modulation frequency, Hz
tau_rep=1/80e6; % laser repetition period, s.
                % The repetition frequency is 80 MHz for TDTR-1, 
                % approximately 74.8 MHz for TDTR-2.
TCR=1e-4; % default coefficient of thermal reflectance
pm = 1.08; % optical power is 1.08x the reading of the model 835 power meter at TDTR-1.
%pm = 1.00; % scaling for readout of TDTR-2 power meter?

%% calparams - calculation parameters
Zdelay = 100; % ps starting time. In automatic fitting, goodness-of-fit
             % is calculated between Zdelay and tdelay_max. In V(in)
             % fitting, Zdelay indicates the time at which to
             % normalize the V(in) signal.

sigfit = 0; % 0: ratio fit, 1: V(in) fit, 2: V(out) fit.

intscheme = 0; % integration scheme: 0 = Legendre-Gauss, 1 = Rombint,
               % 2 = Simpson integration.
nnodes = 35;  % number of nodes for Legendre-Gauss or Simpson integration;
              % affects numerical accuracy. Don't go below 35 nodes
              % unless you know what you're doing.

%% Tparams - parameters for temperature-dependent samples
T0 = -1; % nominal (K). Set T0 = -1 to ignore laser heating and assume room temperature.

% optional parameters:
absC = [295 0.13]; % Transducer absorption coefficient (0.13 for room temperature Al)
perpulse = 1; % TRUE if accounting for per-pulse heating as well as steady-state

%% other operational parameters for MAIN
% Time delay boundaries for thermal modeling AND data presentation.
tdelay_min= 100e-12; % warning: the smaller this is, the longer the
                    % computation time of the thermal model.
tdelay_max= 4000e-9; % approx. max extent of our delay stage.

options = optimset('TolFun',1e-2,'TolX',1e-2); % tolerances for autofitting
                                               % and for error bars.


senseplot = 0; %Generate Sensitivity Plot? This option is available
               %dynamically in MANUALFIT.
              
ebar = 0; %Calculate Errorbars? 0 for no, 1 for yes (takes longer)

importdata = 1; % TRUE if fitting data. FALSE if just running sensitivity
                % plots or error bar calculations.
                
manualfit = 1;  % TRUE if fitting manually, FALSE if auto-fitting.

get_dRdT = 0; % TRUE if extracting dR/dT from the TDTR data.
              % Requires accurate measurement of incident pump
              % and probe powers, as well as photodiode voltage vs. delay
              % time. Keep this off if you don't fully understand it.