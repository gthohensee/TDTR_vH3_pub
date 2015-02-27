
%% Define the directory
directory = '/Users/gregoryhohensee/Documents/MATLAB/TDTR_vH2_1C/TDTR_vH2'; % pwd is the current directory
yymmdd = '150201'; % e.g., 140227 = February 27th, 2014.
savefolder = strcat('/',yymmdd,'_parametric/'); % folder in directory where sensitivities are saved
savedir = strcat(directory, savefolder);

%% Common sysparams - TDTR system settings
f=10e6; % pump laser modulation frequency, Hz
tau_rep=1/80e6; % laser repetition period, s.
                % The repetition frequency is 80 MHz for TDTR-1, 
                % approximately 74.8 MHz for TDTR-2.
TCR=1e-4; % default coefficient of thermal reflectance

%% Beam intensity and spot size conditions
psc = 0;     % TRUE if pump spot changes over delay time.
frac = 0.19; % if psc, fractional change in pump spot change over 
             % USER-SPECIFIED range of delay time.
             
objective = 10; % strength of objective lens

switch objective % YKK values
    case 2, r_pump = 2*15e-6; tc = 0.87; 
    case 5, r_pump = 15e-6;     tc = 0.90; 
    case 10, r_pump = 7.5e-6;   tc = 0.80;
    case 20, r_pump = 3.8e-6;  tc = 0.70;
    %case 50, r_pump = 1.4e-6;  tc = 0.70; % tc = ?? for 50x
    %case 100, r_pump = ?e-6;   tc = 0.??;
end 
r_probe = r_pump;         % probe 1/e^2 radius, meters.

% Laser powers that reach a typical in-air sample.
pm = 1; % deviation of powermeter reading from actual average laser power pre-objective.
A_pump = 10 * 1e-3*pm*tc;
A_probe = 5 * 1e-3*pm*tc;

%% General thermal modeling conditions
Zdelay = 100; % Normalization time delay (ps) for V(in) and V(out).
              % Irrelevant if calculating by ratio or beam offset.
              
BI = 0;         % boolean: TRUE if using bidirectional heat flow model.
n_toplayer = 3; % for BI: # of model layers above the point where heat is deposited.

intscheme = 0; % integration scheme: 0 = Legendre-Gauss, 1 = Rombint,
               % 2 = Simpson integration.
               % Beam offset only works with Legendre-Gauss for now.
nnodes = 35;  % number of nodes for Legendre-Gauss or Simpson integration;
              % affects numerical accuracy. Don't go below 35 nodes
              % unless you know what you're doing.

T0 = -1; % nominal (K). Set T0 = -1 to ignore laser heating and assume room temperature.
perpulse = 0; % TRUE if accounting for per-pulse heating as well as steady-state
% optional parameters:
absC = [295 0.13]; % Transducer absorption coefficient (0.13 for room temperature Al)
                   % (used in place of TCR if self-correcting for
                   % temperature-sensitive thermal parameters).

consider_error = 0; % Boolean; 0 if not performing errorbar calculations.
%% other operational parameters for MAIN
% Time delay boundaries for thermal modeling AND data presentation.
tdelay_min= -10e-12; % warning: the smaller this is, the longer the
                    % computation time of the thermal model.
tdelay_max= -10e-12; % approx. max extent of our delay stage.
N_tdelay = 1; % number of logspace time delay points between min and max.

options = optimset('TolFun',1e-2,'TolX',1e-2); % tolerances for autofitting
                                               % and for error bars.

%% Construct the thermal material properties matrix.
thickness = [1 90 1 1e9]; % thickness in nm of each model layer
stack = {'*','Al', '/', 'Si'}; % identifier for each model layer
Xijc = 0; % LCTE(i,j) indices of the "c" fit parameters; 0 if none.
P0 = -1; % Pressure (GPa); -1 if not high pressure.
T0 = -1; % Temperature (K); -1 if room temperature.
refdir = ''; % directory for T-dependent thermal parameter data files
rho_Al = 4.0e-8; % Al transducer resistivity; -1 for default value.

[LCTE,T_LCTE] = writeLCTE_vH3(stack, thickness,Xijc,P0,T0,refdir,rho_Al);
% could also have written: LCTE = writeLCTE_vH3(stack,thickness).

jabs = 1;           % index of absorption layer (for T-dependence and coupled sensitivities)
jtrans = jabs + 1;  % index of transducer layer (for T-dependence and coupled sensitivities)

% each column of aniso is TRUE if corresponding layer is anisotropic.
% aniso is used in sensitivity calculation to decide how to decouple eta and lambda.
aniso = zeros(1,length(LCTE(1,:)));

%% Choose type of signal: V(out) offset, V(in)(t), V(out)(t), or ratio(t)?
doughnut = 0;   % boolean: TRUE if calculating beam offset sensitivities.
sigfit = 2; % only relevant if not doughnut: 1 == Vin, 2 == Vout, 0 = Ratio

%% Assemble cell params.
matparams = {LCTE aniso BI n_toplayer TCR doughnut};
sysparams = {tau_rep f r_pump r_probe};
Tparams = {T0, T_LCTE, A_pump, A_probe, absC, perpulse, jabs, jtrans};
% calparams assignment needs to wait until Zind is defined.
% datparams assignment needs to wait until tdelay and data are defined.
% If psc == TRUE (varying spot size), sysparams will be updated later,
%   after tdelay is defined.

%% Tweak LCTE
LCTE(1,3) = 0.05; % 50 MW
LCTE(2,4) = 1.66e6;
LCTE(1,4) = 5; % 5 W/m-K
matparams{1} = LCTE;

%% produce model tdelay and/or offset arrays, plus r_pump array (for psc),
% and pack up the rest of the cell params.
TDTR_MAIN_MODEL_vH2; 
%% Define independent variable
xvar = logspace(log10(0.1e6),log10(30e6),10); % range of xvar parameter values in LCTE (SI) units
Xij = 0;
%Xij = [1,4]; % index of LCTE matrix of xvar parameter.
SYS_i = 1; % Scalar. [1,2,3,4] indicates xvar is [f, r_pump, r_probe, tau_rep], respectively.
           % Zero values in Xij or SYS_i indicate that xvar is not in LCTE
           % or a system parameter, respectively.

%% Which sensitivities to calculate?
LCTE_sens_consider = zeros(4,length(LCTE(1,:)));
LCTE_sens_consider(1,4) = 1;
LCTE_sens_consider(3,2) = 1;
sys_consider = [0 0 0 0 1]; %[f, r_pump, r_probe, tau_rep, w0]

%% Execute sensitivity calculations and plotting.
tic
[S_LCTE,S_sys,xvar,SS_LCTE,SS_sys] = ...
    parametric_senseplot_vH3_notN(Xij,SYS_i, xvar, ...
                   datparams,sysparams, calparams, matparams, Tparams, ...
                   LCTE_sens_consider, sys_consider);
toc
