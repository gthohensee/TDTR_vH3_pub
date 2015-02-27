%sense_setup_retro - Template for calling parametric_senseplot_vH3.
% Termed "retro" because the process of assigning the thermal matrix
% parameters is the same as in Joe Feser's original MATLAB code, with
% lambda, C, t, eta arrays, instead of my LCTE system.
%
% Inputs:
%    none
%
% Outputs
%    See outputs of parametric_senseplot_vH3.
%
% Other m-files required: parametric_senseplot_vH3.m
% Subfunctions: none
% MAT-files required: none
%
% See also: sense_setup.m, parametric_senseplot_vH3.m

% Author: Gregory Hohensee
% U. of Illinois Urbana-Champaign
% email: hohense2@illinois.edu
% Website: n/a
% Revision history: 17-Feb-2015 - header comments added.
%------------- BEGIN CODE --------------

%% Define the directory
directory = '/Users/gregoryhohensee/Documents/MATLAB/TDTR_vH3'; % pwd is the current directory
yymmdd = '150201'; % e.g., 140227 = February 27th, 2014.
savefolder = strcat('/',yymmdd,'_parametric/'); % folder in directory where sensitivities are saved
savedir = strcat(directory, savefolder);

%% Old-fashioned manual format for entering system and thermal parameters.
% Inspired by TDTR_MAIN_VV3 (Joe Feser's code).

%Anticipated system properties (initial guess for fitting, if you do fitting/errorbar estimation)
abslayer = 40;
lambda= [150*abslayer   150     0.2         100]; %W/m-K
C=      [2.43*abslayer  2.43    0.1         2]*1e6; %J/m^3-K
t=      [1           80-abslayer  1         1e9]*1e-9; %m 
% ("Al" thickness estimated at tAl=t0+3nm oxide, t0 from picosecond acoustic)
lambda(1)=abslayer*lambda(2); %top 1nm 10X the rest of layer (simulates absorbtion)
C(1)=abslayer*C(2); %top 1nm 10X rest of layer
eta=    ones(1,numel(lambda)); %eta is the anisotropic ratio, eta=kx/ky;

LCTE = vertcat(lambda,C,t,eta); % assembled thermal parameters

f=10e6; %laser Modulation frequency, Hz
r_pump=10e-6; %pump 1/e^2 radius, m
r_probe=10e-6; %probe 1/e^2 radius, m
%Choose range of time delays to fit
tdelay_min=4000e-12; % 100e-12
tdelay_max=4000e-12; % 4e-9
N_tdelay = 1; % number of logspace time delay points between min and max.

tau_rep=1/80e6; %laser repetition period, s
A_pump=20e-3; %laser power (Watts) . . . only used for amplitude est.
A_probe=10e-3; %probe laser power (Watts) . . . only used for amplitude est.
TCR=1e-4; %coefficient of thermal reflectance . . . only used for amplitude est.
tdelay=logspace(tdelay_min,log10(tdelay_max),N_tdelay)'; %vector of time delays (used to generate sensitivity plots)

%% Default values for bells and whistles used in the back-end
Zdelay = 100; % Normalization time delay (ps) for V(in) and V(out).
              % Irrelevant if calculating by ratio or beam offset.

psc = 0; % TRUE if pump spot changes over delay time.
frac = 0.19; % if psc, fractional change in pump spot change over 
             % USER-SPECIFIED range of delay time.

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
                   % temperature-sensitive theansrmal parameters).
                   
jabs = 1;           % index of absorption layer (for T-dependence and coupled sensitivities)
jtrans = jabs + 1;  % index of transducer layer (for T-dependence and coupled sensitivities)

% each column of aniso is TRUE if corresponding layer is anisotropic.
% aniso is used in sensitivity calculation to decide how to decouple eta and lambda.
aniso = zeros(1,length(LCTE(1,:)));
%% Choose type of signal: V(out) offset, V(in)(t), V(out)(t), or ratio(t)?
doughnut = 0;   % boolean: TRUE if calculating beam offset sensitivities.
sigfit = 0; % only relevant if not doughnut: 1 == Vin, 2 == Vout, 0 = Ratio

%% Assemble cell params.
matparams = {LCTE aniso BI n_toplayer TCR doughnut};
sysparams = {tau_rep f r_pump r_probe};
Tparams = {T0, T_LCTE, A_pump, A_probe, absC, perpulse, jabs, jtrans};
% calparams assignment needs to wait until Zind is defined.
% datparams assignment needs to wait until tdelay and data are defined.
% If psc == TRUE (varying spot size), sysparams will be updated later,
%   after tdelay is defined.

%% Tweak LCTE
LCTE(1,3) = 0.2; % prototypical 200 MW
LCTE(2,4) = 2e6;
matparams{1} = LCTE;

%% produce model tdelay and/or offset arrays, plus r_pump array (for psc),
% and pack up the rest of the cell params.
TDTR_MAIN_MODEL_vH2; 
%% Define independent variable
xvar = logspace(-1,3,30); % range of xvar parameter values in LCTE (SI) units
Xij = 0;
Xij = [1,4]; % index of LCTE matrix of xvar parameter.
SYS_i = 0; % Scalar. [1,2,3,4] indicates xvar is [f, r_pump, r_probe, tau_rep], respectively.
           % Zero values in Xij or SYS_i indicate that xvar is not in LCTE
           % or a system parameter, respectively.


%% Which sensitivities to calculate?
LCTE_sens_consider = zeros(4,length(LCTE(1,:)));
LCTE_sens_consider(1,4) = 1;
LCTE_sens_consider(1,3) = 1;
sys_consider = [0 0 0 0 0 0]; %[f, r_pump, r_probe, tau_rep, w0, phase]

%% Execute sensitivity calculations and plotting.
[S_LCTE,S_sys,xvar,SS_LCTE,SS_sys] = ...
    parametric_senseplot_vH3(Xij,SYS_i, xvar, ...
                   datparams,sysparams, calparams, matparams, Tparams, ...
                   LCTE_sens_consider, sys_consider);
