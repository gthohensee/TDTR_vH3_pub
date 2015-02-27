%analyze_yymmdd_template - Template for analyzing typical TDTR data.
% After preparing the data manually or with a script like 
% process_yymmdd_template.m, use this template to efficiently handle
% thermal model fits to a series of TDTR data files.
%
% Basic use of this template:
% (1) "Save As" this MATLAB script (file) with a new name, e.g.
% analyze_140322_[dataseries_ID], such as analyze_140322_AlDAC.
%
% (2) Define the directory and enter the date in the yymmdd variable
%
% (3) Sequentially edit and press "Evaluate Cell and Advance" button, 
%     as indicated by the comments and common sense.
%
% This script can serve as your notebook for how you analyzed this 
% particular data series, your system parameters, and whatever else.
%
% It outputs the fit results to a file and calls scripts to print figures 
% representing the data, fit, and fit parameters for each data set.
%
% Note that this script expects data files to be processed by a
% process_yymmdd_template.m script. Specifically, it expects to find the
% processed data as 5-6 column text files with a "_shifted.txt" suffix
% in a pre-existing folder titled /yymmdd_edit/.
%
% Other m-files required: the TDTR_vH package, specifically TDTR_MAIN_vH.
% Subfunctions: none
% MAT-files required: none
%
% See also: process_yymmdd_template.m

% Author: Gregory T. Hohensee
% University of Illinois Urbana-Champaign
% email: hohense2@illinois.edu
% Website: n/a
% Revision history: 
%       5-August-2014: vH3
%       6-August-2014: MANUALFIT can now fit V(out),
%                      thermal model now accepts pump spot size
%                      as function of delay time, IFF you use the
%                      default Legendre-Gauss integration option
%                                           (see intscheme variable)
%       17-Feb-2015:   vH3, harmonized with various vH3 updates.
%------------- BEGIN CODE --------------
%% Define the directory and filenames
directory = '/Users/gregoryhohensee/Documents/MATLAB/TDTR_vH3'; % pwd is the current directory
yymmdd = '140803'; % e.g., 140227 = February 27th, 2014.
datafolder = strcat('/',yymmdd,'_edit/'); % folder in directory where processed data is kept
savefolder = strcat('/',yymmdd,'_edit/'); % folder in directory where analyzed data is saved
datadir = strcat(directory, datafolder);
savedir = strcat(directory, savefolder);

tagline = 'Al'; % a shared prefix for your series of TDTR data files.

% get filenames in alphanumerical order
filenames = dir(strcat(datadir,tagline,'*'));
filenames.name
nfiles = length(filenames); % number of data files
%% General conditions of modeling
% The regular user may want to look at the init_default scripts in this
% MATLAB package for a condensed way of setting system parameters that
% stay the same from measurement to measurement.

P0 = -1;        % GPa initial. P = -1 assumes atmospheric pressure.

BI = 0;         % boolean: TRUE if using bidirectional heat flow model
n_toplayer = 0; % # of model layers above the point where heat is deposited.

Voutlinfit = 0; % Did you linearize the V(out) for any of the data?
if Voutlinfit
    Voutfitparams = dlmread(strcat(directory, datafolder, yymmdd, ...
                            '_Voutfitparams_', tagline, '.txt'));
end

%% Common sysparams - TDTR system settings
f=9.8e6; % pump laser modulation frequency, Hz
tau_rep=1/80e6; % laser repetition period, s.
                % The repetition frequency is 80 MHz for TDTR-1, 
                % approximately 74.828 MHz for TDTR-2.
TCR=1e-4; % coefficient of thermal reflectance; 1e-4 is roughly TCR(Al) at 785 nm.
          % TCR is only relevant for calculating steady state or per-pulse
          % heating. It affects the amplitude of the model V(out) and
          % V(in), but those aren't absolute calculations in any case.
pm = 1.08; % optical power is 1.08x the reading of the model 835 power meter at TDTR-1.
%pm = 1.00; % scaling for readout of TDTR-2 power meter?

%% calparams - calculation parameters
Zdelay = 100; % ps starting time. In automatic fitting, goodness-of-fit
             % is calculated between Zdelay and tdelay_max. In V(in)
             % fitting, Zdelay indicates the time at which to
             % normalize the V(in) signal, since the thermal model
             % can only produce a normalized V(in).
             % Use Zdelay > 100 ps in case the thermal signal before 
             % Zdelay ps is obscured by heat spreading in the
             % transducer or acoustics.

sigfit = 0; % 1 if fitting by V(in), 2 if fitting by V(out), otherwise fitting by ratio

intscheme = 0; % integration scheme: 0 = Legendre-Gauss, 1 = Rombint,
               % 2 = Simpson integration.
nnodes = 35;  % number of nodes for Legendre-Gauss integration;
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
%% Thickness and resistivity information
% Al/SiO2 - 28.9 ps - 95.7 nm
% Al/Si - 20.9 ps - 70 nm
hAlSiO2 = 95.7;
hAlSi = 70;

rhoAl = 3.90e-8; % Al electrical resistivity; unimportant if not using Al.
                 % 3.90e-8 ==> 186 W/m-K, typical of sputter 1.

%% Arrays for dataset-dependent information
% Thickness information for the multilayer sample. Each row is a different
% sample, the Nth column is the Nth layer for each sample. For a set of
% samples with differing numbers of layers, enter 0 in the rightmost 
% columns of the less layered samples.
% EXAMPLE:
% thickness = [1 hAu 1 100e3;
%     1 hAu 1 100e3;
%     1 hAu 1 100e3;
%     1 hAu 1 100e3;
%     1 hAu 1 100e3;
%     1 hAu 1 100e3];
thickness = [1 hAlSiO2 1 500 500e3;
             1 hAlSiO2 1 500 500e3;
             1 hAlSi   1 500e3 0];

% the "stack": layer IDs for assigning thermal parameters to each model layer.
% Each row is a different sample, the Nth column is the Nth layer for
% each sample. For differing total numbers of layers, enter '' in the
% rightmost columns of the less layered samples. See and add to the
% writeLCTE database function for the string IDs of various materials.
% EXAMPLE:
% stack = {'*', 'Au', 'Ti-MoS2/sapphire',    'sapphire';
%         '*', 'Au', 'Ti-MoS2/sapphire',    'sapphire';
%         '*', 'Au', 'Ti-MoS2/sapphire',    'sapphire';
%         '*', 'Au', 'Ti-MoS2/sapphire',    'sapphire';
%         '*', 'Au', 'Ti-MoS2/sapphire',    'sapphire';
%         '*', 'Au', 'Ti-MoS2/sapphire',    'sapphire'};
stack = {'*', 'Al', '/', 'SiO2','Si';
         '*', 'Al', '/', 'SiO2','Si';
         '*', 'Al', '/', 'Si',''};

% each column of aniso is TRUE if corresponding layer is anisotropic
aniso = zeros(1,length(thickness(1,:))); 

% A 1D array of data sets requires 1D arrays of measurement parameters.
f_list = 9.8 * ones(1,nfiles); % pump modulation frequency
r_list = 5   * ones(1,nfiles); % objective lens ID: 5x, 10x, 20x, ...
jabs_list   = 1 * ones(1,nfiles); % column index of absorption layer
jtrans_list = 2 * ones(1,nfiles); % column index of transducer layer
measured_pump = 20e-3 * ones(1,nfiles);  % raw powermeter reading, W
measured_probe = 10e-3 * ones(1,nfiles); % raw powermeter reading, W

% Xij(m,1:3) = [i j]: assigns LCTE(i,j) as mth fit parameter.
%  Fitting more than two parameters is best done in an exploratory context,
% since TDTR data generally specifies at most two parameters uniquely.
%
%  To fit different parameters for different samples, edit Xij
% and re-run the for loop below for different filename indices.
clear Xij;
Xij(1,1:2) = [1,3]; % [1,3] = thermal conductivity of layer #3
Xij(2,1:2) = [1,4]; % [1,4] = thermal conductivity of layer #4
%Xij(3,1:2) = [3,2]; % [3,2] = thickness of layer #2
%Xij(4,1:2) = [3,4]; % [3,4] = thickness of layer #4

%% Initialize output; careful not to erase your fit results!
output = zeros(nfiles,5);
nf = length(Xij(:,1)); % number of fit parameters
Xguess = zeros(1,nf);
XguessIJ = zeros(nf,3);
%% Call prev_output
caltag = '';

switch sigfit
    case 1
        caltag = strcat(caltag,'_vinfit');
    case 2
        caltag = strcat(caltag,'_voutfit');
    otherwise 
        caltag = strcat(caltag,'_rfit');
end

if exist('frac','var'), caltag = strcat(caltag,'_psc'); end

if manualfit, caltag = strcat(caltag,'_man');
else caltag = strcat(caltag,'_auto'); end

prev_output = dlmread(strcat(datadir, yymmdd, '_solutions_',tagline,caltag,'.txt'));
%% Perform fitting.
for ii = 1:nfiles
   fname = filenames(ii).name;
   datain = strcat(datadir, fname);
   sprintf(datain)
   
   f = f_list(ii) * 1e6;     %laser Modulation frequency, Hz
   jabs = jabs_list(ii);     % index of absorption layer
   jtrans = jtrans_list(ii); % index of transducer layer
   
   % Spot sizes for TDTR-1:
   % According to Xiaojia, as of 1/23/2014:
   % 11.3 um for 5X, 6.1 um for 10X, and 2.7 um for 20X
   % According to me, as of April 9th 2014:
%2x - 28 um
%5x - 11.7 um
%10x - 6.1 um
%20x - 2.65 um
   % Spot sizes for TDTR-2:
   % As of March 24 2014, re: Rich Wilson:
   % 13.3 um, 6.7 um, 3.4 um, and 1.4 um  for 5x, 10x, 20x, and 50x.
   switch r_list(ii) % set pump 1/e^2 intensity focused radius, in meters.
        case 2, r_pump = 2*11.08e-6; tc = 0.87; % "tc": transmission coefficient for this objective
        case 5, r_pump = 11.08e-6; tc = 0.90;
        case 10, r_pump = 6.1e-6; tc = 0.80;
        case 20, r_pump = 2.65e-6; tc = 0.70;
        %case 50, r_pump = 1.4e-6; tc = 0.70; % tc = ?? for 50x
        %case 100, r_pump = ?e-6; tc = 0.??;
   end 
   r_probe = r_pump;         % probe 1/e^2 radius, meters.
   
   % Laser powers that reach a typical in-air sample.
   A_pump = measured_pump(ii)*pm*tc;     
   A_probe = measured_probe(ii)*2*pm*tc; % 2x for optical chopper

   %Construct the thermal material properties matrix.
   %refdir is the path to a folder containing all (interpolated, as
   %necessary) reference information or data for building the T-dependent
   %arrays of thermal parameters in T_LCTE.
   refdir = '';
   
   [LCTE,T_LCTE] = writeLCTE_vH3(stack(ii,:), thickness(ii,:),Xijc,P0,T0,refdir,rho_Al);
   % could also have written: LCTE = writeLCTE_vH3(stack,thickness).
   fprintf('transducer thickness is %f nm\n', thickness(ii,2));
   
   %Which variable(s) are you fitting for?
   for j = 1:nf
       Xguess(j) = LCTE(Xij(j,1),Xij(j,2));
       XguessIJ(j,:) = [Xguess(j), Xij(j,1), Xij(j,2)];
   end
   
   
   % Set initial guess for first two fit parameters; useful for P- or T- 
   % series or re-analyzing data. Uncomment the next three lines and set 
   % the appropriate indices for prev_output and LCTE.
   XguessIJ(1:2,1) = prev_output(ii,2:3);
   LCTE(1,3:4) = prev_output(ii,2:3);
   
   
   % V(out) linear fit functionality is BETA, use at own risk.
   if Voutlinfit
       m1 = Voutfitparams(ii,1); % y = m1*x + m2
       m2 = Voutfitparams(ii,2);
       fitOK = Voutfitparams(ii,3); % boolean, confirms use of this fit
   end
   
   % Convenient re-statement of parameters controlling the fit
   sigfit = 1;
   manualfit = 1;
   psc = 1; % TRUE if pump spot size changes over delay time.
   frac = 0.19; % Fractional change in pump spot size over delay time.
   Vout_offset = 0.1149; % Constant V(out) false "signal" seen
                         % for 5 mW probe, with 5x objective,
                         % while pump beam is blocked. MAIN subtracts
                         % this for you; it's important for high ratio,
                         % small V(out) measurements.
   Vin_offset = 0; % can do the same thing for V(in).
   
   % Call MAIN to begin fit procedure.
   TDTR_MAIN_vH3;  % Now that all system and material parameters are
                  % initialized, TDTR_MAIN_vH begins the thermal modeling,
                  % fitting, printing of figures, saving files, and
                  % everything else for this particular data file.
                  
   % write output row indicating fit results for this sample.
   % Most of these variables were generated in TDTR_MAIN_vH3.
   output(ii,1:length(Xsol)+4) = [r_pump(1)*1e6,Xsol',Z,tdelay_data(Zind)*1e12,frac];
end
%% Write fit results to a text file in the processed data folder.
caltag = '';

switch sigfit
    case 1
        caltag = strcat(caltag,'_vinfit');
    case 2
        caltag = strcat(caltag,'_voutfit');
    otherwise 
        caltag = strcat(caltag,'_rfit');
end

if exist('frac','var'), caltag = strcat(caltag,'_psc'); end

if manualfit, caltag = strcat(caltag,'_man');
else caltag = strcat(caltag,'_auto'); end

dlmwrite(strcat(savedir, yymmdd,'_solutions_',tagline,caltag,'.txt'),output);

%------------- END CODE --------------