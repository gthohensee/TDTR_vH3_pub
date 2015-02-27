%analyze_yymmdd_template - Template for analyzing DAC / TDTR data.
% After preparing the data, use this template to efficiently handle
% thermal model fits to a series of TDTR data files from DAC samples.
%
% Basic use of this template:
%    See analyze_yymmdd_template.m
%
% Other m-files required: the TDTR_vH package.
% Subfunctions: none
% MAT-files required: none
%
% See also: process_yymmdd_template.m

% Author: Gregory T. Hohensee
% University of Illinois Urbana-Champaign
% email: hohense2@illinois.edu
% Website: n/a
% April 2014; Last revision: 2-April-2014
%                            17-Feb-2015: updated from vH1 to vH3. Not
%                                         tested.
%------------- BEGIN CODE --------------
%% Define the directory and filenames
directory = pwd; % pwd is the current directory
yymmdd = '140227'; % e.g., 140227 = February 27th, 2014.
datafolder = strcat('/',yymmdd,'_edit/'); % folder in directory where processed data is kept
savefolder = strcat('/',yymmdd,'_edit/'); % folder in directory where analyzed data is saved
datadir = strcat(directory, datafolder);

tagline = 'AuPdD'; % a shared prefix for your series of TDTR data files.

% get filenames in alphanumerical order
filenames = dir(strcat(datadir,tagline,'*'));
filenames.name
nfiles = length(filenames); % number of data files
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

sigfit = 0; % TRUE if fitting by the V(in) signal instead of ratio.

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

%% Thickness, resistivity, offset information
h0AlD = 150; % (nm) P = 0

rhoAl = 3.90e-8; % Al electrical resistivity; unimportant if not using Al.
                 % 3.90e-8 ==> 186 W/m-K, typical of sputter 1.
 
%Vout_offset = 0.5; % uV, or uV/mW.

%% Pressure data
RDAl = 0.01*[69426 69890 69542 69542 69742 70589]; % ruby fluorescence R1 peak center
PPDAl = [0 13.1 3.3 3.3 8.9 33.9]; % converted pressure in GPa
PP = PPDAl; % PP is the pressure array for analysis.

%% Arrays for dataset-dependent information

% Pressure-volume equation of state for the transducer.
clear vP; clear TDebye; clear DebyeCvAl; clear hAlD;
for i = 1:length(filenames)
    [vP(i),~,TDebye(i)] = Al_EOS(PP(i));
    N0 = 1e29; % placeholder
    DebyeCvAl(i) = DebyeCv(N0 * vP(i), TDebye(i), 300); % 300 K assumed
end

hAlD = h0AlD' ./ vP'; % assume 1D compression (ONLY IF TRANSDUCER IS ON DIAMOND CULET)

thickness = ones(length(filenames),6);
thickness(:,1) = 100e3; % diamond anvil
thickness(:,2) = 1;    % diamond/Al interface
thickness(:,3) = 1;    % absorption layer
thickness(:,4) = hAlD; % transducer
thickness(:,5) = 1;
thickness(:,6) = 30e3; % silicone oil
           
% the "stack": see default analyze template for details.
stack = {'2A','diamond/Al','*', 'Al', 'Al/Soil','Soil_1cst';
         '2A','diamond/Al','*', 'Al', 'Al/Soil','Soil_1cst';
         '2A','diamond/Al','*', 'Al', 'Al/Soil','Soil_1cst';
         '2A','diamond/Al','*', 'Al', 'Al/Soil','Soil_1cst';
         '2A','diamond/Al','*', 'Al', 'Al/Soil','Soil_1cst';
         '2A','diamond/Al','*', 'Al', 'Al/Soil','Soil_1cst'};

% A 1D array of data sets requires 1D arrays of measurement parameters.
f_list = 9.8 * ones(1,nfiles); % pump modulation frequency
r_list = 5   * ones(1,nfiles); % objective lens ID: 5x, 10x, 20x, ...
jabs_list   = 1 * ones(1,nfiles); % column index of absorption layer
jtrans_list = 2 * ones(1,nfiles); % column index of transducer layer
measured_pump = 20e-3 * ones(1,nfiles);  % raw powermeter reading, W
measured_probe = 10e-3 * ones(1,nfiles); % raw powermeter reading, W

% Xij(m,1:3) = [i j]: assigns LCTE(i,j) as mth fit parameter.
clear Xij;
Xij(1,1:2) = [1,3];
Xij(2,1:2) = [1,4];
%Xij(3,1:2) = [3,2];
%Xij(4,1:2) = [3,4];
%%
output = zeros(nfiles,5);
nf = length(Xij(:,1)); % number of fit parameters
Xguess = zeros(1,nf);
XguessIJ = zeros(nf,3);
%% Initialize output and prev_output; careful not to erase your fit results!
caltag = '';

switch sigfit
    case 1
        caltag = strcat(caltag,'_vinfit');
    case 2
        caltag = strcat(caltag,'_voutfit');
    otherwise 
        caltag = strcat(caltag,'_rfit');
end

if psc == 1, caltag = strcat(caltag,'_psc'); end

if manualfit, caltag = strcat(caltag,'_man');
else caltag = strcat(caltag,'_auto'); end

prev_output = dlmread(strcat(datadir, '150215_solutions_',tagline,'_man.txt'));

%% Perform fitting.
j = length(Xij(:,1)); % number of fit parameters
Xguess = zeros(1,j);
XguessIJ = zeros(j,3);

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
   % Spot sizes for TDTR-2:
   % As of March 24 2014, re: Rich Wilson:
   % 13.3 um, 6.7 um, 3.4 um, and 1.4 um  for 5x, 10x, 20x, and 50x.
   switch r_list(ii) % set pump 1/e^2 intensity focused radius, in meters.
        case 2, r_pump = 2*11.3e-6; tc = 0.87; % "tc": transmission coefficient for this objective
        case 5, r_pump = 11.3e-6; tc = 0.90;
        case 10, r_pump = 6.1e-6; tc = 0.80;
        case 20, r_pump = 2.7e-6; tc = 0.70;
        %case 50, r_pump = ?e-6; tc = 0.??;
        %case 100, r_pump = ?e-6; tc = 0.??;
   end 
   r_probe = r_pump;         % probe 1/e^2 radius, meters.
   
   % Laser powers that reach a typical in-air sample.
   A_pump = measured_pump(ii)*pm*tc;     
   A_probe = measured_probe(ii)*2*pm*tc; % 2x for optical chopper

   %Construct the thermal material properties matrix
   P0 = PP(ii);
   refdir = ''; 
   [LCTE,T_LCTE] = writeLCTE_vH3(stack(ii,:), thickness(ii,:),Xijc,P0,T0,refdir,rho_Al);
   % could also have written: LCTE = writeLCTE_vH3(stack,thickness).
   
   fprintf('Pressure is %f GPa\n', P0);
   fprintf('transducer thickness is %f nm\n', thickness(ii,4));
   
   %Which variable(s) are you fitting for?
   for j = 1:nf
       Xguess(j) = LCTE(Xij(j,1),Xij(j,2));
       XguessIJ(j,:) = [Xguess(j), Xij(j,1), Xij(j,2)];
   end
   
   % Set initial guess
   %XguessIJ(1:2,1) = prev_output(ii,2:3);
   %LCTE(1,5:6) = prev_output(ii,2:3);
   
   % V(out) linear fit functionality is BETA, use at own risk.
   if Voutlinfit
       m1 = Voutfitparams(ii,1); % y = m1*x + m2
       m2 = Voutfitparams(ii,2);
       fitOK = Voutfitparams(ii,3); % boolean, confirms use of this fit
   end
   
   % reiterate important modeling parameters for ease of switching
   Zdelay = 100;
   
   sigfit = 1; % 1 == Vin, 2 == Vout, 0 = Ratio
   manualfit = 1;
   psc = 1; % TRUE if pump spot changes over delay time.
   frac = 0.16; % TDTR-2, 16% pump spot size change over delay time.
   doughnut = 0; % no beam offset
   
   % Offsets measured by blocked pump immediately after regular measurement.
   clear Vin_offset; clear Vout_offset;
   
   % Execute MAIN program.
   TDTR_MAIN_vH3;
                  
   output(ii,1:length(Xsol)+3) = [P0,Xsol,Z,tdelay(Zind)*1e12];
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

if psc == 1, caltag = strcat(caltag,'_psc'); end

if manualfit, caltag = strcat(caltag,'_man');
else caltag = strcat(caltag,'_auto'); end

dlmwrite(strcat(savedir, yymmdd,'_solutions_',tagline,caltag,'.txt'),output);

%% Notes
%------------- END CODE --------------