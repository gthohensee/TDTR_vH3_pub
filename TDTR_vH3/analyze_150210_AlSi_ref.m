
%% Define the directory and filenames
directory = '/Users/gregoryhohensee/Documents/MATLAB/TDTR_vH3'; % pwd is the current directory
yymmdd = '150210'; % e.g., 140227 = February 27th, 2014.
datafolder = strcat('/',yymmdd,'_edit/'); % folder in directory where processed data is kept
savefolder = strcat('/',yymmdd,'_edit/'); % folder in directory where analyzed data is saved
datadir = strcat(directory, datafolder);
savedir = strcat(directory, savefolder);

tagline = 'Al'; % a shared prefix for your series of TDTR data files.

% get filenames in alphanumerical order
filenames = dir(strcat(datadir,tagline,'*'));
filenames.name
nfiles = length(filenames); % number of data files
%% General modeling conditions
%init_default_BI_TDTR1; % for TDTR-1
init_default; % not BI
tau_rep = 1/74.82e6; % TDTR-2

P = -1;
P0 = -1;
doughnut = 0;

%% Thickness and resistivity information
                             
rhoAl = 3.90e-8; % Al electrical resistivity; unimportant if not using Al.
                 % 3.90e-8 ==> 186 W/m-K, typical of sputter 1.
                 % 4.90e-8 is more typical of sputter 3.
%% Arrays for dataset-dependent information
hAl = 27.72 * 6.42/2 + 3; % acoustics

thickness = [1 hAl 1 500e3];

stack = {'*', 'Al', '/', 'Si'};

% each column of aniso is TRUE if corresponding layer is anisotropic
aniso = zeros(1,length(thickness(1,:))); 

% A 1D array of data sets requires 1D arrays of measurement parameters.
f_list = [9.1]; % pump modulation frequency
r_list = [10];
jabs_list   = 1 * ones(1,nfiles); % column index of absorption layer
jtrans_list = jabs_list + 1; % column index of transducer layer
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
        case 2, r_pump = 2*11.7e-6; tc = 0.87; % "tc": transmission coefficient for this objective
        % measured 11.5 - 12.5 micron ud - lr spot size, as of Dec. 4th
        % 2014.
        case 5, r_pump = 12e-6; tc = 0.90; 
        case 10, r_pump = 6.7e-6; tc = 0.80; % TDTR-2
        case 20, r_pump = 3.4e-6; tc = 0.70; % TDTR-2
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
   %XguessIJ(1:2,1) = prev_output(ii,2:3);
   %LCTE(1,3:4) = prev_output(ii,2:3);
   %XguessIJ(1,1) = 0.2;
   %LCTE(1,3) = 0.2;
   
   % V(out) linear fit functionality is BETA, use at own risk.
   if Voutlinfit
       m1 = Voutfitparams(ii,1); % y = m1*x + m2
       m2 = Voutfitparams(ii,2);
       fitOK = Voutfitparams(ii,3); % boolean, confirms use of this fit
   end
   
   sigfit = 2; % 1 == Vin, 2 == Vout, 0 = Ratio
   manualfit = 1;
   
   psc = 1; % TRUE if pump spot changes over delay time.
   frac = 0.16;
   clear Vin_offset; clear Vout_offset;
   Vout_offset = 0.04;
   
   TDTR_MAIN_vH3;  % Now that all system and material parameters are
                  % initialized, TDTR_MAIN_vH begins the thermal modeling,
                  % fitting, printing of figures, saving files, and
                  % everything else for this particular data file.
                  
   % write output row indicating fit results for this sample.
   % Most of these variables were generated in TDTR_MAIN_vH3.
   output(ii,1:length(Xsol)+4) = [r_pump(1)*1e6,Xsol',Z,tdelay(Zind)*1e12,frac];
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