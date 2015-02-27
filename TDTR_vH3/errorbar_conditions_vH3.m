% Author: Gregory Hohensee, 4/8/2014
% Usage: called within errorbars_vH1 to initialize various (cell)
% matrices for error bar calculation. Specifies which parameters to
% consider in the uncertainty calculations.
%
% INPUTS: LCTE,XguessIJ
% OUTPUTS: LCTE_consider, r_probe_consider, r_pump_consider, phase_consider;
%          LCTE_err, LCTE_err_temp, r_err, degphase;
% INITIALIZES: LCTE_abs_err, r_probe_err, r_pump_err, phase_err

% Revision history: 7/14/2014 - vH2. No changes.
%                   2/17/2014 - vH3. No changes.
%% Initialize values and uncertainties

% Booleans: TRUE if considering uncertainty from these parameters.
LCTE_consider = ones(4,nl);
r_probe_consider=1;
r_pump_consider=1;
phase_consider=1;

% initialize dimensions
nl = length(LCTE(1,:)); % number of layers
nf = length(XguessIJ(:,1)); % number of fit parameters

% Initialize fractional uncertainty in Xguess due to each parameter.
LCTE_err     = zeros(4,nl);
LCTE_err_temp = zeros(4,nl);
LCTE_abs_err = cell(4,nl);

r_probe_err = zeros(1,nf);
r_pump_err  = zeros(1,nf);
phase_err   = zeros(1,nf);


%% define parameters NOT to consider in error analysis (saves time)
LCTE_consider(3,nl)=0; %last layer is semi-infinite
for i = 1:nl
    if LCTE(2,i) == 1e5 && LCTE(3,i) == 1e-9 % reliable signs of an interface layer
        LCTE_consider(2,i) = 0; %thermal interface layer has no capacitance
        LCTE_consider(3,i) = 0; %thermal interface conductance at fixed t/vary lambda
    end
    
    if ~aniso(i), LCTE_consider(4,i) = 0; end % skip eta for isotropic layers
end

for i = 1:nf
    LCTE_consider(XguessIJ(i,2),XguessIJ(i,3)) = 0; %solving for these!
end

%% define percent uncertainty in each layer/parameter.
% Example LCTE_err for a [*, Al, Al/substrate, substrate]
% model of Al/substrate sample, fitting for G and L(substrate),
% assuming an anisotropic substrate with in-plane conductivity
% uncertain to 5 percent.
% LCTE_err_example = [0.05 0.05 0   0;
%                     0.03 0.03 0   0.03;
%                     0.04 0.04 0   0;
%                     0    0    0   0.05];

% Implement default, generic uncertainties for the remaining layers.
for i = 1:nl
    if LCTE_consider(1,i) ~= 0, LCTE_err(1,i) = 0.05; end
    if LCTE_consider(2,i) ~= 0, LCTE_err(2,i) = 0.03; end
    if LCTE_consider(3,i) ~= 0, LCTE_err(3,i) = 0.04; end
    if LCTE_consider(4,i) ~= 0 && aniso(i), LCTE_err(4,i) = 0.05; end
end
r_err=0.05; % 5% uncertainty in pump/probe beam spot sizes
degphase=0.15;  %phase uncertainty in degrees.
                %You can estimate this by starting with delphase
                %from SetPhase_vH.m, and considering the effect
                %of having many data points to average over.
                %(uncertainty is less than RMS noise)
                