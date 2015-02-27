function [Z,deltaR_model,ratio_model,LCTE,T_adj]=...
    TDTR_FIT_vH3(X, Xij, datparams,sysparams,calparams,matparams,Tparams)
%TDTR_FIT_vH3 - automatically fit thermal model to TDTR ratio data.
%The main program tries to minimize "Z" by optimizing the variable(s) X.
% This version is capable of self-consistent temperature adjustment using
% the Tparams. If you do not wish to do self-consistent T-corrections,
% simply enter the one-element cell {-1} as Tparams.
%
% Syntax:  [Z,deltaR_model,ratio_model,LCTE,T_adj]=...
%    TDTR_FIT_vH3(X, Xij, datparams,sysparams,calparams,matparams,Tparams)
%
% Inputs:
%    X        - Fit parameters. Has to be separate from everything else
%             - for fminsearch functionality.
%    Xij      - Mx2 matrix: The (i,j) indices in LCT of the X fit
%               parameters.
%    datparams - {tdelay ratio_data datadir offset}
%    sysparams - {tau_rep f r_pump r_probe}
%    calparams - {Zind vinfit intscheme nnodes consider_error LCTE_err T0_err}
%    matparams - {LCTE aniso BI n_toplayer TCR doughnut}
%    Tparams   - {T0, T_LCTE, A_pump, A_probe, absC, perpulse, jabs, jtrans}
%[See INITIALIZE_CELLPARAMS_vH3.m for details on the params inputs.]
%
% Outputs:
%    Z           - Goodness-of-fit. In Joe's words: "Typically, this is 
%                  the sum of the squares of the residuals, but you might 
%                  want to weight these by the sensitivity, particularly 
%                  if you don't intend to calculate the errorbars!"
%    deltaR_model - Complex number array. Real part is the model V(in),
%                   imaginary part is the model V(out). Represents change
%                   in reflectance from pump heating.
%    ratio_model  - Ratio signal -V(in)/V(out) from the thermal model.
%    LCTE         - non-fitted material parameters may vary with
%                   self-consistent temperature adjustment from TDEP.
%    T_adj       - T_adj = T0 + dTss + dTpp, the "actual" temperature
%                  in Kelvin adjusted for steady-state (and per-pulse
%                  if perpulse is TRUE) heating.
%
% Example:
%    Xsol = fminsearch(@(X),TDTR_FIT_vH3(X,Xij,...
%                   datparams,sysparams,calparams,matparams,Tparams),...
%                   Xguess,options);
%
% Other m-files required: TDTR_REFL_vH3.m, TDTR_TEMP_vH3.m,
%                         TDTR_TDEP_vH3.m, SS_Heating_vH3.m
% Subfunctions: none
% MAT-files required: none
%
% See also: TDTR_FIT_vH.m, TDTR_MANUALFIT_vH.m, TDTR_FIT_Vin_vH.m

% Author: Gregory Hohensee
% Acknowledgement: built from TDTR_FIT_V4B, my bi-directional tweak to the
% original TDTR_FIT_V4.m, by Joseph P. Feser.
% University of Illinois at Urbana-Champaign
% email: hohense2@illinois.edu
% Website: n/a
% Mar 2014; Last revision: 25-Mar-2014
%                          8-Apr-2014 - harmonized with TDTR_FIT_TTM_vH1.m
%                          14-Apr-2014 - X update can read nm now
%                          14-July-2014 - vH2.
%                          27-Jan-2015 - converted X to nm and J/cm^3-K units,
%                                       so that fminsearch can efficiently
%                                       find Xsol without trying to run to
%                                       6 or 9 significant figures.
%                          27-Jan-2015 - sigfit, beam offset compatibility.
%                          17-Feb-2015 - vH3, no changes.
%------------- BEGIN CODE --------------
%% Check input parameters, assign defaults, errors, warnings as necessary
INITIALIZE_CELLPARAMS_vH3; % unpacks and re-packs cellparams

%% Define the variables to be fit
nrows = length(Xij(:,1));

% convert X units to LCTE SI units
% Xunits = {'W/m-K', 'J/cm^3-K', 'nm', '(Lx/Lz)'};
XtoLCTE = [1 1e6 1e-9 1];

% iterate over fit parameters
for i = 1:nrows, LCTE(Xij(i,1),Xij(i,2)) = X(i) * XtoLCTE(Xij(i,1)); end

%% Perturbations for error bar calculation.
% need to update thermal parameters with perturbations.
if consider_error
    % Note that only one element of LCTE_err should be nonzero at any
    % one call of this function.
    LCTE = LCTE .* (1 + LCTE_err);
    
    if jabs ~= 0 % if there exists an absorption layer...
        % its uncertainty is coupled to that of the transducer layer.
        if sum(LCTE_err(1:2,jabs) == 0) == 0 % double-check that the 
                                             %perturbation is not in what 
                                             %we assume is the absorption 
                                             %layer
            LCTE(1:2,jabs) = LCTE(1:2,jabs) .* (1 + LCTE_err(1:2,jtrans));
        end
    end
    
    % update eta for a perturbation in cross-plane conductivity,
    % if the layer is anisotropic.
    LCTE(4,:) = LCTE(4,:) ./ (1 + aniso.*LCTE_err(1,:));
    
    % measured temperature or pressure uncertainty
    T0 = T0 + T0_err; % perturbation will take effect in TDEP
    %P0 = P0 + P0_err; % perturbation... needs writeLCTE(P).
end

matparams{1} = LCTE; % update matparams after error consideration, before TDEP.

%% Self-consistent steady-state (and optionally per-pulse) heating 
if T0 ~= -1
    % NOTE: TDEP does not consider whether you are modeling
    % beam-offset or conventional TDTR. In either case it
    % will compute the SS and per-pulse heating at zero offset.
    % This is intentional.
    
    % re-use Tparams for the TDEP function, although now Tparams
    % may have some default values in it that were unspecified before.
    Tparams{1} = T0; % update in case of consider_error
    [dTss, dTpp, LCTE] = TDTR_TDEP_vH3(matparams,sysparams,...
                                           Tparams,intscheme,nnodes);
    fprintf('T0 = %0.2f K, dTss = %0.2f, dTpp = %0.2f',T0,dTss,dTpp)
end
matparams{1} = LCTE; % update matparams after TDEP, before fitting.

%% Compute model at zero beam offset (conventional TDTR)
[deltaR_model,ratio_model]=TDTR_REFL_vH3(tdelay,matparams,sysparams,...
                                        A_pump,intscheme,nnodes,0);
switch sigfit
    case 1 % V(in) fit
        if doughnut, error('Code not designed for beam offset fit to V(in).'); end
        % Construct normalized V(in) model and data,
        % relative to its value at Zdelay picoseconds.
        Vin_model = real(deltaR_model);
        Vin_model_Zdelay = Vin_model(Zind) / N;
        nVin_model = Vin_model / Vin_model_Zdelay;

        Vin_data_Zdelay = Vin_data(Zind);
        nVin_data = Vin_data / Vin_data_Zdelay;
    case 2 % V(out) fit
        if doughnut
            Y0 = imag(deltaR_model);  % this from a REFL call with xoffset = 0.
            [~,I] = min(abs(offset));              % get index of r = 0 data point
            Norm = mean(abs(Vout_data(I-1:I+1)));  % take 3pt average of Vout data around r = 0
            nVout_data = abs(Vout_data/Norm);      % remove negative sign, rescale to match DOUGHNUT.
        else
            % Construct normalized V(out) model and data,
            % relative to its mean value near Zdelay picoseconds.
            Vout_model = imag(deltaR_model);
            Vout_model_Zdelay = Vout_model(Zind) / N;
            nVout_model = Vout_model / Vout_model_Zdelay;

            Vout_data_Zdelay = Vout_data(Zind);
            nVout_data = Vout_data / Vout_data_Zdelay;
        end
    otherwise % ratio fit
        % do nothing here
        if doughnut, error('Code not designed for beam offset fit to ratio.'); end
end

%% Do beam offset thermal model.
% input: offset, Nx, REFL params. output: deltaR and ratio vs xvect
if doughnut
    Nx = 10;
    % fyi -- normally no need to model past 3.9x spot size.
    xvect = linspace(0,max(offset),Nx)'; % need to model just one side of the data, by symmetry

    % calculate model from each offset point in xvect.
    [deltaR_modelx, ratio_modelx]=TDTR_REFL_vH3(tdelay,matparams,sysparams,...
                                               A_pump,intscheme,nnodes,xvect);
    % deltaR and ratio may be a matrix (tdelay, offset). Pull out
    % first (only) time delay element, transpose into column.
    Sig_Amp = deltaR_modelx(1,:)'; 
    Sig_ratio = ratio_modelx(1,:)';

    %reconstruct the rest by symmetry
    xvect = [-flipud(xvect(2:end));xvect];
    Sig_Amp = [flipud(Sig_Amp(2:end));Sig_Amp];
    Sig_ratio = [flipud(Sig_ratio(2:end));Sig_ratio];

    Vout_model = interp1(xvect,imag(Sig_Amp),offset,'spline','extrap');
    Vin_model = interp1(xvect,real(Sig_Amp),offset,'spline','extrap');
    deltaR_model = Vin_model + 1i*Vout_model;
    ratio_model = -Vin_model ./ Vout_model;

    nVout_model = abs(Vout_model/Y0);
end
    
%% Uncomment the next few lines to see the non-linear optimization in action!
figure(10)
loglog(tdelay,ratio_data,'ob',tdelay,ratio_model,'r-'); 
xlabel('Time delay (ps)','FontSize',16);

switch sigfit
    case 1 % V(in)
        ylabel('normalized V(in)','FontSize',16);
        axis([1e-10 10e-9 min([nVin_data;nVin_model])/2 max([nVin_data;nVin_model])*2])
    case 2
        ylabel('normalized V(in)','FontSize',16);
        axis([1e-10 10e-9 min([nVout_data;nVin_model])/2 max([nVout_data;nVin_model])*2])
    otherwise % ratio fit
        ylabel('Ratio','FontSize',16); 
        axis([100e-12 10e-9 min([ratio_data;ratio_model])/2 max([ratio_data;ratio_model])*2])
end

set(gca, 'XTick', [1e-10,2e-10,5e-10,1e-9,2e-9,5e-9,1e-8]);
set(gca, 'XTickLabel', [100, 200, 500, 1000, 2000, 5000, 1e4]);
set(gca, 'YTick', [0.5, 1, 2, 5, 10, 20, 50, 100]);
set(gca, 'YTickLabel', [0.5, 1, 2, 5, 10, 20, 50, 100]);
set(gca, 'XMinorTick', 'off');
set(gca, 'FontSize', 16);
pause(0.1) % may be necessary so MATLAB can refresh the figure between iterations

%% Define symbols and units for all possible fit parameters
tag = {'L','C','t','eta'};
units = {'W/m-K', 'J/cm^3-K', 'nm', '(Lx/Lz)'};
LCTEtoX = [1 1e-6 1e9 1]; % from SI to the "units" variable's units.

%% Report current fit parameters
Mat = sprintf('Current material fit parameters: ');

for x = 1:length(X) % iterate through fit parameters
    itag = Xij(x,1);
    Mat = char(Mat,sprintf('%s(%i) = %0.4f %s',...
                  tag{itag},Xij(x,2),...
                  LCTE(Xij(x,1),Xij(x,2))*LCTEtoX(itag),units{itag}));
end
Mat

%% Assign final temperature adjusted for SS and PP heating.
if T0 ~= -1, T_adj = T0 + dTss + dTpp; else T_adj = T0; end

%% goodness-of-fit Z from fractional residuals
    switch sigfit
        % For V(in) and V(out) fits, it may be wiser to weight the 
        % residuals by the amount of time delay away from Zdelay, 
        % where V(in) is pinned to unity.
        case 1
            res=(1-(nVin_model(Zind:length(nVin_model)) ...
                    ./nVin_data(Zind:length(nVin_data))) ).^2;
        case 2
            if doughnut
                res = (nVout_model - nVout_data).^2; % fractional error may overvalue the tails?
            else
                res=(1-(nVout_model(Zind:length(nVout_model)) ...
                        ./nVout_data(Zind:length(nVout_data))) ).^2;
            end
        otherwise % ratiofit
        res=(1-(ratio_model(Zind:length(ratio_model)) ...
                ./ratio_data(Zind:length(ratio_data))) ).^2;
    end 
Z = sum(res)
end