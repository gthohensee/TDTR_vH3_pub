%INITIALIZE_CELLPARAMS_vH3 - Unpacks cell arrays into variables for the
% FIT, MANUALFIT, errorbar, and senseplot functions. If the cell
% arrays are incomplete (not as many elements as expected), then this
% script assigns defaults and issues warnings or errors as necessary.
%
% Inputs:
%    datparams - {tdelay ratio_data datadir offset}
%    sysparams - {tau_rep f r_pump r_probe}
%    calparams - {Zind sigfit intscheme nnodes consider_error LCTE_err T0_err P0_err}
%    matparams - {LCTE aniso BI n_toplayer TCR doughnut}
%    Tparams   - {T0, T_LCTE, A_pump, A_probe, absC, perpulse, jabs, jtrans}
%
% Outputs
%    tdelay  - COLUMN VECTOR of time delay data
%    ratio_data - COLUMN VECTOR of ratio data; can also be "Vin_data",
%                 representing the V(in) signal.
%    datadir - STRING directory path in which (processed) data, fit
%              results, and sensitivity plots are kept and saved.
%    offset - lateral relative offset of pump and probe beams, for beam
%             offset calculation (see doughnut).
%
%    tau_rep - laser repetition time (sec) = 1/f_rep
%    f       - modulation frequency 
%    r_pump  - pump 1/e2 radius (m)
%    r_probe - probe 1/e2 radius (m)
%
%
%    Zind    - Closest index of tdelay such that Zdelay = tdelay(Zind).
%              Zdelay is user-specified time delay (ps) for shortest time
%              delay from which to evaluate goodness-of-fit.
%    sigfit  - 1 == fitting V(in), 2 == fitting V(out), otherwise ratio.
%    intscheme - Numerical integration algorithm.
%                0 = Legendre-Gauss, 1 = rombint, 2 = Simpson.
%    nnodes   - For Legendre-Gauss and Simpson integration,
%               number of nodes or integration points. Default >= 35.
%    consider_error - TRUE if FIT function is called from an errorbars
%                     script. Supports errorbar calculation functionality.
%    LCTE_err -  cell array of uncertainties for each element in LCTE.
%    T0_err - assigned uncertainty (absolute) in the measured temperature.
%    P0_err - assigned uncertainty (absolute) in the measured pressure.
%
%
%    LCTE    - vertcat(lambda,C,t,eta) for the one-channel overlayers
%    aniso   - VECTOR of booleans for each overlayer: aniso(j) is TRUE if
%              jth layer is allowed to have a variable eta between 0 and 1.
%    BI      - TRUE if bidirectional heat flow is present.
%    n_toplayer - Specifies number of layers above the plane where
%                 heat is deposited in the bidirectional heat flow model.
%    TCR     - temperature coefficient of reflectivity  (1/K)...doesn't affect ratio
%    doughnut - TRUE if beam offset. Compatible with BI.
%
%
%    T0      - nominal measured temperature (K)
%    T_LCTE - Cell array, where T_LCTE{i,j} is a two-column array
%             of [T, V(T)] reference information to determine LCTE(i,j)
%             at the self-consistent temperature T0+dT indicated by
%             steady-state and per-pulse laser heating.
%    A_pump  - pump intensity (W)...doesn't effect ratio
%    A_probe - probe intensity (W)...doesn't effect ratio
%    absC - [T, absC(T)] two-column array for temperature dependence of
%           the absorption coefficient of the transducer at the laser
%           wavelength at the specified temperature T.
%    perpulse - TRUE if considering per-pulse heating.
%    jabs - index j for the absorption layer in LCTE.
%    jtrans - index j for the transducer layer in LCTE. jabs and jtrans
%             are necessary for per-pulse heating calculation.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: The TDTR_FIT_vH2.m, TDTR_MANUALFIT_vH2.m

% Author: Gregory Hohensee
% U. of Illinois Urbana-Champaign
% email: hohense2@illinois.edu
% Website: n/a
% Revision history: 27-Mar-2014 - written
%                   8-Apr-2014 - harmonized with TTM version
%                   14-July-2014 - vH2. No changes.
%                   13-Nov-2014 - Include beam offset (doughnut).
%                   Feb-2015 - fixed beam offset related assignments,
%                              retitled vH3.
%------------- BEGIN CODE --------------
%% Check input parameters, assign defaults, errors, warnings as necessary
PL = [length(datparams) length(sysparams) length(calparams) ...
      length(matparams) length(Tparams)];

% pump-probe system parameters
%sysparams = {tau_rep f r_pump r_probe};
if PL(2) < 1, tau_rep = 1/80e6; 
    warning('Defaulting to 80 MHz repetition rate'); else tau_rep = sysparams{1}; end
if PL(2) < 2, f = 9.8e6;
    warning('Defaulting to 9.8 MHz modulation'); else f = sysparams{2}; end
if PL(2) < 3, error('pump spot size not specified.'); else r_pump = sysparams{3}; end
if PL(2) < 4, r_probe = r_pump;
    warning('Defaulting to r_probe = r_pump'); else r_probe = sysparams{4}; end

% material parameters
%matparams = {LCTE aniso BI n_toplayer TCR doughnut};
if PL(4) < 1, error('LCTE not specified'); else LCTE = matparams{1}; end
if PL(4) < 2, aniso = zeros(size(LCTE(4,:))); else aniso = matparams{2}; end
if PL(4) < 3, BI = 0; else BI = matparams{3}; end
if PL(4) < 4, if BI, error('n_toplayer not specified'); else n_toplayer = 0; end 
   else n_toplayer = matparams{4}; end
if PL(4) < 5, TCR = 1e-4; else TCR = matparams{5}; end
if PL(4) < 6, doughnut = 0; else doughnut = matparams{6}; end

% calculation parameters
%calparams = {Zind sigfit intscheme nnodes consider_error LCTE_err T0_err P0_err}
if PL(3) < 1, [~,Zind] = min(abs(tdelay - 100e-12)); else Zind = calparams{1}; end
if PL(3) < 2, sigfit = 0; else sigfit = calparams{2}; end
if PL(3) < 3, intscheme = 0; else intscheme = calparams{3}; end
if PL(3) < 4, nnodes = 35; else nnodes = calparams{4}; end
if PL(3) < 5, consider_error = 0; else consider_error = calparams{5}; end
if PL(3) < 6, LCTE_err = zeros(size(LCTE)); else LCTE_err = calparams{6}; end
if PL(3) < 7, T0_err = 0; else T0_err = calparams{7}; end
%if pl(3) < 8, P0_err = 0; else P0_err = calparams{8}; end

% data parameters
%datparams = {tdelay ratio_data datadir offset} or {tdelay Vin_data datadir offset}
if PL(1) < 1, error('no time vector'); else tdelay = datparams{1}; end
if PL(1) < 2, error('no data vector'); 
else
    switch sigfit
        case 1, Vin_data = datparams{2};
        case 2, Vout_data = datparams{2};
        otherwise ratio_data = datparams{2};
    end
end
if PL(1) < 3, datadir = pwd; else datadir = datparams{3}; end
if PL(1) < 4
    if doughnut, error('no offset vector');
    else offset = 0; end % no pump/probe offset in regular TDTR
else
    offset = datparams{4};
end

% Parameters for self-consistent temperature w/laser heating.
%Tparams = {T0, T_LCTE, A_pump, A_probe, absC, perpulse, jabs, jtrans};
if PL(5) < 1, T0 = -1;     else T0 = Tparams{1}; end
if PL(5) < 2, T_LCTE = -1; else T_LCTE = Tparams{2}; end
if PL(5) < 3, A_pump = 10e-3;  warning('Defaulting to A_pump = 10mW.'); 
   else A_pump = Tparams{3}; end
if PL(5) < 4, A_probe = 10e-3; warning('Defaulting to A_probe = 10mW.'); 
   else A_probe = Tparams{4}; end
if PL(5) < 5, absC = [295 0.13]; else absC = Tparams{5}; end
if PL(5) < 6, perpulse = 0; else perpulse = Tparams{6}; end
if PL(5) < 7, jabs = 0; warning('Defaulting to jabs = 0 null value.'); 
    else jabs = Tparams{7}; end
if PL(5) < 8, jtrans = n_toplayer+2; warning('Defaulting to jtrans = n_toplayer+2.');
    else jtrans = Tparams{8}; end
    
%% Reconstruct all cellparams with any new default values
% if sigfit, f holds 2 values, not directly usable by REFL
sysparams = {tau_rep f r_pump r_probe};

switch sigfit
    case 1, datparams = {tdelay Vin_data datadir offset};
    case 2, datparams = {tdelay Vout_data datadir offset};
    otherwise datparams = {tdelay ratio_data datadir offset};
end

calparams = {Zind sigfit intscheme nnodes consider_error LCTE_err T0_err};
matparams = {LCTE aniso BI n_toplayer TCR doughnut};
Tparams = {T0, T_LCTE, A_pump, A_probe, absC, perpulse, jabs, jtrans};
   