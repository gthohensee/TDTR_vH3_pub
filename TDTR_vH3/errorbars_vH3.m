function [kErr_perc, kErr_abs, ErrSummary_perc, ErrSummary_abs] = ...
    errorbars_vH3(XguessIJ,datparams,sysparams,calparams,matparams,Tparams,options)
%errorbars_vH3 - calculates error bars for the Xijc fit parameters
%based on uncertainties in other thermal/system parameters. Sequentially
%calls TDTR_FIT_vH3.m to calculate how much the fitted Xsol changes
%with a perturbation in another parameter equal to its uncertainty. Sums up
%changes in Xsol in quadrature for all considered perturbations, as
%specified by errorbar_conditions_vH3.m script.
%
% Inputs:
%    Xij      - XguessIJ(:,2:3). For each row, [i j] indicates LCTE(i,j)
%               corresponding to fit parameter Xguess.
%    datparams - {tdelay ratio_data datadir offset}
%    sysparams - {tau_rep f r_pump r_probe}
%    calparams - {Zind rorfit intscheme nnodes consider_error Mcell_err T0_err}
%    matparams - {LCTE aniso BI n_toplayer TCR doughnut}
%    Tparams   - {T0, T_LCTE, A_pump, A_probe, absC, perpulse, jabs, jtrans}
%    options   - fminsearch options, precision tolerances on the fit.
%
% Outputs:
%    kErr_perc - total fractional errorbars from quadrature sum over all
%                    perturbations from uncertainties. Each row corresponds
%                    to the same row in Xguess / the fit parameters.
%    kErr_abs  - Absolute error bars.
%    ErrSummary_perc - individual fractional errors from individual
%                          perturbations in other material / system
%                          parameters.
%    ErrSummary_abs - individual absolute errors in Xsol due to
%                         uncertainty in other parameters.
%
% Other m-files required: TDTR_FIT_vH3.m, TDTR_REFL_vH3.m,
%                         TDTR_TDEP_vH3.m, SS_Heating_vH3.m,
%                         TDTR_TEMP_vH3.m.
% Subfunctions: none
% MAT-files required: none
%
% See also: errorbars_TTM_vH3.m, senseplot_vH3.m

% Author: Gregory Hohensee
% University of Illinois Urbana-Champaign
% email: hohense2@illinois.edu
% Website: n/a
% Revision history: 8-Apr-2014 - converted to function,
%                                harmonized with TTM version.
%                   14-July-2014 - vH2. No changes.
%                   21-Jan-2015 - compatible with sigfit.
%                   28-Jan-2015 - Fixed variable names at ErrSummary stage.
%                                 Compatible with beam offset.
%                   17-Feb-2015 - vH3. No changes.
%------------- BEGIN CODE --------------
%% unpack inputs
INITIALIZE_CELLPARAMS_vH3; % unpacks/cleans the five cellparams

%% Initialize variables and uncertainties
Xij = XguessIJ(:,2:3);
Xguess = XguessIJ(:,1);
errorbar_conditions_vH3; % default uncertainties
                
%% Compute reference fit parameters Xsol
% TolX = 1e-2 reduces strictness of fminsearch
if nargin < 7, options = optimset('TolX',1e-2); end 

fprintf('YErr(n,:) = uncertainty (absolute) in X due to uncertainty in parameter Y(n)\n')
Xsol=fminsearch(@(X) TDTR_FIT_vH3(X,Xij,datparams,sysparams,calparams,...
                                  matparams,Tparams),Xguess,options);
Xguess = Xsol; % these are interchangeable in this script from here onward.
%% Compute change in fit parameters due to independent perturbations
for i = 1:4
    for j=1:nl
        % skip conditions
        if LCTE_consider(i,j) == 0, continue; end
        if j == jabs, continue; end % skip error propagation from absorption layer
                                    % to avoid double-counting; the LCTE_err
                                    % elements for the absorption layer are
                                    % applied simultaneously with the
                                    % LCTE_err elements for the transducer.
                                    % If you have a "1nm transducer", where
                                    % the entire transducer is modeled as
                                    % an absorption layer, you should set
                                    % jabs = 0 in your analyze script.
        
        LCTE_err_temp(i,j) = LCTE_err(i,j);
            
        % couplings between model parameters, e.g. absorption layer and
        % anisotropies.
        %if jabs ~= 0 && j == jabs, LCTE_err_temp(i,jtrans) = LCTE_err(i,jtrans); end
        if jabs ~= 0 && j == jtrans, LCTE_err_temp(i,jabs) = LCTE_err(i,jabs); end
        if i == 1, LCTE_err_temp(4,j) = LCTE_err(i,j); end % Changes in Lz affect eta.
        % Changes in eta do NOT affect Lz. This way, eta represents in-plane L.
        
        % Compute change in X due to uncertainties in LCTE
        matparams{1} = LCTE_err_temp .* LCTE;
        Xsoltemp=fminsearch(@(X) TDTR_FIT_vH3(X,Xij,datparams,sysparams,calparams,...
                              matparams,Tparams),Xguess,options);
        
        LCTE_abs_err{i,j}=abs(Xsoltemp-Xguess); %Errors in X(:) due to variable LCTE(i,j)
        
        LCTE_err_temp = zeros(4,nl); % re-initialize LCTE_err_temp
        sprintf('Done with LCTE(%i,%i) error...',i,j)
    end
end
matparams{1} = LCTE; % re-initialize matparams.

%-------Probe Radius--------------
if r_probe_consider==1
    sysparams{4} =r_probe*(1+r_err);
    Xsoltemp=fminsearch(@(X) TDTR_FIT_vH3(X,Xij,datparams,sysparams,calparams,...
                              matparams,Tparams),Xguess,options);

    for n=1:length(Xguess)
        r_probe_err(1,n)=abs(Xsoltemp(n)-Xguess(n)); %Error in X(n) due to variable r_probe(ii)
    end
    sprintf('Done with r_probe error...')
    r_probe_err
    
    sysparams{4} = r_probe; % return to reference value
end

%-------Pump Radius--------------
if r_pump_consider==1
    sysparams{3} =r_pump*(1+r_err);
    Xsoltemp=fminsearch(@(X) TDTR_FIT_vH3(X,Xij,datparams,sysparams,calparams,...
                              matparams,Tparams),Xguess,options);
    for n=1:length(Xguess)
        r_pump_err(1,n)=abs(Xsoltemp(n)-Xguess(n)); %Error in X(n) due to variable r_pump(ii)
    end
    sprintf('Done with r_pump error...')
    r_pump_err
    
    sysparams{3} = r_pump; % return to reference value
end
%--------Phase Error-------------
if phase_consider==1
    % apply phase shift to data
    radphase=pi/180*degphase;
    Vtemp=(Vin_data+sqrt(-1)*Vout_data)*exp(sqrt(-1)*radphase);
    Vin_phaseshifted=real(Vtemp);
    Vout_phaseshifted=imag(Vtemp);
    ratio_phaseshifted=-Vin_phaseshifted./Vout_phaseshifted;
    
    % update datparams
    if ~doughnut
        switch sigfit 
            case 1
                datparams{2} = Vin_phaseshifted; 
            case 2
                datparams{2} = Vout_phaseshifted;
            otherwise
                datparams{2} = ratio_phaseshifted; 
        end
    else
        datparams{2} = Vout_phaseshifted;
    end
    
    Xsoltemp=fminsearch(@(X) TDTR_FIT_vH3(X,Xij,datparams,sysparams,calparams,...
                              matparams,Tparams),Xguess,options);
    for n=1:length(Xguess)
        phase_err(1,n)=abs(Xsoltemp(n)-Xguess(n));
    end
    sprintf('Done with phase error...')
    phase_err
    
    % return data to reference value
    if ~doughnut
        switch sigfit 
            case 1
                datparams{2} = Vin_data; 
            case 2
                datparams{2} = Vout_data;
            otherwise
                datparams{2} = ratio_data; 
        end
    else
        datparams{2} = Vout_data;
    end
end
%% Assemble error summaries
ErrSummary_abs_LCTE = [];
for i = 1:4
    for j = 1:nl
        % empty cells of LCTE_abs_err{i,j} will have no effect here
        ErrSummary_abs_LCTE = vertcat(ErrSummary_abs_LCTE,LCTE_abs_err{i,j});
    end
end
ErrSummary_abs_sys = vertcat(r_probe_err,r_pump_err,phase_err);

%% Assemble reports for error summaries
ErrSummary_abs = vertcat(ErrSummary_abs_LCTE,...
                         ErrSummary_abs_sys);

repeat_Xsol = ones(length(ErrSummary_abs(:,1)),1)*Xsol; % rows for parameters
ErrSummary_perc = ErrSummary_abs ./ repeat_Xsol         %percent error broken by variable

kErr_perc=sqrt(sum(ErrSummary_perc.^2,1)) %total percent error in each fitted parameter
kErr_abs=kErr_perc.*Xsol                  %total absolute error in each fitted parameter