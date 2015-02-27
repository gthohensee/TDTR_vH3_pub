function [dTss, dTpp, LCTE] = TDTR_TDEP_vH3(matparams,sysparams,...
                                           Tparams,intscheme,nnodes)
% TDTR_TDEP_stack - handles steady state heating, looks at pump transient 
% heating, and adjusts heat capacities of the layers accordingly.
%
% Syntax:
% [dTss, dTpp, LCTE] = TDTR_TDEP_vH3(matparams,sysparams,...
%                                           Tparams,intscheme,nnodes)
%
% Inputs:
%    sysparams - {tau_rep f r_pump r_probe}
%    matparams - {LCTE aniso BI n_toplayer TCR}
%    Tparams   - {T0, T_LCTE, A_pump, A_probe, absC, perpulse, jabs, jtrans}
%[See INITIALIZE_CELLPARAMS_vH3.m for details on the params inputs.]
%
%    intscheme  - 0 = Legendre-Gauss, 1 = rombint, 2 = Simpson. For the
%                 thermal model numerical integration.
%    nnodes     - For Legendre-Gauss integration, number of nodes. Default
%                 should be at least 35.
%
% Outputs:
%    dTss - steady-state heating (T = T0 + dTss + dTpp)
%    dTpp - per-pulse heating
%    LCTE - vertcat(lambda,C,t,eta); sample thermal properties matrix.
%
% Example: 
%    --
%
% Other m-files required: TDTR_TEMP_vH3.m, SS_Heating_vH3.m
% Subfunctions: none
% MAT-files required: none
%
% See also: TDTR_TDEP_TTM_vH3.m

% Author: Gregory Hohensee
% University of Illinois at Urbana-Champaign
% email: hohense2@illinois.edu
% Website: n/a
% Apr 2014; Last revision: 1-Apr-2014
%                          8-Apr-2014 - harmonized with TTM version
%                          14-July-2014 - vH2. Fixed absC perpulse bug.
%                          17-Feb-2015 - vH3. No change.
%------------- BEGIN CODE --------------
%% initialize cell params
PL = [0 length(sysparams) 0 length(matparams) length(Tparams)];

% unpack matparams {Mcell aniso BI n_toplayer TCR};
if PL(4) < 1, error('Mcell not specified in FIT_TTM.'); 
else LCTE = matparams{1}; end
%if PL(4) < 2, aniso = zeros(size(LCTE(4,:))); else aniso = matparams{2}; end
if PL(4) < 3, BI = 0; else BI = matparams{3}; end
if PL(4) < 4, if BI, error('n_toplayer not specified'); else n_toplayer = 0; end 
   else n_toplayer = matparams{4}; end
%if PL(4) < 5, TCR = 1e-4; else TCR = matparams{5}; end


% unpack pump-probe system parameters
%sysparams = {tau_rep f r_pump r_probe};
if PL(2) < 1, tau_rep = 1/80e6; 
    warning('Defaulting to 80 MHz repetition rate'); else tau_rep = sysparams{1}; end
% if PL(2) < 2
%     warning('Defaulting to 9.8 MHz modulation frequency');
%     f = sysparams{2}; 
% end 
if PL(2) < 3, error('pump spot size not specified.'); else r_pump = sysparams{3}; end
if PL(2) < 4, r_probe = r_pump;
    warning('Defaulting to r_probe = r_pump'); else r_probe = sysparams{4}; end


% check Tparams (TDEP check different from INITIALIZE check).
if PL(5) < 1, return; else T0 = Tparams{1}; end
if PL(5) < 2, return; else T_LCTE = Tparams{2}; end
if PL(5) < 3, return; else A_pump = Tparams{3}; end
if PL(5) < 4, A_probe = 0; warning('Defaulting to A_probe = 0.'); 
   else A_probe = Tparams{4}; end
if PL(5) < 5, absC = [295 0.13]; else absC = Tparams{5}; end
if PL(5) < 6, perpulse = 0; else perpulse = Tparams{6}; end
if PL(5) < 7, jabs = 0; else jabs = Tparams{7}; end
if PL(5) < 8, jtrans = n_toplayer+2; else jtrans = Tparams{8}; end

%% initialize remaining variables
% time-averaged laser power at sample surface, before absorption.
A_tot = A_pump + A_probe;
if A_tot == 0, return; end

dTss = 0; dTpp = 0;

[~,iabs] = min(abs(absC(:,1) - T0));
absC_at_T = absC(iabs,2); % should rename absC as T_absC, then use absC for this variable.

%% What is the transducer thickness, given j(abs) and j(transducer)?
if jabs ~= 0 % if there exists an absorption layer...
    tabs = LCTE(2,jabs) / LCTE(2,jtrans);
    ttrans = tabs*1e-9 + LCTE(3,jtrans); % this won't update with temperature.
else % jabs = 0: assume no absorption layer
    ttrans = LCTE(3,jtrans);
end

%% Initialize dT temperature rise(s)
dTss = SS_Heating_vH3(matparams,r_pump,r_probe,absC_at_T,A_tot,intscheme,nnodes);
if perpulse
    % factor of 2 for the 50% duty cycle of square wave pump modulation at
    % EOM. See David's RSI 2004 paper for per-pulse heating equation?
    dTpp = (2*absC_at_T*A_pump*tau_rep) / ...
           (pi * r_pump(1)^2 * ttrans * LCTE(2,jtrans));
end

% Given: T_LCTE(i,j,m,n) = [T(m) V(n)], Vij => LCTE(i,j)
% Steps:
% 1) get individual V(n,i,j) for T(n) = T0 + dTss + dTpp; n can vary
% 2) assign V(nn,i,j) to all elements LCTE(i,j) for the various nn(i,j).
% 3) Recompute dT until self-consistent with thermal parameters.
dTerr = dTss+dTpp; % initialize as deviation from T0
V = LCTE; % temporary values for LCTE
while (dTerr > 0.01) % Loop until dT converges to within 1%. Doesn't take many iterations.
    T = T0 + dTss + dTpp;
    
    % Update V(i,j) from T_LCTE
    for i = 1:4
        for j = 1:length(LCTE(1,:))
            TV = T_LCTE{i,j}; % [T, V(T)] for LCTE(i,j)
            [~,iV] = min(abs(TV(:,1) - T));
            
            % No fit parameter should have any temperature
            % dependence. As such, the writeLCTE(G) databases should be
            % given the (Xijc) indicator, and return a null
            % value (negative) to prevent TDEP from adjusting the fit.
            if TV(iV,2) < 0
                V(i,j) = LCTE(i,j);
            else
                V(i,j) = TV(iV,2);
            end
            
            % % QUESTION: what about anisotropy depending on Lz(T)? % %
            % T-dep of the eta part of LCTE is assumed to be due to changes
            % in Lx/Lz, NOT just the in-plane component. For isotropic layers,
            % T_LCTE(4,j,:,2) should just be 1.
            
            % So, T_L and T_E should be constructed in the writeLCTE database
            % such that they are self-consistent: T_E = T(Lx/Lz), given T_Lx.
            
            % Therefore, I need make no compensation in eta for changes in 
            % Lz with temperature.
            
            % This applies equally to the N-channel substrate; writeLCTEG
            % should contain the necessary information in T_E.
            
            % errorbar perturbations will still have to track coupled
            % variables: eta and absorption layer.
        end
    end
    
    matparamstemp = matparams;
    matparamstemp{1} = V;
    
    % recompute dT_SS
    dTss_new = SS_Heating_vH3(matparamstemp,r_pump,r_probe,absC_at_T,A_tot,intscheme,nnodes);
    if perpulse
        dTpp_new = (2*absC_at_T*A_pump*tau_rep) / ...
                   (pi * r_pump(1)^2 * ttrans * LCTE(2,jtrans));
    else
        dTpp_new = 0;
    end
    dT_new = dTss_new + dTpp_new;
    dT = dTss + dTpp;
    dTerr = (dT - dT_new) / (dT_new); % fractional inconsistency in dT
    
    dTss = dTss_new; % update dTss
    dTpp = dTpp_new; % update dTpp
end
% now V is the self-consistent LCTE, so update LCTE
% Assign final results for thermal properties at temperature T = T0 + dT_SS
LCTE = V;
end
%------------- END CODE --------------