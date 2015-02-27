function dTss = SS_Heating_vH3(matparams,r_pump,r_probe,absC,A_tot,intscheme,nnodes)
%SS_Heating_vH3 - Calculates the steady state heating of an anisotropic 
%multilayer stack. Can handle bidirectional heat flow.
%
% Syntax:  SS_Heating_vH3(matparams,r_pump,r_probe,absC,A_tot,intscheme,nnodes)
%
% Inputs:
%    matparams - {LCTE aniso BI n_toplayer TCR}
%    LCTE    - vertcat(lambda,C,t,eta);
%       lambda  - VECTOR of thermal conductivities (W/mK) (layer 1 is the topmost layer)
%       C       - VECTOR of specific heats (J/m3-K)
%       t       - VECTOR of layer thicknesses (last layer is alway treated semi-inf, but
%                 you still have to enter something).  (m)
%       eta     - VECTOR of anisotropic ratio (kx/ky), use ones(length(lambda)) for
%                 isotropic
%[See INITIALIZE_CELLPARAMS_vH3.m for more details.]
%    r_pump  - pump 1/e2 radius (m)
%    r_probe - probe 1/e2 radius (m)
%    A_tot   - Total laser power (mW) hitting the sample inside of a duty
%              cycle of the ~200Hz probe modulation. This is AFTER
%              correcting the powermeter reading for power meter
%              calibration (1.08x) and the objective's transmission
%              coefficient. See the analyze template for details.
%    absC    - fraction of laser power absorbed by the top layer
%
% Outputs:
%    dTss - steady state temperature rise due to laser heating.
%
% Example: 
%    --
%
% Other m-files required: TDTR_TEMP_vH3.m
% Subfunctions: none
% MAT-files required: none
%
% See also: analyze_yymmdd_template.m

% Author: Gregory Hohensee / built from Joseph Feser's SS_Heating.m
% U. of Illinois Urbana-Champaign
% email: hohense2@illinois.edu
% Website: n/a
% Revision history: 12-Sep-2012 - J. Feser's TDTR_V4 published
%                   26-Mar-2014 - SS_Heating_vH.m
%                   8-Apr-2014 - SS_Heating_vH1.m, harmonized
%                   14-July-2014 - vH2. No change.
%                   15-Nov-2014 - r_pump can now be an array vs. time
%                                 delay.
%                   17-Feb-2015 - vH3. No change.
%------------- BEGIN CODE --------------
%% unpack matparams
PL = [0 0 0 length(matparams) 0];

% unpack matparams {LCTE aniso BI n_toplayer TCR};
if PL(4) < 1, error('LCTE not specified in SS_Heating.'); 
else LCTE = matparams{1}; end
%if PL(4) < 2, aniso = zeros(size(LCTE(4,:))); else aniso = matparams{2}; end
if PL(4) < 3, BI = 0; else BI = matparams{3}; end
if PL(4) < 4, if BI, error('n_toplayer not specified'); else n_toplayer = 0; end 
   else n_toplayer = matparams{4}; end
%if PL(4) < 5, TCR = 1e-4; else TCR = matparams{5}; end

%% assign remaining variables
dTss = 0;

%%
f=0; %laser Modulation frequency, Hz
A_abs = absC*A_tot; %absorbed power
r = sqrt(r_pump(1)^2+r_probe^2);

kmin=1/(10000*r); %smallest wavevector (can't use exactly zero, 
                  %because of a numerical issue -> 0/0 = NaN
                  %when calculating SS heating with f = 0.
                  
if intscheme == 2 % Simpson integration
    % RW preference
    kmax = 2/sqrt(r_pump(1)^2+r_probe^2);
else % go with JPF / prior factors.
    kmax = 1.5/sqrt(r_pump(1)^2+r_probe^2);
end
% RW uses cap2 = 2 instead of 1.5 for his REFL/TEMP scripts,
% which use Simpson integration exclusively. 2x catches higher k-modes,
% is more accurate for fixed node density over k-space.
% For TTM, k(RW) = 4*pi^2*k(JPF), where k(JPF) is in use in TEMP_V4.

switch intscheme
    case 0 % Legendre-Gauss
        [kvect,weights]=lgwt_V4(nnodes,kmin,kmax); %computes weights and node locations...
    case 1 % rombint
        % do nothing %
    case 2 % Simpson (uses nnodes)
        % do nothing %
    otherwise
        error('integration scheme not properly specified. See analyze template.');
end

if BI
    totlayers = length(LCTE(1,:));
    % Split the material parameters into "up" and "down" sections,
    % relative to the heater and thermometer at interface of n_toplayer
    % and n_toplayer+1.
    for j = n_toplayer:-1:1
        LCTEu(:,n_toplayer+1-j) = LCTE(:,j);
    end
    LCTEd = LCTE(:,n_toplayer+1:totlayers);
    
    switch intscheme
        case 0 % Legendre-Gauss
            % calculate the heat flow up, layer n_toplayer up to layer 1
            Iu = TDTR_TEMP_vH3(kvect,f,LCTEu,r_pump,r_probe,A_abs);
            [Nk,Nf,Nt] = size(Iu); % get 3D matrix dimensions (for variable r_pump)
            dTu = reshape(weights' * Iu(:,:), 1, Nf, Nt); %   (for variable r_pump)
            %dTu = weights'*Iu;
            
            % calculate the heat flow down, layer n_toplayer+1 down to bottom
            Id = TDTR_TEMP_vH3(kvect,f,LCTEd,r_pump,r_probe,A_abs);
            [Nk,Nf,Nt] = size(Id); % get 3D matrix dimensions (for variable r_pump)
            dTd = reshape(weights' * Id(:,:), 1, Nf, Nt);   % (for variable r_pump)
            %dTd = weights'*Id;
        case 1 % rombint
            dTu=rombint_VV3(@(kvect) TDTR_TEMP_vH3(kvect,f,LCTEu,r_pump,r_probe,A_abs),kmin,kmax,1);
            dTd=rombint_VV3(@(kvect) TDTR_TEMP_vH3(kvect,f,LCTEd,r_pump,r_probe,A_abs),kmin,kmax,1);
        case 2 % Simpson integration
            kvect = linspace(kmin,kmax,nnodes)';
            Iu = TDTR_TEMP_vH3(kvect,f,LCTEu,r_pump,r_probe,A_abs);
            Id = TDTR_TEMP_vH3(kvect,f,LCTEd,r_pump,r_probe,A_abs);
            dTu=SimpsonInt(kvect,Iu);
            dTd=SimpsonInt(kvect,Id);
    end
    % make the parallel sum of the temperatures for up and down
    dTss = dTu .* dTd ./ (dTu + dTd);
    dTss = mean(dTss);
else % BI = 0    
    switch intscheme
        case 0 % Legendre-Gauss
            I = TDTR_TEMP_vH3(kvect,f,LCTE,r_pump,r_probe,A_abs);
            [Nk,Nf,Nt] = size(I); % get 3D matrix dimensions (for variable r_pump)
            dTss = reshape(weights' * I(:,:), 1, Nf, Nt);   % (for variable r_pump)
            dTss = mean(dTss); % just pull out average dTss across pump spot sizes.
            %dTss = weights'*I;
        case 1 % rombint
            dTss=rombint_VV3(@(kvect) TDTR_TEMP_vH3(kvect,f,LCTE,r_pump,r_probe,A_abs),kmin,kmax,1);
        case 2 % Simpson
            kvect = linspace(kmin,kmax,nnodes)';
            I = TDTR_TEMP_vH3(kvect,f,LCTE,r_pump,r_probe,A_abs);
            dTss=SimpsonInt(kvect,I);
    end
end
end
%----------------- END CODE --------------------