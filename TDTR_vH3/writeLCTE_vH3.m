function [LCTE, T_LCTE] = writeLCTE_vH3(stack, thickness,Xijc,P,T0,refdir,rho_Al)
%writeLCTE - output thermal parameters for a given sample and conditions.
% Database that constructs a Lambda, C, t, eta (LCTE) matrix for thermal 
% modeling given the "stack" of strings representing the layer materials,
% and the thicknesses of each layer.
%
% WARNING: this function is designed as a database for thermal parameters
% at specified (P,T). Some of the values listed for more exotic or less
% well-characterized materials are ARBITRARY STARTING POINTS FOR FITTING
% TDTR DATA.
%
% Self-consistent temperature dependence incorporating steady-state heating
% is handled using writeT_LCTE and the FIT and MANUALFIT functions.
% 
% Syntax:  writeLCTE_vH3(stack, thickness,Xijc,P,T0,refdir,rho_Al)
%
% Inputs:
%    stack     - single-row cell array with strings labeling the model
%                layers of a multilayer sample for TDTR.
%    thickness - thickness of each model layer, in nm.
%    T0        - Nominal (measured) temperature. Set T0 = -1 to ignore
%                temperature dependence and just use room temperature.
%    Xijc      - If T0 is not -1, then Xijc will be necessary for defining
%                the T-independent free parameters in T_LCTE.
%                See analyze scripts for Xijc structure.
%    P         - Pressure in GPa for this particular measurement.
%    refdir - directory path for T-dependent material reference data
%    rho_Al    - [OPTIONAL] aluminum (transducer) measured resistivity.
%                Tends to vary depending on the sputtering chamber used.
%
% Outputs:
%    LCTE - 4xN matrix for sample with N model layers. First row is thermal
%          conductivity, second row is volumetric heat capacity, third row
%          is thickness, fourth is anisotropy eta. All units are SI units.
%          "eta" is Lx/Lz, in-plane divided by cross-plane thermal 
%          conductivity for each model layer.
%    T_LCTE - 4xN cell array representing the T-dependence [T, V(T)] of
%             each element of LCTE.
%
%
% Other m-files required: P-dependence and T-dependence functions.
% Subfunctions: none
% MAT-files required: none
%
% See also: writeLCTEg_vH3.m, analyze_yymmdd_template.m

% Author: Gregory T. Hohensee
% University of Illinois Urbana-Champaign
% email: hohense2@illinois.edu
% Website: n/a
% April 2014; Last revision: 3-April-2014
%                            18-April-2014 - removed a refdir hardcode
%                            11-May-2014 - added Pt
%                            11-Jun-2014 - updated nargin < 6,7 to better
%                                          catch stupid inputs
%                            14-Jul-2014 - vH2: new entries, new
%                                          pressure scaling functions.
%                            Nov-2014    - Dallas DOS for metallic silicon
%                            19-Jan-2015 - Fixed silicone oil heat
%                                          capacity.
%                            1-Feb-2015  - rearranged inputs, made
%                                          all after Xijc optional.
%                            17-Feb-2015 - vH3.
%------------- BEGIN CODE --------------
L0 = 2.445e-8; % constant, Wiedermann-Franz law.
RT = 295; % room temperature (K) default

%% check inputs
if nargin < 3
    error('writeLCTE requires more inputs.');
end

% The user may or may not have specified some or all of these.
if ~exist('Xijc','var') || isempty('Xijc'), Xijc = 0; end
if ~exist('P','var') || isempty('P'), P = 0; end
if ~exist('T0','var') || isempty('T0'), T0 = -1; end
if ~exist('refdir','var') || isempty('refdir'), refdir = ''; end
if ~exist('rho_Al','var') || isempty('rho_Al'), rho_Al = -1; LAl = 150; end
% For Al, sputter 3 gives 150-160 W/m-K. Sputter 2 is similar. Sputter 1
% gives ~180 W/m-K.

if rho_Al < 0 % "ignore rho_Al" default -- gives ~180 W/m-K L(Al).
    rho_Al = 3.90e-8;
end
LAl = L0*RT/rho_Al; % default: 295 K aluminum conductivity given rho_Al
   
if length(Xijc(1,:)) == 2 % user gave Xij instead of Xijc
    Xijc = horzcat(Xijc,ones(length(Xijc(:,1)),1)); % assume all parameters are in LCTE.
end


%% Count number of model layers (discarding filler) and assign thicknesses.
j = length(stack);

for i = 1:length(stack)
    if isempty(stack{i})
        j = j - 1;
    end
end % now j is the actual number of layers for this sample.

% initialize LCTE arrays
lambda = zeros(1,j); 
C = zeros(1,j);
t = zeros(1,j);
eta = ones(1,j);
LCTE = vertcat(lambda,C,t,eta);
T_LCTE = cell(size(LCTE));

%% Go through stack, assign thermal parameters.
for ii = 1:j
    if thickness(ii) == 0, break; end % stop if thickness is invalid
    
    clear lower; % in case the workspace contains a "lower" variable.
    s = lower(stack{ii}); % input layer, convert to lowercase
    
    % If previous layer was an absorption layer, skip next (transducer)
    % layer because it has already been populated.
    if ii > 1
        if lower(stack{ii-1}) == '*', continue; end 
    end
    
    t(ii) = thickness(ii)*1e-9; % assign thickness
    
    if ~isempty(strfind(s,'/')) % handles generic interface layers
        lambda(ii) = 0.1; % 100 MW/m^2-K default
        C(ii) = 0.1e6;    % 0.1 J/cm^3-K over 1 nm (negligible) default
        
    else % Refer to database for thermal parameters.
        switch s % all properties are in SI units
%% TRANSDUCER ABSORPTION LAYERS % % % % % % % % % % % % % % % % % % % % %
            case '*' % absorption layer cases
                ss = lower(stack{ii+1});
                
                switch ss
                    case 'al'
                        abslayer = 40; % nm, absorption layer
                        lambda(ii+1) = LAl; % W/m-K
                        C(ii+1) = 2.43e6;
                        
                        if P ~= -1 % Pressure dependence AT ROOM TEMPERATURE
                            % Scale heat capacity by change in Debye heat
                            % capacity with pressure, according to the Al
                            % pressure-volume Equation-Of-State (EOS).
                            [vP,~,TDebye] = Al_EOS(P, 1e-6);
                            [vP0,~,TDebye0] = Al_EOS(0, 1e-6);
                            N0 = 1e29; % arbitrary placeholder
                            DC = DebyeCv(N0 * vP,  TDebye,  295);
                            DC0 =DebyeCv(N0 * vP0, TDebye0, 295);
                            DC = DC * (2.43e6 / DC0); % normalize N0 to ambient C_Al
                            C(ii+1) = DC;
                        end
                        
                        if T0 ~= -1 % Temperature dependence
                            
                            TL_Al_WF_file = 'Al_WF_interp.txt';
                            TC_Al_file = 'interp_C_Al.txt';
                            
                            [C_T0 L_T0 TC_Al TL_Al] = ...
                                    TDTR_TDEP_general(TC_Al_file, TL_Al_WF_file, T0, refdir);
                                
                            if rho_Al ~= -1 % user specified rho_Al
                                [~, TL_Al] = generate_lambdaAl_WF(rho_Al);
                            end
                            
                            lambda(ii+1) = L_T0;
                            C(ii+1) = C_T0;
                            T_LCTE{1,ii+1} = TL_Al;
                            T_LCTE{2,ii+1} = TC_Al;
                        end
                       
                    case 'pt' % platinum from sputter 3
                        abslayer = 10; % I don't know.
                        lambda(ii+1) = 46; % estimated from 4ptprobe. Ideal is 71.6 W/m-K.
                        C(ii+1) = 2.85e6; % from wolfram-alpha and wiki
                        if P ~= -1 % room temperature assumed.
                            [vP,~,TD] = Pt_EOS(P,1e-6);
                            TD0 = 230; % from wiki
                            rho0 = 21450; % ambient density, Pt, kg/m^3
                            amu = 3.23934612e-25; % atomic mass Pt in kg
                            C0 = DebyeCv(rho0/amu, TD0, 295);
                            c0 = DebyeCv(rho0*vP/amu, TD, 295) / C0;
                            C(ii+1) = c0 * C(ii+1);
                        end
                    case 'aupd' % 95% Au, 5% Pd, as measured by RBS.
                        abslayer = 70; % I don't know, but it's large.
                        lambda(ii+1) = 79; % W/m-K default, typical of sputter 3.
                        C(ii+1) = 2.46e6;
                        if P ~= -1 % room temperature assumed.
                          [rho,~,TD] = Au_EOS(P,1e-6);
                           rho0 = 19300; % ambient density, Au
                           amu = 3.2707062e-25; % atomic mass Au in kg
                           C0 = DebyeCv(rho0/amu, 170, 295);
                           c0 = DebyeCv(rho/amu, TD, 295) / C0;
                           C(ii+1) = c0 * C(ii+1);
                        end
                    case 'au' % gold
                        abslayer = 5; % I don't know, but it's large.
                        lambda(ii+1) = 310; % Jensen et al., "BNL 10200-R, Revised August 1980",
                                            % http://www.bnl.gov/magnets/staff/gupta/cryogenic-data-handbook/Section7.pdf
                        C(ii+1) = 2.46e6;
                        if P ~= -1 % room temperature assumed.
                          [rho,~,TD] = Au_EOS(P,1e-6);
                           rho0 = 19300; % ambient density, Au
                           amu = 3.2707062e-25; % atomic mass Au in kg
                           C0 = DebyeCv(rho0/amu, 170, 295);
                           c0 = DebyeCv(rho/amu, TD, 295) / C0;
                           C(ii+1) = c0 * C(ii+1);
                        end
                    case 'ta' 
                        abslayer = 10; % ??
                        lambda(ii+1) = 57.5; % need to make all these a function of measured resistivity
                        C(ii+1) = 2.30e6;
                    case 'nb' % same as nb_sp3 
                        abslayer = 35; % 35 is best-fit from first range of test samples.
                        lambda(ii+1) = 29; % fits Nb thermalization; 4pt probe gives ~30.
                        C(ii+1) = 2.27e6; % internetting; molar volume 10.84 cm^3/mol, 24.6 J/mol-K
                        if P ~= -1 % Pressure dependence AT ROOM TEMPERATURE
                            % Scale heat capacity by change in Debye heat
                            % capacity with pressure, according to the Al
                            % pressure-volume Equation-Of-State (EOS).
                            [vP,~,TDebye] = Nb_EOS(P, 1e-6);
                            [vP0,~,TDebye0] = Nb_EOS(0, 1e-6);
                            N0 = 1e29; % arbitrary placeholder
                            DC = DebyeCv(N0 * vP,  TDebye,  295);
                            DC0 =DebyeCv(N0 * vP0, TDebye0, 295);
                            DC = DC * (2.27e6 / DC0); % normalize N0 to ambient C_Al
                            C(ii+1) = DC;
                        end
                    case 'nb_sp3' 
                        abslayer = 35; % 35 is best-fit from first range of test samples.
                        lambda(ii+1) = 29; % fits Nb thermalization; 4pt probe gives ~30.
                        C(ii+1) = 2.27e6; % internetting; molar volume 10.84 cm^3/mol, 24.6 J/mol-K
                    case 'nb_sp1' 
                        abslayer = 35; % 35 is best-fit from first range of test samples.
                        lambda(ii+1) = 37; % 37 from 4ptprobe on Sputter-1 Nb film; 53.7 bulk value
                        C(ii+1) = 2.27e6; % internetting; molar volume 10.84 cm^3/mol, 24.6 J/mol-K
                    case 'pb'
                        abslayer = 100; % ??
                        lambda(ii+1) = 32.8;
                        C(ii+1) = 1.47e6;
                        if P ~= -1
                            [rho,~,TD] = Pb_EOS(P,1e-6);
                            rho0 = 11340; % ambient density, kg/m^3, Pb
                            amu = 3.441e-25;
                            C0 = DebyeCv((rho0/amu), 105, 295); % reference
                            c0 = DebyeCv((rho/amu), TD, 295) / C0; % C scaling
                            C(ii+1) = c0 * C(ii+1);
                        end
                    case 'msi' % metallic silicon, P > 13 GPa or so.
                        abslayer = 40; % ??
                        lambda(ii+1) = 200; % roughly
                        C(ii+1) = 1.64e6; % 0 GPa
                        % Next five lines were used in analysis up to November 13th
                        % 2014. They are incorrect in the hex and hcp phases of Si.
                        %if P ~= -1
                        %    [PSi,~,~,~,~,CSi] = empirical_PXK_Si_v2(P); % use new C(Si).
                        %    [~,index] = min(abs(PSi - P));
                        %    C(ii) = CSi(index)*1e6;
                        %end

                        % Added on November 13th, 2014.
                        % See CSi figure with Dallas Trinkle's DOS heat capacities.
                        % These functions are linear interpolations of that.
                        if P ~= -1
                            if P < 11 % diamond cubic
                                y = @(x) 0.0077 * x + 1.65;
                            end
                            if P >= 11 && P < 15 % beta-tin phase
                                y = @(x) 0.2299 * (x-11) + 1.7341;
                            end
                            if P >= 15 && P <= 37 % hexagonal
                                y = @(x) 0.0062 * (x-15) + 2.6766;
                            end
                            if P > 37 && P < 42 % intermediate
                                y = @(x) 0.0896 * (x-37) + 2.8198;
                            end
                            if P >= 42 % hcp (ignore higher pressure phases)
                                y = @(x) 3.1781 + 0.006 * (x-42); % This is a wild guess.
                                % pressure scaling is probably comparable to and
                                % weaker than in primitive hexagonal phase, right?
                            end

                            C(ii+1) = y(P) * 1e6;
                            
                        end
                    otherwise
                        disp('ERROR:writeLCT; absorption layer for what material?\n')
                        abslayer = 1; % minimal default
                end
                
                t_abs = min(abslayer, thickness(ii+1)-1);
                t(ii+1) = (thickness(ii+1) - t_abs)*1e-9;
                
                % set absorption layer properties
                lambda(ii) = lambda(ii+1) * t_abs;
                C(ii) = C(ii+1) * t_abs;
                t(ii) = 1e-9;
                %if T0 ~= -1
                %    T_LCTE{1,ii} = horzcat(TL_Al(:,1),TL_Al(:,2)*t_abs);
                %    T_LCTE{2,ii} = horzcat(TC_Al(:,1),TC_Al(:,2)*t_abs);
                %end
                
%% METALS % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            case 'al'
                lambda(ii) = LAl;
                C(ii) = 2.43e6;
                
                if P ~= -1 % Pressure dependence
                    [vP,~,TDebye] = Al_EOS(P, 1e-6);
                    [vP0,~,TDebye0] = Al_EOS(0, 1e-6);
                    N0 = 1e29; % placeholder
                    DC = DebyeCv(N0 * vP,  TDebye,  295);
                    DC0 =DebyeCv(N0 * vP0, TDebye0, 295);
                    DC = DC * (2.43e6 / DC0); % normalize N0 to ambient C_Al
                    C(ii) = DC;
                end
                
                if T0 ~= -1 % Temperature dependence
                    TL_Al_WF_file = 'Al_WF_interp.txt';
                    TC_Al_file = 'interp_C_Al.txt';

                    if rho_Al ~= -1 % user specified rho_Al
                        [~, TL_Al] = generate_lambdaAl_WF(rho_Al);
                        [C_T0 L_T0 TC_Al ~] = ...
                            TDTR_TDEP_general(TC_Al_file, TL_Al_WF_file, T0, refdir);
                    else
                        [C_T0 L_T0 TC_Al TL_Al] = ...
                            TDTR_TDEP_general(TC_Al_file, TL_Al_WF_file, T0, refdir);
                    end

                    lambda(ii) = L_T0;
                    C(ii) = C_T0;
                    T_LCTE{1,ii} = TL_Al;
                    T_LCTE{2,ii} = TC_Al;
                end
            case 'pt' % from sputter 3
                lambda(ii) = 46; % estimated from 4ptprobe. Ideal is 71.6 W/m-K.
                C(ii) = 2.85e6; % from wolfram-alpha and wiki
                if P ~= -1 % room temperature assumed.
                  [vP,~,TD] = Pt_EOS(P,1e-6);
                  TD0 = 230; % from wiki
                  rho0 = 21450; % ambient density, Pt, kg/m^3
                  amu = 3.23934612e-25; % atomic mass Pt in kg
                  C0 = DebyeCv(rho0/amu, TD0, 295);
                  c0 = DebyeCv(rho0*vP/amu, TD, 295) / C0;
                  C(ii) = c0 * C(ii);
                end              
            case 'aupd'
                lambda(ii) = 79; % W/mK default
                C(ii) = 2.46e6;
                if P ~= -1 % room temperature assumed.
                   [rho,~,TD] = Au_EOS(P,1e-6);
                   rho0 = 19300; % ambient density, Au
                   amu = 3.2707062e-25;
                   C0 = DebyeCv(rho0/amu, 170, 295);
                   c0 = DebyeCv(rho/amu, TD, 295) / C0;
                   C(ii) = c0 * C(ii);
                end
            case 'au'
                lambda(ii) = 310; % Jensen et al., "BNL 10200-R, Revised August 1980)".
                C(ii) = 2.46e6;
                if P ~= -1 % room temperature assumed.
                   [rho,~,TD] = Au_EOS(P,1e-6);
                   rho0 = 19300; % ambient density, Au, kg/m^3
                   amu = 3.2707062e-25;
                   C0 = DebyeCv(rho0/amu, 170, 295);
                   c0 = DebyeCv(rho/amu, TD, 295) / C0;
                   C(ii) = c0 * C(ii);
                end
            case 'ta'
                lambda(ii) = 57.5; % W/mK default
                C(ii) = 2.30e6;
            case 'pb'
                lambda(ii) = 32.8; %32.8 for lowP Pb %24.5 for highP Pb dataset; %35.5;
                C(ii) = 1.47e6;
                if P ~= -1
                    [rho,~,TD] = Pb_EOS(P,1e-6);
                    rho0 = 11340; % ambient density, kg/m^3, Pb
                    amu = 3.441e-25; % atomic mass Pb
                    C0 = DebyeCv((rho0/amu), 105, 295); % reference
                    c0 = DebyeCv((rho/amu), TD, 295) / C0; % C scaling
                    C(ii) = c0 * C(ii);
                end
            case 'nb' % sputter 3
                lambda(ii) = 29; % to be verified by 4ptprobe
                C(ii) = 2.27e6; % internetting; molar volume 10.84 cm^3/mol, 24.6 J/mol-K  
            case 'nb_sp3' % sputter 3
                lambda(ii) = 29; % to be verified by 4ptprobe
                C(ii) = 2.27e6; % internetting; molar volume 10.84 cm^3/mol, 24.6 J/mol-K  
            case 'nb_sp1' 
                lambda(ii) = 37; % from 4ptprobe on Sputter-1 Nb film; 53.7 bulk value
                C(ii) = 2.27e6;
            case 'ti'
                lambda(ii) = 21.9; % wiki
                C(ii) = 2.355e6; % internetting from molar volume
%% OTHER MATERIALS % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            case 'msi'
                lambda(ii) = 200;
                C(ii) = 1.64e6;
                
                % Next five lines were used in analysis up to November 13th
                % 2014. They are incorrect in the hex and hcp phases of Si.
                %if P ~= -1
                %    [PSi,~,~,~,~,CSi] = empirical_PXK_Si_v2(P); % use new C(Si).
                %    [~,index] = min(abs(PSi - P));
                %    C(ii) = CSi(index)*1e6;
                %end
                
                % Added on November 13th, 2014.
                % See CSi figure with Dallas Trinkle's DOS heat capacities.
                % These functions are linear interpolations of that.
                if P ~= -1
                    if P < 11 % diamond cubic
                        y = @(x) 0.0077 * x + 1.65;
                    end
                    if P >= 11 && P < 15 % beta-tin phase
                        y = @(x) 0.2299 * (x-11) + 1.7341;
                    end
                    if P >= 15 && P <= 37 % hexagonal
                        y = @(x) 0.0062 * (x-15) + 2.6766;
                    end
                    if P > 37 && P < 42 % intermediate
                        y = @(x) 0.0896 * (x-37) + 2.8198;
                    end
                    if P >= 42 % hcp (ignore higher pressure phases)
                        y = @(x) 3.1781 + 0.006 * (x-42); % This is a wild guess.
                        % pressure scaling is probably comparable to and
                        % weaker than in primitive hexagonal phase, right?
                    end
                    
                    C(ii) = y(P) * 1e6;
                end
                
            case 'air'
                lambda(ii) = 0.01; % arbitrary low number
                C(ii) = 2e6; % arbitrary ordinary number
            case 'sio2'
                lambda(ii) = 1.30;
                C(ii) = 1.62e6;
            case 'mica'
                lambda(ii) = 0.46;
                C(ii) = 1.62e6; % ??    
            case 'la2cuo4' 
                lambda(ii) = 15; % guess
                C(ii) = 2.6e6; % see Denisov et al., Physics Solid State 2013
                % see email to self on Sept 11 2013;
                % see also latest C_La2CuO4 figure (09/20/2013).
            case 'si99ge01'
                lambda(ii) = 35;
                C(ii) = 1.64e6;
                if P ~= -1
                    [PSi,~,~,~,~,CSi] = empirical_PXK_Si_v2(P); % use new C(Si).
                    [~,index] = min(abs(PSi - P));
                    C(ii) = CSi(index)*1e6;
                end
            case 'si'
                lambda(ii) = 146;
                C(ii) = 1.64e6;
                if P ~= -1
                    [PSi,~,~,~,~,CSi] = empirical_PXK_Si_v2(P); % use new C(Si).
                    [~,index] = min(abs(PSi - P));
                    C(ii) = CSi(index)*1e6;
                end
            case 'sige'
                lambda(ii) = 35;
                C(ii) = 1.64e6;
                if P ~= -1
                    [PSi,~,~,~,~,CSi] = empirical_PXK_Si_v2(P); % use new C(Si).
                    [~,index] = min(abs(PSi - P));
                    C(ii) = CSi(index)*1e6;
                end
            case 'diamondanvil'
                lambda(ii) = 500;
                C(ii) = 1.81e6;
                if P ~= -1
                    [rho0,~,~,~,TD0] = Diamond_EOS(0,1e-6);
                    [rho,~,~,~,TD] = Diamond_EOS(P,1e-6);
                    amu = 3.441e-25; % arbitrary number
                    C0 = DebyeCv((rho0/amu), TD0, 295); % reference
                    c0 = DebyeCv((rho/amu), TD, 295) / C0; % C scaling
                    C(ii) = c0 * C(ii);
                end
            case 'diamond'
                lambda(ii) = 2400;
                C(ii) = 1.81e6;
                if P ~= -1
                    [rho0,~,~,~,TD0] = Diamond_EOS(0,1e-6);
                    [rho,~,~,~,TD] = Diamond_EOS(P,1e-6);
                    amu = 3.441e-25; % arbitrary number
                    C0 = DebyeCv((rho0/amu), TD0, 295); % reference
                    c0 = DebyeCv((rho/amu), TD, 295) / C0; % C scaling
                    C(ii) = c0 * C(ii);
                end
            case '2a'
                lambda(ii) = 2400;
                C(ii) = 1.81e6;
                if P ~= -1
                    [rho0,~,~,~,TD0] = Diamond_EOS(0,1e-6);
                    [rho,~,~,~,TD] = Diamond_EOS(P,1e-6);
                    amu = 3.441e-25; % arbitrary number
                    C0 = DebyeCv((rho0/amu), TD0, 295); % reference
                    c0 = DebyeCv((rho/amu), TD, 295) / C0; % C scaling
                    C(ii) = c0 * C(ii);
                end
            case '1a'
                lambda(ii) = 700;
                C(ii) = 1.81e6;
                if P ~= -1
                    [rho0,~,~,~,TD0] = Diamond_EOS(0,1e-6);
                    [rho,~,~,~,TD] = Diamond_EOS(P,1e-6);
                    amu = 3.441e-25; % arbitrary number
                    C0 = DebyeCv((rho0/amu), TD0, 295); % reference
                    c0 = DebyeCv((rho/amu), TD, 295) / C0; % C scaling
                    C(ii) = c0 * C(ii);
                end
            case 'sic'
                lambda(ii) = 400; % estimate
                C(ii) = 3.21e6;
            case 'ni' % polycrystalline nickel
                lambda(ii) = 78;
                C(ii) = 3.95e6;
            case 'pmma'
                lambda(ii) = 0.18;
                C(ii) = 2e6;
            case 'pdms'
                lambda(ii) = 0.15;
                C(ii) = 1.42e6;
            case 'ar' % Argon
                % K. V. Tretiakov, and S. Scandolo, J. Chem. Phys. 121, 11177 (2004).
                T = 293; % enforce room temperature
                %log_Lambda = 2.6 - 1.31*log(T) + 1.29*log(P);
                lambda(ii) = exp(6.05 - 1.31*log(T) + 1.29*log(P)); % Tretiakov and Scandolo.
                %lambda(ii) = exp(7.61 - 1.4*log(T) + 1.24*log(P)); % Goncharov et al. 2012
                % For some reason, the quoted log_Lambda equation from
                % Tretiakov gives me the wrong answer. The slope is right,
                % but the intersect (2.6) seems wrong. Digitizing Fig 6b
                % from the paper and fitting the intersect to 6.05 seems
                % to work OK.
                
                if P ~= -1 % from Wen-Pin thesis.  
                    % see H. Shimizu, H. Tashiro, T. Kume, and S. Sasaki, Phys. Rev. Lett. 86, 4568 (2001).
                    % based on Wen-Pin Hsieh's thesis page 50 result.
                    C_2GPa = 1.36e6;
                    C(ii) = C_2GPa * Ar_EOS(P) / Ar_EOS(2); 
                    % C(ii) = 3 * n_Ar * kB; % classical limit 3*N*k_Boltzmann
                else
                    C(ii) = 1.2e6; % default value; no meaning, since Ar is a gas at ambient P.
                end
            case 'soil_1cst'
                % Legacy of metal-diamond DAC fits before Jun 14th 2014:
                     %lambda(ii) = 0.095; C(ii) = 1.2e6;
                     % if P ~= -1
                        % from Wen-Pin PRB 2011 paper. All other parameters effectively
                        % constant or irrelevant with pressure, in silicone/AuPd/D case.
                     %    lambda(ii) = (24/38)*(0.15 + 0.11 * sqrt(P)); % 24/38 is the ratio of 1 cSt to 30,000 cSt lambda; see polydimethylsiloxanes1.pdf
                     % end
                     
                % New model, June 14th 2014:
                %C(ii) = 1.227e6; % 1.5 J/g-K x 0.818 g/cc.
                % New new model, January 19th, 2015; the previous was a
                % misreading of specific vs. volumetric heat capacity in
                % the Sundqvist paper. Oops.
                C(ii) = 1.5e6; % 1.8337 J/g-K x 0.818 g/cc;
                % 1.8337 J/g-K chosen to force 1.5e6 J/m^3-K.
                if P ~= -1
                    % See discussion in D_Soil_BrillCalc script.
                    % uses rho0 = 0.818e3 kg/m^3, B0 = 8.15 GPa, B0p = 7.5,
                    % Vinet EOS. See Soil1cSt_EOS as well.
                    [rhoS,KS] = Soil1cSt_EOS(P);
                    C(ii) = (1.8337*rhoS*1e-3)*1e6; % J/m^3-K
                    
                    % lambda:
                    C11_Soil = (3/2) * KS;
                    lambda0 = 0.124; % from sandberg 1982, 1 mm^2/s sample
                    
                    % L0 = 0.0116 coefficient in W/m-K x (cc/g)^1/6 x
                    % GPa^-1/2 units...
                    % July 7th 2014 update: L0 = 0.0186 
                    % so that initial value is ~0.20 W/m-K
                    lambda(ii) = 0.0186*rhoS.^(1/6) .* C11_Soil.^(1/2);
                    % This model L increases by factor of 6 up to 50 GPa.
                    % So, effusivity LC increases by x12 from 0-50 GPa.
                end
            case 'soil_30000cst'
                C(ii) = 1.46e6;
                lambda(ii) = 0.15;
                if P ~= -1
                % from Wen-Pin PRB 2011 paper. All other parameters effectively
                % constant or irrelevant with pressure, in silicone/AuPd/D case.
                    lambda(ii) = 1*(0.15 + 0.11 * sqrt(P)); % see polydimethylsiloxanes1.pdf
                end
            case 'sapphire'
                lambda(ii) = 35;
                C(ii) = 3.10e6; % 3.10 matches David's notes
            case 'indium'
                lambda(ii) = 81;
                C(ii) = 1.68e6;
                
            case 'mgo'
                lambda(ii) = 53.3;
                C(ii) = 3.373e6;
            case 'mgfeo' % ferropericlase Fp
                lambda(ii) = 7; % starting point, for ~10-15% Fe.
                C(ii) = 3.373e6;
            case 'mgzno' % Zinc-doped MgO
                lambda(ii) = 14; % arbitrary starting point
                C(ii) = 3.373e6;
                
            case 'sp86' % Ca8La6Cu24O41 spin ladder
                lambda(ii) = 30;
                eta(ii) = 0.1; % arbitrary starting point
                C(ii) = 2.887e6; % from spinladder lattice Debye heat capacity (T_Debye = 500 K)
            case 'sp95' % Ca9La5Cu24O41 spin ladder
                lambda(ii) = 35;
                C(ii) = 2.887e6; % from spinladder lattice Debye heat capacity (see my paper)
            case 'srcuo' % Sr14Cu24O41 spin ladder
                lambda(ii) = 15;
                eta(ii) = 0.1; % arbitrary starting point
                C(ii) = 2.887e6; % from spinladder lattice Debye heat capacity (T_Debye = 500 K)
            
            %%%%%%%%%%% Jun Liu's Sr14Cu24O41 sample %%%%%%%%%%%
            case 'srcuo' % Ca9La5Cu24O41 spin ladder
                lambda(ii) = 35;
                C(ii) = 2.887e6; % from spinladder lattice Debye heat capacity (see my paper)
                eta(ii)=0.1;
                 if T0 ~= -1 % Temperature dependence

                        refdir='E:\research\research work at UIUC\software\TTM model\spin_ladder_2_one channel\refdir';
                        TC_file = 'interp_C_Ca9La5Cu24O41_lattice.txt';
                        TL2_file = 'Interp_scale_1.00_ka_Shi.txt';
                        TL3_file = 'Interp_scale_1.00_kc_Shi.txt';   % need to change later to the total K

                        [C_T0 bb TC_Al cc] = ...
                                TDTR_TDEP_general(TC_file, TC_file, T0, refdir);     

                        [L2_T0 L2_T0 TL2_Al TL2_Al] = ...
                                TDTR_TDEP_general(TL2_file, TL2_file, T0, refdir); 


                        [L3_T0 L3_T0 TL3_Al TL3_Al] = ...
                                TDTR_TDEP_general(TL3_file, TL3_file, T0, refdir);  


                        C(ii) = C_T0;

                        T_LCTE{2,ii} = TC_Al;

                        for mm=1:length(TL2_Al)
                           TL2_Al(mm,2)=TL2_Al(mm,2)/L3_T0; 
                        end

                        T_LCTE{4,ii} = TL2_Al;     % deal with anisotropy here
                 end
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             
            otherwise
                error('Unable to parse stack of layer IDs, check analyze script or writeLCTE.m\n')
        end % end case structure
    end
end % end iteration over layers
LCTE = vertcat(lambda,C,t,eta); % assemble LCTE matrix.

% Handles T-independent and fit parameters in T_LCTE
for i = 1:4
    for j = 1:length(LCTE(1,:))
        % Assign null values to fit parameters in T_LCTE, to prevent TDEP from
        % "updating" a fitting parameter when handling self-consistent heating.
        if Xijc ~= 0
            for x = 1:length(Xijc(:,1))
                if Xijc(x,1) == i && Xijc(x,2) == j && Xijc(x,3) == 1
                    T_LCTE{i,j} = [RT -1];
                end
            end
        end
        
        % Assign T-independent default values to all remaining voids in
        % T_LCTE
        if isempty(T_LCTE{i,j})
            T_LCTE{i,j} = [RT LCTE(i,j)];
        end
    end
end

end
%------------- END CODE --------------