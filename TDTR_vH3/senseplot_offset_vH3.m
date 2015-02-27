function [S_LCTE,S_sys] = senseplot_offset_vH3(datparams,sysparams, calparams, matparams, Tparams)
%senseplot_offset_vH3 - Calculates sensitivity plots dlogR/dlogX for thermal
%model. The sens_consider booleans can be edited to change which
%sensitivities are calculated. "senseplot_offset" is an expansion of
%"senseplot_vH2" that allows sensitivity calculations for beam
%offset signals, or alternately for the FWHM of the offset signal.
%
% Inputs:
%    datparams - {tdelay ratio_data datadir offset}
%    sysparams - {tau_rep f r_pump r_probe}
%    calparams - {Zind sigfit intscheme nnodes consider_error Mcell_err T0_err P0_err}
%    matparams - {LCTE aniso BI n_toplayer TCR doughnut}
%    Tparams   - {T0, T_LCTE, A_pump, A_probe, absC, perpulse, jabs, jtrans}
%[See INITIALIZE_CELLPARAMS for details.]
%
% Outputs
%    S_LCTE - sensitivities to one-channel overlayer parameters
%    S_sys  - {S_f,S_r_pump,S_r_probe};
%       S_f, S_r_pump, S_r_probe - system parameter sensitivities
%
% Other m-files required: TDTR_REFL_vH3.m
% Subfunctions: none
% MAT-files required: none
%
% See also: errorbars_vH3.m, parametric_senseplot_vH3.m

% Author: Gregory Hohensee
% U. of Illinois Urbana-Champaign
% email: hohense2@illinois.edu
% Website: n/a
% Revision history: 8-Apr-2014 - made into function, comments updated,
%                                harmonized with TTM version.
%                   14-July-2014 - vH2. No changes since June 11th.
%                   24-Nov-2014  - fixed typo in handling eta sensitivity
%                   27-Jan-2015  - beam offset compatibility. S vectors
%                                  are with respect to time or offset.
%                   1-Feb-2015   - retitled vH3, includes FWHM mode.
%                   17-Feb-2015  - updated header comments.
%------------- BEGIN CODE --------------
%%
INITIALIZE_CELLPARAMS_vH3; % unpacks/cleans cellparams (the five inputs)
fprintf('Generating Sensitivity Plot for all Variables...please be patient\n')

if doughnut
    by_FWHM = -1;
    while by_FWHM == -1
        by_FWHM = input('Enter 1 if signal is FWHM of V(out), else 0 if signal is just V(out): ');
        if by_FWHM ~= 1 || by_FWHM ~= 0
            by_FWHM = -1;
            sprintf('Invalid input, try again.\n');
        end
    end
end

%% Generates a reference model based on anticipated parameters
% Normalize V(in) and V(out) reference model, if fitting by these.
% The commented out lines are not used, unless you want to revise
% the sensitivity calculation to use raw data as its reference.
% I find that messy, not intrinsic or helpful, so sensitivities are
% defined as model perturbations relative to the model.

if ~doughnut % conventional TDTR
    [deltaR_model,ratio_model]=TDTR_REFL_vH3(tdelay,matparams,sysparams,A_pump,intscheme,nnodes,0);
    switch sigfit
        case 1 % V(in) fit
            % Construct normalized V(in) model and data,
            % relative to its value at index Zind (from Zdelay picoseconds, or offset).
            Vin_model = real(deltaR_model);
            Vin_model_Zdelay = Vin_model(Zind);
            nVin_model = Vin_model / Vin_model_Zdelay;

            %Vin_data_Zdelay = Vin_data(Zind);
            %nVin_data = Vin_data / Vin_data_Zdelay;
        case 2 % V(out) fit
            % Construct normalized V(out) model and data,
            % relative to its mean value near Zdelay picoseconds.
            Vout_model = imag(deltaR_model);
            Vout_model_Zdelay = mean(Vout_model(Zind:Zind+3));
            nVout_model = Vout_model / Vout_model_Zdelay;

            %Vout_data_Zdelay = mean(Vout_data(Zind:Zind+3));
            %nVout_data = Vout_data / Vout_data_Zdelay;
        otherwise % ratio fit
            % do nothing here
    end
else % beam offset
    if length(tdelay) > 1, error('Senseplot does not handle t-dependent offset calculations. See parametric senseplot instead.'); end
    
    % Thermal model is symmetric around 0 offset, so only use
    % the positive offset side of the offset vector to generate model.
    [~,I] = min(abs(offset));
    xvect = offset(I:length(offset));
    [deltaR_model,ratio_model]=TDTR_REFL_vH3(tdelay,matparams,sysparams,A_pump,intscheme,nnodes,xvect);
    
    % deltaR and ratio may be a matrix (tdelay, offset). Pull out
    % first (only) time delay element, transpose into column.
    Sig_Amp = deltaR_model(1,:)'; 
    Sig_ratio = ratio_model(1,:)';

    %reconstruct the rest by symmetry; no data, so no need to interpolate.
    offset = [-flipud(xvect(2:end));xvect]; % makes sure offset and model line up.
    deltaR_model = [flipud(Sig_Amp(2:end));Sig_Amp];
    ratio_model = [flipud(Sig_ratio(2:end));Sig_ratio];
    Vout_model = imag(deltaR_model);
    Vin_model = real(deltaR_model);
    
    [norm,I] = max(abs(Vout_model)); % get value and index of r = 0 model peak
    nVout_model = sgn(Vout_model(I)) * Vout_model / norm; % normalize to peak, make positive.
    
    if By_FWHM
        [~,J] = min(abs(nVout_model - 0.5)); % gets half-rise index.
        Xguess = 0.4247 * offset(J); % if offset was continuous
                                     % instead of a few points,
                                     % this would be the FWHM.

        if ~exist('options','var'); options = optimset('TolFun',1e-2,'TolX',1e-2); end
        Xsol = fminsearch(@(X) SpotSize_V4(X,nVout_model,offset),Xguess,options);
        FWHM_model = 2 * sqrt(2*log(2)) * Xsol;
    end
end

%% initialize matrices for sensitivities
LCTEtemp = LCTE;

%nt = length(tdelay); 
nl = length(LCTE(1,:));

S_LCTE = cell(4,nl);

deltaR_temp  = cell(4,nl);
ratio_temp   = cell(4,nl);

%% which sensitivities to consider? (saves time)
LCTE_sens_consider = ones(4,nl);

% Skip absorption layer; instead, couple its perturbations into that of
% the transducer layer. If there is no transducer layer, just a "1 nm"
% absorption layer, set jabs = 0 in your analyze script.
if jabs ~= 0, LCTE_sens_consider(:,jabs) = 0; end 

for j = 1:nl
    % Skip thicknesses and heat capacities of all interfaces in LCTE
    if LCTE(2,j) == 1e5 && LCTE(3,j) == 1e-9 % reliable interface indicators
        LCTE_sens_consider(2,j) = 0;
        LCTE_sens_consider(3,j) = 0;
        LCTE_sens_consider(4,j) = 0; % interfaces are not anisotropic.
    end
    
    % skip eta for all isotropic layers
    if ~aniso(j), LCTE_sens_consider(4,j) = 0; end
end
LCTE_sens_consider

%% -----------------Compute sensitivities for LCTE--------------
for i = 1:4
    for j=1:nl
        % skip conditions; includes the aniso variable.
        if LCTE_sens_consider(i,j) == 0, continue; end

        LCTEtemp(i,j) = LCTE(i,j)*1.01;

        % couplings between model parameters, e.g. absorption layers and
        % anisotropies.
        if jabs ~= 0 && j == jabs, LCTEtemp(i,jtrans) = LCTE(i,jtrans)*1.01; end
        if jabs ~= 0 && j == jtrans, LCTEtemp(i,jabs) = LCTE(i,jabs)*1.01; end
        if i == 1 && aniso(j), LCTEtemp(4,j) = LCTE(4,j)/1.01; end 
        % eta and Lz are coupled; need to hold Lx fixed when varying Lz
        
        % Perturbing eta is assumed to be perturbing Lx, not Lz. %

        % Perform sensitivity calculation with LCTEtemp
        matparams{1} = LCTEtemp;
        
        % Thermal model is symmetric around 0 offset, so only use
        % the positive offset side of the offset vector to generate model.
        [~,I] = min(abs(offset));
        xvect = offset(I:length(offset)); % this will be 0 if conventional TDTR
        [deltaR_temp{i,j},ratio_temp{i,j}] = TDTR_REFL_vH3(tdelay,matparams,sysparams,A_pump,intscheme,nnodes,offset);
        
        if ~doughnut
            switch sigfit
                case 1
                    Vin_temp = real(deltaR_temp{i,j}); 
                    norm = Vin_temp(Zind);
                    nVin_temp = Vin_temp / norm;
                    Num=log(nVin_temp)-log(nVin_model);
                case 2
                    Vout_temp = imag(deltaR_temp{i,j}); 
                    norm = Vout_temp(Zind);
                    nVout_temp = Vout_temp / norm;
                    Num=log(nVout_temp)-log(nVout_model);
                otherwise
                    Num=log(ratio_temp{i,j})-log(ratio_model);
            end
        else % beam offset
            % deltaR and ratio may be a matrix (tdelay, offset). Pull out
            % first (only) time delay element, transpose into column.
            temp_deltaR = deltaR_temp{i,j};
            temp_ratio = ratio_temp{i,j};
            
            Sig_Amp = temp_deltaR(1,:)'; 
            Sig_ratio = temp_ratio(1,:)';

            %reconstruct the rest by symmetry; no data, so no need to interpolate.
            offset = [-flipud(xvect(2:end));xvect]; % makes sure offset and model line up.
            temp_deltaR = [flipud(Sig_Amp(2:end));Sig_Amp];
            temp_ratio = [flipud(Sig_ratio(2:end));Sig_ratio];
            
            Vout_temp = imag(temp_deltaR);
            [norm,I] = max(abs(Vout_temp)); % get value and index of r = 0 model peak
            nVout_temp = sgn(Vout_temp(I)) * Vout_temp / norm;
            
            % Alternative mode for depicting beam offset sensitivities:
            % instead of sensitivity to V(out) as a function of offset,
            % condense this information into the scalar change in FWHM
            % of V(out) at specified time delay.
            if By_FWHM && length(offset) > 4
                [~,J] = min(abs(nVout_temp - 0.5)); % gets half-rise index.
                Xguess = 0.4247 * 2*abs(offset(J)); % if offset was continuous
                                             % instead of a few points,
                                             % this would be the FWHM.
                                             
                if ~exist('options','var'); options = optimset('TolFun',1e-2,'TolX',1e-2); end
                Xsol = fminsearch(@(X) SpotSize_V4(X,nVout_temp,offset),Xguess,options);
                FWHM_temp = 2 * sqrt(2*log(2)) * Xsol;
                
                Num = log(FWHM_temp) - log(FWHM_model); % this is scalar
            else
                if length(offset) <= 4
                    by_FWHM = 0;
                    warning('Cannot fit FWHM from so few offset points. Defaulting back to V(out)(x).')
                end
                Num = log(nVout_temp) - log(nVout_model); % this is array of length(offset).
            end
        end % end doughnut
        
        Denom=log(LCTEtemp(i,j))-log(LCTE(i,j));
        S_LCTE{i,j}=Num/Denom;

        LCTEtemp = LCTE;
        matparams{1} = LCTE;
        fprintf('Calculated S_LCTE{%i,%i}...\n',i,j);
    end
end
clear LCTEtemp;

%% -----------------Compute sensitivities for system parameters-----------
sys_consider = [1 2 3];
S_f = []; S_r_pump = []; S_r_probe = [];
for i = sys_consider
    sysbump = 1.01;
    switch i
        case 1, sysparams{2} = f*sysbump;
        case 2, sysparams{3} = r_pump*sysbump;
        case 3, sysparams{4} = r_probe*sysbump;
    end
    
    % Thermal model is symmetric around 0 offset, so only use
    % the positive offset side of the offset vector to generate model.
    [~,I] = min(abs(offset));
    xvect = offset(I:length(offset)); % this will be 0 if conventional TDTR
    [temp_deltaR,temp_ratio]=TDTR_REFL_vH3(tdelay,matparams,sysparams,A_pump,intscheme,nnodes,offset);
    if ~doughnut
        switch sigfit
            case 1
                Vin_temp = real(temp_deltaR); 
                norm = Vin_temp(Zind);
                nVin_temp = Vin_temp / norm;
                Num=log(nVin_temp)-log(nVin_model);
            case 2
                Vout_temp = imag(temp_deltaR); 
                norm = Vout_temp(Zind);
                nVout_temp = Vout_temp / norm;
                Num=log(nVout_temp)-log(nVout_model);
            otherwise
                Num=log(temp_ratio)-log(ratio_model);
        end
    else % beam offset, assumes tdelay has just one element.
        % deltaR and ratio may be a matrix (tdelay, offset). Pull out
        % first (only) time delay element, transpose into column.
        Sig_Amp = temp_deltaR(1,:)'; 
        Sig_ratio = temp_ratio(1,:)';

        %reconstruct the rest by symmetry; no data, so no need to interpolate.
        offset = [-flipud(xvect(2:end));xvect]; % makes sure offset and model line up.
        temp_deltaR = [flipud(Sig_Amp(2:end));Sig_Amp];
        temp_ratio = [flipud(Sig_ratio(2:end));Sig_ratio];
    
        Vout_temp = imag(temp_deltaR);
        [norm,I] = max(abs(Vout_temp)); % get value and index of r = 0 model peak
        nVout_temp = sgn(Vout_temp(I)) * Vout_temp / norm;
        
        if By_FWHM && length(offset) > 4
            [~,J] = min(abs(nVout_temp - 0.5)); % gets half-rise index.
            Xguess = 0.4247 * 2*abs(offset(J)); % if offset was continuous
                                         % instead of a few points,
                                         % this would be the FWHM.

            if ~exist('options','var'); options = optimset('TolFun',1e-2,'TolX',1e-2); end
            Xsol = fminsearch(@(X) SpotSize_V4(X,nVout_temp,offset),Xguess,options);
            FWHM_temp = 2 * sqrt(2*log(2)) * Xsol;

            Num = log(FWHM_temp) - log(FWHM_model); % this is scalar
        else
            if length(offset) <= 4
                by_FWHM = 0;
                warning('Cannot fit FWHM from so few offset points. Defaulting back to V(out)(x).')
            end
            Num = log(nVout_temp) - log(nVout_model); % this is array of length(offset).
        end
        
    end
    Denom=log(sysbump);
    
    switch i
        case 1, S_f=Num/Denom;
        case 2, S_r_pump=Num/Denom;
        case 3, S_r_probe=Num/Denom;
    end
    sysparams = {tau_rep, f, r_pump, r_probe}; % return to reference value
    fprintf('Calculated S_sys #%i...\n',i);
end

%% Plot sensitivities
figure(202)
clf
is_timedelay = ~doughnut || by_FWHM;

if is_timedelay, axes('XScale','log'); end
set(gca,'Box','on');
hold on;
ColorOrder = get(gcf,'DefaultAxesColorOrder');

%% label LCTE sensitivities
LCTElegend = [];
LCTElab = cell(size(LCTE));
LCTEmarker = {'o','*','x','+'};

for i = 1:4
    for j = 1:nl
        switch i
            case 1, LCTElab{i,j} = sprintf('L%i',j);
            case 2, LCTElab{i,j} = sprintf('C%i',j);
            case 3, LCTElab{i,j} = sprintf('t%i',j);
            case 4, LCTElab{i,j} = sprintf('e%i',j);
        end
        
        if isempty(S_LCTE{i,j}), continue;
        else
            if ~doughnut || by_FWHM,
                semilogx(tdelay,S_LCTE{i,j},LCTEmarker{i},'Color',ColorOrder(j,:));
            else
                plot(xvect,S_LCTE{i,j},LCTEmarker{i},'Color',ColorOrder(j,:));
            end
            LCTElegend = [LCTElegend;LCTElab{i,j}];
        end
    end
end

%% plot and label system sensitivities
syslegend = [];
for i = sys_consider
    switch i
        case 1
            if is_timedelay,
                semilogx(tdelay,S_f,'-','Color',ColorOrder(i,:));
            else
                plot(xvect,S_f,'-','Color',ColorOrder(i,:));
            end
            syslegend = [syslegend;'ff'];
        case 2
            if is_timedelay,
                semilogx(tdelay,S_r_pump,'-','Color',ColorOrder(i,:));
            else
                plot(xvect,S_r_pump,'-','Color',ColorOrder(i,:));
            end
            syslegend = [syslegend;'Rp'];
        case 3
            if is_timedelay,
                semilogx(tdelay,S_r_probe,'-','Color',ColorOrder(i,:));
            else
                plot(xvect,S_r_probe,'-','Color',ColorOrder(i,:));
            end
            syslegend = [syslegend;'Rb'];
    end
end

%% Other plot details
figure(202);
set(gca,'FontSize',16);
set(gca, 'TickLength' , [.02 .02]);

legend([LCTElegend;syslegend])
if is_timedelay,
    xlabel('time delay (ps)','FontSize',16);
else
    xlabel('offset (microns)','FontSize',16);
end

if is_timedelay
    if ~doughnut
        switch sigfit
            case 1, ylabel(sprintf('Sensitivity:  dlog[nV(in) @ %0.0f ps]/dlogX',tdelay(Zind)*1e12), 'FontSize',16);
            case 2, ylabel(sprintf('Sensitivity:  dlog[nV(out) @ %0.0f ps]/dlogX',tdelay(Zind)*1e12), 'FontSize',16);
            otherwise, ylabel('Sensitivity:  dlogR/dlogX', 'FontSize',16);
        end
    else % by_FWHM = TRUE
        ylabel(sprintf('Sensitivity:  dlog[V(out) FWHM]/dlogX'), 'FontSize',16);
    end
    
    set(gca, 'XTick', [1e-10,2e-10,5e-10,1e-9,2e-9,5e-9,1e-8]);
    set(gca, 'XTickLabel', [100, 200, 500, 1000, 2000, 5000, 1e4]);
    set(gca, 'XMinorTick', 'off');
    axis([100e-12 10e-9 -2 2]);
else % beam offset
    ylabel(sprintf('Sensitivity:  dlog[nV(out) @ %0.0f ps]/dlogX',tdelay(1)*1e12), 'FontSize',16);
    axis([min(xvect) max(xvect) -2 2]);
end

%% export sensitivities
%S_LCTE = S_LCTE;
S_sys = {S_f,S_r_pump,S_r_probe};
end
