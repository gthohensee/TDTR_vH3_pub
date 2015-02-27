function [S_LCTE,S_sys,xvar,SS_LCTE,SS_sys] = parametric_senseplot_vH3(Xij,...
    SYS_i, xvar,datparams,sysparams, calparams, matparams, Tparams,LCTE_sens_consider,sys_consider)
%parametric_senseplot_vH3 - Calculates sensitivity plots dlogR/dlogX for 
%thermal model. The sens_consider booleans can be edited to change which
%sensitivities are calculated. The difference here is that dlogR/dlogX
%is calculated at a series of parametric points in LCTE and system
%parameter space, rather than the time or spatial axis.
%
% Inputs:
%    Xij       - BOOLEAN 2x1 (i,j) index of LCTE elements
%                marking which LCTE variable is the independent variable 
%                over which sensivities are calculated. Zero = none.
%    SYS_i     - BOOLEAN scalar key with value of 1, 2, 3, or 0 
%                indicating which of f, r_pump, r_probe, or none are the
%                independent variable over which sensitivities are
%                calculated.
%    xvar      - 1D scalar array for x-axis (independent variable). 
%                LCTE_ij and SYS_i serve to mark which physical
%                or system parameter it represents.
%    datparams - {tdelay ratio_data datadir offset}
%    sysparams - {tau_rep f r_pump r_probe}
%    calparams - {Zind sigfit intscheme nnodes consider_error LCTE_err T0_err P0_err}
%    matparams - {LCTE aniso BI n_toplayer TCR doughnut}
%    Tparams   - {T0, T_LCTE, A_pump, A_probe, absC, perpulse, jabs, jtrans}
%[See INITIALIZE_CELLPARAMS_vH3 for details.]
%    LCTE_sens_consider - Boolean 4xN matrix indicating which elements
%                         of LCTE are to have their sensitivities
%                         calculated.
%    sys_consider - Boolean 4x1 matrix indicating which elements
%                   [f, r_pump, r_probe, tau_rep] of system parameters
%                   to have their sensitivities calculated.
%
% Outputs
%    S_LCTE - sensitivities to one-channel overlayer parameters
%    S_sys  - {S_f,S_r_pump,S_r_probe,S_tau_rep};
%       S_f, S_r_pump, S_r_probe - system parameter sensitivities
%    xvar   - the x-axis independent variable over which the selected 
%             sensitivities S_LCTE and S_sys are calculated.
%    SS_LCTE - S_LCTE cell, except the elements are arrays in time
%              or offset and there is an extra cell dimension of
%              length(xvar).
%    SS_sys  - same idea as SS_LCTE and S_sys.
%
% Other m-files required: TDTR_REFL_vH3.m
% Subfunctions: none
% MAT-files required: none
%
% See also: errorbars_vH3.m, senseplot_offset_vH3.m

% Author: Gregory Hohensee
% U. of Illinois Urbana-Champaign
% email: hohense2@illinois.edu
% Website: GitHub
% Revision history: 28-Jan-2015 - written.
%                   17-Feb-2015 - vH3.
%------------- BEGIN CODE --------------
%%
INITIALIZE_CELLPARAMS_vH3; % unpacks/cleans cellparams (the five inputs)
fprintf('Generating Arbitrary Sensitivity Plot for All Variables...please be patient\n')

if doughnut
    by_FWHM = -1;
    while by_FWHM == -1
        by_FWHM = input('Enter 1 if signal is FWHM of V(out), else 0 if signal is just V(out): ');
        if by_FWHM ~= 1 || by_FWHM ~= 0
            by_FWHM = -1;
            sprintf('Invalid input, try again.\n');
        end
    end
else
    by_FWHM = 0; % default; not used for conventional TDTR
end

%% initialize matrices for sensitivities

%nt = length(tdelay); 
nl = length(LCTE(1,:));
nx = length(xvar);
SS_LCTE = cell(4,nl,nx);
SS_sys = cell(4,nx);

S_LCTE = cell(4,nl);
S_sys = cell(4);

deltaR_temp  = cell(4,nl);
ratio_temp   = cell(4,nl);

%% which sensitivities to consider? (saves time)
if ~exist('sys_consider','var')
    sys_consider = [1 1 1 0]; % f, r_pump, r_probe, tau_rep
end

if ~exist('LCTE_sens_consider','var')
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
end

%% Begin loop
for k = 1:length(xvar)
    fprintf('Iteration %i out of %i...\n',k,length(xvar));
    %% Update independent variable in thermal or system parameters
    if Xij(1)
        LCTE(Xij(1), Xij(2)) = xvar(k);
    else
        switch SYS_i
            case 1
                f = xvar(k);
            case 2
                r_probe = xvar(k);
            case 3
                r_pump = xvar(k);
            case 4
                tau_rep = xvar(k);
            otherwise % includes case 0
                error('You have not properly specified an independent variable for this sensitivity plot.');
        end
    end
    %% Refresh reference model for new independent variable value
     matparams{1} = LCTE; 
     sysparams = {tau_rep f r_pump r_probe};
     
    [nVin_model, nVout_model, deltaR_model,ratio_model, FWHM_model] = ...
    TDTR_SENS_MODEL(doughnut, sigfit, by_FWHM, tdelay, offset, ...
                    matparams, sysparams, A_pump, Zind, intscheme, nnodes);
%     return;
    %% Now that the LCTE(i,j), f, r_probe, or r_pump have been changed to
    % match the next element of xvar, we compute the chosen sensitivities
    % (according to LCTE_sens_consider and sys_consider) of the chosen
    % signal (Vout offset, ratio, Vin, Vout) at a single delay time or
    % offset position, about this new reference point set by xvar.
    % It is compatible with N>1 or M>1 time or space elements, but
    % this script isn't set up for 3D plotting.
    %
    % Everything from here proceeds identically as in senseplot_offset_vH3,
    % except the S_LCTE and S_sys are built of k elements.
    %% -----------------Compute sensitivities for LCTE--------------
    LCTEtemp = LCTE;
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
            [deltaR_temp{i,j},ratio_temp{i,j}] = TDTR_REFL_vH3(tdelay,matparams,sysparams,A_pump,intscheme,nnodes,xvect);
            
            % subfunction for calculating Num, depending on the signal.
            Num = sense_getNum(doughnut,sigfit,by_FWHM,tdelay,xvect,...
                               deltaR_temp{i,j},ratio_temp{i,j},...
                               nVin_model,nVout_model,ratio_model);
            
            Denom=log(LCTEtemp(i,j))-log(LCTE(i,j));
            SS_LCTE{i,j,k}=Num/Denom;

            LCTEtemp = LCTE;
            matparams{1} = LCTE;
        end
    end
    clear LCTEtemp;

    %% -----------------Compute sensitivities for system parameters-----------
    for i = 1:6
        if sys_consider(i) == 0, continue; end
        
        sysbump = 1.01;
        switch i
            case 1, sysparams{2} = f*sysbump;
            case 2, sysparams{3} = r_pump*sysbump;
            case 3, sysparams{4} = r_probe*sysbump;
            case 4, sysparams{1} = tau_rep*sysbump;
            case 5 % "w0", overall spot size is perturbed
                sysparams{3} = r_pump*sysbump;
                sysparams{4} = r_probe*sysbump;
            case 6
                if sigfit == 0
                    SS_sys{i,k} = ratio_model + ratio_model.^-1; % Definition of phase sensitivity.
                else
                    error('Sorry, phase sensitivity is only programmed for ratio signals as of vH3.')
                end
        end
        
        if i ~= 6 % non-phase sensitivities
            % Thermal model is symmetric around 0 offset, so only use
            % the positive offset side of the offset vector to generate model.
            [~,I] = min(abs(offset));
            xvect = offset(I:length(offset)); % this will be 0 if conventional TDTR
            [temp_deltaR,temp_ratio]=TDTR_REFL_vH3(tdelay,matparams,sysparams,A_pump,intscheme,nnodes,xvect);
            
            Num = sense_getNum(doughnut,sigfit,by_FWHM,tdelay,xvect,...
                               temp_deltaR,temp_ratio,...
                               nVin_model,nVout_model,ratio_model);
                           
            Denom=log(sysbump);
            SS_sys{i,k}=Num/Denom;
        end
        
        sysparams = {tau_rep, f, r_pump, r_probe}; % return to reference value
    end
    % S_sys is now a 4 by k cell, with each element being the sensitivity
    % for a given time and space range (preferably singular), and for
    % the kth element of xvar, against which we will plot sensitivities.
end
fprintf('Calculated SS_LCTE and SS_sys...\n',i);

%% Compressing sensitivities into xvar space, collapsing time and offset.
Ix = 1; % default index for compressing SS_LCTE, SS_sys into xvar arrays.
It = 1; % same, for time axis
if sum(sum(LCTE_sens_consider)) > 0
    test_index = find(~cellfun(@isempty,SS_LCTE),1);
    test_sens = SS_LCTE{test_index};
else
    test_index = find(~cellfun(@isempty,SS_sys),1);
    SS_sys{test_index};
    %test_sens = SS_sys{test_index};
end
[Mt,Mx] = size(test_sens);

if Mt > 1
    toff = input('You have an array of time delays. Please enter just one, in picoseconds: ');
    [~, It] = min(abs(tdelay*1e12 - toff));
    tval = tdelay(It)*1e12;
    fprintf('The nearest computed time delay is: %i ps\n',tval);
else
    tval = tdelay*1e12;
end

if Mx > 1
    xoff = input('You have an array of offset positions. Please enter just one, in microns: ');
    [~, Ix] = min(abs(offset - xoff));
    offval = offset(Ix);
    fprintf('The nearest computed beam offset position is: %i microns\n',offval);
else
    offval = offset;
end

% Rebuild SS_LCTE and SS_sys as functions of xvar only. Pick the (It,Ix) 
% index of the first dimension of each element of SS_LCTE{i,j,k} and 
% SS_sys[i,k}. Make an array out of the k elements. Produce 
% S_LCTE{i,j} and S_sys{i}, analogous to output of senseplot_offest_vH3.

% SS_LCTE to S_LCTE...
[ni, nj, nk] = size(SS_LCTE);
for i = 1:ni
    for j = 1:nj
        array_temp = [];
        for k = 1:nk
            if ~isempty(SS_LCTE{i,j,k})
                temp = SS_LCTE{i,j,k};
                array_temp(k) = temp(It,Ix);
            end
        end
        S_LCTE{i,j} = array_temp;
    end
end

% SS_sys to S_sys...
[ni, nk] = size(SS_sys);
for i = 1:ni
    array_temp = [];
    for k = 1:nk
        if ~isempty(SS_sys{i,k})
            temp = SS_sys{i,k};
            array_temp(k) = temp(It,Ix);
        end
    end
    S_sys{i} = array_temp;
end
S_f = S_sys{1}; 
S_r_pump = S_sys{2}; 
S_r_probe = S_sys{3}; 
S_tau_rep = S_sys{4};
S_w0 = S_sys{5};
S_phase = S_sys{6};
% Now S_LCTE and S_sys contain elements in xvar parameter space.
%% Plot sensitivities
fignum = 404;
figure(fignum)
clf

set(gca,'Box','on');
hold on;
ColorOrder = get(gcf,'DefaultAxesColorOrder');

%% label LCTE sensitivities
LCTElegend = [];
LCTElab = cell(size(LCTE));
LCTEmarker = {'o-','*-','x-','+-'};

% define rescaling factor for x-axis
if Xij(1)
    LCTEtoX = [1, 1e-6, 1e9, LCTE(1,Xij(2))]; 
    scale = LCTEtoX(Xij(1));
else
    SYStoX = [1e-6, 1e6, 1e6, 1e-6, 1e6]; % note: phase shift is not a valid independent variable in this program.
    scale = SYStoX(SYS_i); 
end
    
for j = 1:nl
    for i = 1:4
        switch i
            case 1, LCTElab{i,j} = sprintf('L%i',j);
            case 2, LCTElab{i,j} = sprintf('C%i',j);
            case 3, LCTElab{i,j} = sprintf('t%i',j);
            case 4, LCTElab{i,j} = sprintf('e%i',j);
        end
        
        if isempty(S_LCTE{i,j}), continue;
        else
            loglog(xvar*scale,abs(S_LCTE{i,j}),LCTEmarker{i},...
                'MarkerSize',8.0,'LineWidth',1,'Color',ColorOrder(j,:));
            LCTElegend = [LCTElegend;LCTElab{i,j}];
        end
    end
    
    if isempty(S_LCTE{1,j}) || isempty(S_LCTE{4,j}), continue;
    else
        loglog(xvar*scale,abs(S_LCTE{1,j} - S_LCTE{4,j}),'s-',...
            'MarkerSize',8.0,'LineWidth',1,'Color',ColorOrder(5,:));
        LCTElegend = [LCTElegend;sprintf('Z%i',j)];
    end
end



%% plot and label system sensitivities
syslegend = [];
for i = 1:length(sys_consider)
    switch i
        case 1
            S_dat = S_f;
            legendtxt = 'ff';
        case 2
            S_dat = S_r_pump;
            legendtxt = 'Rp';
        case 3
            S_dat = S_r_probe;
            legendtxt = 'Rb';
        case 4
            S_dat = S_tau_rep;
            legendtxt = 'TR';
        case 5
            S_dat = S_w0;
            legendtxt = 'w0';
        case 6
            S_dat = S_phase;
            legendtxt = 'ph';
    end
    if sys_consider(i)
        loglog(xvar*scale,abs(S_dat),'--',...
            'MarkerSize',8.0,'LineWidth',1,'Color',ColorOrder(i,:));
        syslegend = [syslegend;legendtxt];
    end
end

%% Other plot details
figure(fignum);
set(gca,'FontSize',18);
set(gca, 'TickLength' , [.02 .02]);

legend([LCTElegend;syslegend])

% determine which parameter xvar is, again, to specify its xlabel.
xlabelstr = '';
if Xij(1)
    switch Xij(1)
        case 1, xlabelstr = 'Thermal conductivity #%i (W m^{-1} K^{-1})';
        case 2, xlabelstr = 'Heat capacity #%i (J cm^{-3} K^{-1})';
        case 3, xlabelstr = 'Thickness #%i (nm)';
        case 4, xlabelstr = 'In-plane conductivity #%i (W m^{-1} K^{-1})';
    end
    xlabelstr = sprintf(xlabelstr,Xij(2));
    xmin = xvar(1) * LCTEtoX(Xij(1)); % xvar is in SI units
    xmax = max(2*xvar(1),xvar(end)) * LCTEtoX(Xij(1));
else
    switch SYS_i
        case 1, xlabelstr = 'Pump modulation frequency (MHz)';
        case 2, xlabelstr = 'Probe spot size (microns)';
        case 3, xlabelstr = 'Pump spot size (microns)';
        case 4, xlabelstr = 'Pulse repetition rate (MHz)';
        case 5, xlabelstr = 'Pump and probe spot size (microns)';
    end
    xlabelstr = sprintf(xlabelstr);
    xmin = xvar(1) * scale;
    xmax = max(2*xvar(1),xvar(end)) * scale;
end
    
xlabel(xlabelstr,'FontSize',18); % assign xlabel

% generate ylabel
if ~doughnut
    switch sigfit
        case 1, ylabel(sprintf('Sensitivity:  dlog[nV(in) @ %0.0f ps]/dlogX',tval), 'FontSize',18);
        case 2, ylabel(sprintf('Sensitivity:  dlog[nV(out) @ %0.0f ps]/dlogX',tval), 'FontSize',18);
        otherwise, ylabel('Sensitivity:  dlogR/dlogX', 'FontSize',18);
    end
else % beam offset signal
    if by_FWHM
        ylabel(sprintf('Sensitivity:  dlog[V(out) FWHM]/dlogX'), 'FontSize',18);
    else
        ylabel(sprintf('Sensitivity:  dlog[nV_{out}(%0.0f ps, %0.0f microns)]/dlogX',tval,offval), 'FontSize',18);
    end
end

axis([xmin xmax 0.01 2]); % set x-y axis bounds
set(gca,'XScale','Log');
set(gca,'YScale','Log');
%% export sensitivities
end % End of function

%% Subfunction: Compute model signal for sensitivities
function [nVin_model, nVout_model, deltaR_model,ratio_model, FWHM_model] = ...
    TDTR_SENS_MODEL(doughnut, sigfit, by_FWHM, tdelay, offset, ...
                    matparams, sysparams, A_pump, Zind, intscheme, nnodes)
%% Initialize reference models (outputs)
nVin_model = [];
nVout_model = [];
deltaR_model = [];
ratio_model = [];
FWHM_model = [];

%% Calculate desired reference model signal
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
            Vout_model_Zdelay = mean(Vout_model(Zind));
            nVout_model = Vout_model / Vout_model_Zdelay;

            %Vout_data_Zdelay = mean(Vout_data(Zind));
            %nVout_data = Vout_data / Vout_data_Zdelay;
        otherwise % ratio fit
            % do nothing here
    end
else % beam offset
    % Thermal model is symmetric around 0 offset, so only use
    % the positive offset side of the offset vector to generate model.
    [~,I] = min(abs(offset));
    xvect = offset(I:length(offset));
    [deltaR_model,ratio_model]=TDTR_REFL_vH3(tdelay,matparams,sysparams,A_pump,intscheme,nnodes,xvect);
    
    % deltaR and ratio may be a matrix (tdelay, offset).
    for t_index = 1:length(tdelay)
        Sig_Amp = deltaR_model(t_index,:)'; 
        Sig_ratio = ratio_model(t_index,:)';

        %reconstruct the rest by symmetry; no data, so no need to interpolate.
        offset = [-flipud(xvect(2:end));xvect]; % makes sure offset and model line up.
        deltaR_model = [flipud(Sig_Amp(2:end));Sig_Amp];
        ratio_model = [flipud(Sig_ratio(2:end));Sig_ratio];
        Vout_model = imag(deltaR_model);
        Vin_model = real(deltaR_model);

        % get value and index of r = 0 model peak
        [norm,I] = max(abs(Vout_model)); 
        
        % normalize to peak, make positive.
        nVout_model(t_index,:) = sgn(Vout_model(I)) * Vout_model' / norm; 
        
        if by_FWHM
            [~,J] = min(abs(nVout_model - 0.5)); % gets half-rise index.
            Xguess = 0.4247 * 2*abs(offset(J)); % if offset was continuous
                                                % instead of a few points,
                                                % this would be the FWHM.

            if ~exist('options','var'); options = optimset('TolFun',1e-2,'TolX',1e-2); end
            Xsol = fminsearch(@(X) SpotSize_V4(X,nVout_model,offset),Xguess,options);
            FWHM_model(t_index,1) = 2 * sqrt(2*log(2)) * Xsol; % column vector over tdelay
        end
    end
end

end

%% subfunction for calculating Num based on signal type.
function Num = sense_getNum(doughnut, sigfit, by_FWHM, tdelay, xvect, temp_deltaR,temp_ratio,...
                            nVin_model,nVout_model,ratio_model)
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
else % doughnut means beam offset.
    for t_index = 1:length(tdelay) % for loop is crude hack for non-scalar tdelay.
        % deltaR and ratio may be a matrix (tdelay, offset).
        Sig_Amp = temp_deltaR(t_index,:)'; 
        Sig_ratio = temp_ratio(t_index,:)';

        %reconstruct the rest by symmetry; no data, so no need to interpolate.
        offset = [-flipud(xvect(2:end));xvect]; % makes sure offset and model line up.
        temp_deltaR = [flipud(Sig_Amp(2:end));Sig_Amp];
        temp_ratio = [flipud(Sig_ratio(2:end));Sig_ratio];

        Vout_temp = imag(temp_deltaR);
        [norm,I] = max(abs(Vout_temp)); % get value and index of r = 0 model peak
        nVout_temp = sgn(Vout_temp(I)) * Vout_temp / norm;

        % FWHM mode for depicting beam offset sensitivities.
        if by_FWHM && length(offset) > 4
            [~,J] = min(abs(nVout_temp - 0.5)); % gets half-rise index.
            Xguess = 0.4247 * 2*abs(offset(J)); % if offset was continuous
                                         % instead of a few points,
                                         % this would be the FWHM.

            if ~exist('options','var'); options = optimset('TolFun',1e-2,'TolX',1e-2); end
            Xsol = fminsearch(@(X) SpotSize_V4(X,nVout_temp,offset),Xguess,options);
            FWHM_temp = 2 * sqrt(2*log(2)) * Xsol;

            Num(t_index,1) = log(FWHM_temp) - log(FWHM_model(t_index));
        else
            if length(offset) <= 4
                by_FWHM = 0;
                warning('Cannot fit FWHM from so few offset points. Defaulting back to V(out)(x).')
            end
            % this will be a matrix of size [length(tdelay),length(offset).
            Num(t_index,:) = log(nVout_temp)' - log(nVout_model)'; 
        end
    end
    % now Num is a vertical vector in tdelay and a horizontal
    % vector in offset, unless by_FWHM = TRUE, in which case
    % Num is a scalar in the offset dimension.
end
end