%% Subfunction of TDTR_MAIN, for importing processed TDTR data.
% Input: tdelay, datadir, datain, mid, sigfit, 
%        psc, doughnut, r_pump scalar, Zdelay,
%        calparams, Vout_offset*, detV_raw*, 
%        Voutlinfit*, fitOK*, (*:optional).
%
% Output: signal_data and offset arrays,
%         r_pump array, datparams, calparams.

DM1 = dlmread(datain);
if ~doughnut % conventional TDTR
    offset = 0; % pump and probe beams are overlapped.

    tdelay_raw=DM1(:,2)*1e-12; %imported in picoseconds.  Change to seconds.
    Vin_raw=DM1(:,3); 
    Vout_raw=DM1(:,4);
    ratio_raw=DM1(:,5);

    [~,ratio_data]          =extract_interior(tdelay_raw,ratio_raw, tdelay_min,tdelay_max);
    [~,Vin_data]            =extract_interior(tdelay_raw,Vin_raw,   tdelay_min,tdelay_max);
    [tdelay,Vout_data] =extract_interior(tdelay_raw,Vout_raw,  tdelay_min,tdelay_max);
    if psc % psc: is pump spot changing over time delay? Relies on
           % "frac", the fractional change, assigned in analyze script.
        r_pump = r_pump*(1 + frac*tdelay/tdelay(end));
        sysparams = {tau_rep f r_pump r_probe}; % update sysparams
    end

    % Define Zind: Zdelay = tdelay(Zind), approximately.
    [~,Zind] = min(abs(tdelay - Zdelay*1e-12));

    % If the analysis comes with a reasonable linear fit to V_out, 
    % use that to smooth ratio_data. [NOT TESTED]
    if exist('Voutlinfit','var') && exist('fitOK','var')
        if Voutlinfit && fitOK
            Vout_lin = m1*1e12 * tdelay + m2; % [m1 uV/ps]*[1e12 ps/s]*[tdelay seconds] + uV
            r_lin = -Vin_data ./ Vout_lin;
            ratio_data = r_lin;
        end
    end
else % beam offset TDTR

    offset = DM1(:,1);
    Vout_data = DM1(:,8);
    Vin_data = DM1(:,7);
    ratio_data = -Vin_data./Vout_data;

    offset = offset-mid; % mid is defined in the process or analyze script
end

% shifting V(out) up or down to account for false signal component
if exist('Vout_offset','var')
    Vout_data = Vout_data - Vout_offset;
    ratio_data = -Vin_data ./ Vout_data;
end

% shifting V(in) up or down to account for false signal component
if exist('Vin_offset','var')
    Vin_data = Vin_data - Vin_offset;
    ratio_data = -Vin_data ./ Vout_data;
end % (added Feb. 7th, 2015)

% introducing detector voltage for dR/dT and such.
if length(DM1(1,:)) > 5 && ~doughnut
    detV_raw = DM1(:,6); % correct for beam offset data, too.
end 
if exist('detV_raw','var') && ~doughnut
    [~,detV_data] =extract_interior(tdelay_raw,detV_raw,tdelay_min,tdelay_max);
end

%% Compose remaining parameter cells, 
% now that Zind, tdelay, and data are defined.
switch sigfit
    case 1, datparams = {tdelay Vin_data datadir offset};
    case 2, datparams = {tdelay Vout_data datadir offset};
    otherwise datparams = {tdelay ratio_data datadir offset}; 
end

calparams = {Zind sigfit intscheme nnodes consider_error};