function [t0] = set_t0(data,t_window, norm)
%PickEchoes - User picks time delay of acoustic travel time(s) in thin films.
%This code does not edit the data matrix. The zero of delay time for Al
%transducers is typically half-rise. The zero of delay time of other
%transducers may be different, and should be referenced to that of an
%Al-coated sample on your particular pump-probe system.
%
%
% Syntax:  [t0] = set_t0(data, t_window)
%
% Inputs:
%    data - 5-6 column matrix of raw or processed TDTR data.
%           Expected format is [stageposition, time_delay, V(in), V(out),
%           ratio, detV].
%    t_window - optional; default is [-10 10] ps. Domain for viewing signal near t0.
%    norm - optional; indicates whether to use raw V(in) and ratio, or
%           V(in) normalized to square of probe power drift, as measured
%           by the detector voltage, column 6 of data).
% Outputs:
%    t0 - Scalar, units of picoseconds. Zero of delay time for input TDTR data. 
% 
% Example: 
%    [t0] = set_t0(data, [-10 10])
%
% Other m-files required: extract_interior.m
% Subfunctions: none
% MAT-files required: none
%
% See also: SetPhase_vH, extract_interior.m

% Author: Gregory T. Hohensee
% University of Illinois Urbana-Champaign
% email: hohense2@illinois.edu
% Website: n/a
% March 2014; Last revision: 22-March-2014

%------------- BEGIN CODE --------------

% Interpret data
%Encodercount=data(:,1);
tdelay_raw=data(:,2); %imported in picoseconds.
Vin_raw=data(:,3); %Units (uV)
%Vout_raw=data(:,4);
ratio_raw=data(:,5);

% if there's a sixth column in the data, assume it is Detector Voltage (mV)
if length(data(1,:)) > 5 
    detV_raw = data(:,6);
end

% handle optional inputs
if nargin < 3
    norm = 0;
    if nargin < 2 || length(t_window) ~= 2
        t_window = [-10 10]; % default time window
    end
end

if norm % normalize V(in) if indicated.
    Vin_raw = Vin_raw .* (detV_raw ./ detV_raw(1)).^2;
    %Vout_raw = Vout_raw .* (detV_raw ./ detV_raw(1)).^2;
end

%Choose range of time delays
tdelay_min   = min(t_window,-20);  % min time delay
tdelay_max   = max(t_window,80);   % max time delay

[~,ratio_data]=extract_interior(tdelay_raw,ratio_raw,tdelay_min,tdelay_max);
[tdelay_data,Vin_data]=extract_interior(tdelay_raw,Vin_raw,tdelay_min,tdelay_max);
%[~,Vout_data]=extract_interior(tdelay_raw,Vout_raw,tdelay_min,tdelay_max);
%[~,detV_data]=extract_interior(tdelay_raw,detV_raw,tdelay_min,tdelay_max);

% Normalize the rise in V(in) and ratio across zero time delay
% for ease of picking the half-rise point for t0.
[~,ipeak] = max(abs(Vin_data));
scaled_ratio_data = (ratio_data - min(ratio_data)) ./ ...
                  (max(ratio_data - min(ratio_data)));
if Vin_data(ipeak) > 0
    scaled_Vin_data = (Vin_data - min(Vin_data)) ./ ...
                  (max(Vin_data - min(Vin_data)));
else
    temp = -1 * Vin_data;
    scaled_Vin_data = (temp - min(temp)) ./ ...
                  (max(temp - min(temp)));
end

% Plot and have user pick zero of time delay
figure(100);
plot(tdelay_data,scaled_Vin_data,'-bo',...
     tdelay_data,scaled_ratio_data,'-ro')
axis([t_window(1) t_window(2) -0.2 1.2001])
set(gca, 'TickLength' , [.02 .02]);
set(gca,'FontSize',20);
xlabel('Time delay (ps)');
ylabel('Signal');
legend('V(in)', 'Ratio');

fprintf('Select t0 by clicking on the plot...\n')
[t0,~]=ginput(1);
fprintf('You chose t0 = %0.2f ps\n\n',t0)

figure(100);
close;
end

