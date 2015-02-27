function [Bguess,t0] = PickEchoesB(data,tmax,t0,t_window)
%PickEchoesB - Given a data set with periodic oscillations, this script
%              assists user in making a sinusoidal fit to the data
%              after subtracting the thermal (exponential) background.
%              Good for extracting Brillouin backscattering data.
%
% Syntax:  [Bguess,t0] = PickEchoesB(data,tmax,t0,t_window)
%
% Inputs:
%    data - 5-6 column matrix of raw or processed TDTR data.
%           Format must be [~, time_delay, V(in), ~, ratio, detV],
%           where ~ is unimportant.
%    tmax - optional: maximum time delay in ps, beyond which no
%           acoustic or other features of interest exist in the data.
%           Default is 200 ps.
%    t0 - optional: zero of time delay in ps. If specified, program does
%         not call set_t0.
%    t_window - optional 1x2 array: window for picking t0 in set_t0.m.
%               default is [-10 10].
% Outputs:
%    echotime - 1xN array of time delays corresponding to acoustic echoes
%               or travel times in films, or any times of interest.
%    t0 - Zero of delay time, as determined with set_t0.m
%
%
% Other m-files required: extract_interior.m, set_t0.m
% Subfunctions: none
% MAT-files required: none
%
% See also: SetPhase_vH2, PickEchoes.m

% Author: Gregory T. Hohensee
% University of Illinois Urbana-Champaign
% email: hohense2@illinois.edu
% Website: n/a
% March 2014; Last revision: 22-March-2014
%                            14-July-2014: vH2. More comments in header.

%------------- BEGIN CODE --------------
% Check input
if nargin < 1, error('At least input the TDTR data file!'); end
if nargin < 2, tmax = 200; end
if nargin < 4, t_window = [-10 10]; end
if nargin < 3, [t0] = set_t0(data, t_window); end % set t0
Bguess = 0;

% Interpret data
Encodercount=data(:,1);
tdelay_raw=data(:,2) - t0; %imported in picoseconds. Adjust t0.
Vin_raw=data(:,3); %Units (uV ?)
Vout_raw=data(:,4);
ratio_raw=data(:,5);

%Choose range of time delays to show
tdelay_min= -20;
tdelay_max= tmax;

%Extract Points in that time range
[~,ratio_data]=extract_interior(tdelay_raw,ratio_raw,tdelay_min,tdelay_max);
[tdelay_data,Vin_data]=extract_interior(tdelay_raw,Vin_raw,tdelay_min,tdelay_max);

scrsz = get(0,'ScreenSize');
h1 = figure(114);
h2 = figure(115);
% left bottom width height
set(h1,'Position',[scrsz(3)/5 scrsz(4)/2 scrsz(3)/4 scrsz(4)/4])
set(h2,'Position',[scrsz(3)/2 scrsz(4)/2 scrsz(3)/4 scrsz(4)/4])

%% present the figures
figure(h1)
plot(tdelay_data,Vin_data,'bo-')
set(gca, 'TickLength' , [.02 .02]);
set(gca,'FontSize',20);
axis([tdelay_min tdelay_max min(Vin_data) max(Vin_data)]);
xlabel('Time delay (ps)')
ylabel('V_{in}')

figure(h2)
plot(tdelay_data,ratio_data,'ro-')
set(gca, 'TickLength' , [.02 .02]);
set(gca,'FontSize',20);
axis([tdelay_min tdelay_max min(ratio_data) max(ratio_data)]);
xlabel('Time delay (ps)')
ylabel('Ratio')

% User decides which plot is less noisy, has clearer echoes
%pickV = input('Enter 1 or 0 to pick echoes from V(in) or ratio, respectively.\n');
pickV = 1;
if pickV
    figure(h1);
    close 115;
    [stdelay_data, sVin_data] = Subtract_Background(tdelay_data,Vin_data);
    signal = sVin_data;
else
    figure(h2);
    close 114;
    [stdelay_data, sratio_data] = Subtract_Background(tdelay_data,ratio_data);
    signal = sratio_data;
end
if isempty(signal), return; end
if pickV
    figure(h1);
    plot(stdelay_data,signal,'bo-')
    axis([tdelay_min tdelay_max min(signal) max(signal)]);
    xlabel('Time delay (ps)')
    ylabel('V_{in}')
else
    figure(h2);
    plot(stdelay_data,signal,'ro-')
    axis([tdelay_min tdelay_max min(signal) max(signal)]);
    xlabel('Time delay (ps)')
    ylabel('Ratio')
end
%%
sinefit = @(b,x)  b(1).*(sin(2*pi*x./b(2) + (2*pi/180)*b(3))) + b(4);    % Function to fit

Bguess = [1 30 0 0]; % 30 ps, amplitude 1 sine wave
tag = {'Amplitude', 'Period (ps)', 'Shift (deg)', 'Offset'};

figure(88)
plot(stdelay_data,signal,'-bo',stdelay_data,sinefit(Bguess,stdelay_data),'-r')

done = 0;
while ~done
    tag
    Bguess
    for i = 1:4 % iterate through fit parameters
        temp = input(sprintf('Adjust parameter %s: ',tag{i}));
        if ~isempty(temp), Bguess(i) = temp; end
    end % end iteration through fit parameters

    % Show fit
    figure(88)
    plot(stdelay_data,signal,'-bo',stdelay_data,sinefit(Bguess,stdelay_data),'-r')
    done = input('If done, enter 1: ');
    if isempty(done) || done ~= 1, done = 0; end
end
if pickV, close 114; else close 115; end
end

