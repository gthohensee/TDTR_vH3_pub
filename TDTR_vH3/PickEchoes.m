function [echotime,t0] = PickEchoes(data,tmax,t0,t_window)
%PickEchoes - User picks time delay of acoustic travel time(s) in thin films.
%
% Syntax:  [echotime,t0] = PickEchoes(data,tmax,t0,t_window)
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
% Example: 
%    echotime = PickEchoes(data,200,1.32)
%    [echotime,t0] = PickEchoes(data,200)
%
% Other m-files required: extract_interior.m, set_t0.m
% Subfunctions: none
% MAT-files required: none
%
% See also: SetPhase_vH, extract_interior.m, set_t0.m

% Author: Gregory T. Hohensee
% University of Illinois Urbana-Champaign
% email: hohense2@illinois.edu
% Website: n/a
% March 2014; Last revision: 22-March-2014
%                            14-July-2014: vH2.
%------------- BEGIN CODE --------------
% Check input
if nargin < 1, error('At least input the TDTR data file!'); end
if nargin < 2, tmax = 200; end
if nargin < 4, t_window = [-10 10]; end
if nargin < 3, [t0] = set_t0(data, t_window); end % set t0

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
pickV = input('Enter 1 or 0 to pick echoes from V(in) or ratio, respectively.\n');
if pickV
    figure(h1);
    close 115;
else
    figure(h2);
    close 114;
end

% loop for user input on echo times
n = 1;
test=1;
while test
    if pickV
        figure(h1);
        plot(tdelay_data,Vin_data,'bo-')
        axis([tdelay_min tdelay_max min(Vin_data) max(Vin_data)]);
        xlabel('Time delay (ps)')
        ylabel('V_{in}')
    else
        figure(h2);
        plot(tdelay_data,ratio_data,'ro-')
        axis([tdelay_min tdelay_max min(ratio_data) max(ratio_data)]);
        xlabel('Time delay (ps)')
        ylabel('Ratio')
    end
    
    input('Examine the plots, zoom in as desired, and hit RETURN when ready.');
    fprintf('Select t(echo)...the press RETURN to continue\n')
    
    if pickV
        figure(h1);
        hpick = h1;
    else
        figure(h2);
        hpick = h2;
    end
    
    [echotime(n),~] = ginput(1);
    sprintf('You chose t(%i) = %0.2f ps.\n',n,echotime(n))
    n = n + 1;
    
    test = input('If done, enter 0: ');
    if isempty(test) || test ~= 0, test = 1; end
end
if pickV, close 114; else close 115; end
end

