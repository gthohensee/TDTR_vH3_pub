function [Psol,fitOK] = VoutLinearFit_vH3(DATAMATRIX)
%VoutLinearFit_vH3 - Assists user in specifying a linear fit to the V(out)
% data, in case V(out) is very small and one wishes to smooth
% the ratio from the disproportionate V(out) noise. Relies on the
% assumption that V(out) data is strictly linear with random noise.
%
% [Psol,fitOK] = VoutLinearFit_vH3(DATAMATRIX)
%
% Inputs:
%    DATAMATRIX - matrix of TDTR data in the standard format
%
% Outputs:
%    Psol - [p1 p2], for the curve V(out) = p1*tdelay + p2, 
%           where p1 is in uV/ps and p2 is in uV.
%    fitOK - TRUE if user trusts the fit.
%
%
% Other m-files required: extract_interior.m
% Subfunctions: none
% MAT-files required: none
%
% See also: SetPhase_vH3.m, extract_interior.m

% Author: Gregory T. Hohensee
% University of Illinois Urbana-Champaign
% email: hohense2@illinois.edu
% Website: n/a
% April 2014; Last revision: 3-April-2014
%                            14-July-2014: vH2. Header comments.
%                            17-Feb-2015: vH3. No changes.
% BEGIN CODE %

%Encodercount=DATAMATRIX(:,1);
tdelay_raw=DATAMATRIX(:,2); % imported in picoseconds
Vin_raw=DATAMATRIX(:,3); %Units (uV ?)
Vout_raw=DATAMATRIX(:,4);
ratio_raw=DATAMATRIX(:,5);

%Choose range of time delays to fit
tdelay_min=min(tdelay_raw);
tdelay_max=max(tdelay_raw);

% extract data points
[~,ratio_data]=extract_interior(tdelay_raw,ratio_raw,tdelay_min,tdelay_max);
[~,Vin_data]=extract_interior(tdelay_raw,Vin_raw,tdelay_min,tdelay_max);
[tdelay_data,Vout_data]=extract_interior(tdelay_raw,Vout_raw,tdelay_min,tdelay_max);

%% Initialize V(out) linear fit
%[~,ind] = min(abs(tdelay_data - 100)); % restrict fitting to start at 100 ps
%[FO, G] = fit(tdelay_data(ind:length(tdelay_data)),Vout_data(ind:length(Vout_data)),'poly1')
[FO, G] = fit(tdelay_data,Vout_data,'poly1')

p1 = FO.p1;
p2 = FO.p2;
Pguess = [p1 p2];

% Optimize the V(out) linear fit such that the ratio signal is minimally
% changed (linearizing V(out) should only reduce noise in the ratio).
Psol = fminsearch(@(M) VoutLinearFit_ratiofit(M,tdelay_data,Vin_data,ratio_data),Pguess);
p1 = Psol(1);
p2 = Psol(2);

%% Compute V(out) and ratio fits
Vout_lin = p1 * tdelay_data + p2;

r_lin = -Vin_data ./ Vout_lin;
if mean(r_lin) < 0
    r_lin = r_lin*-1;
    ratio_data = ratio_data*-1;
end

%% Position the V(out) and ratio figures

% % Initialize the two V(out) figures % %
set(0,'Units','pixels') 
scnsize = get(0,'ScreenSize');

figv = 111;
figr = 112;
if ishandle(figv), close 111; end
if ishandle(figr), close 112; end
h1 = figure(figv);
h2 = figure(figr);

% The figure Position property only includes the drawable extent of the window, 
%exclusive of the window borders. Obtain the entire window's size from the 
%OuterPosition property, and compare the two:
position = get(h1,'Position');
outerpos = get(h1,'OuterPosition');
borders = outerpos - position;
%The left, right, and bottom borders are each 4 pixels wide. The top border, 
%which contains a menu bar and a figure toolbar is 75-4, or 71 pixels wide.

%Define the desired size and location of the figures. Leave a space equal 
%to their border width between them:
edge = -borders(1)/2; % half a border's width
pos1 = [edge,... % left
        scnsize(4) * (1/3),... % fig1 will be in middle 1/3 of screen
        scnsize(3)/3 - edge,... % width
        scnsize(4)/3]; % height
pos2 = [scnsize(3)/3 + edge,... % right of fig1
        pos1(2),... % same bottom
        pos1(3),... % same width
        pos1(4)]; % same height
    
%Reposition the two figures by changing both of their OuterPosition properties:
set(h1,'OuterPosition',pos1)
set(h2,'OuterPosition',pos2)
%% Plot V(out) and the linear fit together
figure(figv)
plot(tdelay_data, Vout_data, '-bo', tdelay_data, Vout_lin, '-ro');
vmin = min([Vout_data;Vout_lin]);
vmax = max([Vout_data;Vout_lin]);
axis([tdelay_min 4000 vmin vmax]);

set(gca, 'TickLength' , [.02 .02]);
set(gca,'FontSize',20);
xlabel('Time delay (ps)');
ylabel('V_{out} (uV)');
legend('data','linearized');

%% Plot the ratio and the ratio fit together
[~,i10] = min(abs(tdelay_data - 10));
nd = length(tdelay_data);

figure(figr)
loglog(tdelay_data(i10:nd),ratio_data(i10:nd),'-bo',...
       tdelay_data(i10:nd),r_lin(i10:nd),'-ro');
rmin = (1/2)*min([ratio_data(i10:nd);r_lin(i10:nd)]);
rmax = 2*max([ratio_data(i10:nd);r_lin(i10:nd)]);
axis([10 5000 rmin rmax]);
set(gca, 'XTick', [1 2 5 10 20 50 100, 200, 500, 1000, 2000, 5000, 1e4]);
set(gca, 'XMinorTick', 'off');
set(gca, 'YTick', [0.1 0.2 0.5 1 2 5 10 20 30 50 100 200]);
set(gca, 'TickLength' , [.02 .02]);
set(gca,'FontSize',20);
xlabel('Time delay (ps)');
ylabel('Ratio');
legend('data','linearized');

%% Let user decide on the best fit now
fitOK = input('Hit RETURN if this fit is OK, hit 0 to change it: ');
if isempty(fitOK), fitOK = 1;
else fitOK = 0; end

if ~fitOK
    done = 0;
    while done ~= 1
        Vout_lin = p1 * tdelay_data + p2;
        
        r_lin = -Vin_data ./ Vout_lin;
        if mean(r_lin) < 0
            r_lin = r_lin*-1;
            ratio_data = ratio_data*-1;
        end
        
        % refresh figures
        figure(figv)
        plot(tdelay_data, Vout_data, '-bo', tdelay_data, Vout_lin, '-ro');
        vmin = min([Vout_data;Vout_lin]);
        vmax = max([Vout_data;Vout_lin]);
        axis([tdelay_min 4000 vmin vmax]);
        xlabel('Time delay (ps)');
        ylabel('V_{out} (uV)');
        legend('data','linearized');
        
        figure(figr)
        loglog(tdelay_data(i10:nd),ratio_data(i10:nd),'-bo',...
               tdelay_data(i10:nd),r_lin(i10:nd),'-ro');
        rmin = (1/2)*min([ratio_data(i10:nd);r_lin(i10:nd)]);
        rmax = 2*max([ratio_data(i10:nd);r_lin(i10:nd)]);
        axis([10 5000 rmin rmax]);
        set(gca, 'XTick', [1 2 5 10 20 50 100, 200, 500, 1000, 2000, 5000, 1e4]);
        set(gca, 'XMinorTick', 'off');
        set(gca, 'YTick', [0.1 0.2 0.5 1 2 5 10 20 30 50 100 200]);
        xlabel('Time delay (ps)');
        ylabel('Ratio');
        legend('data','linearized');
        
        % Inform user of current fit parameters
        Mat = 'Current linear fit parameters V(out) = p1*x + p2: ';
        Mat = char(Mat,sprintf('p1 = %0.3d uV/ps',p1));
        Mat = char(Mat,sprintf('p2 = %0.3f uV',p2));
        Mat
        
        % get and clean input
        done = input('Enter 1 if done, else hit "Enter": ');
        if isempty(done) 
            done = 0;
        else
            if done ~= 1
                done = 0;
            end
        end

        if ~done
            tp1 = input(sprintf('Adjust parameter p1: '));
            tp2 = input(sprintf('Adjust parameter p2: '));
            if ~isempty(tp1)
                p1 = tp1;
            end
            if ~isempty(tp2)
                p2 = tp2;
            end
        end
    end
    Psol = [p1 p2];
end

end % end main function

%% This function checks how closely the linearized ratio follows the original ratio
function [Z,P] = VoutLinearFit_ratiofit(P, tdelay_data, Vin_data, ratio_data)

p1 = P(1);
p2 = P(2);
Vout_lin = p1 * tdelay_data + p2;
r_lin = -Vin_data ./ Vout_lin;

res = (r_lin - ratio_data).^2;
Z = sum(res);
end

% ------- END CODE ------- %
