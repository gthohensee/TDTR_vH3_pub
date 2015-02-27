function [LambdaAl_roomtemp, TL] = generate_lambdaAl_WF(rho_measured, outputdir)
%generate_lambdaAl_WF - Computes temperature-dependent thermal
%conductivity of Aluminum by referring to the measured four-point probe
%resistivity, the Wiedermann-Franz Law, and literature data for the
%intrinsic resistivity of Al at low temperatures.
%
% Syntax:  [LambdaAl_roomtemp,TL] = generate_lambdaAl_WF(rho_measured)
%
% Inputs:
%    rho_measured - Number: measured resistivity of Al sample, SI units.
%    outputdir - String: directory to which to save output
%                temperature-dependent calculated thermal conductivity 
%                for this Al. If unspecified, file is not saved.
%
% Outputs:
%    LambdaAl_roomtemp - Number: Al room temperature thermal conductivity 
%                        predicted from Wiedermann-Franz law using measured
%                        resitivity.
%    LT - Two-column matrix of temperatures (K) and thermal conductivities
%         (W/m-K) for Al of the given resistivity from 20 K to 923 K.
%         Literature data is in increments of 20 K; this code interpolates
%         100x, for increments of 0.2 K.
%
% Example: 
%    [LAl, ~] = generate_lambdaAl_WF(3.90e-8)
%    result: LAl = 192 W/m-K, no file saved.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author: Gregory Hohensee
% University of Illinois Urbana-Champaign
% email: hohense2@illinois.edu
% Website: n/a
% Apr 2014; Last revision: 1-April-2014

%------------- BEGIN CODE --------------
%% Set standard initialization parameters
L0 = 2.445e-8; % constant, Wiedermann-Franz law.
rho_i_Al_295K = 2.678e-8; % intrinsic resistivity of Al at 295 K

% Check input parameters
if nargin < 2, outputdir = pwd; output = 0; end % Don't save output.
if nargin < 1
    % If no input given, calculate Lambda by specifying the extrinsic
    % resistivity, rather than the unspecified measured resistivity.
    rho_e_Al = 1.37e-8; % yields 178 W/m-K at RT. Used for spinladders #10 and #12
    % rho_e_Al = 1.22e-8 for spmini's 192 W/m-K Al from sputter 1, from
    % rho_e_Al = rho_measured - rho_intrinsic = 3.90e-8 - 2.678e-8 = 1.22e-8.
    ID = strcat(num2str(rho_e_Al),'e');
else % nargin = 1 or 2
    rho_e_Al = rho_measured - rho_i_Al_295K;
%   % Next two lines back-calculate rho_e_Al from lambda0_Al, for the lazy.
    %lambda0_Al = 180.5744; % alloy check
    %rho_e_Al = (L0 * 295 / lambda0_Al) - 2.678e-8
    ID = strcat(num2str(rho_measured), 'm');
end

% output filename indicates the used resistivity with an "e" or "m"
% appended to indicate extrinsic or measured resistivity.
outputfile = strcat(outputdir,'/Al_WF_interpolated_',ID,'.txt');

%%
% Established data regarding intrinsic Al resistivity, from Landolt and
% Bornstein, section 1.2.2, page 15, tables 1 and 3, extracts from a book
% whose author's name is "Bass". I received scans from Rich on 5/2/2012.
% Units are Kelvins and microOhm-cm
rho_i_Al = [20	0.0007
30	0.0046
40	0.0180
50	0.0474
60	0.0957
70	0.1626
80	0.2449
90	0.3395
100	0.4420
120	0.6627
140	0.8930
160	1.127
180	1.361
200	1.593
220	1.824
240	2.053
260	2.280
273.2	2.430
293	2.654
295	2.678
373	3.556
473	4.692
573	5.846
673	7.038
773	8.280
823	8.942
873	9.642
898	10.01
923	10.39];

rho_i_Al(:,2) = rho_i_Al(:,2) * 1e-8; % convert units to Ohm-m.

%% calculate Wiedermann-Franz predicted thermal conductivity Lambda 
% from total resistivity.
rho_Al = rho_e_Al + rho_i_Al(:,2);
Lambda_WF = ((rho_Al)./(L0*rho_i_Al(:,1))).^-1;

% Using the rho_i_Al(:,1) temperature data and Lambda_WF, interpolate by
% 100x to produce a two-column array predicting L(T) for Al with this
% room-temperature resistivity.
T_interp = interp(rho_i_Al(:,1),100);
Lambda_WF_interp = interp(Lambda_WF,100);
[~,index] = min(abs(T_interp - 295));
LambdaAl_roomtemp = Lambda_WF_interp(index);
TL = [T_interp,Lambda_WF_interp];
% Save to specified filename and directory, if given.
if output, dlmwrite(outputfile, TL); end

end
    
    
    