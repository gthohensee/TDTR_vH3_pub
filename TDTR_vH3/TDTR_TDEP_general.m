function [C_T0 L_T0 CC LL] = TDTR_TDEP_general(fnameC, fnameL, T0, directory)
%TDTR_TDEP_general - Import reference (C,L) data, report C(T0) and L(T0)
%and return the entire imported [T,C] and [T,L] arrays. Can be used as a
%one-off call for C and L at a specific temperature, or to return the full
%C(T) and L(T) variables for self-consistent laser-heating thermal
%parameter assignment.
%
% For N-heat-channel materials, like the two-temperature model for the
% magnon and phonon heat channels in Ca9La5Cu24O41 spin ladders, or for
% electrons and phonons in metals, this function can be called N times with
% different fnameC and fnameL strings for the different channels.
%
% This program will NOT interpolate sparse T-dependent data. You should
% instead perform your own separate interpolation, save that to a file, and
% import that file with this function.
% 
% This program will automatically convert from J/cm^3-K to J/m^3-K heat
% capacity data, but it will not respond correctly if the heat capacity
% data is in any other units or representations, like J/g-K or J/mol-K.
% The program will not convert thermal conductivity values.
%
% Syntax:  [C_T L_T CC LL] = TDTR_TDEP_general(T0, directory, fnameC, fnameL)
%
% Inputs:
%    fnameC    - STRING name of file containing [T,C] data in SI units.
%    fnameL    - STRING name of file containing [T,L] data in SI units.
%    T0        - Designated temperature in Kelvin.
%    directory - STRING directory where files are kept.
%
% Outputs:
%    C_T0 - Volumetric heat capacity at temperature T in SI units.
%    L_T0 - Thermal conductivity at temperature T in SI units.
%    CC   - Two-column [T,C] vector in SI units.
%    LL   - Two-column [T,L] vector in SI units.
%
% Example: 
%    [C_T0 L_T0] = TDTR_TDEP_general(fnameC, fnameL, 80);
%      % gets C(80K) and L(80K) assuming reference files are in current
%      % directory.
%    [~,~,CC_SiO2,LL_SiO2] = TDTR_TDEP_general('SiO2_C.txt','SiO2_L.txt');
%      % imports SiO2 reference files from current directory.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author: Gregory Hohensee
% University of Illinois Urbana-Champaign
% email: hohense2@illinois.edu
% Website: n/a
% Mar 2014; Last revision: 26-Mar-2014
%                          14-July-2014: vH2. No change.
%                          17-Feb-2015: vH3. No change.
%------------- BEGIN CODE --------------
if nargin < 4, directory = pwd; end
if nargin < 3, T0 = 295; end % default temperature.
if nargin < 2, error('Reference values for L(T) and/or C(T) not specified.'); end

pathC = strcat(directory, '/', fnameC);
pathL = strcat(directory, '/', fnameL);

% import (interpolated) reference files for C and L.
% Format MUST be two columns of real numbers: temperature in K,
% and heat capacity in SI units.
CC = dlmread(pathC);
LL = dlmread(pathL);

if max(CC(:,2)) > 0 && max(CC(:,2)) < 10  % the only conceivable mistake that would imply a volumetric
                              % heat capacity of less than "10 units" is if someone wrote 
                              % their heat capacity data in J/cm^3-K instead of J/m^3-K.
                              % I assume nothing has a heat capacity of >10 J/cm^3-K.
                              % I assume nothing has a heat capacity of <10 J/m^3-K.
      CC(:,2) = CC(:,2)*1e6; % Convert to SI units.
      fprintf('WARNING: TDEP has converted a heat capacity reference from cm to SI units.\n');
end

% get indices at temperature T
[~,iC] = min(abs(CC(:,1) - T0));
[~,iL] = min(abs(LL(:,1) - T0));

% get values at temperature T
C_T0 = CC(iC,2);
L_T0 = LL(iL,2);

end