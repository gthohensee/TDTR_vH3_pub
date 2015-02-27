function [LCTEtext] = genLCTEtext(LCTEsol,XsolIJ,stack)
%genLCTEtext - Formats LCTE solution XsolIJ with (stack) layer IDs into a
%text block with a format similar to the standard .par file for the
%original FORTRAN TDTR thermal models.
%
% Syntax:  [LCTEtext] = genLCTEtext(LCTEsol,XsolIJ,stack)
%
% Inputs:
%    LCTEsol - 4xN matrix of thermal parameters: thermal conductivity, heat
%              capacity, thicknesses, and "eta" of each model layer. SI units.
%    XsolIJ  - Mx3 matrix. First column is a fit parameter. Second and third
%              columns are the (i,j) position in LCTsol of that fit
%              parameter. Each row is for a different fit parameter; there 
%              are usually M = 2 rows.
%    stack   - 1xN string array of model layer identifiers. See
%              analyze_yymmdd_template.m and writeLCT.m.
%
% Outputs:
%    LCTEtext - Array of text formatted for printing onto a MATLAB figure
%               text box. Intended for a fit results figure, showing the
%               data, the fitted curve, and the thermal parameters.
%
% Example: 
%    [LCTEtext] = genLCTEtext(LCTEsol,XsolIJ,stack)
%      LCTEstr = 
% 
%         'L_x/L_z(eta),[W/cm-K, J/cm^3-K, nm]...layerID'
%         '1,[7 1.81 1e+05]...Diamond'
%         '1,[0.004 0.1 1]...Pb/Diamond'
%         '1,[32.8 152 1]...*'
%         '1,[0.328 1.52 73.4]...Pb'
%         '1,[0.0018 1.2 3e+04]...Soil_1cSt'
%         'Fit to L(Diamond)=7 W/cm-K and G(Pb/Diamond)=400 MW/m^2-K'
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: TDTR_MAIN_vH1.m, analyze_yymmdd_template.m, writeLCTE_vH1.m

% Author: Gregory Hohensee
% University of Illinois Urbana-Champaign
% email: hohense2@illinois.edu
% Website: n/a
% Revision history:  25-Mar-2014
%                    8-Apr-2014
% TO DO: incorporate system parameters (f, r_pump, tau_rep, TDTR-1/2,...)

%% ------------- BEGIN CODE --------------
fitstr = 'Fit to';
for jj = 1:length(XsolIJ(:,1)) % iterate through fit parameters
    if jj > 1
        fitstr = strcat(fitstr, ' and');
    end

    LCTEsol(XsolIJ(jj,2),XsolIJ(jj,3)) = XsolIJ(jj,1);

    switch XsolIJ(jj,2)
        case 1 % thermal conductivity row of LCT
            if ~isempty(strfind(stack{XsolIJ(jj,3)},'/')) % interface
                fitstr1 = ' G(';
                units = ' MW/m^2-K';
                unitfactor = 1e3;
            else
                fitstr1 = ' L(';
                units = ' W/cm-K';
                unitfactor = 1e-2;
            end
        case 2 % heat capacity row of LCT
            fitstr1 = ' C(';
            units = ' J/m^3K';
            unitfactor = 1e-6;
        case 3 % thickness row of LCT
            fitstr1 = ' t(';
            units = ' nm';
            unitfactor = 1e9;
        case 4
            fitstr1 = ' eta(';
            units = 'L_x/L_z';
            unitfactor = 1;
    end
    fitstr2 = stack{XsolIJ(jj,3)}; % get layer ID

    fitstr = strcat(fitstr, fitstr1, fitstr2, ')=', num2str(XsolIJ(jj,1)*unitfactor,3),units);
end

% Switch to standard .par format for Fortran TDTR. For visibility, use
% SI-cm units for L and C, and nm units for T.
LCTEvis = LCTEsol'; % standard format puts the layers into rows, not columns.
LCTEvis(:,1) = LCTEvis(:,1) * 1e-2; % to W/cm K
LCTEvis(:,2) = LCTEvis(:,2) * 1e-6; % from J/m^3K to J/cm^3K
LCTEvis(:,3) = LCTEvis(:,3) * 1e9;  % from m to nm

% Create a cell of strings from LCTvis
LCTEtext = cell(length(LCTEvis(:,1))+2,1);
LCTEtext{1} = strcat('[L_x/L_z, W/cm-K, J/cm^3-K, nm]...layerID');
for kk = 1:length(LCTEvis(:,1)) % iterating through the layers
    temp = strcat('[',mat2str(LCTEvis(kk,4),3), ', ', ...
                      mat2str(LCTEvis(kk,1),3), ', ',...
                      mat2str(LCTEvis(kk,2),3), ', ',...
                      mat2str(LCTEvis(kk,3),3), ']...', stack(kk));
    LCTEtext{kk+1} = temp{1};
end
%LCTEstr{length(LCTvis(:,1))+1+1} = '';
LCTEtext{length(LCTEvis(:,1))+1+1} = fitstr;
end

