function [deltaR,ratio]=TDTR_REFL_vH3(tdelay,matparams,sysparams,A_pump,...
                                     intscheme,nnodes,offset)
%TDTR_REFL_vH - Calculates the Reflectivity Signal and Ratio
%In order to speed up the code, it is parallelized...the convention is...
%tdelay is COLUMN vector of desired delay times
%Mvect (the fourier components) are ROW vectors
%Matrices have size, length(tdelay) x length(Mvect)
%
% Syntax:  [deltaR,ratio]=TDTR_REFL_vH3(tdelay,matparams,sysparams,A_pump,...
%                                     intscheme,nnodes,Fortran)
%
% Inputs (scalars unless told otherwise):  
%    tdelay  - COLUMN VECTOR of time delay data
%    A_pump  - pump intensity (W)...doesn't effect ratio
% SYSPARAMS %
%    tau_rep - laser repetition time (sec) = 1/f_rep
%    f       - modulation frequency 
%    r_pump  - pump 1/e2 radius (m)
%    r_probe - probe 1/e2 radius (m)
% MATPARAMS %
%    LCTE    - vertcat(lambda,C,t,eta);
%     lambda  - VECTOR of thermal conductivities (W/mK) (layer 1 is the topmost layer)
%     C       - VECTOR of specific heats (J/m3-K)
%     t       - VECTOR of layer thicknesses (last layer is alway treated semi-inf, but
%               you still have to enter something).  (m)
%     eta     - VECTOR of anisotropic ratio (kx/ky), use ones(length(lambda)) for
%               isotropic
%    BI      - TRUE if bidirectional heat flow is present.
%    n_toplayer - Specifies number of layers above the plane where
%                 heat is deposited in the bidirectional heat flow model.
%    TCR     - temperature coefficient of reflectivity  (1/K)...doesn't affect ratio
%    doughnut - TRUE if this is a beam offset calculation. Can combine with
%               BI for bidirectional beam offset.
% OTHER %
%    intscheme - Numerical integration algorithm.
%                0 = Legendre-Gauss, 1 = rombint, 2 = Simpson.
%    nnodes   - For Legendre-Gauss integration, number of nodes. Default
%               should be at least 35.
%    Fortran  - if TRUE, calls FORTRAN executable instead of MATLAB thermal
%               model. May run faster or slower.
%    offset   - 1D vector of pump-probe relative offset distances (MICRONS)
%
% Outputs:
%    deltaR - Complex number VECTOR. Real part is model V(in), imaginary
%             part is model V(out).
%    ratio  - model -V(in)/V(out).
%
% Other m-files required: TDTR_TEMP_vH3.m
% Subfunctions: none
% MAT-files required: none
%
% See also: The TDTR_V4 package.

% Author: Gregory Hohensee / built from Joseph Feser's TDTR_REFL_V4.
% U. of Illinois Urbana-Champaign
% email: hohense2@illinois.edu
% Website: n/a
% Revision history: 12-Sep-2012 - J. Feser's REFL_V4 published
%                   2013 - added bidirectional heat flow, REFL_V4B.
%                   25-Mar-2014 - standardized comments into TDTR_vH package.
%                   26-Mar-2014 - use LCTE instead of lambda,C,t,eta. 
%                                 Allow various integration schemes.
%                   8-Apr-2014 - harmonized with TTM version
%                   14-July-2014 - vH2. No change.
%                   13-Nov-2014 - Merged with beam offset code.
%                   21-Jan-2015 - accepts offset as vector.
%                   28-Jan-2015 - error messages now block rombint and
%                                 Simpson int from choking on a variable
%                                 pump spot size. Only works with L-G int.
%                   17-Feb-2015 - vH3.
%------------- BEGIN CODE --------------
%% check input
PL = [0 length(sysparams) 0 length(matparams) 0];

% pump-probe system parameters
%sysparams = {tau_rep f r_pump r_probe};
if PL(2) < 1, tau_rep = 1/80e6; 
    warning('Defaulting to 80 MHz rep. rate'); else tau_rep = sysparams{1}; end
if PL(2) < 2, f = 9.8e6;
    warning('Defaulting to 9.8 MHz modulation'); else f = sysparams{2}; end
if PL(2) < 3, error('pump spot size not specified.'); else r_pump = sysparams{3}; end
if PL(2) < 4, r_probe = r_pump;
    warning('Defaulting to r_probe = r_pump'); else r_probe = sysparams{4}; end

% material parameters
%matparams = {LCTE aniso BI n_toplayer TCR};
if PL(4) < 1, error('LCTE not specified'); else LCTE = matparams{1}; end
%if pl(4) < 2, aniso = zeros(size(LCTE(4,:))); else aniso = matparams{2}; end
if PL(4) < 3, BI = 0; else BI = matparams{3}; end
if PL(4) < 4, if BI, error('n_toplayer not specified'); else n_toplayer = 0; end 
   else n_toplayer = matparams{4}; end
if PL(4) < 5, TCR = 1e-4; else TCR = matparams{5}; end
if PL(4) < 6, doughnut = 0; else doughnut = matparams{6}; end


if nargin < 7
    offset = 0; 
    if doughnut, warning('Assuming zero offset in REFL_vH3.'); end 
end
offset = offset*1e-6; % correct units for TEMP_DOUGHNUT

%if nargin < 7, Fortran = 0; end
if nargin < 6, nnodes = 35; end
if nargin < 5, intscheme = 0; end
if nargin < 4, A_pump = 10e-3; end % arbitrary positive value
if nargin < 3, error('Insufficient parameters for TDTR_REFL_vH3.'); end


% %% FORTRAN thermal model - TDTR_Bi is ISOTROPIC
% if BI && Fortran % change this and next if-statement if you change FORTRAN to non-bidirectional.
%     %%
%     directory= pwd; % Directory that contains the FORTRAN executable.
%     
%     % Change next five lines of code if using different FORTRAN thermal
%     % model executable. E.g., non-bidirectional, or with a repetition
%     % rate that is NOT 80.00 MHz.
%     filein=strcat(directory,'bifdv5.par');
%     fileTDTR=strcat(directory,'TDTR_Bi_NoText_80MHz_opt');
%     fileout=strcat(directory,'bifdv5.dat');
%     
%     lambda = LCTE(1,:);
%     C = LCTE(2,:);
%     t = LCTE(3,:);
%     M=[r_pump,f,n_toplayer;lambda',C',t'];
%     dlmwrite(filein,M,' '); %write param file
%     %%
%     setenv('DYLD_LIBRARY_PATH', '/opt/local/lib'); % Necessary for Mac OS X MATLAB to properly call the Mac-compiled Fortran executable.
%     unix(fileTDTR); % Calls the Fortran executable.
%     RAWDATA=dlmread(fileout); % reads the model output.
% 
%     dRin=interp1q(RAWDATA(:,1),RAWDATA(:,2),tdelay);
%     dRout=interp1q(RAWDATA(:,1),RAWDATA(:,3),tdelay);
%     deltaR=dRin+sqrt(-1)*dRout;
%     ratio=interp1q(RAWDATA(:,1),RAWDATA(:,4),tdelay);
% end

%% MATLAB port, bi-directional anisotropic (!) thermal model.
% That is, I see no reason why this can't also handle anisotropic cases.

    if intscheme == 2 % Simpson integration
        % implement identically to RW scripts, which use Simpson with
        % nnodes = 35.
        cap1 = 15;
        cap2 = 2;
    else % go with JPF / prior factors.
        cap1 = 10;
        cap2 = 1.5;
    end
    
    % RW uses cap1 = 15 in his REFL/TEMP scripts with Simpson
    % integration. More accuracy in integration at the cost of speed.
    fmax=cap1/min(abs(tdelay)); %maximum frequency considered (see Cahill 2004 RSI paper)
    ii=sqrt(-1);
    
    M=cap1*ceil(tau_rep/min(abs(tdelay))); %Highest Fourier component considered
    mvect=-M:M; %Range of Fourier components to consider (Vector)
    fudge1=exp(-pi*((mvect/tau_rep+f)/fmax).^2);%artificial decay (see RSI paper)
    fudge2=exp(-pi*((mvect/tau_rep-f)/fmax).^2);
    
    kmin = 0;
    kmax=cap2/sqrt(min(r_pump)^2+r_probe^2);
    % RW uses cap2 = 2 instead of 1.5 for his REFL/TEMP scripts,
    % which use Simpson integration exclusively. 2x catches higher k-modes,
    % is more accurate assuming comparable node density over k-space.
    % For TTM, k(RW) = 4*pi^2*k(JPF), where k(JPF) is in use in TEMP_V4.
    
    if BI
        dT1u = zeros(1,length(mvect))';
        dT2u = zeros(1,length(mvect))';

        dT1d = zeros(1,length(mvect))';
        dT2d = zeros(1,length(mvect))';
    end

    dT1=zeros(1,length(mvect))';
    dT2=zeros(1,length(mvect))';
    
% loop through pump-probe offset positions, because TDTR_TEMP_DOUGHNUT and
% associated Hankel calculations aren't parallelized in offset. Not sure
% how to do that, and we're only looping ~10 elements in offset anyway.
deltaR = zeros(length(tdelay),length(offset));
ratio = deltaR;
for jj = 1:length(offset)
    xo = offset(jj); % still fine if offset = 0 (time-delay TDTR)
    
    if BI
        totlayers = length(LCTE(1,:));

        % Split the material parameters into "up" and "down" sections,
        % relative to the heater and thermometer at interface of n_toplayer
        % and n_toplayer+1.
        for j = n_toplayer:-1:1
            LCTEu(:,n_toplayer+1-j) = LCTE(:,j);
        end
        LCTEd = LCTE(:,n_toplayer+1:totlayers);
        
        switch intscheme
            case 0 % Legendre-Gauss, allows variable pump spot size.
                [kvect,weights]=lgwt_V4(nnodes,kmin,kmax); %computes weights and node locations...
                % calculate the heat flow up, layer n_toplayer up to layer 1
                if doughnut
                    I1u = TDTR_TEMP_DOUGHNUT_vH3(kvect,mvect/tau_rep+f,LCTEu,r_pump,r_probe,A_pump,xo);
                    I2u = TDTR_TEMP_DOUGHNUT_vH3(kvect,mvect/tau_rep-f,LCTEu,r_pump,r_probe,A_pump,xo);
                else
                    I1u = TDTR_TEMP_vH3(kvect,mvect/tau_rep+f,LCTEu,r_pump,r_probe,A_pump);
                    I2u = TDTR_TEMP_vH3(kvect,mvect/tau_rep-f,LCTEu,r_pump,r_probe,A_pump);
                end
                
                [Nk,Nf,Nt] = size(I1u); % get 3D matrix dimensions
                
                dT1u = reshape(weights' * I1u(:,:), 1, Nf, Nt);
                dT2u = reshape(weights' * I2u(:,:), 1, Nf, Nt);
                
                % calculate the heat flow down, layer n_toplayer+1 down to bottom
                if doughnut
                    I1d = TDTR_TEMP_DOUGHNUT_vH3(kvect,mvect/tau_rep+f,LCTEd,r_pump,r_probe,A_pump,xo);
                    I2d = TDTR_TEMP_DOUGHNUT_vH3(kvect,mvect/tau_rep-f,LCTEd,r_pump,r_probe,A_pump,xo);
                else
                    I1d = TDTR_TEMP_vH3(kvect,mvect/tau_rep+f,LCTEd,r_pump,r_probe,A_pump);
                    I2d = TDTR_TEMP_vH3(kvect,mvect/tau_rep-f,LCTEd,r_pump,r_probe,A_pump);
                end
                
                dT1d = reshape(weights' * I1d(:,:), 1, Nf, Nt);
                dT2d = reshape(weights' * I2d(:,:), 1, Nf, Nt);
            case 1 % rombint [had a +/-f typo here before Nov. 13, 2014.]
                if length(r_pump)>1, error('Sorry, variable pump spot size is not yet compatible with rombint.'); end
                if doughnut
                    dT1u=rombint_VV3(@(kvect) TDTR_TEMP_DOUGHNUT_vH3(kvect,mvect/tau_rep+f,LCTEu,r_pump,r_probe,A_pump,offset),0,kmax,length(mvect));
                    dT2u=rombint_VV3(@(kvect) TDTR_TEMP_DOUGHNUT_vH3(kvect,mvect/tau_rep-f,LCTEu,r_pump,r_probe,A_pump,offset),0,kmax,length(mvect));
                    dT1d=rombint_VV3(@(kvect) TDTR_TEMP_DOUGHNUT_vH3(kvect,mvect/tau_rep+f,LCTEd,r_pump,r_probe,A_pump,offset),0,kmax,length(mvect));
                    dT2d=rombint_VV3(@(kvect) TDTR_TEMP_DOUGHNUT_vH3(kvect,mvect/tau_rep-f,LCTEd,r_pump,r_probe,A_pump,offset),0,kmax,length(mvect));
                else
                    dT1u=rombint_VV3(@(kvect) TDTR_TEMP_vH3(kvect,mvect/tau_rep+f,LCTEu,r_pump,r_probe,A_pump),0,kmax,length(mvect));
                    dT2u=rombint_VV3(@(kvect) TDTR_TEMP_vH3(kvect,mvect/tau_rep-f,LCTEu,r_pump,r_probe,A_pump),0,kmax,length(mvect));
                    dT1d=rombint_VV3(@(kvect) TDTR_TEMP_vH3(kvect,mvect/tau_rep+f,LCTEd,r_pump,r_probe,A_pump),0,kmax,length(mvect));
                    dT2d=rombint_VV3(@(kvect) TDTR_TEMP_vH3(kvect,mvect/tau_rep-f,LCTEd,r_pump,r_probe,A_pump),0,kmax,length(mvect));
                end
            case 2 % Simpson integration
                if length(r_pump)>1, error('Sorry, variable pump spot size is not yet compatible with Simpson int.'); end
                kvect = linspace(0,kmax,nnodes)';
                if doughnut
                    I1u = TDTR_TEMP_DOUGHNUT_vH3(kvect,mvect/tau_rep+f,LCTEu,r_pump,r_probe,A_pump,xo);
                    I2u = TDTR_TEMP_DOUGHNUT_vH3(kvect,mvect/tau_rep-f,LCTEu,r_pump,r_probe,A_pump,xo);
                    I1d = TDTR_TEMP_DOUGHNUT_vH3(kvect,mvect/tau_rep+f,LCTEd,r_pump,r_probe,A_pump,xo);
                    I2d = TDTR_TEMP_DOUGHNUT_vH3(kvect,mvect/tau_rep-f,LCTEd,r_pump,r_probe,A_pump,xo);
                else
                    I1u = TDTR_TEMP_vH3(kvect,mvect/tau_rep+f,LCTEu,r_pump,r_probe,A_pump);
                    I2u = TDTR_TEMP_vH3(kvect,mvect/tau_rep-f,LCTEu,r_pump,r_probe,A_pump);
                    I1d = TDTR_TEMP_vH3(kvect,mvect/tau_rep+f,LCTEd,r_pump,r_probe,A_pump);
                    I2d = TDTR_TEMP_vH3(kvect,mvect/tau_rep-f,LCTEd,r_pump,r_probe,A_pump);
                end
                dT1u=SimpsonInt(kvect,I1u);
                dT2u=SimpsonInt(kvect,I2u);
                dT1d=SimpsonInt(kvect,I1d);
                dT2d=SimpsonInt(kvect,I2d);
            otherwise
                error('integration scheme not properly specified. See analyze template.');
        end
        % make the parallel sum of the temperatures for up and down
        dT1 = dT1u .* dT1d ./ (dT1u + dT1d);
        dT2 = dT2u .* dT2d ./ (dT2u + dT2d);
    else % BI = 0, not bidirectional.    
        switch intscheme
            case 0 % Legendre-Gauss [compatible with time delay dependent spot size]
                [kvect,weights]=lgwt_V4(nnodes,kmin,kmax); %computes weights and node locations...
                
                if doughnut % TRUE if this is a beam offset measurement
                    I1 = TDTR_TEMP_DOUGHNUT_vH3(kvect,mvect/tau_rep+f,LCTE,r_pump,r_probe,A_pump,xo);
                    I2 = TDTR_TEMP_DOUGHNUT_vH3(kvect,mvect/tau_rep-f,LCTE,r_pump,r_probe,A_pump,xo);
                else
                    I1 = TDTR_TEMP_vH3(kvect,mvect/tau_rep+f,LCTE,r_pump,r_probe,A_pump);
                    I2 = TDTR_TEMP_vH3(kvect,mvect/tau_rep-f,LCTE,r_pump,r_probe,A_pump);
                end
        
                %dT1 = weights'*I1;
                %dT2 = weights'*I2;
                [Nk,Nf,Nt] = size(I1); % get 3D matrix dimensions
                
                dT1 = reshape(weights' * I1(:,:), 1, Nf, Nt);
                dT2 = reshape(weights' * I2(:,:), 1, Nf, Nt);
            case 1 % rombint
                if length(r_pump)>1, error('Sorry, variable pump spot size is not yet compatible with rombint.'); end
                if doughnut
                    dT1=rombint_VV3(@(kvect) TDTR_TEMP_DOUGHNUT_vH3(kvect,mvect/tau_rep+f,LCTE,r_pump,r_probe,A_pump,offset),0,kmax,length(mvect));
                    dT2=rombint_VV3(@(kvect) TDTR_TEMP_DOUGHNUT_vH3(kvect,mvect/tau_rep+f,LCTE,r_pump,r_probe,A_pump,offset),0,kmax,length(mvect));
                else
                    dT1=rombint_VV3(@(kvect) TDTR_TEMP_vH3(kvect,mvect/tau_rep+f,LCTE,r_pump,r_probe,A_pump),0,kmax,length(mvect));
                    dT2=rombint_VV3(@(kvect) TDTR_TEMP_vH3(kvect,mvect/tau_rep+f,LCTE,r_pump,r_probe,A_pump),0,kmax,length(mvect));
                end
            case 2 % Simpson
                if length(r_pump)>1, error('Sorry, variable pump spot size is not yet compatible with Simpson int.'); end
                kvect = linspace(0,kmax,nnodes)';
                if doughnut
                    I1 = TDTR_TEMP_DOUGHNUT_vH3(kvect,mvect/tau_rep+f,LCTE,r_pump,r_probe,A_pump,xo);
                    I2 = TDTR_TEMP_DOUGHNUT_vH3(kvect,mvect/tau_rep-f,LCTE,r_pump,r_probe,A_pump,xo);
                else
                    I1 = TDTR_TEMP_vH3(kvect,mvect/tau_rep+f,LCTE,r_pump,r_probe,A_pump);
                    I2 = TDTR_TEMP_vH3(kvect,mvect/tau_rep-f,LCTE,r_pump,r_probe,A_pump);
                end
                dT1=SimpsonInt(kvect,I1);
                dT2=SimpsonInt(kvect,I2);
            otherwise
                error('integration scheme not properly specified. See analyze template.');
        end
    end
    % now, Legendre-Gauss option creates dT1, dT2 of size [1,Nf,Nt],
    % where, without pump spot size variation over time delay,
    % dT1, dT2 were [1,Nf]. Nt is the length of r_pump, which
    % should be a scalar or the length of tdelay.
    
    expterm=exp(ii*2*pi/tau_rep*(tdelay*mvect)); % has size [Nt Nf]
    
    NNt = ones(length(tdelay),1);
    if length(r_pump) > 1
        % For simplicity let's convert [1 Nf Nt] to [Nt Nf],
        % as was the original format for Retemp and Imtemp.
        dT1 = permute(squeeze(dT1),[2,1]);
        dT2 = permute(squeeze(dT2),[2,1]);
        
        % NNt(Nt,1) * fudge(1,Nf) has size [Nt Nf].
        Retemp =     (dT1.*(NNt*fudge1)+dT2.*(NNt*fudge2)).*expterm;
        Imtemp = -ii*(dT1        -dT2        ).*expterm;
        
    else % r_pump is scalar (for beam offset, this is always true)
        % These are the original Retemp and Imtemp, for scalar r_pump
        Retemp=    (NNt*(dT1.*fudge1+dT2.*fudge2)).*expterm;
        Imtemp=-ii*(NNt*(dT1        -dT2        )).*expterm;
    end
    
    Resum=sum(Retemp,2); %Sum over all Fourier series components
    Imsum=sum(Imtemp,2);
    

    deltaRm=TCR*(Resum+ii*Imsum); %
    
    % columns are offset, rows are time delay.
    deltaR(:,jj) =deltaRm.*exp(ii*2*pi*f*tdelay); %Reflectance Fluxation (Complex)
    ratio(:,jj) =-real(deltaR)./imag(deltaR);
    
end % end for loop for offset positions
end