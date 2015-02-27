function Integrand=TDTR_TEMP_DOUGHNUT_vH3(kvectin,freq,LCTE,r_pump,r_probe,A_pump,xoffset)
%TDTR_TEMP_DOUGHNUT_vH3 - For beam offset measurements, computes frequency 
% domain average temperature response to periodic gaussian pump beam, 
% probed by another gaussian beam.
%
% This program is vectorized/optimized to handle all frequencies 
% simultaneously (f is a ROW vector)
%
% Syntax:  [Integrand,G]=TDTR_TEMP_vH3(kvectin,freq,LCTE,r_pump,r_probe,A_pump)
%
% Inputs:
%    kvectin - vector of wavenumber (m^-1)
%    freq    - excitation frequency (Hz), ROW vector
%    LCTE    - vertcat(lambda,C,t,eta)
%    lambda  - vector of thermal conductivities, 
%              lambda(1)=top surface,(W/m-K)
%    C       - vector of volumetric specific heat (J/m3-K)
%    t       - thicknesses of each layer (layer N will NOT be used, semiinfinite)
%    eta     - anisotropy of each layer.
%    r_pump  - Pump spot size (m)
%    r_probe - Probe spot size (m)
%    A_pump  - Pump power (W), used to ESTIMATE amplitude (not used for fitting)
%    xoffset - lateral relative offset of pump and probe beams.
%
% Outputs:
%    Integrand - G.*Kernal; see David Cahill's 2004 paper.
%    G         - The layer G(k), where k is spatial wavenumber
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: TDTR_REFL_vH3.m

% Author: Joseph P. Feser
% --
% email: --
% Website: TDTR_V4 package: http://users.mrl.illinois.edu/cahill/tcdata/tcdata.html
% September 2012; 
% Revision history: 13-November-2014: Merged with vH2 package. Renamed
%                                     EXTRA to xoffset.
%                   17-Feb-2014: vH3. No change.
%------------- BEGIN CODE --------------
%% Unpacked from my other data structures so I don't need to edit the code.
lambda = LCTE(1,:);
C = LCTE(2,:);
t = LCTE(3,:);
eta = LCTE(4,:);

%% Joe's code from here onward.
Nfreq=length(freq);
kvect=kvectin(:)*ones(1,Nfreq);
Nlayers=length(lambda); %# of layers
Nint=length(kvectin); %# of different frequencies to calculate for

%k is a COLUMN vector (actually a matrix that changes down the rows)
%f is a ROW vector

%kmax=1/sqrt(r_pump^2+r_probe^2)*1.5; %cutoff wavevector for integration
ii=sqrt(-1);
alpha=lambda./C;
omega=2*pi*freq;
q2=ones(Nint,1)*(ii*omega./alpha(Nlayers));
kvect2=kvect.^2;

un=sqrt(4*pi^2*eta(Nlayers)*kvect2+q2);
gamman=lambda(Nlayers)*un;
Bplus=zeros(Nint,Nfreq);
Bminus=ones(Nint,Nfreq);
kterm2=4*pi^2*kvect2;
if Nlayers~=1
    for n=Nlayers:-1:2
        q2=ones(Nint,1)*(ii*omega./alpha(n-1));
        unminus=sqrt(eta(n-1)*kterm2+q2);
        gammanminus=lambda(n-1)*unminus;
        AA=gammanminus+gamman;
        BB=gammanminus-gamman;
        temp1=AA.*Bplus+BB.*Bminus;
        temp2=BB.*Bplus+AA.*Bminus;
        expterm=exp(unminus*t(n-1));
        Bplus=(0.5./(gammanminus.*expterm)).*temp1;
        Bminus=0.5./(gammanminus).*expterm.*temp2;
        % These next 3 lines fix a numerical stability issue if one of the
        % layers is very thick or resistive;
        penetration_logic=logical(t(n-1)*abs(unminus)>100);  %if pentration is smaller than layer...set to semi-inf
        Bplus(penetration_logic)=0;
        Bminus(penetration_logic)=1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        un=unminus;
        gamman=gammanminus;
    end
end

G=(Bplus+Bminus)./(Bminus-Bplus)./gamman; %The layer G(k)
P=A_pump*exp(-pi^2*r_pump^2/2*kvect.^2);

%S=besselj(0,2*pi*r_probe*kvect); %Dirac ring
%xo = EXTRA;
S = GetHankel_OffsetBeam(kvect,xoffset,r_probe,1);
%S=exp(-pi^2*r_probe^2/2.*kvect.^2); %Normal Gaussian Probe
%S=(1-pi^2.*(kvect.^2)*r_probe^2).*exp(-pi^2*r_probe^2*kvect.^2); %The Gaussian Doughnut.  Yum
%S=1/(pi*(a^2-b^2))*(a*BesselJ(1,2*pi*a*kvect)./kvect-b*BesselJ(1,2*pi*b*kvect)./kvect);

Integrand=G.*P.*S.*kvect; %The rest of the integrand

        
