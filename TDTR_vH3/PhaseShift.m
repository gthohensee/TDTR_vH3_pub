%% Simple subfunction: execute given phase shift on V(out), V(in).
function [res, delphase, Vin_shifted_A, Vout_shifted_A, ...
                         Vin_shifted_B, Vout_shifted_B] ...
                = PhaseShift(phase,ishort,Vin_shifted_A, Vout_shifted_A, ...
                                          Vin_shifted_B, Vout_shifted_B)
    radphase=pi/180*phase; % in radians
    
    
    % Complex representation of phase shift
    VB=(Vin_shifted_B+sqrt(-1)*Vout_shifted_B)*exp(sqrt(-1)*radphase);
    VA=(Vin_shifted_A+sqrt(-1)*Vout_shifted_A)*exp(sqrt(-1)*radphase);
    
    % update V(in), V(out) for phase shift
    Vin_shifted_A=real(VA);
    Vout_shifted_A=imag(VA);
    Vin_shifted_B=real(VB);
    Vout_shifted_B=imag(VB);
    
    % Calculate stats for the shifted V(out)
    res=abs(mean(Vout_shifted_A(1:ishort))-mean(Vout_shifted_B));
    noiseY=std([Vout_shifted_A(1:ishort);Vout_shifted_B]);
    delX=max(Vin_shifted_A(1:ishort))-min(Vin_shifted_B);
    delphase=2*noiseY/delX*180/pi;
end