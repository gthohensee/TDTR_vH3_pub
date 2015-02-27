function [QQ] = polyadd(P,S)
%Adds polynomial P&S: makes same size as maximum Size
NP = length(P);
NS = length(S);

N = max(NP,NS);

PP = zeros(1,N);
SS = zeros(1,N);

PP(end-NP+1:end)=P;
SS(end-NS+1:end)=S;

QQ = PP+SS;

end

