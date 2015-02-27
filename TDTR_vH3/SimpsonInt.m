%Simpson integration
%int(f(x),x=a..b)
%int=
function Integral=SimpsonInt(x,f)

N=length(x)-1;
Nints=length(f(1,:));
delta=(x(N+1)-x(1))/N;
%generate weighting factors
w=zeros(N+1,1);
w(1:2:N+1)=2;
w(2:2:N)=4;
w(1)=1;
w(N+1)=1;
warray=w*ones(1,Nints);
Integral=sum(warray.*f);
Integral=delta/3*Integral;
