function [sol,nfinal]=rombint_VV3(f,a,b,NN)

relerr=1;
abserr=1;
limit=1e-5;
nmax=25;
mmax=25;
R=zeros(NN,nmax+1,mmax+1);

n=0;
m=0;
R(1:NN,n+1,m+1)=0.5*(b-a)*(f(a)+f(b));
while max(relerr)>limit
    n=n+1;
    h=(b-a)/(2^(n));
    k=1:2^(n-1);
    x=a*ones(1,length(k))+h*(2*k-1);
    feval=(f(x))';
    R(1:NN,n+1,1)=0.5*R(1:NN,n,1)+h.*sum(feval,2);
    for m=1:n
       fourm=4^m;
       R(1:NN,n+1,m+1)=1/(fourm-1)*(fourm*R(1:NN,n+1,m)-R(1:NN,n,m)); 
    end
    R(1:NN,n+1,m+1);
    relerr=abs(1-R(1:NN,n,m)./R(1:NN,n+1,m+1));
    abserr=abs(R(1:NN,n+1,m+1)-R(1:NN,n,m));
    if n==nmax
        fprintf('max iterations reached')
        sol=R(1:NN,n+1,m+1);
        nfinal=n;
        return
    end
end
sol=R(1:NN,n+1,m+1);
nfinal=n;
end