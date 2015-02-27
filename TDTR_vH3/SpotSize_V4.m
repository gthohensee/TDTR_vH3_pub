%The main program tries to minimize "Z" by optimizing the variable(s) X
%This program Lets you:
%   1) Define the vector X: example, if lambda(3) is what you want to solve for,
%   then set X=lambda(3)....if you whish to simulatenous solve for more than one
%   variable, you can just define multiple variables (eg. X(1)=lambda(3), X(2)=lambda(4))
%
%   2) Define the fit.  Typically, this is the sum of the squares of the
%   residuals, but you might want to weight these by the sensitivity,
%   particularly if you don't intend to calculate the errorbars!

function [Z,A,mid,w0]=SpotSize_V4(X,y,x,A,mid,w0)

%Define the variables to be fit
w0 = X(1);
A = X(2);
mid= X(3);

y2=gaussian(A,mid,x,w0);
%Uncomment the next three lines to see the non-linear optimization in
%action!

figure(10)
plot(x,y,'ok',x,y2,'r');
pause(0.1)

% X
res=(y-y2).^2;
Z=sum(res);