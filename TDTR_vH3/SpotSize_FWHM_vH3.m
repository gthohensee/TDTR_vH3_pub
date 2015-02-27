%The main program tries to minimize "Z" by optimizing the variable(s) X
%This program Lets you:
%   1) Define the vector X: example, if lambda(3) is what you want to solve for,
%   then set X=lambda(3)....if you whish to simulatenous solve for more than one
%   variable, you can just define multiple variables (eg. X(1)=lambda(3), X(2)=lambda(4))
%
%   2) Define the fit.  Typically, this is the sum of the squares of the
%   residuals, but you might want to weight these by the sensitivity,
%   particularly if you don't intend to calculate the errorbars!

% SpotSize_FWHM_vH3 is a variation on SpotSize_V4, where the Gaussian
% fit is done on a sanitized input (normalized to 1, centered at 0).
function [Z,sigma]=SpotSize_FWHM_vH3(X,y,x)

%Define the variables to be fit
sigma = X;
A = 1;   % assume y(x) is normalized
mid = 0; % assume y(x) is simulated or centered

y2=Gaussian_vH3(A,mid,x,sigma);
%Uncomment the next three lines to see the non-linear optimization in
%action!

%figure(10)
%plot(x,y,'ok',x,y2,'r');
%pause(0.1)

% X
res=(y-y2).^2;
Z=sum(res);
end

% gaussian function can be a subfunction.
function ret = Gaussian_vH3(A,mid,z,sigma)
% the original "gaussian.m" script with a lower case "g" fitted for
% the correlated pump-probe spot size w0^2 = sigma^2+sigma^2, where sigma
% is the actual "width" or variance of the Gaussian distribution.
ret = A.*exp(-(z-mid).^2/(2*sigma^2));
end