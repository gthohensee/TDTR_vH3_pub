function ret = gaussian(A,mid,z,w0);
ret = A.*exp(-(z-mid).^2/w0^2);