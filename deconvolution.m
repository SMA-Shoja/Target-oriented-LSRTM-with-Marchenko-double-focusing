% Wavelet deconvolution function
% Author: Aydin Shoja Email: se.mo.ay.sh@gmail.com; s.m.a.shoja@tudelft.nl

function [out,Winv] = deconvolution(in,wavelet,eps)

[t r s z] = size(in);

Wf = fft(wavelet).'./length(wavelet);
Wf2 = conj(Wf).*(Wf);
m = max(Wf2);
Winv = repmat(conj(Wf),1,r,s,z)./(repmat(Wf2,1,r,s,z)+m*eps); 

in_decon = fft(in).* Winv;
out = ifft(in_decon,[],1,'symmetric');