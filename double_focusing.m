% This function computes the double-focused wavefields
% perc = taper precentage
% Author: Aydin Shoja Email: se.mo.ay.sh@gmail.com; s.m.a.shoja@tudelft.nl

function [Gamma_r,Source,Gmf1p,df1p] = double_focusing(f1p,Gplus,Gmin,Gd,NR,nztar,NS,NFocus,dxR,dt,dxmodel,nt,acquisitiondepth,focusingdepth,density,velocity,pad,perc)
% velocity(1,1) = velocity of the actual source location, velocity(2,1) =
% velocity of the virtual source location. Same is true for the density.
% Deconvelve the G+ and Transmission with the source signature before give
% them as input

[df1p,Fi,iKZ] = MonoToDi(f1p,NFocus,nztar,NR,dxR,dt,dxmodel,nt,acquisitiondepth,density(1,1),velocity(1,1),pad,'rcv',1e-5,1);

Source = MDC(Gplus,df1p,dt,dxR,100,0,'none','notaper',perc);


Gdinv = real(ifft(conj(fft(Gd))));
Gamma_r = MDC(Gdinv,Gd,dt,dxR,100,0,'both','taper',perc);
Gamma_r = ifft(conj(fft(Gamma_r,[],1)),[],'symmetric');

Gmf1p = MDC(Gmin,df1p,dt,dxR,100,0,'none','taper',0.4);


