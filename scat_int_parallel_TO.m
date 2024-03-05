% Computes forward and inverse Lipmann-Schewinger integral
% Based on code by Peter van den berg book's "Forward and Inverse Scattering Algorithms based on
% Contrast Source Integral Equations" codes.
% flag = 'notransp' is forward modeling, flag = 'transp' is backward
% propagation
% Modified by: Aydin Shoja Email: se.mo.ay.sh@gmail.com; s.m.a.shoja@tudelft.nl

function [yu,y, pertinc] = scat_int_parallel_TO(Gsu,Gr,x,flag,input)

df = input.df;
c_0 = input.c_0;
NS = input.NS;
dx = input.dx;
dxR = input.dxR;
src = input.src;
scl = (dxR)/input.den(2,1);
w = input.Wavelfrq;
fsamples = input.fsamples;

if strcmp(flag,'notransp')
    srcmulti = input.srcmulti;
    x = real(x);
    x = repmat(x,1,input.NS);
    yu = zeros(input.NR,input.NS,input.Nfft);
    pertinc = zeros(input.N1*input.N2,input.NS,input.Nfft);


    parfor fl = 2 : fsamples


        freq = (fl-1) * df;
        s = 1i*2*pi*freq;
        gamma_0 = s./c_0;
        gamma_0 = repmat(gamma_0,1,NS);

        if srcmulti == 0
            % compute contribution
            yu(:,:,fl) =  squeeze(Gr(:,:,fl)) *((gamma_0*dx).^2.*x.*squeeze(Gsu(:,:,fl)).*squeeze(w(1,fl)));

        elseif srcmulti == 1
            yu(:,:,fl) =  squeeze(Gr(:,:,fl)) *((gamma_0*dx).^2.*x.*(squeeze(Gsu(:,:,fl))*squeeze(src(:,:,fl))))*scl;
            pertinc(:,:,fl) = (squeeze(Gsu(:,:,fl))*squeeze(src(:,:,fl)))*scl;


        end
    end

    % inverse fourier transform
    yu = real(ifft((yu),[],3,'symmetric'));
    pertinc = real(ifft((pertinc),[],3,'symmetric'));
    y = 0;

% adjoint mode y = A'*x
elseif strcmp(flag,'transp')
    % reshape input
    x = reshape(x,[input.NR,input.NS,input.nt]);
    x = fft(x,[],3);
    Gamma = input.Gamma;
    srcmulti = input.srcmulti;
    % initiaze output
    y = zeros(input.N1,input.N2);
    y = y(:);
    % forward fourier transform


    parfor fl = 2 : fsamples

        % set parameters
        freq = (fl-1) * df;
        s = 1i*2*pi*freq;
        gamma_0 = s./c_0;

        % migration
        if srcmulti == 0

            Pinc =  squeeze(Gsu(:,:,fl)).*squeeze(w(1,fl));

            y = y + ( (gamma_0.* dx).^2 .* sum(((squeeze(Gr(:,:,fl))).' * conj(squeeze(x(:,:,fl)))) .* Pinc,2));

        elseif srcmulti == 1
            Pinc =  squeeze(Gsu(:,:,fl))*squeeze(src(:,:,fl))*scl;
            %             Gr(:,:,fl) = squeeze(Gamma(fl,:,:))*Gr(:,:,fl);

            y = y + ( ( gamma_0*dx).^2 .* sum(((squeeze(Gr(:,:,fl))).' * conj(squeeze(x(:,:,fl)))) .* Pinc,2));

        end
    end
    % end frequency loop

    % take real part of contrast
    y = real(y(:));
    yu = 0;


end
end

