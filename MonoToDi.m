% Tranform a monopole (pressure) to dipole (vertical particel velocity)
%diploc = 'src' changes source to dipole, diploc = 'rcv' changes receiver
%vertical particle velocity
% Author: Aydin Shoja Email: se.mo.ay.sh@gmail.com; s.m.a.shoja@tudelft.nl

function [R,Fi,iKZ] = MonoToDi(G,nr,nztar,ns,dx,dt,dxmodel,nt,recieverdepth,density,velocity,pad,diploc,eps,TempInt)
%% Initialization

nrext= nr + 200;
ntpad = nt + pad;

Nyq_fq = (2*pi)/(2*dt);                                     % Nyquist frequency
Nyq_kx = (2*pi)/(2*dx);                                   % Nyquist wavenumber
df      = (2*pi)/(ntpad*dt);                         %Angular frequency sampling
dkx     = (2*pi)/(nrext*dx);                         %Wavenumber sampling in x-direction


if mod(ntpad,2) == 0 % Nt is even
    fvec = -Nyq_fq : df : Nyq_fq-df;
else            % Nt is odd
    fvec = [sort(-1*(df:df:Nyq_fq)) (0:df:Nyq_fq)];
end
if mod(nrext,2) == 0 % nr is even
    kvec = -Nyq_kx : dkx : Nyq_kx-dkx;
else            % nr is odd
    kvec = [sort(-1*(dkx:dkx:Nyq_kx)) (0:dkx:Nyq_kx)];
end

[~,mx] = size(velocity);
if mx == 1
    cp = velocity;
elseif recieverdepth <= dxmodel
    density = density(1,mx/2);
    cp = velocity(1,mx/2);
else
    cp = velocity(recieverdepth/dxmodel,mx/2);
    density = density(recieverdepth/dxmodel,mx/2);
end

[KX,omega]   = meshgrid(kvec,fvec);

% Create K^2 grid ((omega/c)^2)
K2 = (omega./cp).^2;
% Create KX^2 grid
KX2 = KX.^2;
% Get KZ^2 by using K^2=KX^2+kZ^2 | KZ^2=K^2-KX^2
KZ2 = K2 - KX2;
KZ2(KZ2<0) = 0;
KZ  = sqrt(KZ2);
iKZ = 1i.*KZ;

if TempInt == 1
denominator = (1i.*omega.*density).*conj(1i.*omega.*density);
Fi = 2*iKZ.*conj(1i.*omega.*density)./(denominator +eps*max(max(denominator)));
% Fi(Fi<0) = -1*Fi(Fi<0);
elseif TempInt == 2
denominator = (1i.*omega).*conj(1i.*omega);
Fi = conj(1i.*omega)./(denominator +eps*max(max(denominator)));
else
Fi = 2*iKZ./density;
% Fi(Fi<0) = -1*Fi(Fi<0);
end
% Fi(Fi<0) = -Fi(Fi<0);


if strcmp(diploc,'rcv')
    Fi = repmat(Fi,1,1);
    %     Fi = permute(Fi,[1 3 2]);
    %     Fi = [zeros(nt,100),Fi,zeros(nt,100)];
    %     G = [zeros(nt,100,ns),G,zeros(nt,100,ns)];
    R = zeros(size(G));
    for i = 1:nr
        G2 = squeeze(G(:,:,i));
        G2 = [zeros(nt,100),G2,zeros(nt,100)];
        Gfk = fftshift(fft2(G2));
        tmp = Gfk.*Fi;
        tmp = ifftshift(tmp);
        tmp = ifft(tmp,[],2);
        tmp = ifft(tmp,[],1,'symmetric');
        tmp = tmp(:,101:end-100,:);
        R(:,:,i) = tmp;
    end
    %     R = R(:,101:end-100,:);
elseif strcmp(diploc,'src')
    if size(G,2) > nr
        Fi = repmat(Fi,1,1);
        G = reshape(G,nt,nr,nztar,ns);
        %     Fi = permute(Fi,[1 3 2]);
        %     Fi = [zeros(nt,100),Fi,zeros(nt,100)];
        %     G = [zeros(nt,100,ns),G,zeros(nt,100,ns)];
        R = zeros(size(G));
       for j = 1:nztar
            temp = squeeze(G(:,:,j,:));
            
            for i = 1:ns
                G2 = squeeze(temp(:,:,i));
                G2 = [zeros(nt,100),G2,zeros(nt,100)];
                Gfk = fftshift(fft2(G2));
                tmp = Gfk.*Fi;
                tmp = ifftshift(tmp);
                tmp = ifft(tmp,[],2);
                tmp = ifft(tmp,[],1,'symmetric');
                tmp = tmp(:,101:end-100,:);
                R(:,:,j,i) = tmp;
            end
        end
    else
        Fi = repmat(Fi,1,1);
        R = zeros(size(G));
        for i = 1:ns
            G2 = squeeze(G(:,:,i));
            G2 = [zeros(nt,100),G2,zeros(nt,100)];
            Gfk = fftshift(fft2(G2));
            tmp = Gfk.*Fi;
            tmp = ifftshift(tmp);
            tmp = ifft(tmp,[],2);
            tmp = ifft(tmp,[],1,'symmetric');
            tmp = tmp(:,101:end-100,:);
            R(:,:,i) = tmp;
        end
    end
end
