% Multi-diemnsional convolution
% Author: Aydin Shoja Email: se.mo.ay.sh@gmail.com; s.m.a.shoja@tudelft.nl
% corr = 1 conjugates gather1, corr=2 conjugates gather2, else no
% conjugation
% trans = 'none' no transposition, trans = 'A' transposes gather1, trans =
% 'B' transposes gather2, trans='both' tranposes both
% flag = 'taper' apllies ataper, perc = taper precentage

function d = MDC(gather1,gather2,dt,dx,f,corr,trans,flag,perc)


[nt,nr,ns] = size(gather1);
[~,nr2,ns2] = size(gather2);
scl = dx;
nfft = nt;
Nyq_fq = 1/(2*dt);                                     % Nyquist frequency
df      = 1/(nfft*dt);                         %Angular frequency sampling


if mod(nfft,2) == 0 % Nt is even
    fvec = -Nyq_fq : df : Nyq_fq-df;
else            % Nt is odd
    fvec = [sort(-1*(df:df:Nyq_fq)) (0:df:Nyq_fq)];
end

fvec = fftshift(fvec);
[~,index] = min(abs(fvec-f));
nf = index;

if strcmp(flag,'taper')
    tmp = tukeywin(nr,perc);
    win1 = zeros(nr,ns);
    for i = -(ns-1)/2:(ns-1)/2
        k = 1+i+(ns-1)/2;
        win1(:,k) = circshift(tmp,i);
        if i <=0
            win1(ns+i:end,k) =0;
        elseif i>0
            win1(1:i,k) =0;
            
        end
    end
    win1 = repmat(win1,1,1,nt);
    win1 = permute(win1,[3 1 2]);
    gather1 = win1.*gather1;
    tmp = tukeywin(nr2,perc);
    win2 = zeros(nr2,ns2);
    for i = -(ns2-1)/2:(ns2-1)/2
        k = 1+i+(ns2-1)/2;
        win2(:,k) = circshift(tmp,i);
        if i <=0
            win2(ns2+i:end,k) =0;
        elseif i>0
            win2(1:i,k) =0;
            
        end
    end
    win2 = repmat(win2,1,1,nt);
    win2 = permute(win2,[3 1 2]);
    gather2 = win2.*gather2;
end

if corr==1
    gather1f = conj(fft(gather1,nfft,1));
    gather2f = fft(gather2,nfft,1);
elseif corr==2
    gather2f = conj(fft(gather2,nfft,1));
    gather1f = fft(gather1,nfft,1);
else
    gather1f = fft(gather1,nfft,1);
    gather2f = fft(gather2,nfft,1);
end

if strcmp(trans,'none')
    convD = zeros(nf,nr,ns2);
elseif strcmp(trans,'A')
    convD = zeros(nf,ns,ns2);
elseif strcmp(trans,'B')
    convD = zeros(nf,nr,nr2);
elseif strcmp(trans,'both')
    convD = zeros(nf,ns,nr2);
end
for i = 1:nf
    
    
    if strcmp(trans,'none')
        convD(i,:,:) = squeeze(gather1f(i,:,:))*squeeze(gather2f(i,:,:))*(scl);
    elseif strcmp(trans,'A')
        convD(i,:,:) = squeeze(gather1f(i,:,:)).'*squeeze(gather2f(i,:,:))*(scl);
    elseif strcmp(trans,'B')
        convD(i,:,:) = squeeze(gather1f(i,:,:))*squeeze(gather2f(i,:,:)).'*(scl);
    elseif strcmp(trans,'both')
        convD(i,:,:) = squeeze(gather1f(i,:,:)).'*squeeze(gather2f(i,:,:)).'*(scl);
    end
end

d = ifft(convD,nt,1,'symmetric');