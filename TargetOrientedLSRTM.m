% Target-oriented LSRTM with Marchenko double-focusing
% Author: Aydin Shoja Email: se.mo.ay.sh@gmail.com; s.m.a.shoja@tudelft.nl
function [numWorkers,time] = TargetOrientedLSRTM() % Dummy function for submitting to a cluster
% cleaning and clearing
close all; clear; clc;

cd /scratch/ashoja/TargetEnclosed/
% set paths (Change to your own paths)
addpath(genpath('/home/ashoja/seismic'));
addpath(genpath('/scratch/ashoja/TargetEnclosed'));
addpath(genpath('/home/ashoja/Redatumed_LSmigration/'));
addpath(genpath('/home/ashoja/Redatumed_LSmigration/Functions/'));
addpath(genpath('/home/ashoja/Redatumed_LSmigration/Functions/Born/scripts_peter'));
addpath(genpath('/home/ashoja/Redatumed_LSmigration/Functions/Born/scripts_joost'));
addpath(genpath('/home/ashoja/Redatumed_LSmigration/Functions/Born/scripts_extern'));
addpath(genpath('/scratch/ashoja/TargetEnclosed/SAGA_deep_2100_2600/'));

% set options
opt.plotconf    = 1; % plot configuration
opt.plotwavelet = 0; % plot wavelet
opt.scth        = 1; % scatter thickness
opt.specshot    = 161; % shot id to be plotted
opt.plotimage   = 1; % plot image
opt.srctype     = 1; % source type, 0 = monopole, 1 = dipole
opt.rectype     = 0; % receiver type,  0 = monopole, 1 = dipole
CG = 1;               % CG = 1 is Conjugate-Gradient algorithme otherwsie Steepest decent

% define chapters, 1 means the chapter is active
chapters = zeros(10,1);
chapters(1) = 1; % initialization
chapters(2) = 1; % 1 = loading data
chapters(3) = 1; % double focusing and making PSF
chapters(4) = 1; % make input list
chapters(5) = 1; % Green's function for migration with Born approximation of Lipmann-Schwinger integral
chapters(6) = 1; % LSM
chapters(7) = 1; % save data


% start chapter 1: spatial parameter setting
if chapters(1) == 1
    
    disp(['1. initialization']);   

% target computational grid
acquisitiondepth = 0;
TargetLateralOrigin = 7000; % Targets grid origin in m
focusingdeptha = 2100;      % targets upper boundary depth in m
focusingdepthb = 2600;      % targets lower boundary depth in m
nz1 = 1601;                 % original full model's number of depth sampels
nx1 = 6401;                 % original full model's number of lateral sampels
nxtar1 = 961;               % original target's number of latral sampels
dx1 = 3.125;                % original spatial sampling rate in m
nztar1 =  focusingdepthb/dx1 - focusingdeptha/dx1 +1;  % original target's number of depth sampels
dxmodel = 5;                                           % Reduced spatial sampling rate
mratio = dxmodel/dx1;                                  % Ratio of original and reduced spatial sampling rate
nz = ceil(nz1/mratio);                                 % Reduced full model's number of depth sampels
nx = ceil(nx1/mratio);                                 % Reduced full model's number of lateral sampels
input.dx = dxmodel;
nxtar=ceil(nxtar1/mratio);                             % Reduced target's number of latral sampels
nztar= focusingdepthb/dxmodel - focusingdeptha/dxmodel +1;  % Reduced target's number of depth sampels
perc = 0.2;

x1 = (0:nxtar)*dxmodel + TargetLateralOrigin;     
x2 = (0:nztar-1)*dxmodel + focusingdeptha;


% define your receiver geometry
 NR = 399;                  % number of receivers
 zR = 0;                    % Depth of receivers
 dxR = 25;                  % receiver spacing
 xR = zeros(2,NR);
 xR(1,1:NR) = (0:NR-1)*dxR;
 xR(2,1:NR) = zR; 
 
% define your source geometry  
 NS = 399;                  % number of sources
 zS = 0;                    % Depth of sources
 dxS = 25;                  % source spacing
 xS = zeros(2,NS);
 xS(1,1:NR) = (0:NR-1)*dxR;
 xS(2,1:NS) = zS;
 
% define target area grid
NS_tar = 242;              % number of virtual sources
dxF = 12.5;                % virtual sources spacing
x1_tar = (0:nxtar-1)*dxmodel + TargetLateralOrigin;     
x2_tar = (0:nztar-1)*dxmodel + focusingdeptha;  %-(nx+1)*dx/2+(1:nx)*dx;
[X1_tar,X2_tar] = ndgrid(x1_tar,x2_tar);

% define virtual source geometry
 zS_tar = focusingdeptha; 
 xS_tar = zeros(2,NS_tar);
 xS_tar(1,1:NS_tar) = (0:NS_tar-1)*dxF;
 xS_tar(2,1:NS_tar) = zS_tar*ones(1,NS_tar);  
 
% define virtual receiver geometry  
 NFocus = 242;             % number of virtual receivers
 zR_tar = focusingdeptha;
 xR_tar = zeros(2,NFocus);
 xR_tar(1,1:NFocus) = (0:NFocus-1)*dxF;
 xR_tar(2,1:NFocus) = zR_tar;

pad=0;
end

%% loading
if chapters(2) == 1  % Load your model and wavelet. The loading and reading functions and paths are just an example here.
    
    disp('Modle and temporal parameters')
    
    [AB,~] = ReadSumax('/scratch/ashoja/TargetEnclosed/SAGA_deep_2100_2600/saga_cp.su','cwp');  % Full smooth model
    AB = cell2mat(AB);
    AB = reshape(AB,nz1,nx1);
    AB = AB(1:2:end,1:2:end);
    
    [b,~] = ReadSumax('/scratch/ashoja/TargetEnclosed/SAGA_deep_2100_2600/target3_cp.su','cwp'); % Target smooth model
    b = cell2mat(b);
    b = reshape(b,nztar1,nxtar1);
    b = b(1:mratio:end,1:mratio:end);
    b2 = b';
     
    b_true = b;                                                                                  % True model if available
    CHI_mig = 1 - (b./b).^2;                                                                     % Initial perturbation model of the target
    s = (1./b).';                                                                                % Target's slowness
    CHI_tar = 1 - (b./b_true).^2;                                                                % True perturbation of the target if available
    CHI_tar = CHI_tar.';
    CHI_mig = CHI_mig.';
    velocity(1,1) = AB(1,1);
    velocity(2,1) = AB(focusingdeptha/dxmodel,1);
    velocity(3,1) = AB(focusingdepthb/dxmodel,1);
    
    
    [AB_den,~] = ReadSumax('/scratch/ashoja/TargetEnclosed/SAGA_deep_2100_2600/saga_ro.su','cwp');  % density of the entire medium
    AB_den = cell2mat(AB_den);
    AB_den = reshape(AB_den,nz1,nx1);
    AB_den = AB_den(1:mratio:end,1:mratio:end);
    
    
    density(1,1) = AB_den(1,1);
    density(2,1) = AB_den(focusingdeptha/dxmodel,1);
    density(3,1) = AB_den(focusingdepthb/dxmodel,1);
    
    [Wavelt,~] = ReadSumax('/scratch/ashoja/TargetEnclosed/SAGA_deep_2100_2600/Ricker30_zero_512.su','cwp'); % Modeling wavelet
    Wavelt = cell2mat(Wavelt).';
    Wavelt = [Wavelt(1:256),zeros(1,2002-512),Wavelt(257:end)];
    
    % temporal and frequency components
    nt = length(Wavelt);
    nt2 = 256;
    dt = 0.004;
    nfh = ceil(nt/2);
    Nyq_fq = (2*pi)/(2*dt);
    dw      = (2*pi)/(nt*dt);
    
    
    if mod(nt,2) == 0 % Nt is even
        fvec = -Nyq_fq : dw : Nyq_fq-dw;
    else            % Nt is odd
        fvec = [sort(-1*(dw:dw:Nyq_fq)) (0:dw:Nyq_fq)];
    end
    
   dWavelt = ifft(1i*fftshift(fvec).*fft(Wavelt),'symmetric')./nt; % Temporal derivative of the wavelet
   Wavelf = fft(Wavelt);
   dWavelf = fft(dWavelt)./nt;

   
disp('2. loading');  
% loading Marchenko Green's functions and data sets. Pay attention to the dimensions

[Gplus,~] = ReadSumax('/scratch/ashoja/TargetEnclosed/SAGA_deep_2100_2600/Marchenko/gplus/gplus_tot_2100.su','cwp'); % Marchenko downgoing Green's function at targets upper boundary
Gplus = reshape(cell2mat(Gplus),nt,NS,NFocus);
Gplus = permute(Gplus,[1 3 2]); % Sources ay surface should be the third dimension,
                                % and the focusing points should be in the second dimension, so I used the permute function.
Gplus = Gplus(:,1:end-1,:);

[Gmin,~] = ReadSumax('/scratch/ashoja/TargetEnclosed/SAGA_deep_2100_2600/Marchenko/gmin/gmin_tot_2100.su','cwp');    % Marchenko upgoing Green's function at targets upper boundary
Gmin = reshape(cell2mat(Gmin),nt,NS,NFocus);
Gmin = permute(Gmin,[1 3 2]);
Gmin = Gmin(:,1:end-1,:);

end

%% Duoble focusing
disp('Duoble focusing')
NR_surf = 399;
[Gd,~] = ReadSumax('/scratch/ashoja/TargetEnclosed/SAGA_deep_2100_2600/Trans_mono/Trans_direct.su','cwp');  % Direct arrival from surface to the target boundary
Gd = reshape(cell2mat(Gd),nt,[],NFocus);

[f1p,~] = ReadSumax('/scratch/ashoja/TargetEnclosed/SAGA_deep_2100_2600/Marchenko/f1plus/f1plus_tot1.su','cwp'); % Marchenko downgoing focusing function
f1p = reshape(cell2mat(f1p),nt,NS,NFocus);
f1p = permute(f1p,[1 3 2]);
Gd = Gd(:,:,1:end-1);  % I had an extra focusing point, hence the -1 one. You may not need it.
f1p = f1p(:,1:end-1,:);

%% double focusing and making PSF

[Gamma_r,Gpf1p,Gmf1p,df1p] = double_focusing(f1p(:,:,:),squeeze(Gplus(:,:,:)),squeeze(Gmin(:,:,:)),Gd,NR_surf,nztar,NS,NFocus-1,dxR,dt,dxmodel,nt,acquisitiondepth,focusingdeptha,density,velocity,pad,perc);
Gpf1p = permute(Gpf1p,[1 2 4 3]);

Gamma_r = deconvolution(Gamma_r,Wavelt,1e-3);

% Normalizing Gamma
for i = 1:NFocus-1
    
    m(:,i) = max(squeeze(abs(fft2(squeeze(Gamma_r(1,:,i))))));
    Gammafk = fft2(squeeze(Gamma_r(:,:,i)))./m(:,i);
    Gamma_r(:,:,i) = real(ifft(ifft(Gammafk,[],2),[],1));
end
Gamma_r = [Gamma_r(1:128,:,:);Gamma_r(1875:end,:,:)]; % Removing zeros from the middle of the gathers to reduce the number of the time sample.

% Virtual source computation
Source = squeeze(Gpf1p(:,:,:));
Source = deconvolution(squeeze(Source(:,:,:)),Wavelt,1e-3);
Source = [Source(1:128,:,:);Source(1875:end,:,:)];

Source= MonoToDi(squeeze(Source(:,:,:)),NFocus-1,nztar,NFocus-1,dxmodel,dt,dxmodel,nt2,acquisitiondepth,density(2,1),velocity(2,1),pad,'src',1e-5,1);
Source = permute(Source,[2 3 1]);
% Virtual observed data computation 
Gmf1p_w = deconvolution(squeeze(Gmf1p(:,:,:,1)),Wavelt,1e-3);
R_obs = squeeze(Gmf1p_w(:,:,:,1));
R_obs = permute(R_obs,[2 3 1]);
Gmf1p = permute(Gmf1p,[1 2 4 3]);
R_obs = double(R_obs);
R_obs = R_obs(:,:,1:256);

clear Gmin Gmin_decon Gplus Gplus_decon Gd Gd_decon f1p f1p_decon      % Deleting some variables to release memory

if chapters(4) == 1
disp(['4. make input list']); 
           

    % make input list for forward and inverse modeling function
    input.c_0 = b2(:);
    input.c_0_mig = b2(:);
    input.NR = NFocus-1;
    input.NS = NS_tar-1;
    input.xR = xR_tar;
    input.dxR = dxF;
    input.xS = xS_tar;
    input.dxS = dxF;
    input.N1 = nxtar;
    input.N2 = nztar;
    input.dx = dxmodel;
    input.X1 = X1_tar;
    input.X2 = X2_tar;
    input.fsamples  = nfh2;
    input.Wavelfrq = density(2)*2*Wavelf(1,1:nfh);
    input.Nfft = nf2;
    input.f = 120;
    input.df = 1/((nt2)*dt);
    input.nt = nt2;
    input.dt = dt;
    input.srctype = 0;
    input.rectype = 0;
    input.trans = 0;
    input.src = fft(Source,[],3);
    input.srcmulti = 1;
    input.iter = 20;
    input.Gpp = Source;
    input.Gamma = 0;
    input.psf = 0;
    input2=input;
    input2.c_0 = b(1,1);
    input.den = density;
    
    



end

if chapters(5) == 1
    disp(['5. Green function']);
    tic
    
        % Load Green's function of the entire target computed by the Finite difference algorithm

        [Gsu,~] = ReadSumax('/scratch/ashoja/TargetEnclosed/SAGA_deep_2100_2600/Green_tar/Green_tar_0.su','cwp');
        Gsu = reshape(cell2mat(Gsu),256,nxtar*(nztar),NFocus-1);
        size(Gsu)
        Gr = permute(Gsu,[3 2 1]);      % Source side Green's function
        Gsu = permute(Gsu,[2 3 1]);     % Using source-receiver reciprocity to create receiver side Green's functions
        Gsu = (fft(Gsu,[],3))*dt;
        Gr = (fft(Gr,[],3))*dt;
        
                
    toc
end



if chapters(6) == 1
disp(['6. Least squares imaging']); 
iter = 35;
CHI_mig = zeros(nxtar,nztar);
e = zeros(1,iter);
im_old = zeros((nztar)*nxtar,1);
p_old_vec = zeros((nztar)*nxtar,1);
itvec = 1:iter;

for it = 1:iter
    
    if it == 1
        res = -R_obs;
        R_pred = zeros(NFocus-1,NS_tar-1,nt2);
    else
        
        % compute predicted data
        tic
        [R_predu,~] = scat_int_parallel_TO(Gsu,Gr,CHI_mig(:),'notransp',input);
        toc
        R_predu = reshape(R_predu,[NFocus-1,NS_tar-1,nt2]);
        R_predu = permute(R_predu,[3 1 2]);
        R_predu = MDC(Gamma_r,R_predu,dt,dxR,100,0,'none','notaper',0.4);
        R_pred = R_predu;
        R_pred = permute(R_pred,[2 3 1]);
         
        res = R_pred - R_obs;
        
    end
    
    e(1,it) = res(:)'*res(:);
    
    figure(4);
    plot(itvec(1:end),e(1:end)./(R_obs(:)'*R_obs(:)));
    pause(0.1)
    print('-f4','costMO_2100','-dpng');
    
    
    % imaging
  
    input2 = input;
    input2.srcmulti =1;
    input2.Wavelfrq = density(2)*2*Wavelf(1,1:nfh);
    input2.src = fft(Source,[],3);
    input2.Gamma = fft(Gamma_r);
    tic
    [~,grad] = scat_int_parallel_TO(Gsu,Gr,res(:),'transp',input2);
    toc

    grad = -reshape(grad,[nxtar,nztar]);

    
    disp('CG step')
    if it == 1 || CG == 0
        beta = 0;
        
    else
        beta = (grad(:)'*(grad(:) - im_old(:)))/((grad(:) - im_old(:))'*p_old_vec); % Hestenes - Steifel
        % beta = (norm(grad_vec).^2)/(norm(grad_old_vec).^2);                       % Fletcher - Reeves
        % beta = (grad_vec'*(grad_vec - grad_old_vec))/(norm(grad_vec).^2);         % Polak - Rebier74
        % beta = (norm(grad_vec)^2) / (p_old_vec'*(grad_vec - grad_old_vec));       % Dai - Yuan
    end
    
    
    p_vec = grad(:) + beta*p_old_vec;
    p = reshape(p_vec,[nxtar,nztar]);
    
    im_old(:) = grad(:);
    p_old_vec = p_vec;
    
    
    
    disp(' Step length ')
    
    % step length
    alpha_t = 1e5;
    
    if it==1
        while max(max(abs(alpha_t*p(:)))) > max(max(abs(grad(:))))/10
            alpha_t = 0.9 * alpha_t;
        end
    else
        while max(max(abs(alpha_t*p(:)))) > max(max(abs(CHI_mig(:))))/10
            alpha_t = 0.9 * alpha_t;
        end
    end
    
    
    vv = CHI_mig + alpha_t*p;
    
    tic
    [d_cal_newu,~] = scat_int_parallel_TO(Gsu,Gr,vv(:),'notransp',input);
    toc
    d_cal_newu = reshape(d_cal_newu,[NFocus-1,NS_tar-1,nt2]);
    d_cal_newu = permute(d_cal_newu,[3 1 2]);
    d_cal_newu  = MDC(Gamma_r,d_cal_newu,dt,dxR,100,0,'none','notaper',0.4);
    d_cal_new = d_cal_newu;
    d_cal_new = permute(d_cal_new,[2 3 1]);
    r_d_cal =  -d_cal_new(:) + R_pred(:);
    alpha = real(((r_d_cal'*res(:))/(r_d_cal'*r_d_cal))*alpha_t)
    
    Alp(1,it) = alpha;
    P(:,:,it) = alpha*p;
    
    CHI_mig = CHI_mig + alpha*p;
    
    
    

    % start plot image
        if opt.plotimage == 1
            fig_image = figure(5);
            
            
            subplot(3,1,1);
            imagesc(x1,x2,imresize(real(P(:,:,1)),16).'); hold on;
            xlabel('x_1 (m)','fontsize',14);
            ylabel('x_2 (m)','fontsize',14);
            set(gca,'fontsize',14); axis ij;
            axis([min(x1(:)),max(x1(:)),min(x2(:)),max(x2(:))]);
            caxis([-0.5*max(max(P(:,:,1))),0.5*max(max(P(:,:,1)))]);
            colormap(flipud(colormap('gray')));
            colorbar
            
            subplot(3,1,2);
            imagesc(x1,x2,imresize(real(P(:,:,it)),16).'); hold on;
            xlabel('x_1 (m)','fontsize',14);
            ylabel('x_2 (m)','fontsize',14);
            set(gca,'fontsize',14); axis ij;
            axis([min(x1(:)),max(x1(:)),min(x2(:)),max(x2(:))]);
            caxis([-max(max(P(:,:,it))),max(max(P(:,:,it)))]);
            colormap(flipud(colormap('gray')));
            colorbar
                        
            subplot(3,1,3);
            imagesc(x1,x2,imresize(real(CHI_mig),16).'); hold on;
            xlabel('x_1 (m)','fontsize',14);
            ylabel('x_2 (m)','fontsize',14);
            set(gca,'fontsize',14); axis ij;
            axis([min(x1(:)),max(x1(:)),min(x2(:)),max(x2(:))]);
            caxis([-max(abs(CHI_mig(:))),max(abs(CHI_mig(:)))]);
            colormap(flipud(colormap('gray')));
            colorbar
            pause(0.1)
            print('-f5','imageMO_2100','-dpng');
            
            
           figure(6)
           subplot(1,3,1)
           imagesc(xR_tar(1,:),[0:(nt2)-1]*dt,imresize(squeeze(R_obs(:,opt.specshot,:)),16).'); hold on;
           xlabel('x_1 (m)','fontsize',14);
           ylabel('t (s)','fontsize',14);
           set(gca,'fontsize',14); axis ij;
           caxis([-max(max(R_obs(:,opt.specshot,:))),max(max(R_obs(:,opt.specshot,:)))]);
           colormap(flipud(colormap('promax')));
           title('Observed');
            
           subplot(1,3,2)
           imagesc(xR_tar(1,:),[0:(nt2)-1]*dt,imresize(squeeze(R_pred(:,opt.specshot,:)),16).'); hold on;
           xlabel('x_1 (m)','fontsize',14);
           ylabel('t (s)','fontsize',14);
           set(gca,'fontsize',14); axis ij;
           caxis([-max(max(R_obs(:,opt.specshot,:))),max(max(R_obs(:,opt.specshot,:)))]);
           colormap(flipud(colormap('promax')));
           title('Predicted');
            
           subplot(1,3,3)
           imagesc(xR_tar(1,:),[0:(nt2)-1]*dt,imresize(squeeze(res(:,opt.specshot,:)),16).'); hold on;
           xlabel('x_1 (m)','fontsize',14);
           ylabel('t (s)','fontsize',14);
           set(gca,'fontsize',14); axis ij;
           caxis([-max(max(R_obs(:,opt.specshot,:))),max(max(R_obs(:,opt.specshot,:)))]);
           colormap(flipud(colormap('promax')));
           title('Residual');
           pause(0.1)
           print('-f6','dataMO_2100','-dpng');
        end
    
end
end
if chapters(7)==1
      save('SAGA_TargetOriented_Marchenko_7000_2100_50iter','R_obs','R_pred','CHI_mig','P','e','input','Source');
end
