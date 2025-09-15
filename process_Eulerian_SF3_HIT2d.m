

%%%%%%%%%%%%%%%%%%%y direct

clear
addpath(genpath('D:\colorbar\'));
aaaa=dir(['E:\hitd\2dturbulence_w_alp01\HIT2D_t_*',])
tic

% load('two_dim_Parameters.mat')
addpath('E:\hitd\2dturbulence_w_alp01')
load HIT2D_Parameters

Nt=199;
rescale=2.*pi./512;
% kk(kk==0)=1e-10;
% mm(mm==0)=1e-10;

%%%%%%%%%% quantities used for angle average

xc=x-max(x)/2;
zc=z-max(z)/2;

dr=1*(xc(2)-xc(1));
rmax=max(xc);

RR=((xx-max(x)/2).^2+(zz-max(z)/2).^2).^0.5;
N=512;


for j=1:Nt

    load(['HIT2D_t_' num2str(4000+j-1) '.mat'])
    hpsi=Diag.*hq;
    hu1=-1i.*mm.*hpsi;
    hw1=1i.*kk.*hpsi;

    u=real(ifft2(hu1));
    v=real(ifft2(hw1));
    % [S3L_Roy1(j,:,:),S3T_Roy1(j,:,:)]=calc_SF3(u,w,N);
    [X, Y] = meshgrid(-N/2:N/2-1, -N/2:N/2-1);

    % for direct=1:size(u,1)
    %     u1=u(direct,:);
    %     v1=v(direct,:);
    % 
    %     u2=u(direct,:).^2;v2=v(direct,:).^2;uv=u(direct,:).*v(direct,:);

        % uh = fft(u1); %
        % vh = fft(v1);
        % uuh= fft(u2);
        % vvh= fft(v2);
        % uvh= fft(uv);
        % 
        % 
        % Cu_uu=fftshift(ifft(uh.*conj(uuh)./N.^2));
        % Cuu_u=fftshift(ifft(uuh.*conj(uh)./N.^2));
        % Cv_vv=fftshift(ifft(vh.*conj(vvh)./N.^2));
        % Cvv_v=fftshift(ifft(vvh.*conj(vh)./N.^2));
        % 
        % Cuv_u=fftshift(ifft(uvh.*conj(uh)./N.^2));
        % Cuv_v=fftshift(ifft(uvh.*conj(vh)./N.^2));
        % Cu_uv=fftshift(ifft(uh.*conj(uvh)./N.^2));
        % Cv_uv=fftshift(ifft(vh.*conj(uvh)./N.^2));
        % 
        % Cuu_v=fftshift(ifft(uuh.*conj(vh)./N.^2));
        % Cvv_u=fftshift(ifft(vvh.*conj(uh)./N.^2));
        % Cv_uu=fftshift(ifft(vh.*conj(uuh)./N.^2));
        % Cu_vv=fftshift(ifft(uh.*conj(vvh)./N.^2));
        % %
        % 
        % %%%%%%%%%%%%%% x direct
        % SLx(j,direct,:)=3.*Cu_uu-3.*Cuu_u;
        % STx(j,direct,:)=Cu_vv-Cvv_u-2.*Cuv_v+2.*Cv_uv;
    [x1,y1,S3L(j,:,:),S3T(j,:,:)]=test4_calc_SF3(u,v,N);
    % end
    disp(j)
end

[~, S3L1_iso] = calc_ispec2(x1(1,:), y1(:,1), nanmean(S3L,1), 2);
[rbin, S3T1_iso] = calc_ispec2(x1(1,:), y1(:,1), nanmean(S3T,1), 2);

rbin1=diag(rbin);
SF3=(S3L1_iso+S3T1_iso);

loglog(rbin1,abs(SF3));hold on
loglog(rbin1,-(SF3),'Marker','+');
loglog(rbin1,(SF3),'Marker','o');