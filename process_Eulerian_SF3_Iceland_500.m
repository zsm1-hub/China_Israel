clear all;close all
clc
%%%%%%%%%%%%%%%%calc Eulerian SF3 for Iceland%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('/meddy/simingzhang/Analysis/matlab/Parcels_SF/')
addpath('/meddy/simingzhang/Data/Parcels_data')
addpath(genpath('/meddy/simingzhang/Data/iceland_500m'))
input_dir='/meddy/simingzhang/Data/iceland_500m/';
Case='nowave_500';
nworkers=4;

if strcmpi(Case, 'wave_500')
    fname=[input_dir,'iceland_wave/','afterspinup_z_his_depth.0002.nc'];
end

if strcmpi(Case, 'nowave_500')
    fname=[input_dir,'iceland_nowave/','afterspinup_z_his_depth.0002.nc'];
end

gname=[input_dir,'sample_niskin_500m_grd.nc'];
dx=mean(mean(1./ncread(gname,'pm')));
dy=mean(mean(1./ncread(gname,'pn')));

xscale=mean([dx,dy]);
disp(['******the scale of grid mean value is ',num2str(xscale),'m ********']);

u=(squeeze(ncread(fname,'u')));
u1=u2rho_3d(permute(u,[3,2,1]));
v=(squeeze(ncread(fname,'v')));
v1=v2rho_3d(permute(v,[3,2,1]));
N=size(u1,2);
clear u;clear v;


if isempty(gcp('nocreate'))
    parpool('local', 4);
end
tic
parfor t = 1:size(u1,1)
    u2=squeeze(u1(t,:,:));
    v2=squeeze(v1(t,:,:));
    %
    % [S3L(t,:,:),S3T(t,:,:)]=test4_calc_SF3(u2,v2,N);
    % rescale=2;xscale1=xscale;
    % [~,~,S3L1_alltime(t,:),S3T1_alltime(t,:)]=calc_radial(S3L(t,:,:),S3T(t,:,:),N,xscale1);
    % disp(t)

    [S3L11,S3T11]=test4_calc_SF3(u2,v2,N);
    xscale1=xscale;
    [r,~,S3L1_alltime{t}(:),S3T1_alltime{t}(:)]=calc_radial(S3L11,S3T11,N,xscale1);
    disp(t)
end
delete(gcp('nocreate'));
toc

t=1
u2=squeeze(u1(t,:,:));
v2=squeeze(v1(t,:,:));

[S3L11,S3T11]=test4_calc_SF3(u2,v2,N);
xscale1=xscale;
[r,~,~,~]=calc_radial(S3L11,S3T11,N,xscale1);
% [~, S3L1_iso] = calc_ispec2(x1(1,:), y1(:,1), nanmean(S3L,1), 2);
% [rbin, S3T1_iso] = calc_ispec2(x1(1,:), y1(:,1), nanmean(S3T,1), 2);
% rbin1=diag(rbin);
% clear rbin
% dist_axis=rbin1(2:end).*xscale;
% SF3=(S3L1_iso(2:end)+S3T1_iso(2:end));
% outputname=[Case,'_Eulerian_SF3.mat']
% save(outputname,'SF3','dist_axis');
% 
% S3L1=squeeze(nanmean(S3L,1));
% S3T1=squeeze(nanmean(S3T,1));
outputname=[Case,'_Eulerian_SF3.mat']
save(outputname,'S3L1_alltime','S3T1_alltime','r');

disp('***************done**************************')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% check quick %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;close all;clc
str1=1;
end1=500e3;
cgname='s2sflux_spec_hf_500m_kaiser1.0002.nc';
filtscale=ncread(cgname,'filtscale');
Thm_hf=ncread(cgname,'Thm');
filtscale=filtscale(2:end);
Thm_hf=Thm_hf(2:end);

cgname='s2sflux_spec_smooth_500m_kaiser1.0002.nc';
filtscale=ncread(cgname,'filtscale');
Thm_sm=ncread(cgname,'Thm');
filtscale=filtscale(2:end);
Thm_sm=Thm_sm(2:end);

load wave_500_Eulerian_SF3.mat
% double(S3T1_alltime)
S3L1=vertcat(S3L1_alltime{:});
S3T1=vertcat(S3T1_alltime{:});
r=r(2:end);
SF3=(S3T1(:,2:end)+S3L1(:,2:end))';
SF3_mean=nanmean(SF3,2);SF3_mean_hf=SF3_mean;
% lambda=[10,1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6, ...
%     1e-7,1e-8,1e-9,1e-10,1e-11,1e-12];
% lambda=[1e-7,1e-8,1e-9,1e-10,1e-11];
lambda=[1e-7,1e-8,1e-9,1e-10,1e-11];

clear SpecFlux_E;clear kf_E;
kf1=fliplr(1./filtscale'.*2.*pi);
[SpecFlux_E_hf,Vt_E_hf,ebs_E_hf,kf_E,lf_E]=Fk_fitting_SF3_Lcurve(SF3_mean, ...
    r,str1,end1,'fuc','RLS',lambda,kf1,1);

load nowave_500_Eulerian_SF3.mat
% double(S3T1_alltime)
S3L1=vertcat(S3L1_alltime{:});
S3T1=vertcat(S3T1_alltime{:});
r=r(2:end);
SF3=(S3T1(:,2:end)+S3L1(:,2:end))';
SF3_mean=nanmean(SF3,2);SF3_mean_sm=SF3_mean;
% lambda=[1e-1,1e-2,1e-3,1e-4,1e-5,1e-6, ...
%     1e-7,1e-8,1e-9,1e-10,1e-11,1e-12,1e-13];
clear SpecFlux_E;clear kf_E;
kf1=fliplr(1./filtscale'.*2.*pi);
[SpecFlux_E_sm,Vt_E_sm,ebs_E_sm,kf_E,lf_E]=Fk_fitting_SF3_Lcurve(SF3_mean, ...
    r,str1,end1,'fuc','RLS',lambda,kf1,1);

figure(2)
subplot(2,2,1)
semilogx(r./1e3,SF3_mean_hf./r');hold on
semilogx(lf_E./1e3,Vt_E_hf./lf_E','Marker','+','LineStyle','none')
grid on


subplot(2,2,2)

semilogx(filtscale./1e3,Thm_hf,'color','b','Linewidth',1.5)
hold on
semilogx(filtscale./1e3,Thm_sm,'color','r','Linewidth',1.5)

semilogx(1./kf_E.*2.*pi./1e3,SpecFlux_E_hf,'color','b', ...
    'Linewidth',1.5,'LineStyle','--');
semilogx(1./kf_E.*2.*pi./1e3,SpecFlux_E_sm,'color','r', ...
    'Linewidth',1.5,'Linewidth',1.5,'LineStyle','--');

grid on

subplot(2,2,3)
semilogx(filtscale./1e3,Thm_hf);hold on
semilogx(filtscale./1e3,Thm_sm);
grid on

subplot(2,2,4)
semilogx(r./1e3,SF3_mean_sm./r');hold on
semilogx(lf_E./1e3,Vt_E_sm./lf_E','Marker','+','LineStyle','none')
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;close all;clc
str1=1;
end1=500e3;
cgname='s2sflux_spec_hf_500m_kaiser1.0002.nc';
filtscale=ncread(cgname,'filtscale');
Thm_hf=ncread(cgname,'Thm');
filtscale=filtscale(2:end);
Thm_hf=Thm_hf(2:end);

kc=flipud(1./filtscale);
Thmc_hf=flipud(Thm_hf);
kc_mid_s=0.5.*(kc(2:end)+kc(1:end-1));
dPidk_hf_s=(Thmc_hf(2:end)-Thmc_hf(1:end-1));

cgname='s2sflux_spec_smooth_500m_kaiser1.0002.nc';
filtscale=ncread(cgname,'filtscale');
Thm_sm=ncread(cgname,'Thm');
filtscale=filtscale(2:end);
Thm_sm=Thm_sm(2:end);
Thmc_sm=flipud(Thm_sm);
dPidk_sm_s=(Thmc_sm(2:end)-Thmc_sm(1:end-1));

load wave_500_Eulerian_SF3.mat
S3L1=vertcat(S3L1_alltime{:});
S3T1=vertcat(S3T1_alltime{:});
r=r(2:end);
SF3=(S3T1(:,2:end)+S3L1(:,2:end))';
SF3_mean=nanmean(SF3,2);SF3_mean_hf=SF3_mean;
lambda=[1e-7,1e-8,1e-9,1e-10,1e-11];

clear SpecFlux_E;clear kf_E;
kf1=fliplr(1./filtscale'.*2.*pi);
[SpecFlux_E_hf,Vt_E_hf,ebs_E_hf,kf_E,lf_E]=Fk_fitting_SF3_Lcurve(SF3_mean, ...
    r,str1,end1,'fuc','RLS',lambda,kf1,0);
dk_E=v2rho_2d(abs(diff(kf_E)));


load nowave_500_Eulerian_SF3.mat
S3L1=vertcat(S3L1_alltime{:});
S3T1=vertcat(S3T1_alltime{:});
r=r(2:end);
SF3=(S3T1(:,2:end)+S3L1(:,2:end))';
SF3_mean=nanmean(SF3,2);SF3_mean_sm=SF3_mean;

kf1=fliplr(1./filtscale'.*2.*pi);
[SpecFlux_E_sm,Vt_E_sm,ebs_E_sm,kf_E,lf_E]=Fk_fitting_SF3_Lcurve(SF3_mean, ...
    r,str1,end1,'fuc','RLS',lambda,kf1,0);

%
colors={'#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30','#4DBEEE'};


c4=semilogx(1./kf_E(1:end-2)./1e3.*2.*pi,ebs_E_sm(1:end-3).*dk_E(1:end-2),...
    'LineWidth', 1.5,'Color',colors{5});
hold on
c1=semilogx(1./kc_mid_s./1e3,dPidk_sm_s,'LineWidth',1.5,'Color',colors{1});

