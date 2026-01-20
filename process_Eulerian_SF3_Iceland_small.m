clear all;close all
clc
%%%%%%%%%%%%%%%%calc Eulerian SF3 for Iceland%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('/meddy/simingzhang/Analysis/matlab/Parcels_SF/')
addpath('/meddy/simingzhang/Data/Parcels_data')
addpath('/meddy/simingzhang/Data/RB_iceland_data')
addpath('/meddy/simingzhang/Data/PYC');
input_dir='/meddy/simingzhang/Data/RB_iceland_data/';
% Case='wave';
Case='wave';

if strcmpi(Case, 'wave')
    fname=[input_dir,'z_niskin2km_his_hf_depth_500m_grd.0002.nc'];
    gname=[input_dir,'niskin2km_500m_grd.nc'];
    N=287;
end

if strcmpi(Case, 'nowave')
    fname=[input_dir,'z_niskin2km_his_smooth_depth_500m_grd.0002.nc'];
    gname=[input_dir,'niskin2km_500m_grd.nc'];
    N=287;
end


if strcmpi(Case, 'SWC2km_tide')
    input_dir='/meddy/simingzhang/Data/PYC/'
    fname=[input_dir,'z_sample_SWC2km_winter.0002.nc'];
    gname=[input_dir,'sample_grid.nc'];
    N=241;

end


dx=mean(mean(1./ncread(gname,'pm')));
dy=mean(mean(1./ncread(gname,'pn')));

xscale=mean([dx,dy]);
disp(['******the scale of grid mean value is ',num2str(xscale),'m ********']);

u=(squeeze(ncread(fname,'u')));
u1=u2rho_3d(permute(u,[3,2,1]));
v=(squeeze(ncread(fname,'v')));
v1=v2rho_3d(permute(v,[3,2,1]));

% u1=u1(:,85:195,85:195);
% v1=v1(:,85:195,85:195);
N=size(u1,2);

for t = 1:size(u1,1)
   u2=squeeze(u1(t,:,:));
    v2=squeeze(v1(t,:,:));
    % [X, Y] = meshgrid(-N/2:N/2-1, -N/2:N/2-1);
    [S3L(t,:,:),S3T(t,:,:)]=test4_calc_SF3(u2,v2,N);
    rescale=2;xscale1=xscale;
    
    [~,~,S3L1_alltime(t,:),S3T1_alltime(t,:)]=calc_radial(S3L(t,:,:),S3T(t,:,:),N,xscale1);
    disp(t)
end

% [~, S3L1_iso] = calc_ispec2(x1(1,:), y1(:,1), nanmean(S3L,1), 2);
% [rbin, S3T1_iso] = calc_ispec2(x1(1,:), y1(:,1), nanmean(S3T,1), 2);
% rbin1=diag(rbin);
% clear rbin
% dist_axis=rbin1(2:end).*xscale;
% SF3=(S3L1_iso(2:end)+S3T1_iso(2:end));
% outputname=[Case,'_Eulerian_SF3.mat']
% save(outputname,'SF3','dist_axis');

S3L1=squeeze(nanmean(S3L,1));
S3T1=squeeze(nanmean(S3T,1));
outputname=[Case,'small_Eulerian_SF3_small.mat']
save(outputname,'S3L1','S3T1','xscale','S3L1_alltime','S3T1_alltime');

disp('***************done**************************')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
% load wave_Eulerian_SF3_small.mat
load wavesmall_Eulerian_SF3_small.mat
% load SWC2km_tide_Eulerian_SF3.mat
% magn=8e-4;

xscale1=xscale;
N=111;
[r,~,~,~]=calc_radial(S3L1,S3T1,N,xscale1);
SF3=(S3T1_alltime+S3L1_alltime)';
SF3_mean=nanmean(SF3,2);
% lambda=[1e-5,1e-6, ...
%     1e-7,1e-8,1e-9,1e-10,1e-11,1e-12,1e-13,1e-14];
lambda=[1e-7,1e-8,1e-9,1e-10,1e-11];
% lambda=[1e-7,1e-8,1e-9,1e-10,1e-11];

fname='s2sflux_spec_hf_kaiser_corr.0002.nc '
% fname='s2sflux_spec_SWC2km_tide.0002_kaiser.nc'

filtscale=ncread(fname,'filtscale');
Thm=ncread(fname,'Thm');


clear SpecFlux_E;clear kf_E;
kf1=fliplr(1./filtscale'.*2.*pi);

ns = find(r >= 1.5e3, 1);
ne = find(r <= 200e3, 1, 'last');

R = r(ns:ne); 
% kf1 = logspace(log10(1/max(R))-log10(2.*pi), log10(1/min(R)), 10).*2.*pi;
kf1 = logspace(log10(1/max(R)), log10(1/min(R)), 20).*2.*pi;
% kf1 = logspace(log10(1/max(R)), log10(1/min(R))+log10(2.*pi), 20);

[SpecFlux_E_hf,Vt_E_hf,ebs_E_hf,kf_E,lf_E]=Fk_fitting_SF3_Lcurve(SF3_mean, ...
    r,2.5e3,200e3,'log','RLS',lambda,kf1,1);
rescale=2.*pi;
% rescale=1;
figure(1)
subplot(2,2,1)
semilogx(filtscale,Thm,'Marker','+')
hold on
semilogx(1./kf_E.*rescale,SpecFlux_E_hf)

subplot(2,2,2)
semilogx(r,SF3_mean./r');hold on
semilogx(lf_E,Vt_E_hf./lf_E');

subplot(2,2,3)
semilogx(1./kf_E.*rescale./1e3,ebs_E_hf(2:end).*kf_E');hold on

%%%%%%%%%%%%%%%%
clear
% load wave_Eulerian_SF3_small.mat
load wavesmall_Eulerian_SF3_small.mat
% load SWC2km_tide_Eulerian_SF3.mat
% magn=8e-4;

xscale1=xscale;
N=111;
[r1,~,~,~]=calc_radial(S3L1,S3T1,N,xscale1);
SF3=(S3T1_alltime+S3L1_alltime)';
SF3_mean1=nanmean(SF3,2);

load wave_Eulerian_SF3.mat
xscale1=xscale;
N=287;
[r,~,~,~]=calc_radial(S3L1,S3T1,N,xscale1);
SF3=(S3T1_alltime+S3L1_alltime)';
SF3_mean=nanmean(SF3,2);


semilogx(r,SF3_mean./r');hold on
semilogx(r1,SF3_mean1./r1');
ylim([-1e-7,1e-7])
load wave_pars_P289T1940_roughsmallbootstrap.mat
semilogx(dist_axis,SF3_mean./dist_axis');
