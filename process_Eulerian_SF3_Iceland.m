clear all;close all
clc
%%%%%%%%%%%%%%%%calc Eulerian SF3 for Iceland%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('/meddy/simingzhang/Analysis/matlab/Parcels_SF/')
addpath('/meddy/simingzhang/Data/Parcels_data')
addpath('/meddy/simingzhang/Data/RB_iceland_data')
input_dir='/meddy/simingzhang/Data/RB_iceland_data/';
Case='nowave';

if strcmpi(Case, 'wave')
    fname=[input_dir,'z_niskin2km_his_hf_depth_500m_grd.0002.nc'];
end

if strcmpi(Case, 'nowave')
    fname=[input_dir,'z_niskin2km_his_smooth_depth_500m_grd.0002.nc'];
end

gname=[input_dir,'niskin2km_500m_grd.nc'];
dx=mean(mean(1./ncread(gname,'pm')));
dy=mean(mean(1./ncread(gname,'pn')));

xscale=mean([dx,dy]);
disp(['******the scale of grid mean value is ',num2str(xscale),'m ********']);

u=(squeeze(ncread(fname,'u')));
u1=u2rho_3d(permute(u,[3,2,1]));
v=(squeeze(ncread(fname,'v')));
v1=v2rho_3d(permute(v,[3,2,1]));
N=size(u1,2);

for t = 1:size(u1,1)
    u2=squeeze(u1(t,:,:));
    v2=squeeze(v1(t,:,:));
    % [X, Y] = meshgrid(-N/2:N/2-1, -N/2:N/2-1);
    [x1,y1,S3L(t,:,:),S3T(t,:,:)]=test4_calc_SF3(u2,v2,N);
    disp(t)
end

[~, S3L1_iso] = calc_ispec2(x1(1,:), y1(:,1), nanmean(S3L,1), 2);
[rbin, S3T1_iso] = calc_ispec2(x1(1,:), y1(:,1), nanmean(S3T,1), 2);
rbin1=diag(rbin);
clear rbin
dist_axis=rbin1(2:end);
SF3=(S3L1_iso(2:end)+S3T1_iso(2:end));


outputname=[Case,'_Eulerian_SF3.mat']
save(outputname,'SF3','dist_axis');
disp('***************done**************************')

