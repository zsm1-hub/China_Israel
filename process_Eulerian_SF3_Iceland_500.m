clear all;close all
clc
%%%%%%%%%%%%%%%%calc Eulerian SF3 for Iceland%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('/meddy/simingzhang/Analysis/matlab/Parcels_SF/')
addpath('/meddy/simingzhang/Data/Parcels_data')
addpath(genpath('/meddy/simingzhang/Data/iceland_500m'))
input_dir='/meddy/simingzhang/Data/iceland_500m/';
Case='wave_500';
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
    [~,~,S3L1_alltime{t}(:),S3T1_alltime{t}(:)]=calc_radial(S3L11,S3T11,N,xscale1);
    disp(t)
end
delete(gcp('nocreate'));

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
save(outputname,'S3L1_alltime','S3T1_alltime');

disp('***************done**************************')

