%   data source: Iceland (wave and no wave case)
%   utility: import parcels data to calc SF2,SF3
%   doesn't use bootstrap to resample,insteadly, calc time-mean SF2 and SF3 directly
%   code writer: zsm, modified from Balwada 2022 sciadv supplyment
clear all;close all;clc
% addpath('D:\LIN2023\model\RoyBarkan\LLC4320/')
% addpath('D:\LIN2023\crocotools\Preprocessingtools') % add function "spheric_dist.m"
% 
addpath('/meddy/simingzhang/Analysis/matlab/Parcels_SF/')
addpath('/meddy/simingzhang/Data/Parcels_data')
addpath('/meddy/simingzhang/Data/RB_iceland_data')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          1. Basic setup and read data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Case='nowave'; % wave
nparticles=625; % numbers of particles
days=89.5;  % days
dt=3600; % s  Advection_RK4 delta_t drift时间间隔
% input_dir='D:\LIN2023\model\RoyBarkan\LLC4320/'; % drift所在文件夹
input_dir='/meddy/simingzhang/Data/Parcels_data/';
grd='/meddy/simingzhang/Data/RB_iceland_data/niskin2km_500m_grd.nc'
% timerange=24*10:24*11-6; % 计算结构函数用的时间范围
timerange=1:2140;

if strcmpi(Case, 'wave')
    fname=['wave_pars_P',num2str(nparticles),'T',num2str(days),'days.nc'];
end

if strcmpi(Case, 'nowave')
    fname=['nowave_pars_P',num2str(nparticles),'T',num2str(days),'days.nc'];
end

lon=ncread(fname,'lon');
lat=ncread(fname,'lat');

lon_g=ncread(grd,'lon_rho');
lat_g=ncread(grd,'lat_rho');
% save([fname(1:end-3),'traj.mat'])
clear
clf
load nowave_pars_P625T89.5daystraj.mat

for ii=1:size(lon,2)
    plot(lon(:,ii),lat(:,ii),'LineWidth',1.2);
    hold on
end
plot(lon_g(1,:),lat_g(1,:),'Color','r','LineWidth',1,'LineStyle','--');
plot(lon_g(end,:),lat_g(end,:),'Color','r','LineWidth',1,'LineStyle','--');
plot(lon_g(:,1),lat_g(:,1),'Color','r','LineWidth',1,'LineStyle','--');
plot(lon_g(:,end),lat_g(:,end),'Color','r','LineWidth',1,'LineStyle','--');
xlabel('Longitude')
ylabel('Latitude')
title([Case,': P',num2str(nparticles),' trajectories'])
set(gca,'fontsize',14,'FontWeight','b')
saveas(gcf,[fname(1:end-3),'traj.png'],'png')
disp('done')

%%
colors={'#0072BD','#D95319','#EDB120','#7E2F8E'};

jj=4;
load wave_pars_P15376T89.5daystraj.mat
a1=scatter(lon(1,:),lat(1,:),'LineWidth',1.2,'MarkerFaceColor',colors{jj}, ...
    'MarkerEdgeColor','none');hold on

jj=jj-1;
load wave_pars_P2500T89.5daystraj.mat
a2=scatter(lon(1,:),lat(1,:),'LineWidth',1.2,'MarkerFaceColor',colors{jj}, ...
    'MarkerEdgeColor','none');

jj=jj-1;
load wave_pars_P625T89.5daystraj.mat
a3=scatter(lon(1,:),lat(1,:),'LineWidth',1.2,'MarkerFaceColor',colors{jj}, ...
    'MarkerEdgeColor','none');

jj=jj-1;
load wave_pars_P289T89.5daystraj.mat
a4=scatter(lon(1,:),lat(1,:),'LineWidth',1.2,'MarkerFaceColor',colors{jj}, ...
    'MarkerEdgeColor','none');

plot(lon_g(1,:),lat_g(1,:),'Color','r','LineWidth',1,'LineStyle','--');
plot(lon_g(end,:),lat_g(end,:),'Color','r','LineWidth',1,'LineStyle','--');
plot(lon_g(:,1),lat_g(:,1),'Color','r','LineWidth',1,'LineStyle','--');
plot(lon_g(:,end),lat_g(:,end),'Color','r','LineWidth',1,'LineStyle','--');
xlabel('Longitude')
ylabel('Latitude')
title([Case,': P',num2str(nparticles),' initial position'])
legend([a1,a2,a3,a4],{'15376','2500','625','289'})
set(gca,'fontsize',14,'FontWeight','b')