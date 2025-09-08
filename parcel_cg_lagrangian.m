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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          1. Basic setup and read data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Case='HIT2d'; % wave
nparticles=2500; % numbers of particles
days=89.5;  % days
seconds=0.2;
dt=3600; % s  Advection_RK4 delta_t drift时间间隔
% input_dir='D:\LIN2023\model\RoyBarkan\LLC4320/'; % drift所在文件夹
% input_dir='/meddy/simingzhang/Data/Parcels_data/';
% input_dir='/meddy/simingzhang/Data/Parcels_data/onetime_spectukey/';
input_dir='/meddy/simingzhang/Data/Parcels_data/onetime_tukey/';

% timerange=24*10:24*11-6; % 计算结构函数用的时间范围
timerange=1:2140;

if strcmpi(Case, 'wave')
    fname=[input_dir,'wave_pars_P',num2str(nparticles),'T',num2str(days),'days.nc'];
end

if strcmpi(Case, 'nowave')
    fname=[input_dir,'nowave_pars_P',num2str(nparticles),'T',num2str(days),'days.nc'];
end

if strcmpi(Case, 'HIT2d')
    fname=[input_dir,Case,'_pars_P',num2str(nparticles),'T',num2str(seconds),'seconds.nc'];
end

lon=ncread(fname,'lon');
lat=ncread(fname,'lat');

% ue=ncread(fname,'ue').*1852.*60.*cos(lat.*pi./180);
% ve=ncread(fname,'ve').*1852.*60;

ue=ncread(fname,'ue');
ve=ncread(fname,'ve');

% read coarse-graining
xscale=[1:18,21:3:48,54:6:114];
PI=zeros(1,length(xscale));
for iii=1:length(xscale)
    pistr=['th',num2str(xscale(iii))];
    eval(['th',num2str(xscale(iii)),'=ncread(fname,','''',pistr,'''',');'])
    eval(['Th(',num2str(iii),')=nanmean(','th',num2str(xscale(iii)),'(:));'])
end
% semilogx(xscale.*2e3,Th)
save([fname(1:end-3),'CG_Lag_tukey.mat'],'Th');

%% plot
% hf
clear
fname='s2sflux_spec_hf.0002.nc';
ncdisp(fname)
Thm_Eulerian=ncread(fname,'Thm');
filtscale=ncread(fname,'filtscale');
% wave_pars_P625T89.5daysCG_Lag_tukey.mat

jj=1;
colors={[.7,.7,.7],'#0072BD','#D95319','#EDB120','#7E2F8E'};
a1=semilogx(filtscale./1e3,Thm_Eulerian,'LineWidth',1.5,'Color',colors{jj});
hold on

jj=jj+1;
load wave_pars_P289T89.5daysCG_Lag_tukey.mat
a2=semilogx(filtscale./1e3,Th','LineWidth',1.5,'Color',colors{jj});

jj=jj+1;
load wave_pars_P625T89.5daysCG_Lag_tukey.mat
a3=semilogx(filtscale./1e3,Th','LineWidth',1.5,'Color',colors{jj});

jj=jj+1;
load wave_pars_P2500T89.5daysCG_Lag_tukey.mat
a4=semilogx(filtscale./1e3,Th','LineWidth',1.5,'Color',colors{jj});

jj=jj+1;
load wave_pars_P15376T89.5daysCG_Lag_tukey.mat
a5=semilogx(filtscale./1e3,Th','LineWidth',1.5,'Color',colors{jj});
grid on
legend(['a1','a2','a3','a4','a5'],{'Eulerian','287','625','2500','15376'}, ...
    'Location', 'northeast')
ylim([-3e-8,3e-8]);
xlim([1e3,1e6]./1e3);
xlabel('km')
ylabel('m^{2}/s^{3}')
title('HF case: Coarsegraining flux')
set(gca,'fontsize',16,'FontWeight','b')

%%
clear
fname='s2sflux_spec_smooth.0002.nc';
ncdisp(fname)
Thm_Eulerian=ncread(fname,'Thm');
filtscale=ncread(fname,'filtscale');
% wave_pars_P625T89.5daysCG_Lag_tukey.mat

jj=1;
colors={[.7,.7,.7],'#0072BD','#D95319','#EDB120','#7E2F8E'};
a1=semilogx(filtscale./1e3,Thm_Eulerian,'LineWidth',1.5,'Color',colors{jj});
hold on

jj=jj+1;
load nowave_pars_P289T89.5daysCG_Lag_tukey.mat
a2=semilogx(filtscale./1e3,Th','LineWidth',1.5,'Color',colors{jj});

jj=jj+1;
load nowave_pars_P625T89.5daysCG_Lag_tukey.mat
a3=semilogx(filtscale./1e3,Th','LineWidth',1.5,'Color',colors{jj});

jj=jj+1;
load nowave_pars_P2500T89.5daysCG_Lag_tukey.mat
a4=semilogx(filtscale./1e3,Th','LineWidth',1.5,'Color',colors{jj});

jj=jj+1;
load nowave_pars_P15376T89.5daysCG_Lag_tukey.mat
a5=semilogx(filtscale./1e3,Th','LineWidth',1.5,'Color',colors{jj});
grid on
legend(['a1','a2','a3','a4','a5'],{'Eulerian','287','625','2500','15376'}, ...
    'Location', 'northeast')
ylim([-3e-8,3e-8]);
xlim([1e3,1e6]./1e3);
xlabel('km')
ylabel('m^{2}/s^{3}')
title('Smooth case: cg flux')
set(gca,'fontsize',16,'FontWeight','b')

%% hit2d
clear
fname='s2sflux_spec_hit_tukey.0002.nc';
ncdisp(fname)
Thm_Eulerian=ncread(fname,'Thm');
filtscale=ncread(fname,'filtscale');
% wave_pars_P625T89.5daysCG_Lag_tukey.mat

jj=1;
colors={[.7,.7,.7],'#0072BD','#D95319','#EDB120','#7E2F8E'};
a1=semilogx(filtscale,Thm_Eulerian,'LineWidth',1.5,'Color',colors{jj});
hold on

jj=jj+1;
load HIT2d_pars_P2500T0.2secondsCG_Lag_tukey.mat
a2=semilogx(filtscale,Th','LineWidth',1.5,'Color',colors{jj});
% 
% jj=jj+1;
% load wave_pars_P15376T89.5daysCG_Lag_tukey.mat
% a5=semilogx(filtscale./1e3,Th','LineWidth',1.5,'Color',colors{jj});

grid on
% legend(['a1','a2','a3','a4','a5'],{'Eulerian','287','625','2500','15376'}, ...
%     'Location', 'northeast')
ylim([-3e-8,3e-8]);
xlim([1e3,1e6]./1e3);
xlabel('km')
ylabel('m^{2}/s^{3}')
title('HF case: Coarsegraining flux')
set(gca,'fontsize',16,'FontWeight','b')


