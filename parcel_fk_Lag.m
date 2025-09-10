clear all;close all;clc
% addpath('D:\LIN2023\model\RoyBarkan\LLC4320/')
% addpath('D:\LIN2023\crocotools\Preprocessingtools') % add function "spheric_dist.m"
% 
addpath('/meddy/simingzhang/Analysis/matlab/Parcels_SF/')
addpath('/meddy/simingzhang/Data/Parcels_data')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          1. Basic setup and read data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Case='wave'; % wave
nparticles=[289,625,2500,15376]; % numbers of particles
% nparticles=[289]; % numbers of particles

days=89.5;  % days
dt=3600; % s  Advection_RK4 delta_t drift时间间隔
% input_dir='D:\LIN2023\model\RoyBarkan\LLC4320/'; % drift所在文件夹
input_dir='/meddy/simingzhang/Data/Parcels_data/';
% timerange=24*10:24*11-6; % 计算结构函数用的时间范围
% timerange=1:2140;
inv_style='RLS';
lambda=1e-10;

% if strcmpi(Case, 'wave')
%     fname=['wave_pars_P',num2str(nparticles),'T',num2str(days),'days.nc'];
% end
% 
% if strcmpi(Case, 'nowave')
%     fname=['nowave_pars_P',num2str(nparticles),'T',num2str(days),'days.nc'];
% end
colors={'#0072BD','#D95319','#EDB120','#7E2F8E'};
for ii=1:length(nparticles)
    fname{ii,:}=[Case,'_pars_P',num2str(nparticles(ii)),'T',num2str(days),'days.nc'];
    eval(['load ',fname{ii}(1:end-3),'SF123.mat']);
    if ii>2
        SF3=(SF3lll_time+SF3ltt_time)';
    else
        SF3=(SF3lll+SF3ltt)';
    end
    [SpecFlux,Vt,ebs,kf,lf]=Fk_fitting_SF3(SF3,dist_axis,2,300e3,'log','RLS',lambda);
    
    a{ii}=semilogx(1./kf./1e3,SpecFlux,'Marker','x','Color',colors{ii},'LineWidth',1.5);hold on
    grid on
    
end
fname='s2sflux_spec_smooth.0002.nc';
ncdisp(fname)
Thm_Eulerian=ncread(fname,'Thm');
filtscale=ncread(fname,'filtscale');
a{5}=semilogx(filtscale./1e3,Thm_Eulerian,'LineWidth',1.5,'Color',[.7,.7,.7]);

legend([a{1},a{2},a{3},a{4},a{5}],{'287','625','2500','15376','Eulerian'}, ...
        'Location', 'northeast')
ylim([-2e-7,2e-7]);
xlim([1e3,1e6]./1e3);
xlabel('km')
ylabel('m^{2}/s^{3}')
title([Case,': cg flux'])
set(gca,'fontsize',16,'FontWeight','b')
%% Lagrangian cg flux vs sf3-fitting fk



clear all;close all;clc
% addpath('D:\LIN2023\model\RoyBarkan\LLC4320/')
% addpath('D:\LIN2023\crocotools\Preprocessingtools') % add function "spheric_dist.m"
% 
addpath('/meddy/simingzhang/Analysis/matlab/Parcels_SF/')
addpath('/meddy/simingzhang/Data/Parcels_data')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          1. Basic setup and read data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Case='wave'; % wave
nparticles=[289,625,2500,15376]; % numbers of particles
% nparticles=[289]; % numbers of particles

days=89.5;  % days
dt=3600; % s  Advection_RK4 delta_t drift时间间隔
% input_dir='D:\LIN2023\model\RoyBarkan\LLC4320/'; % drift所在文件夹
input_dir='/meddy/simingzhang/Data/Parcels_data/';
% timerange=24*10:24*11-6; % 计算结构函数用的时间范围
% timerange=1:2140;
inv_style='RLS';
lambda=1e-10;

% if strcmpi(Case, 'wave')
%     fname=['wave_pars_P',num2str(nparticles),'T',num2str(days),'days.nc'];
% end
% 
% if strcmpi(Case, 'nowave')
%     fname=['nowave_pars_P',num2str(nparticles),'T',num2str(days),'days.nc'];
% end
if strcmpi(Case, 'wave')
    cgname='s2sflux_spec_hf.0002.nc';
end

if strcmpi(Case, 'nowave')
    cgname='s2sflux_spec_smooth.0002.nc';
end

% cgname='s2sflux_spec_smooth.0002.nc';
ncdisp(cgname)
Thm_Eulerian=ncread(cgname,'Thm');
filtscale=ncread(cgname,'filtscale');


colors={'#0072BD','#D95319','#EDB120','#7E2F8E'};
A={'a1','a2','a3','a4'}
B={'b1','b2','b3','b4'}
C={'c1','c2','c3','c4'}
D={'d1','d2','d3','d4'}
E={'e1','e2','e3','e4'}


for ii=1:length(nparticles)
    fname{ii,:}=[Case,'_pars_P',num2str(nparticles(ii)),'T',num2str(days),'days.nc'];
    eval(['load ',fname{ii}(1:end-3),'SF123.mat']);
    if ii>2
        SF3=(SF3lll_time+SF3ltt_time)';
        SF1=(SF1l_time+SF1t_time)';
        SF1L=(SF1l_time)';
        SF1T=(SF1t_time)';
    else
        SF3=(SF3lll+SF3ltt)';
        SF1=(SF1l+SF1t)';
        SF1L=(SF1l)';
        SF1T=(SF1t)';

    end
    [SpecFlux,Vt,ebs,kf,lf]=Fk_fitting_SF3(SF3,dist_axis,2,300e3,'log','RLS',lambda);

    figure(1)
    subplot(2,2,ii)
    B{ii}=semilogx(1./kf./1e3,SpecFlux,'Marker','x','Color',colors{ii},'LineWidth',1.5);
    hold on
    eval(['load ',Case,'_pars_P',num2str(nparticles(ii)),'T89.5daysCG_Lag_tukey.mat']);
    C{ii}=semilogx(filtscale./1e3,Th','LineWidth',1.5,'Color',colors{ii},'LineStyle','--');
    % semilogx(filtscale./1e3,Thm_Eulerian','LineWidth',1.5,'Color',[.7, .7, .7],'LineStyle','-');
    grid on
    legend([B{ii},C{ii}],{'SF3-fitting fk','Lag cg fk'},'Location','southeast');

    ylim([-1e-7,1.5e-7]);
    xlim([1e3,1e6]./1e3);
    xlabel('km')
    ylabel('m^{2}/s^{3}')
    title([Case,'P',num2str(nparticles(ii)),': cg v.s. sf3-fitting'])
    set(gca,'fontsize',16,'FontWeight','b')

    figure(2)
    semilogx(1./kf./1e3,ebs(1:end-1),'Marker','x','Color',colors{ii},'LineWidth',1.5)
    hold on
    grid on
    xlim([1e3,1e6]./1e3);
    ylim([-2.5e-3,2.5e-3]);
    xlabel('km')
    ylabel('m^{2}/s^{3}')
    title([Case,': energy injection'])
    set(gca,'fontsize',16,'FontWeight','b')


    figure(3)
    semilogx(dist_axis./1e3,SF1,'Marker','x','Color',colors{ii},'LineWidth',1.5)
    % semilogx(dist_axis./1e3,SF1L,'Marker','x','Color',colors{ii},'LineWidth',1.5,'LineStyle','--')
    % semilogx(dist_axis./1e3,SF1T,'Marker','x','Color',colors{ii},'LineWidth',1.5,'LineStyle','-.')
    hold on
    grid on
    xlim([1e3,1e6]./1e3);
    ylim([0,1]);
    xlabel('km')
    ylabel('m/s')
    title([Case,': SF1'])
    set(gca,'fontsize',16,'FontWeight','b')


    figure(4)
    A{ii}=loglog(dist_axis./1e3,abs(SF3),'Color',colors{ii},'LineWidth',1.5)
    % semilogx(dist_axis./1e3,SF1L,'Marker','x','Color',colors{ii},'LineWidth',1.5,'LineStyle','--')
    % semilogx(dist_axis./1e3,SF1T,'Marker','x','Color',colors{ii},'LineWidth',1.5,'LineStyle','-.')
    hold on
    semilogx(dist_axis./1e3,-(SF3),'Color',colors{ii},'LineWidth',1.5,'Marker','+')
    semilogx(dist_axis./1e3,(SF3),'Color',colors{ii},'LineWidth',1.5,'Marker','o')
    grid on
    
    xlim([1e3,1e6]./1e3);
    ylim([1e-5,1e-2]);
    xlabel('km')
    ylabel('m^{3}/s^{3}')
    title([Case,': SF3'])
    set(gca,'fontsize',16,'FontWeight','b')

    figure(5)
    subplot(2,2,ii)
    D{ii}=semilogx(dist_axis./1e3,(SF3),'Color',colors{ii},'LineWidth',1.5)
    hold on
    E{ii}=semilogx(lf./1e3,(Vt),'Color',colors{ii},'LineWidth',1.5,'LineStyle','--')
    legend([D{ii},E{ii}],{'SF3','SF3 reconstruct'},'Location','southeast');
    grid on
    xlim([1e3,1e6]./1e3);
    ylim([-2e-3,10e-3]);
    xlabel('km')
    ylabel('m^{3}/s^{3}')
    title([Case,': SF3 v.s. SF3 reconstruct'])
    set(gca,'fontsize',16,'FontWeight','b')
end
legend([A{1},A{2},A{3},A{4}],{'289','625','2500','15376'},'Location','southeast');

