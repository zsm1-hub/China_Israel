
%% work on meddy 
%%%%%%%%%%%%%%% get traj %%%%%%%%%%%%%%%%%%%%%
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
Case='hit'; % wave
% nparticles=625; % numbers of particles
second=0.05;  % seconds
dt=3600; % s  Advection_RK4 delta_t drift时间间隔
% input_dir='D:\LIN2023\model\RoyBarkan\LLC4320/'; % drift所在文件夹
ini='_rough'
if strcmpi(ini, '_grid')
    input_dir='/meddy/simingzhang/Data/Parcels_data/HIT2d_grid/';
end
if strcmpi(ini, '_rough')
    input_dir='/meddy/simingzhang/Data/Parcels_data/HIT2d_rough/';
end
% input_dir='/meddy/simingzhang/Data/Parcels_data/tranV_onetime_spectukey/';
grd='/meddy/simingzhang/Data/RB_iceland_data/HIT2d_new.0002.nc'
% timerange=24*10:24*11-6; % 计算结构函数用的时间范围
timerange=1:180;
np1=[324,676,2704];
% np1=[289,625,2500];

for ii=1:length(np1)
    nparticles=np1(ii);
    fname=[input_dir,'HIT2d_pars_P',num2str(nparticles),'T',num2str(second),'seconds.nc'];

    disp(fname)

    lon=ncread(fname,'lon');
    lat=ncread(fname,'lat');

    ue=ncread(fname,'ue');
    ve=ncread(fname,'ve');

    lon_g=ncread(grd,'lon_rho');
    lat_g=ncread(grd,'lat_rho');
    
    ue=ncread(fname,'ue');
    ve=ncread(fname,'ve');
    
    % read coarse-graining
    xscale=[2:18,21:3:48,54:6:120,132:12:240,264:24:504];
    PI=zeros(1,length(xscale));
    for iii=1:length(xscale)
        pistr=['th',num2str(xscale(iii))];
        eval(['th',num2str(xscale(iii)),'=ncread(fname,','''',pistr,'''',');'])
        eval(['Th(',num2str(iii),')=nanmean(','th',num2str(xscale(iii)),'(:));'])
        eval(['Th_all(',num2str(iii),',:)=nanmean(','th',num2str(xscale(iii)),'(1:183,:),2);'])
    end
    % semilogx(xscale.*2e3,Th)
    
    % Th_all=Th_all;
    % save([fname(1:end-3),'traj.mat'])
    save([fname(1:end-3),'traj',ini,'.mat'])

    
end


%%%%%%%%%%%%%%% get snap vor th4 %%%%%%%%%%%%%%%%%%%%%
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
Case='hit'; % wave
% nparticles=625; % numbers of particles
second=0.05;  % seconds
dt=3600; % s  Advection_RK4 delta_t drift时间间隔
% input_dir='D:\LIN2023\model\RoyBarkan\LLC4320/'; % drift所在文件夹
ini='_rough'
if strcmpi(ini, '_grid')
    input_dir='/meddy/simingzhang/Data/Parcels_data/HIT2d_grid/';
end
if strcmpi(ini, '_rough')
    input_dir='/meddy/simingzhang/Data/Parcels_data/HIT2d_rough/';
end
% input_dir='/meddy/simingzhang/Data/Parcels_data/tranV_onetime_spectukey/';
grd='/meddy/simingzhang/Data/Parcels_data/HIT2dcase_modified_cg_tukey.nc'
% timerange=24*10:24*11-6; % 计算结构函数用的时间范围
tin=100;

fname=[input_dir,'HIT2d_new.0002.nc'];

lon_g=ncread(grd,'lon_rho');
lat_g=ncread(grd,'lat_rho');
u=ncread(grd,'u');
v=ncread(grd,'v');
Th4=ncread(grd,'Th4');
th4=squeeze(Th4(:,:,:,tin));
pm=ncread(grd,'pm');
pn=ncread(grd,'pn');
u1=(u2rho_2d([squeeze(u(:,:,:,tin))]'))';
v1=(v2rho_2d([squeeze(v(:,:,:,tin))]'))';
vx=(v1(2:end,:)-v1(1:end-1,:)).*0.5.*(pm(2:end,:)+pm(1:end-1,:));
uy=(u1(:,2:end)-u1(:,1:end-1)).*0.5.*(pn(:,2:end)+pn(:,1:end-1));
vor=(u2rho_2d(vx')-v2rho_2d(uy'))';

% save([fname(1:end-3),'traj.mat'])
save([fname(1:end-3),'snap',ini,'.mat'],'vor','lat_g',"lon_g",'th4');

    
% end

%% work on laptop


%%%%%%%%%%%%%%% plot ini %%%%%%%%%%%%%%%%%%%%%
clear all;close all;clc
% addpath('D:\LIN2023\model\RoyBarkan\LLC4320/')
% addpath('D:\LIN2023\crocotools\Preprocessingtools') % add function "spheric_dist.m"
%
addpath('/meddy/simingzhang/Analysis/matlab/Parcels_SF/')
addpath('/meddy/simingzhang/Data/Parcels_data')
addpath('/meddy/simingzhang/Data/RB_iceland_data')

Case='hit'; % wave
% nparticles=625; % numbers of particles
second=0.05;  % seconds
dt=3600; % s  Advection_RK4 delta_t drift时间间隔
% input_dir='D:\LIN2023\model\RoyBarkan\LLC4320/'; % drift所在文件夹
ini='_rough'
if strcmpi(ini, '_grid')
    input_dir='/meddy/simingzhang/Data/Parcels_data/HIT2d_grid/';
end
if strcmpi(ini, '_rough')
    input_dir='/meddy/simingzhang/Data/Parcels_data/HIT2d_rough/';
end
% input_dir='/meddy/simingzhang/Data/Parcels_data/tranV_onetime_spectukey/';
grd='/meddy/simingzhang/Data/RB_iceland_data/HIT2d_new.0002.nc'
% timerange=24*10:24*11-6; % 计算结构函数用的时间范围
timerange=1:180;
np1=[2704,676,324];
% np1=[289,625,2500];
colors={'#0072BD','#D95319','#EDB120'};

A1={'a1','a2','a3'};
for ii=1:length(np1)
    np1=[2704,676,324];
    nparticles=np1(ii);
    fname=['HIT2d_pars_P',num2str(nparticles),'T',num2str(second),'seconds.nc'];

    disp(fname)
    
    eval(['load ',fname(1:end-3),'traj',ini,'.mat'])

    
    figure(1)
    
    hold on
    plot(lon_g(1,:),lat_g(1,:),'Color','r','LineWidth',1,'LineStyle','-');
    plot(lon_g(end,:),lat_g(end,:),'Color','r','LineWidth',1,'LineStyle','-');
    plot(lon_g(:,1),lat_g(:,1),'Color','r','LineWidth',1,'LineStyle','-');
    plot(lon_g(:,end),lat_g(:,end),'Color','r','LineWidth',1,'LineStyle','-');
    pio2=pi/2
    pi3o2=pi*3/2
    plot([pio2,pio2],[pio2,pi3o2],color='k')
    plot([pi3o2,pi3o2],[pio2,pi3o2],color='k')
    plot([pio2,pi3o2],[pio2,pio2],color='k')
    plot([pio2,pi3o2],[pi3o2,pi3o2],color='k')
    xlabel('x [m]')
    ylabel('y [m]')
    xtick([0,pi/2,pi,3*pi/2,2*pi])
    xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
    ytick([0,pi/2,pi,3*pi/2,2*pi])
    yticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
    xlim([0 2*pi])
    ylim([0 2*pi])
    A1{ii}=scatter(lon(1,:),lat(1,:),5,'LineWidth',1.2,'MarkerFaceColor',colors{ii}, ...
    'MarkerEdgeColor','none');
    title(['HIT2d: Particles initial position'],'Interpreter','latex')

    set(gca,'fontsize',14,'FontWeight','b')

end
legend([A1{1},A1{2},A1{3}],{'324','676','2704'})


%%%%%%%%%%%%%%% plot traj %%%%%%%%%%%%%%%%%%%%%
clear all;close all;clc
% addpath('D:\LIN2023\model\RoyBarkan\LLC4320/')
% addpath('D:\LIN2023\crocotools\Preprocessingtools') % add function "spheric_dist.m"
%
addpath('/meddy/simingzhang/Analysis/matlab/Parcels_SF/')
addpath('/meddy/simingzhang/Data/Parcels_data')
addpath('/meddy/simingzhang/Data/RB_iceland_data')

Case='hit'; % wave
% nparticles=625; % numbers of particles
second=0.05;  % seconds
dt=3600; % s  Advection_RK4 delta_t drift时间间隔
% input_dir='D:\LIN2023\model\RoyBarkan\LLC4320/'; % drift所在文件夹
ini='_rough'
if strcmpi(ini, '_grid')
    input_dir='/meddy/simingzhang/Data/Parcels_data/HIT2d_grid/';
end
if strcmpi(ini, '_rough')
    input_dir='/meddy/simingzhang/Data/Parcels_data/HIT2d_rough/';
end
% input_dir='/meddy/simingzhang/Data/Parcels_data/tranV_onetime_spectukey/';
grd='/meddy/simingzhang/Data/RB_iceland_data/HIT2d_new.0002.nc'
% timerange=24*10:24*11-6; % 计算结构函数用的时间范围
timerange=1:180;
np1=[2704,676,324];
% np1=[289,625,2500];
colors={'#0072BD','#D95319','#EDB120'};

A1={'a1','a2','a3'};
for ii=1:length(np1)
    np1=[2704,676,324];
    nparticles=np1(ii);
    fname=['HIT2d_pars_P',num2str(nparticles),'T',num2str(second),'seconds.nc'];

    disp(fname)
    
    eval(['load ',fname(1:end-3),'traj',ini,'.mat'])

    
    figure(1)
    
    hold on
    plot(lon_g(1,:),lat_g(1,:),'Color','r','LineWidth',1,'LineStyle','-');
    plot(lon_g(end,:),lat_g(end,:),'Color','r','LineWidth',1,'LineStyle','-');
    plot(lon_g(:,1),lat_g(:,1),'Color','r','LineWidth',1,'LineStyle','-');
    plot(lon_g(:,end),lat_g(:,end),'Color','r','LineWidth',1,'LineStyle','-');
    pio2=pi/2
    pi3o2=pi*3/2
    plot([pio2,pio2],[pio2,pi3o2],color='k')
    plot([pi3o2,pi3o2],[pio2,pi3o2],color='k')
    plot([pio2,pi3o2],[pio2,pio2],color='k')
    plot([pio2,pi3o2],[pi3o2,pi3o2],color='k')
    xlabel('x [m]')
    ylabel('y [m]')
    xtick([0,pi/2,pi,3*pi/2,2*pi])
    xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
    ytick([0,pi/2,pi,3*pi/2,2*pi])
    yticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
    xlim([0 2*pi])
    ylim([0 2*pi])
    for jj=1:size(lon,2)
        plot(lon(:,jj),lat(:,jj),'LineWidth',1.2);
        hold on
    end
    % A1{ii}=scatter(lon(1,:),lat(1,:),5,'LineWidth',1.2,'MarkerFaceColor',colors{ii}, ...
    % 'MarkerEdgeColor','none');
    title(['HIT2d: P', num2str(nparticles),' trajectories'],'Interpreter','latex')

    set(gca,'fontsize',14,'FontWeight','b')
    saveas(gcf,['HIT2d_P',num2str(nparticles),'traj'],'png')
    clf
end

%%%%%%%%%%%%%%% plot snap %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;close all;clc
% addpath('D:\LIN2023\model\RoyBarkan\LLC4320/')
% addpath('D:\LIN2023\crocotools\Preprocessingtools') % add function "spheric_dist.m"
%
addpath('/meddy/simingzhang/Analysis/matlab/Parcels_SF/')
addpath('/meddy/simingzhang/Data/Parcels_data')
addpath('/meddy/simingzhang/Data/RB_iceland_data')
addpath(genpath('D:\LIN2023\crocotools'))

load HIT2d_new.0002snap_rough.mat

figure(1)
colorbar1=textread('BlWhRe.txt');

pcolor(lon_g,lat_g,vor./1e3);shading interp;hold on
colormap(colorbar1)
colorbar
caxis([-1.5,1.5])

plot(lon_g(1,:),lat_g(1,:),'Color','r','LineWidth',1,'LineStyle','-');
plot(lon_g(end,:),lat_g(end,:),'Color','r','LineWidth',1,'LineStyle','-');
plot(lon_g(:,1),lat_g(:,1),'Color','r','LineWidth',1,'LineStyle','-');
plot(lon_g(:,end),lat_g(:,end),'Color','r','LineWidth',1,'LineStyle','-');
pio2=pi/2
pi3o2=pi*3/2
plot([pio2,pio2],[pio2,pi3o2],color='k')
plot([pi3o2,pi3o2],[pio2,pi3o2],color='k')
plot([pio2,pi3o2],[pio2,pio2],color='k')
plot([pio2,pi3o2],[pi3o2,pi3o2],color='k')
xlabel('x [m]')
ylabel('y [m]')
xtick([0,pi/2,pi,3*pi/2,2*pi])
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
ytick([0,pi/2,pi,3*pi/2,2*pi])
yticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
xlim([0 2*pi])
ylim([0 2*pi])
text(pi/2,3.2/2*pi,'Launch particles','fontsize',14,'FontWeight','b')
title(['HIT2d: tindex=100 vorticity/1000'],'Interpreter','latex')
set(gca,'fontsize',14,'FontWeight','b')


figure(2)
colorbar1=textread('BlWhRe.txt');

pcolor(lon_g,lat_g,th4);shading interp;hold on
colormap(colorbar1)
cb = colorbar
caxis([-1000,1000])
title(cb, 'm^{2}/s^{3}', 'FontSize', 12, 'FontWeight', 'bold');

plot(lon_g(1,:),lat_g(1,:),'Color','r','LineWidth',1,'LineStyle','-');
plot(lon_g(end,:),lat_g(end,:),'Color','r','LineWidth',1,'LineStyle','-');
plot(lon_g(:,1),lat_g(:,1),'Color','r','LineWidth',1,'LineStyle','-');
plot(lon_g(:,end),lat_g(:,end),'Color','r','LineWidth',1,'LineStyle','-');
pio2=pi/2
pi3o2=pi*3/2
plot([pio2,pio2],[pio2,pi3o2],color='k')
plot([pi3o2,pi3o2],[pio2,pi3o2],color='k')
plot([pio2,pi3o2],[pio2,pio2],color='k')
plot([pio2,pi3o2],[pi3o2,pi3o2],color='k')
xlabel('x [m]')
ylabel('y [m]')
xtick([0,pi/2,pi,3*pi/2,2*pi])
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
ytick([0,pi/2,pi,3*pi/2,2*pi])
yticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
xlim([0 2*pi])
ylim([0 2*pi])
% text(pi/2,3.2/2*pi,'Launch particles','fontsize',14,'FontWeight','b')
title(['HIT2d: tindex=100 CG flux $\Pi^{4}_{r}$'],'Interpreter','latex')
set(gca,'fontsize',14,'FontWeight','b')

%%%%%%%%
%%%%%%%%%%%%%%% plot Lag cg %%%%%%%%%%%%%%%%%%%%%
clear all;close all;clc
% addpath('D:\LIN2023\model\RoyBarkan\LLC4320/')
% addpath('D:\LIN2023\crocotools\Preprocessingtools') % add function "spheric_dist.m"
%
addpath('/meddy/simingzhang/Analysis/matlab/Parcels_SF/')
addpath('/meddy/simingzhang/Data/Parcels_data')
addpath('/meddy/simingzhang/Data/RB_iceland_data')

Case='hit'; % wave
% nparticles=625; % numbers of particles
second=0.05;  % seconds
dt=3600; % s  Advection_RK4 delta_t drift时间间隔
% input_dir='D:\LIN2023\model\RoyBarkan\LLC4320/'; % drift所在文件夹
ini='_rough'
if strcmpi(ini, '_grid')
    input_dir='/meddy/simingzhang/Data/Parcels_data/HIT2d_grid/';
end
if strcmpi(ini, '_rough')
    input_dir='/meddy/simingzhang/Data/Parcels_data/HIT2d_rough/';
end
% input_dir='/meddy/simingzhang/Data/Parcels_data/tranV_onetime_spectukey/';
grd='/meddy/simingzhang/Data/RB_iceland_data/HIT2d_new.0002.nc'
% timerange=24*10:24*11-6; % 计算结构函数用的时间范围
timerange=1:180;
% np1=[2704,676,324];
np1=[324,676,2704];

% np1=[289,625,2500];
colors={'#EDB120','#D95319','#0072BD'};

fname='s2sflux_spec_hit_tukey.0002.nc';
ncdisp(fname)
% Thm_Eulerian=ncread(fname,'Thm');
filtscale=ncread(fname,'filtscale');
filtscale=filtscale(1:end-1);

A1={'a1','a2','a3'};
for ii=1:length(np1)
    np1=[324,676,2704];
    nparticles=np1(ii);
    fname=['HIT2d_pars_P',num2str(nparticles),'T',num2str(second),'seconds.nc'];

    disp(fname)
    
    eval(['load ',fname(1:end-3),'traj',ini,'.mat'])

    
    figure(1)
    A{ii}=semilogx(filtscale,nanmean(Th_all,2),'LineWidth',1.5,'LineStyle','-')
    
    hold on
    grid on
    ylim([-10,5]);
    xlabel('r [m]')
    ylabel('\Pi [m^{2}/s^{3}]')
    % legend([a1,a2,a3,a4,a5],{'cg spec','cg uni','fft-specflux','RLS','RLS (angular k)'})
    % legend([a1,a2,a3,a4],{'FFT specflux','CG (sfilt)','CG (unifilt)','RLS (angular k)'})
    % legend([a1,a2,a3,a4], ...
    %     {'FFT specflux','CG (sfilt)','CG (unifilt)','RLS (angular k)'}, ...
    %     "location",'northeast')
    title(['HIT2d: CG vs RLS cross-scale energy flux'],'Interpreter','latex')
    set(gca,'fontsize',14,'FontWeight','b')
    % plot(lon_g(1,:),lat_g(1,:),'Color','r','LineWidth',1,'LineStyle','-');
    % plot(lon_g(end,:),lat_g(end,:),'Color','r','LineWidth',1,'LineStyle','-');
    % plot(lon_g(:,1),lat_g(:,1),'Color','r','LineWidth',1,'LineStyle','-');
    % plot(lon_g(:,end),lat_g(:,end),'Color','r','LineWidth',1,'LineStyle','-');
    % pio2=pi/2
    % pi3o2=pi*3/2
    % plot([pio2,pio2],[pio2,pi3o2],color='k')
    % plot([pi3o2,pi3o2],[pio2,pi3o2],color='k')
    % plot([pio2,pi3o2],[pio2,pio2],color='k')
    % plot([pio2,pi3o2],[pi3o2,pi3o2],color='k')
    % xlabel('x [m]')
    % ylabel('y [m]')
    % xtick([0,pi/2,pi,3*pi/2,2*pi])
    % xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
    % ytick([0,pi/2,pi,3*pi/2,2*pi])
    % yticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
    % xlim([0 2*pi])
    % ylim([0 2*pi])
    % for jj=1:size(lon,2)
    %     plot(lon(:,jj),lat(:,jj),'LineWidth',1.2);
    %     hold on
    % end
    % % A1{ii}=scatter(lon(1,:),lat(1,:),5,'LineWidth',1.2,'MarkerFaceColor',colors{ii}, ...
    % % 'MarkerEdgeColor','none');
    % title(['HIT2d: P', num2str(nparticles),' trajectories'],'Interpreter','latex')
    % 
    % set(gca,'fontsize',14,'FontWeight','b')
    % saveas(gcf,['HIT2d_P',num2str(nparticles),'traj'],'png')
    % clf
end
legend([A{1},A{2},A{3}],{'P324','P676','P2704'})
set(gca,'fontsize',14,'FontWeight','b')

%% check lagrangian SF3 and fitting
clear
% load test2500_1.mat
% load test5000_2.mat; % 65536 randoms 5000
load test.mat; % 
kf1=1./dist_axis'.*2.*pi

str1=0.0;
en1=1.7;
figure(3)
a1=semilogx(dist_axis,SF3_mean./dist_axis','LineWidth',1.5)
hold on
load HIT2d_Eulerian_SF3.mat
lambda=8e-1;
dot=362;
a2=semilogx(r,SF3./r,'LineWidth',1.5)

grid on
xlabel('r [m]')
ylabel('SF3/r [m^{2}/s^{3}]')
% legend([a1,a2,a3,a4,a5],{'cg spec','cg uni','fft-specflux','RLS','RLS (angular k)'})
% legend([a1,a2,a3,a4],{'FFT specflux','CG (sfilt)','CG (unifilt)','RLS (angular k)'})
legend([a1,a2], ...
    {'Lagrangian SF3/r (P5000)','Eulerian SF3/r (fft)'}, ...
    "location",'southwest')
title(['HIT2d: The convergence of Third-order structure function'],'Interpreter','latex')
set(gca,'fontsize',14,'FontWeight','b')

ylim([-5,20])
% [SpecFlux,Vt,ebs,kf,lf]=Fk_fitting_SF3(SF3_mean,dist_axis,0.009,6.32, ...
%     'log','RLS',1e-10);

[residual_norms, solution_norms,...
lambda_opt_idx]=Fk_fitting_SF3_Lcurve(SF3_mean(1:end),dist_axis(1:end),str1,en1, ...
    'fuc','RLS',[10,1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6, ...
    1e-7,1e-8,1e-9,1e-10,1e-11,1e-12,1e-13,1e-14,1e-15],kf1);


clf
kf1=1./dist_axis'.*2.*pi
[SpecFlux_L,Vt_L,ebs_L,kf_L,lf_L]=Fk_fitting_SF3(SF3_mean(1:end),dist_axis(1:end),str1,en1, ...
    'fuc','RLS',1e-3,kf1);

figure(2)
a1=semilogx(dist_axis,SF3_mean./dist_axis','LineWidth',1.5);hold on
a2=semilogx(lf_L,Vt_L./lf_L','LineStyle','none','Marker','+','LineWidth',1.5);
grid on
xlabel('r [m]')
ylabel('SF3/r [m^{2}/s^{3}]')
legend([a1,a2], ...
    {'Lagrangian SF3/r (P5000)','RLS-fitting'}, ...
    "location",'southwest')
title(['HIT2d: The convergence of Third-order structure function'],'Interpreter','latex')
set(gca,'fontsize',14,'FontWeight','b')

ylim([-5,20])


figure(1)
load HIT2d_fftfk.mat
a1=semilogx(1./K1D,specFlux_mean,'LineWidth',1.5, ...
'Color','k','LineStyle','-'); 
hold on


fname='s2sflux_spec_hit_tukey.0002.nc';
ncdisp(fname)
Thm_Eulerian=ncread(fname,'Thm');
filtscale=ncread(fname,'filtscale');
a2=semilogx(filtscale,Thm_Eulerian,'LineWidth',1.5,'Color', ...
    'b');

kf1=K1D.*2.*pi;
% kf1=1./filtscale.*2.*pi;
[SpecFlux,Vt,ebs,kf,lf]=Fk_fitting_SF3(SF3(1:dot)',r(1:dot),0.009,6.32, ...
    'fuc','RLS',lambda,kf1);
% [residual_norms, solution_norms,...
% lambda_opt_idx]=Fk_fitting_SF3_Lcurve(SF3(1:dot)',r(1:dot),0.009,6.32, ...
%     'fuc','RLS',[10,1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6, ...
%     1e-7,1e-8,1e-9,1e-10,1e-11,1e-12],kf1);
a3=semilogx(1./kf.*(2.*pi),SpecFlux,'LineWidth',1.5, ...
'Color','r');

a4=semilogx(1./kf_L.*2.*pi,SpecFlux_L,'LineWidth',1.5,'Color',[0.95, 0.8, 0.9]);

load HIT2d_pars_P5000T0.05secondstraj.mat
a5=semilogx(filtscale(2:end),nanmean(Th_all,2),'LineWidth',1.5,'Color', ...
    [0.95, 0.85, 0.7]);

grid on
ylim([-12,8])
xlabel('r [m]')
ylabel('\Pi [m^{2}/s^{3}]')
% legend([a1,a2,a3,a4,a5],{'cg spec','cg uni','fft-specflux','RLS','RLS (angular k)'})
% legend([a1,a2,a3,a4],{'FFT specflux','CG (sfilt)','CG (unifilt)','RLS (angular k)'})
legend([a1,a2,a3,a4,a5], ...
    {'FFT specflux','CG (sfilt)','RLS (angular k)','Lag P5000 RLS','Lag P5000 CG'}, ...
    "location",'northeast')
title(['HIT2d: CG vs RLS cross-scale energy flux'],'Interpreter','latex')
set(gca,'fontsize',14,'FontWeight','b')

figure(8)
% dkk=mean(abs(K1D(2:end)-K1D(1:end-1)));
semilogx(K1D.*2.*pi,-divFlux_mean,'LineWidth',1.5, ...
'Color','k')
hold on
grid on
% dk=-(kf(2:end)-kf(1:end-1)).*2.*pi
semilogx(kf_L,ebs_L(1:end-1).*v2rho_2d(abs(diff(kf_L))),'LineWidth',1.5, ...
'Color','r');

% size(kf_L)
