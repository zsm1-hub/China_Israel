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
Case='wave'; % wave
% nparticles=625; % numbers of particles
days=89.5;  % days
dt=3600; % s  Advection_RK4 delta_t drift时间间隔
% input_dir='D:\LIN2023\model\RoyBarkan\LLC4320/'; % drift所在文件夹
input_dir='/meddy/simingzhang/Data/Parcels_data/tranV_onetime_spectukey/';
grd='/meddy/simingzhang/Data/RB_iceland_data/niskin2km_500m_grd.nc'
% timerange=24*10:24*11-6; % 计算结构函数用的时间范围
timerange=1:2140;
np1=[289,625,2500,15376];
ini='_rough';%grid, rough, repeat
if strcmpi(ini, '_grid')
    input_dir='/meddy/simingzhang/Data/Parcels_data/tranV_onetime_spectukey/';
end

if strcmpi(ini, '_rough')
    input_dir='/meddy/simingzhang/Data/Parcels_data/tranV_onetime_roughdistr_tukey/';
end

for ii=1:4
    nparticles=np1(ii);
    if strcmpi(Case, 'wave')
        fname=[input_dir,'wave_pars_P',num2str(nparticles),'T',num2str(days),'days.nc'];
    end

    if strcmpi(Case, 'nowave')
        fname=[input_dir,'nowave_pars_P',num2str(nparticles),'T',num2str(days),'days.nc'];
    end
    disp(fname)

    lon=ncread(fname,'lon');
    lat=ncread(fname,'lat');

    lon_g=ncread(grd,'lon_rho');
    lat_g=ncread(grd,'lat_rho');
    
    save([fname(1:end-3),'traj.mat'])
    save([fname(1:end-3),'traj.mat'])
    % clear
    % clf
    % load nowave_pars_P625T89.5daystraj.mat

    % for ii=1:size(lon,2)
    %     plot(lon(:,ii),lat(:,ii),'LineWidth',1.2);
    %     hold on
    % end
    % plot(lon_g(1,:),lat_g(1,:),'Color','r','LineWidth',1,'LineStyle','--');
    % plot(lon_g(end,:),lat_g(end,:),'Color','r','LineWidth',1,'LineStyle','--');
    % plot(lon_g(:,1),lat_g(:,1),'Color','r','LineWidth',1,'LineStyle','--');
    % plot(lon_g(:,end),lat_g(:,end),'Color','r','LineWidth',1,'LineStyle','--');
    % xlabel('Longitude')
    % ylabel('Latitude')
    % title([Case,': P',num2str(nparticles),' trajectories'])
    % set(gca,'fontsize',14,'FontWeight','b')
    % saveas(gcf,[fname(1:end-3),'traj.png'],'png')
    % clf
end
disp('done')

%% initial
colors={'#0072BD','#D95319','#EDB120','#7E2F8E'};

jj=4;
load wave_pars_P15376T89.5daystraj_grid.mat
a1=scatter(lon(1,:),lat(1,:),'LineWidth',1.2,'MarkerFaceColor',colors{jj}, ...
    'MarkerEdgeColor','none');hold on

jj=jj-1;
load wave_pars_P2500T89.5daystraj_grid.mat
a2=scatter(lon(1,:),lat(1,:),'LineWidth',1.2,'MarkerFaceColor',colors{jj}, ...
    'MarkerEdgeColor','none');

jj=jj-1;
load wave_pars_P625T89.5daystraj_grid.mat
a3=scatter(lon(1,:),lat(1,:),'LineWidth',1.2,'MarkerFaceColor',colors{jj}, ...
    'MarkerEdgeColor','none');

jj=jj-1;
load wave_pars_P289T89.5daystraj_grid.mat
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
%% traj

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
% nparticles=625; % numbers of particles
days=89.5;  % days
dt=3600; % s  Advection_RK4 delta_t drift时间间隔
% input_dir='D:\LIN2023\model\RoyBarkan\LLC4320/'; % drift所在文件夹
ini='_grid'
if strcmpi(ini, '_grid')
    input_dir='/meddy/simingzhang/Data/Parcels_data/tranV_onetime_spectukey/';
end
if strcmpi(ini, '_rough')
    input_dir='/meddy/simingzhang/Data/Parcels_data/tranV_onetime_roughdistr_tukey/';
end
% input_dir='/meddy/simingzhang/Data/Parcels_data/tranV_onetime_spectukey/';
grd='/meddy/simingzhang/Data/RB_iceland_data/niskin2km_500m_grd.nc'
% timerange=24*10:24*11-6; % 计算结构函数用的时间范围
timerange=1:2140;
np1=[289,625,2500,15376];

for ii=1:4
    nparticles=np1(ii);
    if strcmpi(Case, 'wave')
        fname=['wave_pars_P',num2str(nparticles),'T',num2str(days),'days.nc'];
    end

    if strcmpi(Case, 'nowave')
        fname=['nowave_pars_P',num2str(nparticles),'T',num2str(days),'days.nc'];
    end
    disp(fname)

    % lon=ncread(fname,'lon');
    % lat=ncread(fname,'lat');
    % 
    % lon_g=ncread(grd,'lon_rho');
    % lat_g=ncread(grd,'lat_rho');
    % save([fname(1:end-3),'traj.mat'])
    % clear
    % clf
    % load nowave_pars_P625T89.5daystraj.mat
    eval(['load ',fname(1:end-3),'traj.mat'])
    clear fname
    if strcmpi(Case, 'wave')
        fname=['wave_pars_P',num2str(nparticles),'T',num2str(days),'days.nc'];
    end

    if strcmpi(Case, 'nowave')
        fname=['nowave_pars_P',num2str(nparticles),'T',num2str(days),'days.nc'];
    end
    disp(fname)

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
    saveas(gcf,[fname(1:end-3),'traj',ini,'.png'],'png')
    clf
end
disp('done')
%% lagragian spec

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
% nparticles=625; % numbers of particles
days=89.5;  % days
dt=3600; % s  Advection_RK4 delta_t drift时间间隔
% input_dir='D:\LIN2023\model\RoyBarkan\LLC4320/'; % drift所在文件夹
ini='_grid'
if strcmpi(ini, '_grid')
    input_dir='/meddy/simingzhang/Data/Parcels_data/tranV_onetime_spectukey/';
end
if strcmpi(ini, '_rough')
    input_dir='/meddy/simingzhang/Data/Parcels_data/tranV_onetime_roughdistr_tukey/';
end
% input_dir='/meddy/simingzhang/Data/Parcels_data/tranV_onetime_spectukey/';
grd='/meddy/simingzhang/Data/RB_iceland_data/niskin2km_500m_grd.nc'
% timerange=24*10:24*11-6; % 计算结构函数用的时间范围
timerange=1:2140;
np1=[289,625,2500,15376];
colors={'#0072BD','#D95319','#EDB120','#7E2F8E'};
colors_rgb = {...
    [0, 114/255, 189/255], ...   % #0072BD
    [217/255, 83/255, 25/255], ... % #D95319
    [237/255, 177/255, 32/255], ... % #EDB120
    [126/255, 47/255, 142/255] ... % #7E2F8E
};

for ii=1:4
    nparticles=np1(ii);
    if strcmpi(Case, 'wave')
        fname=['wave_pars_P',num2str(nparticles),'T',num2str(days),'days.nc'];
    end

    if strcmpi(Case, 'nowave')
        fname=['nowave_pars_P',num2str(nparticles),'T',num2str(days),'days.nc'];
    end
    disp(fname)

    eval(['load ',fname(1:end-3),'traj',ini,'.mat'])
    clear fname
    if strcmpi(Case, 'wave')
        fname=['wave_pars_P',num2str(nparticles),'T',num2str(days),'days.nc'];
    end

    if strcmpi(Case, 'nowave')
        fname=['nowave_pars_P',num2str(nparticles),'T',num2str(days),'days.nc'];
    end
    disp(fname)
    lon=lon(1:2148,:);
    lat=lat(1:2148,:);
    for t=1:size(lon,1)-1
        U(t,:)=(spheric_dist(lat(t,:),lat(t,:),lon(t,:),lon(t+1,:)))./dt.*...
            sign(lon(t+1,:)-lon(t,:));
        V(t,:)=(spheric_dist(lat(t+1,:),lat(t,:),lon(t,:),lon(t,:)))./dt.*...
            sign(lat(t+1,:)-lat(t,:));
    end
    u1=zeros(2148,nparticles);v1=zeros(2148,nparticles);
    u1(1:end-1,:)=U;v1(1:end-1,:)=V;
    u1(end,:)=u1(end-1,:);v1(end,:)=v1(end-1,:);

    % find whole time exist particles
    J=find(sum(~isnan(lat),1)==2148);
    u2=u1(:,J);
    v2=v1(:,J);
    omega1=[-2148/2:2148/2-1]./2148*1/1;
    om=omega1((end+1)/2+1:end);
    
    f=4*pi*sin(pi*60/180)*366.25/(24*3600*365.25)*3600/2/pi;
    for iii=1:size(J,2)
        ke=abs(fftshift(fft(u2(:,iii).*hann(2148))))+...
            abs(fftshift(fft(v2(:,iii).*hann(2148))));
        ke1(:,iii)=ke((end+1)/2+1:end);
    end
    disp('spec done');
    % std1 = std(ke1, 0, 2, 'omitnan');
    % x_fill = [om;fliplr(om)];
    % y_fill = [nanmean(ke1,2)+std1,fliplr(nanmean(ke1,2)-std1)]';
    for i = 1:size(ke1,1)    
        CI_ke1(:,i) = prctile(ke1(i,:), [95,5]);
    end
    x_fill = [om, fliplr(om)];
    y_fill = [CI_ke1(1,:), fliplr(CI_ke1(2,:))];
    subplot(2,2,ii)
    loglog(om,nanmean(ke1,2),'LineWidth',1.5,'Color',colors{ii});hold on
    % loglog(om,nanmean(ke1,2)+std1,'LineWidth',1.5);
    % loglog(om,nanmean(ke1,2)-std1,'LineWidth',1.5);
    fill(x_fill, y_fill,colors_rgb{ii}, ...
        'FaceAlpha', 0.3, 'EdgeColor', 'none');    
    clear ke1;clear U;clear V;clear u1;clear v1;clear u2;clear v2
    a11=plot([1/12,1/12],[1e-2,1e2],'LineWidth',1.5)
    a22=plot([1/6,1/6],[1e-2,1e2],'LineWidth',1.5)
    a33=plot([1/24,1/24],[1e-2,1e2],'LineWidth',1.5)
    a44=plot([f,f],[1e-2,1e2],'color',[.7 .7 .7],'linestyle','--','LineWidth',1.5)
    % text(1/12,1e-1,'12 h')
    % text(1/6,1e-1,'6 h')
    % text(1/24,1e-1,'24 h')
    % text(f,1e-1,'1/f')
    xlabel('1/hour')
    ylabel('m^{2}/s^{2}')
    legend([a11,a22,a33,a44],{'12 h','6 h','24 h','1/f'},'Location','southwest')
    grid on
    title(['P',num2str(nparticles)]);
    set(gca,'fontsize',12,'FontWeight','bold')

end
sgtitle([Case,': Lagrangian KE spec']);
set(gca,'fontsize',12,'FontWeight','bold')
saveas(gcf,[Case,'Lagspec.png'],'png')
