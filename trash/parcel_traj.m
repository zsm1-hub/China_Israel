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
% nparticles=625; % numbers of particles
days=89.5;  % days
dt=3600; % s  Advection_RK4 delta_t drift时间间隔
% input_dir='D:\LIN2023\model\RoyBarkan\LLC4320/'; % drift所在文件夹
input_dir='/meddy/simingzhang/Data/Parcels_data/tranV_onetime_spectukey/';
grd='/meddy/simingzhang/Data/RB_iceland_data/niskin2km_500m_grd.nc'
% timerange=24*10:24*11-6; % 计算结构函数用的时间范围
timerange=1:2140;
np1=[289,625,2500,15376];
% np1=[289,625,2500];

ini='_grid';%grid, rough, repeat
if strcmpi(ini, '_grid')
    input_dir='/meddy/simingzhang/Data/Parcels_data/tranV_onetime_spectukey/';
end

if strcmpi(ini, '_rough')
    input_dir='/meddy/simingzhang/Data/Parcels_data/tranV_onetime_roughdistr_tukey/';
end

for ii=1:length(np1)
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

    ue=ncread(fname,'ue');
    ve=ncread(fname,'ve');

    lon_g=ncread(grd,'lon_rho');
    lat_g=ncread(grd,'lat_rho');
    
    % save([fname(1:end-3),'traj.mat'])
    save([fname(1:end-3),'traj',ini,'.mat'])
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
ini='_grid'
jj=4;
load wave_pars_P15376T89.5daystraj_grid.mat
a1=scatter(lon(1,:),lat(1,:),'LineWidth',1.2,'MarkerFaceColor',colors{jj}, ...
    'MarkerEdgeColor','none');hold on

jj=jj-1;
load wave_pars_P2500T89.5daystraj_grid.mat
a2=scatter(lon(1,:),lat(1,:),'LineWidth',1.2,'MarkerFaceColor',colors{jj}, ...
    'MarkerEdgeColor','none');hold on

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
title([ini(2:end),'  initial position'])
legend([a1,a2,a3,a4],{'15376','2500','625','289'})
% legend([a2,a3,a4],{'2500','625','289'})
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
Case='wave'; % wave
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
% np1=[289,625,2500];

for ii=1:length(np1)
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
    eval(['load ',fname(1:end-3),'traj',ini,'.mat'])
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
%% Eulerian spec
clear all;close all;clc
% addpath('D:\LIN2023\model\RoyBarkan\LLC4320/')
% addpath('D:\LIN2023\crocotools\Preprocessingtools') % add function "spheric_dist.m"
%
addpath('/meddy/simingzhang/Analysis/matlab/Parcels_SF/')
addpath('/meddy/simingzhang/Data/Parcels_data')
addpath('/meddy/simingzhang/Data/RB_iceland_data')
%
tic
Case='nowave';
if strcmpi(Case, 'wave')
    fname='/meddy/simingzhang/Data/RB_iceland_data/z_niskin2km_his_hf_depth_500m_grd.0002.nc'
end
if strcmpi(Case, 'nowave')
    fname='/meddy/simingzhang/Data/RB_iceland_data/z_niskin2km_his_smooth_depth_500m_grd.0002.nc'
end

ncdisp(fname)
u=u2rho_3d(permute((squeeze(ncread(fname,'u'))),[3,2,1]));
v=v2rho_3d(permute((squeeze(ncread(fname,'v'))),[3,2,1]));

for i=1:287
    for j=1:287
        ke(:,j,i) = abs(fftshift(fft(u(:,j,i).*hann(2148)))+...
            fftshift(fft(v(:,j,i).*hann(2148))));

    end
end

[T, J, I] = size(ke);

% 将第二维和第三维展平为一个维度
ke1 = reshape(ke, T, J*I);
ke1_eul=ke1((end+1)/2+1:end,:);

for i = 1:size(ke1,1)    
        CI_ke1(:,i) = prctile(ke1_eul(i,:), [95,5]);
end
ke1_eul_mean=nanmean(ke1_eul,2);

save([Case,'_freqspec_eul.mat'],'ke1_eul_mean','CI_ke1');
toc

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
Case='wave'; % wave
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
% np1=[289,625,2500];

colors={'#0072BD','#D95319','#EDB120','#7E2F8E'};
colors_rgb = {...
    [0, 114/255, 189/255], ...   % #0072BD
    [217/255, 83/255, 25/255], ... % #D95319
    [237/255, 177/255, 32/255], ... % #EDB120
    [126/255, 47/255, 142/255] ... % #7E2F8E
};
	
% color_eul=[0.4660, 0.6740, 0.1880];
color_eul=[0.7, 0.7, 0.7];
tic
for ii=1:length(np1)
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
    % ue=ue(1:2148,:);
    % ve=ve(1:2148,:);
    for t=1:size(lon,1)-1
        U(t,:)=(spheric_dist(lat(t,:),lat(t,:),lon(t,:),lon(t+1,:)))./dt.*...
            sign(lon(t+1,:)-lon(t,:));
        V(t,:)=(spheric_dist(lat(t+1,:),lat(t,:),lon(t,:),lon(t,:)))./dt.*...
            sign(lat(t+1,:)-lat(t,:));
    end
    u1=zeros(2148,nparticles);v1=zeros(2148,nparticles);
    u1(1:end-1,:)=U;v1(1:end-1,:)=V;
    u1(end,:)=u1(end-1,:);v1(end,:)=v1(end-1,:);
    % u1=ue;v1=ve;
   

    % find whole time exist particles
    % J=find(sum(~isnan(lat),1)==2148);
    J=1:nparticles;
    u2=u1(:,J);
    v2=v1(:,J);
    J=1:nparticles;
    % u2=u1(:,J);
    % v2=v1(:,J);
    % u2(isnan(u2))=0;
    % v2(isnan(v2))=0;
    omega1=[-2148/2:2148/2-1]./2148*1/1;
    om=omega1((end+1)/2+1:end);
    
    f=4*pi*sin(pi*60/180)*366.25/(24*3600*365.25)*3600/2/pi;
    for iii=1:size(J,2)
        ke=abs(fftshift(fft(u2(:,iii).*hann(2148))))+...
            abs(fftshift(fft(v2(:,iii).*hann(2148))));
        ke1(:,iii)=ke((end+1)/2+1:end);
    end
    disp('spec done');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Lag %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:size(ke1,1)    
        CI_ke1(:,i) = prctile(ke1(i,:), [95,5]);
    end
    x_fill = [om, fliplr(om)];
    y_fill = [CI_ke1(1,:), fliplr(CI_ke1(2,:))];
    subplot(2,2,ii)
    b1=loglog(om,nanmean(ke1,2),'LineWidth',1.5,'Color',colors{ii});hold on
    % loglog(om,nanmean(ke1,2)+std1,'LineWidth',1.5);
    % loglog(om,nanmean(ke1,2)-std1,'LineWidth',1.5);
    fill(x_fill, y_fill,colors_rgb{ii}, ...
        'FaceAlpha', 0.3, 'EdgeColor', 'none');    
    clear ke1;clear U;clear V;clear u1;clear v1;clear u2;clear v2
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Eul %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    eval(['load ',Case,'_freqspec_eul.mat']);
    x_fill = [om, fliplr(om)];
    y_fill = [CI_ke1(1,:), fliplr(CI_ke1(2,:))];
    b2=loglog(om,ke1_eul_mean,'LineWidth',1.5,'Color',color_eul);
    fill(x_fill, y_fill,color_eul, ...
        'FaceAlpha', 0.3, 'EdgeColor', 'none');    

    a11=plot([1/12,1/12],[1e-2,1e2],'LineWidth',1.5,'Color',[0.4660, 0.6740, 0.1880])
    a22=plot([1/6,1/6],[1e-2,1e2],'LineWidth',1.5,'Color',[0.3010, 0.7450, 0.9330])
    a33=plot([1/24,1/24],[1e-2,1e2],'LineWidth',1.5,'Color',[0.6350, 0.0780, 0.1840])
    a44=plot([f,f],[1e-2,1e2],'color',[.7 .7 .7],'linestyle','--', ...
        'LineWidth',1.5,'color','k')
    % text(1/12,1e-1,'12 h')
    % text(1/6,1e-1,'6 h')
    % text(1/24,1e-1,'24 h')
    % text(f,1e-1,'1/f')
    xlabel('1/hour')
    ylabel('m^{2}/s^{2}')
    legend([b1,b2,a11,a22,a33,a44],{'Lag','Eul','12 h','6 h','24 h','1/f'}, ...
        'Location','southwest')
    grid on
    title(['P',num2str(nparticles)]);
    set(gca,'fontsize',12,'FontWeight','bold')

end
sgtitle([Case,': Lagrangian KE spec']);
set(gca,'fontsize',12,'FontWeight','bold')
toc

saveas(gcf,[Case,'Lagspec.png'],'png')

% saveas(gcf,[Case,'testLagspec.png'],'png')


%% snapshot of Ro and particles
clear all;close all;clc
addpath(genpath('D:\colorbar'))
addpath(genpath('E:\DATA\RoyBarkan\RB_iceland_data\RB_iceland_data'))
addpath(genpath('D:\LIN2023\crocotools'))

colorbar1=textread('BlWhRe.txt');
tindex=250;

Case='nowave'; % wave
% nparticles=625; % numbers of particles
days=89.5;  % days
dt=3600; % s  Advection_RK4 delta_t drift时间间隔
nparticles=15376;
% input_dir='D:\LIN2023\model\RoyBarkan\LLC4320/'; % drift所在文件夹
ini='_grid'


grdname='niskin2km_500m_grd.nc'
fnamew='z_niskin2km_his_hf_depth_500m_grd.0002.nc'

if strcmpi(Case, 'wave')
    fnamew='z_niskin2km_his_hf_depth_500m_grd.0002.nc'
end
if strcmpi(Case, 'nowave')
    fnamew='z_niskin2km_his_smooth_depth_500m_grd.0002.nc'
end

ncdisp(fnamew)
ncdisp(grdname)

lon_r=ncread(grdname,'lon_rho');
lat_r=ncread(grdname,'lat_rho');
ocean_time=ncread(fnamew,'ocean_time');
dx=1./ncread(grdname,'pm');
dy=1./ncread(grdname,'pn');
f=ncread(grdname,'f');

u=ncread(fnamew,'u');
v=ncread(fnamew,'v');

u1=(u2rho_2d((squeeze(u(:,:,:,tindex)))'))';
v1=(v2rho_2d((squeeze(v(:,:,:,tindex)))'))';
clear u;clear v;

vx=v2rho_2d((v1(2:end,:)-v1(1:end-1,:))./(0.5.*(dx(2:end,:)+dx(1:end-1,:))));
uy=u2rho_2d((u1(:,2:end)-u1(:,1:end-1))./(0.5.*(dy(:,2:end)+dy(:,1:end-1))));
Ro=(vx-uy)./f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(ini, '_grid')
    input_dir='/meddy/simingzhang/Data/Parcels_data/tranV_onetime_spectukey/';
end
if strcmpi(ini, '_rough')
    input_dir='/meddy/simingzhang/Data/Parcels_data/tranV_onetime_roughdistr_tukey/';
end

if strcmpi(Case, 'wave')
    fname=['wave_pars_P',num2str(nparticles),'T',num2str(days),'days.nc'];
end

if strcmpi(Case, 'nowave')
    fname=['nowave_pars_P',num2str(nparticles),'T',num2str(days),'days.nc'];
end
disp(fname)

eval(['load ',fname(1:end-3),'traj',ini,'.mat'])
lonp=lon(tindex,:);
latp=lat(tindex,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf
figure(1)

pcolor(lon_r,lat_r,Ro);shading interp
hold on
scatter(lonp,latp,.5,'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
colormap(colorbar1)
colorbar
caxis([-2,2])
set(gca, 'Color', [.7,.7,.7]);
title([Case,'_P',num2str(nparticles),' snapshot: ',num2str(tindex),' h'], ...
    'Interpreter','none')
set(gca, 'fontsize',14,'FontWeight','b');

Ro_E=Ro;
eval(['load ',Case,'_P',num2str(nparticles),'_Rodivof',ini,'.mat'])
Ro_L15376=Ro(tindex,:);
eval(['load ',Case,'_P',num2str(2500),'_Rodivof',ini,'.mat'])
Ro_L2500=Ro(tindex,:);
eval(['load ',Case,'_P',num2str(625),'_Rodivof',ini,'.mat'])
Ro_L625=Ro(tindex,:);
eval(['load ',Case,'_P',num2str(289),'_Rodivof',ini,'.mat'])
Ro_L289=Ro(tindex,:);

x_min=-5;x_max=5;
xi = linspace(x_min, x_max, 700); % 500个点用于平滑曲线
dxi=xi(2)-xi(1);
% [pdf_L, ~] = ksdensity(Ro_L, xi, 'BoundaryCorrection', 'reflection', ...
%     'Support', [x_min, x_max]);
% [pdf_E, ~] = ksdensity(reshape(Ro_E,[1,287*287]), xi, 'BoundaryCorrection', 'reflection', ...
%     'Support', [x_min, x_max]);

[pdf_L15376, ~] = ksdensity(Ro_L15376, xi, ...
    'Support', [x_min, x_max]);
[pdf_L2500, ~] = ksdensity(Ro_L2500, xi, ...
    'Support', [x_min, x_max]);
[pdf_L625, ~] = ksdensity(Ro_L625, xi, ...
    'Support', [x_min, x_max]);
[pdf_L289, ~] = ksdensity(Ro_L289, xi, ...
    'Support', [x_min, x_max]);
[pdf_E, ~] = ksdensity(reshape(Ro_E,[1,287*287]), xi, ...
    'Support', [x_min, x_max]);
colors={'#0072BD','#D95319','#EDB120','#7E2F8E'};
figure(2)
semilogy(xi, pdf_E,'LineWidth', 2,'LineStyle','--','Color',[.7,.7,.7]);hold on
semilogy(xi, pdf_L289,'LineWidth', 2,'Color',colors{1});
semilogy(xi, pdf_L625,'LineWidth', 2,'Color',colors{2});
semilogy(xi, pdf_L2500,'LineWidth', 2,'Color',colors{3});
semilogy(xi, pdf_L15376,'LineWidth', 2,'Color',colors{4});
% plot(xi, pdf_E,'LineWidth', 2,'LineStyle','--','Color',[.7,.7,.7]);hold on
% plot(xi, pdf_L289,'LineWidth', 2,'Color',colors{1});
% plot(xi, pdf_L625,'LineWidth', 2,'Color',colors{2});
% plot(xi, pdf_L2500,'LineWidth', 2,'Color',colors{3});
% plot(xi, pdf_L15376,'LineWidth', 2,'Color',colors{4});

% plot(xi, pdf_L,'LineWidth', 2);hold on
% plot(xi, pdf_E,'LineWidth', 2,'LineStyle','--');
grid on
% ylim([1e-3,1e1])
xlim([-5,5])

% sum(dxi.*pdf_L)
%% hot pot
datestr(double(ocean_time./86400)+datenum(2000,1,1))