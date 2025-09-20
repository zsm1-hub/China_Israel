


clear all;close all;clc
% addpath('D:\LIN2023\model\RoyBarkan\LLC4320/')
% addpath('D:\LIN2023\crocotools\Preprocessingtools') % add function "spheric_dist.m"
% 
addpath('/meddy/simingzhang/Analysis/matlab/Parcels_SF/')
addpath('/meddy/simingzhang/Data/Parcels_data')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          1. Basic setup and read data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Case='nowave'; % wave
nparticles=[289,625,2500,15376]; % numbers of particles
% nparticles=[289,625,2500]; % numbers of particles

% nparticles=[289]; % numbers of particles

days=89.5;  % days
dt=3600; % s  Advection_RK4 delta_t drift时间间隔
% input_dir='D:\LIN2023\model\RoyBarkan\LLC4320/'; % drift所在文件夹
input_dir='/meddy/simingzhang/Data/Parcels_data/';
% timerange=24*10:24*11-6; % 计算结构函数用的时间范围
% timerange=1:2140;
inv_style='RLS';
lambda=1e-10;
ini='_grid'

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
    eval(['load ',fname{ii}(1:end-3),'SF123',ini,'.mat']);
    if ii>2
        % SF3=(SF3lll_time+SF3ltt_time)';
        SF3=(SF3lll_time+SF3ltt_time);
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
    dstr=14;
    % dstr=1;
    rescale=1;
    figure(1)
    subplot(2,2,ii)
    B{ii}=semilogx(1./kf(dstr:end)./1e3.*rescale,SpecFlux(dstr:end), ...
        'Marker','x','Color',colors{ii},'LineWidth',1.5);
    hold on
    
    eval(['load ',Case,'_pars_P',num2str(nparticles(ii)),'T89.5daysCG_Lag_spectukey', ...
        ini,'.mat']);
    C{ii}=semilogx(filtscale./1e3,Th','LineWidth',1.5,'Color',colors{ii},'LineStyle','--');
    % eval(['load ',Case,'_pars_P',num2str(nparticles(ii)),'T89.5daysCG_Lag_uni.mat']);
    % F{ii}=semilogx(filtscale./1e3,Th','LineWidth',1.5,'Color',colors{ii},'LineStyle',':');
    % eval(['load ',Case,'_Eulerian_SF3.mat'])
    % N=287;dstr=14;dot=202;
    % [r,SF3,S3L1,S3T1]=calc_radial(S3L1,S3T1,N,xscale);
    % [SpecFlux,Vt,ebs,kf,lf]=Fk_fitting_SF3(SF3(1:dot)',r(1:dot),1,200e3,'log','RLS',lambda);
    % F{ii}=semilogx(1./kf(dstr:end)./1e3,SpecFlux(dstr:end),'LineWidth',1.5, ...
    % 'Color','k');
    % semilogx(filtscale./1e3,Thm_Eulerian','LineWidth',1.5,'Color',[.7, .7, .7],'LineStyle','-');
    grid on
    legend([B{ii},C{ii}],{'SF3-fitting fk','Lag cg-spec fk'},'Location','southeast');
    % if ii>2
    %     ylim([-0.3e-7,0.3e-7]);
    % else
        ylim([-1e-7,1.5e-7]);
        
    % end
    xlim([1e3,1e6]./1e3);
    xlabel('km')
    ylabel('m^{2}/s^{3}')
    title([Case,'P',num2str(nparticles(ii)),': cg v.s. sf3-fitting'])
    set(gca,'fontsize',16,'FontWeight','b')

    figure(2)
    % semilogx(1./kf./1e3,ebs(1:end-1).*kf','Marker','x','Color',colors{ii},'LineWidth',1.5)
    semilogx(1./kf./1e3,ebs(1:end-1),'Marker','x','Color',colors{ii},'LineWidth',1.5)
    hold on
    grid on
    xlim([1e3,1e6]./1e3);
    % ylim([-8e-8,8e-8]);
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

% saveas(gcf,[Case,'_fk_Lag.png'],'png')

%% all E and L
clear;close all
Case='wave';
ini='_grid';
if strcmpi(Case, 'wave')
    fname='s2sflux_spec_hf.0002.nc';
end

if strcmpi(Case, 'nowave')
    fname='s2sflux_spec_smooth.0002.nc';
end

ncdisp(fname)
Thm_Eulerian=ncread(fname,'Thm');
filtscale=ncread(fname,'filtscale');
lambda=1e-10;
ii=4
nparticles=[289,625,2500,15376];

% wave_pars_P625T89.5daysCG_Lag_tukey.mat

jj=1;
colors={[.7,.7,.7],'#0072BD','#D95319','#EDB120','#7E2F8E'};
a1=semilogx(filtscale./1e3,Thm_Eulerian,'LineWidth',1.5,'Color',colors{jj});
hold on


eval(['load ',Case,'_Eulerian_SF3.mat'])
N=287;dstr=6;dot=202;rescale=2;xscale=xscale.*sqrt(2);

[r,SF3,S3L1,S3T1]=calc_radial(S3L1,S3T1,N,xscale);
[SpecFlux,Vt,ebs,kf,lf]=Fk_fitting_SF3(SF3(1:dot)',r(1:dot),1,200e3,'log','RLS',lambda);
a6=semilogx(1./kf(dstr:end)./1e3,SpecFlux(dstr:end),'LineWidth',1.5, ...
'Color','k');
% a7=semilogx(1./kf(dstr:end)./1e3.*rescale,SpecFlux(dstr:end),'LineWidth',1.5, ...
% 'Color','k','LineStyle','--');

% nparticles=[15376];
clear fname
fname{ii,:}=[Case,'_pars_P',num2str(nparticles(ii)),'T',num2str(89.5),'days.nc'];
eval(['load ',fname{ii}(1:end-3),'SF123',ini,'.mat']);
if ii>2
    SF3=(SF3lll_time+SF3ltt_time);
    SF1=(SF1l_time+SF1t_time)';
    SF1L=(SF1l_time)';
    SF1T=(SF1t_time)';
else
    SF3=(SF3lll+SF3ltt)';
    SF1=(SF1l+SF1t)';
    SF1L=(SF1l)';
    SF1T=(SF1t)';

end
[SpecFlux,Vt,ebs,kf,lf]=Fk_fitting_SF3(SF3,dist_axis,1,300e3,'log','RLS',lambda);
dstr=14;
% dstr=1;
a4=semilogx(1./kf(dstr:end)./1e3,SpecFlux(dstr:end),'Marker','x','Color',colors{ii+1},'LineWidth',1.5);
eval(['load ',Case,'_pars_P',num2str(nparticles(ii)),'T89.5daysCG_Lag_spectukey',ini,'.mat']);
a2=semilogx(filtscale./1e3,Th','LineWidth',1.5,'Color',colors{ii+1},'LineStyle','--');

grid on
legend(['a1','a2','a4','a6'],{'Eulerian cg fk','Eulerian-SF3-fitting fk', ...
    'P15376 Lag SF3-fitting fk','P15376 Lag cg fk'}, ...
    'Location', 'northeast')
ylim([-0.2e-7,0.5e-7]);
xlim([1e3,1e6]./1e3);
xlabel('km')
ylabel('m^{2}/s^{3}')
title('Smooth case')
set(gca,'fontsize',16,'FontWeight','b')
%% test
clear; close all
lambda=1e-10;
dstr=14;
rescale=1
colors={'#0072BD','#D95319','#EDB120','#7E2F8E'};

fname='s2sflux_spec_hf.0002.nc';
ncdisp(fname)
filtscale=ncread(fname,'filtscale');

figure(1)
ii=1
load wave_pars_P15376T89.5daysSF123_grid.mat
SF3=(SF3lll_time+SF3ltt_time);
[SpecFlux,Vt,ebs,kf,lf]=Fk_fitting_SF3(SF3,dist_axis,1,300e3,'log','RLS',lambda);
a4=semilogx(1./kf(dstr:end)./1e3.*rescale,SpecFlux(dstr:end),'Marker','x','LineWidth',1.5, ...
    'Color',colors{ii});
hold on

load wave_pars_P15376T89.5daysCG_Lag_spectukey_grid.mat
a2=semilogx(filtscale./1e3,Th','LineWidth',1.5,'Color',colors{ii},'LineStyle','--');

%%%%%%%%%%%%%%%%%%%
ii=2
load nowave_pars_P15376T89.5daysSF123_grid.mat
SF3=(SF3lll_time+SF3ltt_time);
[SpecFlux,Vt,ebs,kf,lf]=Fk_fitting_SF3(SF3,dist_axis,1,300e3,'log','RLS',lambda);

load nowave_pars_P15376T89.5daysCG_Lag_spectukey_grid.mat
a3=semilogx(filtscale./1e3,Th','LineWidth',1.5,'Color',colors{ii},'LineStyle','--');

a5=semilogx(1./kf(dstr:end)./1e3.*rescale,SpecFlux(dstr:end),'Marker','x','LineWidth', ...
    1.5,'Color',colors{ii});
grid on
ylim([-0.2e-7,0.5e-7]);
xlim([1e3,1e6]./1e3);
xlabel('km')
ylabel('m^{2}/s^{3}')
% title('Smooth case')
set(gca,'fontsize',16,'FontWeight','b')
%% check hit2d
clear;close all
fname='s2sflux_spec_hit_tukey.0002.nc';
ncdisp(fname)
Thm_Eulerian=ncread(fname,'Thm');
filtscale=ncread(fname,'filtscale');
lambda=1e-10;
ii=4
jj=1
rescale=1;
colors={[.7,.7,.7],'#0072BD','#D95319','#EDB120','#7E2F8E'};
colors_rgb = {...
    [.7,.7,.7],...
    [0, 114/255, 189/255], ...   % #0072BD
    [217/255, 83/255, 25/255], ... % #D95319
    [237/255, 177/255, 32/255], ... % #EDB120
    [126/255, 47/255, 142/255] ... % #7E2F8E
};
figure(1)
semilogx(filtscale,Thm_Eulerian,'LineWidth',1.5,'Color',colors{jj});hold on
fname='s2sflux_spec_hit_uni.0002.nc';
ncdisp(fname)
Thm_Eulerian=ncread(fname,'Thm');
filtscale=ncread(fname,'filtscale');
semilogx(filtscale,Thm_Eulerian,'LineWidth',1.5,'Color','k');

grid on

load HIT2d_Eulerian_SF3.mat
lambda=8e-1;
dot=362;
% [SpecFlux,Vt,ebs,kf,lf]=Fk_fitting_SF3(SF3(1:dot)',r(1:dot),0.009,6.32,'log','RLS',lambda);

[SpecFlux,Vt,ebs,kf,lf]=Fk_fitting_SF3(SF3(1:dot)',r(1:dot),0.009,6.32,'log','RLS',lambda);


a6=semilogx(1./kf.*rescale,SpecFlux,'LineWidth',1.5, ...
'Color',colors{jj+1});


figure(2)
loglog(r,abs(SF3),'Color',colors{ii},'LineWidth',1.5)
% semilogx(dist_axis./1e3,SF1L,'Marker','x','Color',colors{ii},'LineWidth',1.5,'LineStyle','--')
% semilogx(dist_axis./1e3,SF1T,'Marker','x','Color',colors{ii},'LineWidth',1.5,'LineStyle','-.')
hold on
loglog(r,-(SF3),'Color',colors{ii},'LineWidth',1.5,'Marker','+')
loglog(r,(SF3),'Color',colors{ii},'LineWidth',1.5,'Marker','o')
[x32,y32]=get_line_loglog(3/2,1e-2,1e-2,-2,0);
loglog(x32,y32,'LineWidth',1.5,'color','k')
[x45,y45]=get_line_loglog(4/5,1e-2,1e-2,-2,0);
loglog(x45,y45,'LineWidth',1.5,'color','b')
grid on

figure(3)
%div
SF3or=(SF3(2:end)-SF3(1:end-1))./(r(2:end)-r(1:end-1));
r2=0.5.*(r(2:end)+r(1:end-1));
loglog(r2,SF3or,'LineWidth',1.5,'color','b')

figure(4)
semilogx(r,(SF3),'Color',colors{ii},'LineWidth',1.5);hold on
semilogx(r,nanmean(SF3_time,1),'Color',colors{ii+1},'LineWidth',1.5);

figure(5)
for j=1:199
    [SpecFlux(:,j),Vt(:,j),ebs(:,j),kf,lf]=Fk_fitting_SF3(SF3_time(j,1:dot)', ...
        r(1:dot),0.009,6.32,'log','RLS',lambda);
end
std1=std(SpecFlux, 0, 2, 'omitnan');
x_fill = [1./kf, fliplr(1./kf)];
% y_fill = [CI_SpecFlux(1,:), fliplr(CI_SpecFlux(2,:))];
y_fill = [(nanmean(SpecFlux,2)+std1)', fliplr((nanmean(SpecFlux,2)-std1)')];
semilogx(1./kf,nanmean(SpecFlux,2),'Color',colors{ii},'LineWidth',1.5);hold on
fill(x_fill, y_fill, colors_rgb{ii}, 'FaceAlpha', 0.3, 'EdgeColor', 'none')

semilogx(filtscale,Thm_Eulerian,'LineWidth',1.5,'Color',colors{jj});hold on
fname='s2sflux_spec_hit_uni.0002.nc';
ncdisp(fname)
Thm_Eulerian=ncread(fname,'Thm');
filtscale=ncread(fname,'filtscale');
semilogx(filtscale,Thm_Eulerian,'LineWidth',1.5,'Color','k');

fname='s2sflux_spec_hit_tukey.0002.nc';
ncdisp(fname)
Thm_Eulerian=ncread(fname,'Thm');
filtscale=ncread(fname,'filtscale');
semilogx(filtscale,Thm_Eulerian,'LineWidth',1.5,'Color',[.7,.7,.7]);
grid on
ylim([-20,20])

% load HIT2d_pars_P2500T0.2secondsCG_Lag_uni_grid.mat
% a3=semilogx(filtscale,Th','LineWidth',1.5,'Color',colors{ii+1},'LineStyle','--');

% a7=semilogx(1./kf(dstr:end)./1e3.*rescale,SpecFlux(dstr:end),'LineWidth',1.5, ...
% 'Color','k','LineStyle','--');

% nparticles=[15376];
% clear fname
% fname{ii,:}=[Case,'_pars_P',num2str(nparticles(ii)),'T',num2str(89.5),'days.nc'];
% eval(['load ',fname{ii}(1:end-3),'SF123',ini,'.mat']);
% if ii>2
%     SF3=(SF3lll_time+SF3ltt_time);
%     SF1=(SF1l_time+SF1t_time)';
%     SF1L=(SF1l_time)';
%     SF1T=(SF1t_time)';
% else
%     SF3=(SF3lll+SF3ltt)';
%     SF1=(SF1l+SF1t)';
%     SF1L=(SF1l)';
%     SF1T=(SF1t)';
% 
% end
% [SpecFlux,Vt,ebs,kf,lf]=Fk_fitting_SF3(SF3,dist_axis,1,300e3,'log','RLS',lambda);
% dstr=14;
% % dstr=1;
% a4=semilogx(1./kf(dstr:end)./1e3,SpecFlux(dstr:end),'Marker','x','Color',colors{ii+1},'LineWidth',1.5);
% eval(['load ',Case,'_pars_P',num2str(nparticles(ii)),'T89.5daysCG_Lag_spectukey',ini,'.mat']);
% a2=semilogx(filtscale./1e3,Th','LineWidth',1.5,'Color',colors{ii+1},'LineStyle','--');
% 
% grid on
% legend(['a1','a2','a4','a6'],{'Eulerian cg fk','Eulerian-SF3-fitting fk', ...
%     'P15376 Lag SF3-fitting fk','P15376 Lag cg fk'}, ...
%     'Location', 'northeast')
% ylim([-0.2e-7,0.5e-7]);
% xlim([1e3,1e6]./1e3);
% xlabel('km')
% ylabel('m^{2}/s^{3}')
% title('Smooth case')
% set(gca,'fontsize',16,'FontWeight','b')