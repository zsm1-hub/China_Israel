clear
load wave_pars_P289T89.5daysSF123_c.mat
SF3_c=SF3lll+SF3ltt;

load wave_pars_P289T89.5daysSF123_alltime.mat
SF3_p=nanmean(SF3lll_time+SF3ltt_time,2);

% load wave289test.mat
load wave_pars_P289T89.5daysSF123_p.mat
SF3_p4=nanmean(SF3lll_time+SF3ltt_time,2);

load wave_pars_P625T89.5daysSF123_alltime.mat
SF3_625_p=nanmean(SF3lll_time+SF3ltt_time,2);


a1=semilogx(dist_axis,SF3_c,'LineWidth',1.5);hold on
a2=semilogx(dist_axis,SF3_p,'LineWidth',1.5);
a3=semilogx(dist_axis,SF3_p4,'LineWidth',1.5);
a4=semilogx(dist_axis,SF3_625_p,'Marker','+','LineWidth',1.5);
legend([a1,a2,a3,a4],{'chuanxing','bingx289','bingx289P4','bingx625'})


loglog(dist_axis,abs(SF3_p));hold on
loglog(dist_axis,-(SF3_p),'Marker','+');
loglog(dist_axis,(SF3_p),'Marker','o');
loglog(dist_axis,abs(SF3_625_p));
loglog(dist_axis,-(SF3_625_p),'Marker','+');
loglog(dist_axis,(SF3_625_p),'Marker','o');


loglog(dist_axis,abs(SF3_p));hold on
loglog(dist_axis,-(SF3_p),'Marker','+');
loglog(dist_axis,(SF3_p),'Marker','o');

%%
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
colors_rgb = {...
    [0, 114/255, 189/255], ...   % #0072BD
    [217/255, 83/255, 25/255], ... % #D95319
    [237/255, 177/255, 32/255], ... % #EDB120
    [126/255, 47/255, 142/255] ... % #7E2F8E
};
B={'b1','b2','b3','b4'}
C={'c1','c2','c3','c4'}

for ii=1:length(nparticles)
    fname{ii,:}=[Case,'_pars_P',num2str(nparticles(ii)),'T',num2str(days),'days.nc'];
    eval(['load ',fname{ii}(1:end-3),'SF123_alltime.mat']);
    % eval(['load ',fname{ii}(1:end-3),'SF123.mat']);
    % if ii~=2
    %     SF3=(SF3lll_time+SF3ltt_time);
    % else
    %     SF3=(SF3lll+SF3ltt)';
    % end
    SF3=(SF3lll_time+SF3ltt_time);
    dstr=14;
    for tt=1:size(SF3,2)
         [SpecFlux(:,tt),Vt(:,tt),ebs(:,tt),kf,lf]=Fk_fitting_SF3(SF3(dstr:end,tt),dist_axis(dstr:end), ...
             1,300e3,'log','RLS',lambda);
    end

    for i = 1:size(SpecFlux,1)    
        CI_SpecFlux(:,i) = prctile(SpecFlux(i,:), [95,5]);
    end
    % [SpecFlux,Vt,ebs,kf,lf]=Fk_fitting_SF3(nanmean(SF3,2),dist_axis,2,300e3,'log','RLS',lambda);
    % std_over_time = std(data, 'omitnan', 2);
    % shadedErrorBar_semilogx((1./kf)/1e3, mean_ebs(1:end-1), CI_ebs(:,1:end-1) ...
    %             ,{'o-','linewidth',2,'color', colors(1,:)}, 0.6)

    % [SpecFlux,Vt,ebs,kf,lf]=Fk_fitting_SF3(SF3,dist_axis,2,300e3,'log','RLS',lambda);
    figure(1)
    subplot(2,2,ii)
    B{ii}=semilogx(1./kf./1e3,nanmean(SpecFlux,2),'Marker','x','Color',colors{ii}, ...
        'LineWidth',1.5);
    % semilogx(1./kf./1e3,(SpecFlux),'Marker','x','Color',colors{ii},'LineWidth',1.5);

    hold on
    % clear x_fill;clear y_fill
    x_fill = [1./kf./1e3, fliplr(1./kf./1e3)];
    y_fill = [CI_SpecFlux(1,:), fliplr(CI_SpecFlux(2,:))];
    fill(x_fill, y_fill, colors_rgb{ii}, 'FaceAlpha', 0.3, 'EdgeColor', 'none');

    eval(['load ',Case,'_pars_P',num2str(nparticles(ii)),'T89.5daysCG_Lag_tukey.mat']);
    C{ii}=semilogx(filtscale./1e3,Th','LineWidth',1.5,'Color',colors{ii},'LineStyle','--');
    % semilogx(filtscale./1e3,Thm_Eulerian','LineWidth',1.5,'Color',[.7, .7, .7],'LineStyle','-');
    grid on
    ylim([-5e-7,5e-7]);
    xlim([1e3,1e6]./1e3);
    xlabel('km')
    ylabel('m^{2}/s^{3}')
    title([Case,': cg flux v.s. sf3-fitting fk'])
    legend([B{ii},C{ii}],{'SF3-fitting fk','Lag cg-spec fk'},'Location','southeast');
    set(gca,'fontsize',16,'FontWeight','b')

    SF3t=nanmean(SF3,2);
    figure(4)
    A{ii}=loglog(dist_axis./1e3,abs(SF3t),'Color',colors{ii},'LineWidth',1.5)
    % semilogx(dist_axis./1e3,SF1L,'Marker','x','Color',colors{ii},'LineWidth',1.5,'LineStyle','--')
    % semilogx(dist_axis./1e3,SF1T,'Marker','x','Color',colors{ii},'LineWidth',1.5,'LineStyle','-.')
    hold on
    semilogx(dist_axis./1e3,-(SF3t),'Color',colors{ii},'LineWidth',1.5,'Marker','+')
    semilogx(dist_axis./1e3,(SF3t),'Color',colors{ii},'LineWidth',1.5,'Marker','o')
    grid on
    
    xlim([1e3,1e6]./1e3);
    ylim([1e-5,1e-2]);
    xlabel('km')
    ylabel('m^{3}/s^{3}')
    title([Case,': SF3'])
    set(gca,'fontsize',16,'FontWeight','b')
end


    figure(2)
    semilogx(1./kf./1e3,ebs(1:end-1),'Marker','x','Color',colors{ii},'LineWidth',1.5)
    hold on
    grid on
    xlim([1e3,1e6]./1e3);
    ylim([-2.5e-3,2.5e-3]);
    xlabel('km')
    ylabel('m^{2}/s^{3}')
    title([Case,': ebs'])
    set(gca,'fontsize',16,'FontWeight','b')
    end

% std(SpecFlux,0, 'omitnan', 2);
% std(SpecFlux, 0, 'omitnan', 2);
std1=std(SpecFlux, 0, 2, 'omitnan');


for i = 1:size(SpecFlux,1)    
    CI_SpecFlux(:,i) = prctile(SpecFlux(i,:), [95,5]);
end
% shadedErrorBar_semilogx((1./kf)/1e3, nanmean(SpecFlux,2), CI_SpecFlux...
%                 ,{'o-','linewidth',2,'color', colors{ii}}, 0.6)
% semilogx(1./kf./1e3,nanmean(SpecFlux,2));hold on
% semilogx(1./kf./1e3,CI_SpecFlux(1,:));
% semilogx(1./kf./1e3,CI_SpecFlux(2,:));

semilogx(1./kf./1e3,nanmean(SpecFlux,2));hold on
fill([1./kf./1e3, fliplr(1./kf./1e3)], [CI_SpecFlux(1,:), fliplr(CI_SpecFlux(2,:))], ...
    'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([1./kf./1e3, fliplr(1./kf./1e3)], [[nanmean(SpecFlux,2)+std1]', ...
    [fliplr(nanmean(SpecFlux,2)-std1)]'], ...
    'g', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
eval(['load ',Case,'_pars_P',num2str(nparticles(ii)),'T89.5daysCG_Lag_tukey.mat']);
    semilogx(filtscale./1e3,Th','LineWidth',1.5,'Color',colors{ii},'LineStyle','--');
    % semilogx(filtscale./1e3,Thm_Eulerian','LineWidth',1.5,'Color',[.7, .7, .7],'LineStyle','-');
    grid on
    ylim([-0.5e-7,0.5e-7]);
    xlim([1e3,1e6]./1e3);
    xlabel('km')
    ylabel('m^{2}/s^{3}')
    title([Case,': cg flux v.s. sf3-fitting fk'])
    set(gca,'fontsize',16,'FontWeight','b')

%%

clear
fname1='wave_pars_P289T89.5days.nc'
fname2='wave_pars_P289T89.5days_extratime.nc'

lon1=ncread(fname1,'lon');lat1=ncread(fname1,'lat');
u1=ncread(fname1,'ue');v1=ncread(fname1,'ve');

lon2=ncread(fname2,'lon');lat2=ncread(fname2,'lat');

figure(1)
hold on
traj=250
ts=100
tt=2000
% for ii=1:size(lat1,2)
    plot(lon1(ts:tt,traj),lat1(ts:tt,traj));
% end

for ii=1:size(lat1,2)
    plot(lon1(ts:tt,ii),lat1(ts:tt,ii));
end