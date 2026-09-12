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

gamma = 1.5;
dist_bin = gamma.^(0:100) * 10; % 初始距离箱
id_stop = find(dist_bin > 1000e3, 1); % 找到第一个超过1000km的箱子
if ~isempty(id_stop)
    dist_bin = dist_bin(1:id_stop);
end
dist_bin = [0, dist_bin]; % 在开头插入0
dist_axis = 0.5 * (dist_bin(1:end-1) + dist_bin(2:end));

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
    eval(['load ',fname{ii}(1:end-3),'SF123_alltime.mat']);
    nvalid=nansum(nvaild_time,2)./nansum(nansum(nvaild_time,2));
    % nvalid2=zeros(16,1);
    % nvalid2(1,1)=nansum(nvalid(1:15));
    % for i=2:16
    %   nvaild2(i)=nvalid(i+14);
    % end
    nvalid2=[nansum(nvalid(1:15));nvalid(16:25);nansum(nvalid(26:end))];
    % semilogx(dist_axis./1e3,nvalid,'Color',colors{ii},'LineWidth',1.5);
    semilogx(dist_axis(15:26)./1e3,nvalid2,'Color',colors{ii},'LineWidth',1.5, ...
        'Marker','+');
    hold on
    grid on
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%log
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
if strcmpi(Case, 'wave')
    Case1='HF'
end

if strcmpi(Case, 'nowave')
    Case1='Smooth'
end

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

gamma = 1.5;
dist_bin = gamma.^(0:100) * 10; % 初始距离箱
id_stop = find(dist_bin > 1000e3, 1); % 找到第一个超过1000km的箱子
if ~isempty(id_stop)
    dist_bin = dist_bin(1:id_stop);
end
dist_bin = [0, dist_bin]; % 在开头插入0
dist_axis = 0.5 * (dist_bin(1:end-1) + dist_bin(2:end));

% if strcmpi(Case, 'wave')
%     fname=['wave_pars_P',num2str(nparticles),'T',num2str(days),'days.nc'];
% end
% 
% if strcmpi(Case, 'nowave')
%     fname=['nowave_pars_P',num2str(nparticles),'T',num2str(days),'days.nc'];
% end
colors={'#0072BD','#D95319','#EDB120','#7E2F8E'};
A={'a1','a2','a3','a4'}

figure(2)
for ii=1:length(nparticles)
    fname{ii,:}=[Case,'_pars_P',num2str(nparticles(ii)),'T',num2str(days),'days.nc'];
    eval(['load ',fname{ii}(1:end-3),'SF123_alltime.mat']);
    nvalid=nansum(nvaild_time,2)./nansum(nansum(nvaild_time,2));
    nvalid2=[nansum(nvalid(1:15));nvalid(16:25);nansum(nvalid(26:end))];
    % 1. 计算箱宽度（对数空间）
    log_bin_edges = log(dist_bin(15:27)); % 对数边界
    log_bin_widths = diff(log_bin_edges); % 对数宽度
    linear_bin_widths = diff(dist_bin(15:27)); % 线性宽度
    
    % 2. 计算概率密度函数
    total_pairs = sum(nvalid2);
    prob = nvalid2 / total_pairs;
    pdf_values = prob ./ linear_bin_widths'; % 关键：概率密度 = 概率 / 箱宽
    
    % 3. 创建专业对数坐标图
    % figure('Position', [100, 100, 900, 700], 'Color', 'w');
    % ax = axes('NextPlot', 'add', 'Box', 'on', 'FontSize', 12);
    
    % 4. 绘制PDF曲线
    A{ii}=loglog(dist_axis(15:26)./1e3, pdf_values, 'o-', ...
        'MarkerSize', 4, ...
        'MarkerFaceColor', colors{ii}, ...
        'LineWidth', 1.5, ...
        'Color', colors{ii});hold on
end   
grid on
set(gca, 'XScale', 'log', 'YScale', 'log', ...
    'XMinorGrid', 'on', 'YMinorGrid', 'on', ...
    'MinorGridLineStyle', '-', 'MinorGridAlpha', 0.1, ...
     'FontSize', 14,'fontweight','b');
xlabel('r (km)', 'FontSize', 14);
ylabel('p(r) (m^{-1})', 'FontSize', 14);
legend(['a1','a2','a3','a4'],{'289','625','2500','15376'})
title([Case1,': PDF of Particle Separations'], 'FontSize', 16);
xlim([1e3, 1e6]./1e3)

% xlim([min(dist_axis(dist_axis>0))./1e3, max(dist_axis)./1e3])
   