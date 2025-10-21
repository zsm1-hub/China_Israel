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
win='kaiser'
ini='rough'
gamma = 1.5;

dist_bin(1) = 10; % in m
dist_bin = gamma.^[0:100]*dist_bin(1);
id = find(dist_bin>1000*10^3,1);
dist_bin = dist_bin(1:id);
dist_bin(2:end+1) = dist_bin(1:end);
dist_bin(1) = 0;
dist_axis = 0.5*(dist_bin(1:end-1) + dist_bin(2:end));


timerange=1:2148;
param_Eul_cg_sf3fk
lambda=[10,1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6, ...
    1e-7,1e-8,1e-9,1e-10,1e-11,1e-12];
dend=9;
% all time
% load wave_pars_P289T89.5daysSF123_alltime_rough.mat
% SF3_Lag=nanmean(SF3lll_time+SF3ltt_time,2);
% SF3_Lag_all=(SF3lll_time+SF3ltt_time);
% 
% for jj=1:size(SF3lll_time,2)
% [SpecFlux_Lag(:,jj),Vt_Lag,ebs_Lag,kf_Lag,...
%     lf_Lag]=Fk_fitting_SF3_Lcurve(SF3_Lag_all(:,jj), ...
%     dist_axis,2e3,500e3,'log','RLS',lambda,kf1,0);
% end

% 480 h bootstrap
timerange=1:720;

[SpecFlux_Lag289,Vt_Lag289,ebs_Lag289,~,...
    ~,CI_ebs289,CI_Vt289,...
    CI_SpecFlux289,Th_Lag289]=get_param_Lag_sf3fk_bootstrap(Case,289,lambda,timerange);
fname=[Case,'_pars_P289','T',num2str(timerange(end)),'bootstrap.mat'];
load(fname);
SF3_Lag289=nanmean(SF3,2);

[SpecFlux_Lag625,Vt_Lag625,ebs_Lag625,~,...
    ~,CI_ebs625,CI_Vt625,...
    CI_SpecFlux625,Th_Lag625]=get_param_Lag_sf3fk_bootstrap(Case,625,lambda,timerange);
fname=[Case,'_pars_P625','T',num2str(timerange(end)),'bootstrap.mat'];
load(fname);
SF3_Lag625=nanmean(SF3,2);

[SpecFlux_Lag2500,Vt_Lag2500,ebs_Lag2500,kf_Lag,...
    lf_Lag,CI_ebs2500,CI_Vt2500,...
    CI_SpecFlux2500,Th_Lag2500]=get_param_Lag_sf3fk_bootstrap(Case,2500,lambda,timerange);
fname=[Case,'_pars_P2500','T',num2str(timerange(end)),'bootstrap.mat'];
load(fname);
SF3_Lag2500=nanmean(SF3,2);

x_fill = [1./kf_Lag(1:dend)./1e3.*2.*pi, fliplr(1./kf_Lag(1:dend)./1e3.*2.*pi)];

% if timerange(end)<500
%     if strcmpi(Case, 'wave')
%         load wave_CG_Eul_T480_rough.mat
% 
%     end
%     if strcmpi(Case, 'nowave')
%         load nowave_CG_Eul_T480_rough.mat
%     end
%     Thm_Eulerian=nanmean(Th_all480,2);
% end


% b1=semilogx(filtscale./1e3,Thm_Eulerian,'LineWidth',1.5,'Color',colors{1});
% r_cg = findXatYZero(filtscale./1e3, Thm_Eulerian);

% b2=semilogx(1./kf_E(dstr:end)./1e3.*2.*pi,SpecFlux_E(dstr:end), ...
%     'LineWidth',1.5, 'Color','k');
% r_fit = findXatYZero(1./kf_E(dstr:end)./1e3.*2.*pi, SpecFlux_E(dstr:end));

x_coords = [4, 200, 200, 4];  % x坐标：左边界4，右边界200
y_coords = [-4e-8, -4e-8, 4e-8, 4e-8]; % y坐标：下边界-2e-8，上边界4e-8
% 
% patch(x_coords, y_coords, [0.9 0.9 0.9], 'EdgeColor', 'none'); % 中灰色，无边线
screenSize = get(0, 'ScreenSize');
% 创建figure并设置位置
figure('Position', [0, 0, screenSize(3), screenSize(4)]);
subplot(2,3,1)

patch(x_coords, y_coords, [0.9 0.9 0.9],'FaceAlpha', 0.3, ...
    'EdgeColor', 'k'); 
set(gca,'xscale','log')
hold on

b1=semilogx(dist_axis./1e3,SF3_Lag289./dist_axis', ...
    'LineWidth',1.5,'Color',colors{4});
b2=semilogx(lf_Lag./1e3,Vt_Lag289./lf_Lag', ...
    'LineWidth',1.5,'Color',colors{1},'LineStyle','none','Marker','+');

grid on
ylim([-4e-8,4e-8]);
xlim([1e3,1e6]./1e3);
xtick([1,4,1e1,1e2,200,1e3])
% xticklabels({'10^{0}','4','10^{1}','10^{2}','200','10^{3}'})

xlabel('$$\mathbf{r \ [km]}$$','Interpreter','latex')
ylabel('$$\mathbf{D3(r)/r \ [m^{2}/s^{3}]}$$','Interpreter','latex')
% % title('Hf: CG vs sf3-RLS flux')
% legend([b1,b2],{['CG (',win,')'],'SF3-fitting (RLS)'})
legend([b1,b2],{['Lag D3(r)/r'],'RLS-fitting'})
text(0.03, 0.95, 'a) Hf P289', 'Units', 'normalized', ...
     'FontSize', 12, 'FontWeight', 'bold', ...
     'BackgroundColor', [1, 1, 0.8, 0.6], ... % 半透明背景
     'EdgeColor', [.7,.7,.7], ... % 边框颜色
     'Margin', 3, ... % 边距
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')

set(gca,'fontsize',12,'fontweight','bold')

subplot(2,3,2)
patch(x_coords, y_coords, [0.9 0.9 0.9],'FaceAlpha', 0.3, ...
    'EdgeColor', 'k'); 
set(gca,'xscale','log')
hold on

b1=semilogx(dist_axis./1e3,SF3_Lag625./dist_axis', ...
    'LineWidth',1.5,'Color',colors{5});
b2=semilogx(lf_Lag./1e3,Vt_Lag625./lf_Lag', ...
    'LineWidth',1.5,'Color',colors{1},'LineStyle','none','Marker','+');

grid on
ylim([-4e-8,4e-8]);
xlim([1e3,1e6]./1e3);
xtick([1,4,1e1,1e2,200,1e3])
xlabel('$$\mathbf{r \ [km]}$$','Interpreter','latex')
ylabel('$$\mathbf{D3(r)/r \ [m^{2}/s^{3}]}$$','Interpreter','latex')
% % title('Hf: CG vs sf3-RLS flux')
% legend([b1,b2],{['CG (',win,')'],'SF3-fitting (RLS)'})
legend([b1,b2],{['Lag D3(r)/r'],'RLS-fitting'})
text(0.03, 0.95, 'b) Hf P625', 'Units', 'normalized', ...
     'FontSize', 12, 'FontWeight', 'bold', ...
     'BackgroundColor', [1, 1, 0.8, 0.6], ... % 半透明背景
     'EdgeColor', [.7,.7,.7], ... % 边框颜色
     'Margin', 3, ... % 边距
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')
set(gca,'fontsize',12,'fontweight','bold')

subplot(2,3,3)
patch(x_coords, y_coords, [0.9 0.9 0.9],'FaceAlpha', 0.3, ...
    'EdgeColor', 'k'); 
set(gca,'xscale','log')
hold on

b1=semilogx(dist_axis./1e3,SF3_Lag625./dist_axis', ...
    'LineWidth',1.5,'Color',colors{6});
b2=semilogx(lf_Lag./1e3,Vt_Lag625./lf_Lag', ...
    'LineWidth',1.5,'Color',colors{1},'LineStyle','none','Marker','+');

grid on
ylim([-4e-8,4e-8]);
xlim([1e3,1e6]./1e3);
xtick([1,4,1e1,1e2,200,1e3])
xlabel('$$\mathbf{r \ [km]}$$','Interpreter','latex')
ylabel('$$\mathbf{D3(r)/r \ [m^{2}/s^{3}]}$$','Interpreter','latex')


legend([b1,b2],{['Lag D3(r)/r'],'RLS-fitting'})
text(0.03, 0.95, 'c) Hf P2500', 'Units', 'normalized', ...
     'FontSize', 12, 'FontWeight', 'bold', ...
     'BackgroundColor', [1, 1, 0.8, 0.6], ... % 半透明背景
     'EdgeColor', [.7,.7,.7], ... % 边框颜色
     'Margin', 3, ... % 边距
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')

set(gca,'fontsize',12,'fontweight','bold')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% nowave 
Case='nowave'
timerange=1:2148;
param_Eul_cg_sf3fk
lambda=[10,1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6, ...
    1e-7,1e-8,1e-9,1e-10,1e-11,1e-12];
dend=9;
% all time
% load wave_pars_P289T89.5daysSF123_alltime_rough.mat
% SF3_Lag=nanmean(SF3lll_time+SF3ltt_time,2);
% SF3_Lag_all=(SF3lll_time+SF3ltt_time);
% 
% for jj=1:size(SF3lll_time,2)
% [SpecFlux_Lag(:,jj),Vt_Lag,ebs_Lag,kf_Lag,...
%     lf_Lag]=Fk_fitting_SF3_Lcurve(SF3_Lag_all(:,jj), ...
%     dist_axis,2e3,500e3,'log','RLS',lambda,kf1,0);
% end

% 480 h bootstrap
timerange=1:720;

[SpecFlux_Lag289,Vt_Lag289,ebs_Lag289,~,...
    ~,CI_ebs289,CI_Vt289,...
    CI_SpecFlux289,Th_Lag289]=get_param_Lag_sf3fk_bootstrap(Case,289,lambda,timerange);
fname=[Case,'_pars_P289','T',num2str(timerange(end)),'bootstrap.mat'];
load(fname);
SF3_Lag289=nanmean(SF3,2);

[SpecFlux_Lag625,Vt_Lag625,ebs_Lag625,~,...
    ~,CI_ebs625,CI_Vt625,...
    CI_SpecFlux625,Th_Lag625]=get_param_Lag_sf3fk_bootstrap(Case,625,lambda,timerange);
fname=[Case,'_pars_P625','T',num2str(timerange(end)),'bootstrap.mat'];
load(fname);
SF3_Lag625=nanmean(SF3,2);

[SpecFlux_Lag2500,Vt_Lag2500,ebs_Lag2500,kf_Lag,...
    lf_Lag,CI_ebs2500,CI_Vt2500,...
    CI_SpecFlux2500,Th_Lag2500]=get_param_Lag_sf3fk_bootstrap(Case,2500,lambda,timerange);
fname=[Case,'_pars_P2500','T',num2str(timerange(end)),'bootstrap.mat'];
load(fname);
SF3_Lag2500=nanmean(SF3,2);



subplot(2,3,4)
patch(x_coords, y_coords, [0.9 0.9 0.9],'FaceAlpha', 0.3, ...
    'EdgeColor', 'k'); 
set(gca,'xscale','log')
hold on
b1=semilogx(dist_axis./1e3,SF3_Lag289./dist_axis', ...
    'LineWidth',1.5,'Color',colors{4});
b2=semilogx(lf_Lag./1e3,Vt_Lag289./lf_Lag', ...
    'LineWidth',1.5,'Color',colors{2},'LineStyle','none','Marker','+');

grid on
ylim([-4e-8,4e-8]);
xlim([1e3,1e6]./1e3);
xtick([1,4,1e1,1e2,200,1e3])

xlabel('$$\mathbf{r \ [km]}$$','Interpreter','latex')
ylabel('$$\mathbf{D3(r)/r \ [m^{2}/s^{3}]}$$','Interpreter','latex')
% % title('Hf: CG vs sf3-RLS flux')
% legend([b1,b2],{['CG (',win,')'],'SF3-fitting (RLS)'})
legend([b1,b2],{['Lag D3(r)/r'],'RLS-fitting'})
text(0.03, 0.95, 'd) Sm P289', 'Units', 'normalized', ...
     'FontSize', 12, 'FontWeight', 'bold', ...
     'BackgroundColor', [1, 1, 0.8, 0.6], ... % 半透明背景
     'EdgeColor', [.7,.7,.7], ... % 边框颜色
     'Margin', 3, ... % 边距
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')
set(gca,'fontsize',12,'fontweight','bold')

subplot(2,3,5)
patch(x_coords, y_coords, [0.9 0.9 0.9],'FaceAlpha', 0.3, ...
    'EdgeColor', 'k'); 
set(gca,'xscale','log')
hold on
b1=semilogx(dist_axis./1e3,SF3_Lag625./dist_axis', ...
    'LineWidth',1.5,'Color',colors{5});
b2=semilogx(lf_Lag./1e3,Vt_Lag625./lf_Lag', ...
    'LineWidth',1.5,'Color',colors{2},'LineStyle','none','Marker','+');

grid on
ylim([-4e-8,4e-8]);
xlim([1e3,1e6]./1e3);
xtick([1,4,1e1,1e2,200,1e3])

xlabel('$$\mathbf{r \ [km]}$$','Interpreter','latex')
ylabel('$$\mathbf{D3(r)/r \ [m^{2}/s^{3}]}$$','Interpreter','latex')

legend([b1,b2],{['Lag D3(r)/r'],'RLS-fitting'})
text(0.03, 0.95, 'e) Sm P625', 'Units', 'normalized', ...
     'FontSize', 12, 'FontWeight', 'bold', ...
     'BackgroundColor', [1, 1, 0.8, 0.6], ... % 半透明背景
     'EdgeColor', [.7,.7,.7], ... % 边框颜色
     'Margin', 3, ... % 边距
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')

set(gca,'fontsize',12,'fontweight','bold')

subplot(2,3,6)
patch(x_coords, y_coords, [0.9 0.9 0.9],'FaceAlpha', 0.3, ...
    'EdgeColor', 'k'); 
set(gca,'xscale','log')
hold on
b1=semilogx(dist_axis./1e3,SF3_Lag2500./dist_axis', ...
    'LineWidth',1.5,'Color',colors{6});
b2=semilogx(lf_Lag./1e3,Vt_Lag2500./lf_Lag', ...
    'LineWidth',1.5,'Color',colors{2},'LineStyle','none','Marker','+');

grid on
ylim([-4e-8,4e-8]);
xlim([1e3,1e6]./1e3);
xtick([1,4,1e1,1e2,200,1e3])

xlabel('$$\mathbf{r \ [km]}$$','Interpreter','latex')
ylabel('$$\mathbf{D3(r)/r \ [m^{2}/s^{3}]}$$','Interpreter','latex')

legend([b1,b2],{['Lag D3(r)/r'],'RLS-fitting'})
text(0.03, 0.95, 'f) Sm P2500', 'Units', 'normalized', ...
     'FontSize', 12, 'FontWeight', 'bold', ...
     'BackgroundColor', [1, 1, 0.8, 0.6], ... % 半透明背景
     'EdgeColor', [.7,.7,.7], ... % 边框颜色
     'Margin', 3, ... % 边距
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')

set(gca,'fontsize',12,'fontweight','bold')
saveas(gcf,['Eul_CG',win,'_vs_RLS_SF3or'],'png')
