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
timerange=1:720;
[SpecFlux_Lag625,Vt_Lag625,ebs_Lag625,~,...
    ~,CI_ebs625,CI_Vt625,...
    CI_SpecFlux625,Th_Lag625]=get_param_Lag_sf3fk_bootstrap(Case,625,lambda,timerange);
timerange=1:720;

[SpecFlux_Lag2500,Vt_Lag2500,ebs_Lag2500,kf_Lag,...
    lf_Lag,CI_ebs2500,CI_Vt2500,...
    CI_SpecFlux2500,Th_Lag2500]=get_param_Lag_sf3fk_bootstrap(Case,2500,lambda,timerange);
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
y_coords = [-2e-8, -2e-8, 4e-8, 4e-8]; % y坐标：下边界-2e-8，上边界4e-8
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
b1=semilogx(filtscale./1e3,Thm_Eulerian,'LineWidth',1.5,'Color',colors{1});
hold on
r_cg_eul = findXatYZero(filtscale./1e3, Thm_Eulerian);

b2=semilogx(1./kf_E(dstr:end)./1e3.*2.*pi,SpecFlux_E(dstr:end), ...
    'LineWidth',1.5, ...
'Color','k');
r_RLS_eul = findXatYZero(1./kf_E(dstr:end)./1e3.*2.*pi,SpecFlux_E(dstr:end));

B1=semilogx(filtscale./1e3,Th_Lag289,'LineWidth',1.5,'Color',colors{4}, ...
    'LineStyle','--')
% hold on
r_cg_Lag = findXatYZero(filtscale./1e3,Th_Lag289);

B2=semilogx(1./kf_Lag(1:dend).*2.*pi./1e3,SpecFlux_Lag289(1:dend),'LineWidth',1.5,'Color',colors{4})
r_RLS_Lag = findXatYZero(1./kf_Lag(1:dend).*2.*pi./1e3,SpecFlux_Lag289(1:dend));

fill(x_fill, [CI_SpecFlux289(1,1:dend), fliplr(CI_SpecFlux289(2,1:dend))], ...
    colors_rgb{4}, 'FaceAlpha', 0.1, 'EdgeColor', 'none');


grid on
ylim([-2e-8,4e-8]);
xlim([1e3,1e6]./1e3);
xtick([1,4,1e1,1e2,200,1e3])
% xticklabels({'10^{0}','4','10^{1}','10^{2}','200','10^{3}'})

xlabel('$$\mathbf{r \ [km]}$$','Interpreter','latex')
ylabel('$$\mathbf{F(r) \ [m^{2}/s^{3}]}$$','Interpreter','latex')
% % title('Hf: CG vs sf3-RLS flux')
% legend([b1,b2],{['CG (',win,')'],'SF3-fitting (RLS)'})
legend([b1,b2,B1,B2],{['Eul CG (',win,')'],'SF3-fitting (RLS)',['Lag CG (',win,')'],'Lag SF3-fitting (RLS)'})
text(0.03, 0.95, 'a) Hf P289', 'Units', 'normalized', ...
     'FontSize', 12, 'FontWeight', 'bold', ...
     'BackgroundColor', [1, 1, 0.8, 0.6], ... % 半透明背景
     'EdgeColor', [.7,.7,.7], ... % 边框颜色
     'Margin', 3, ... % 边距
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')
text(10^0.1,-0.6e-8,['$\mathbf{\ r_{CG \ Eul}= }$',num2str(r_cg_eul)],'Interpreter','latex', ...
    'Color',colors{1})
text(10^0.1,-0.9e-8,['$\mathbf{\ r_{RLS \ Eul}= }$',num2str(r_RLS_eul)],'Interpreter','latex')
text(10^0.1,-1.2e-8,['$\mathbf{\ r_{CG \ Lag}= }$',num2str(r_cg_Lag)],'Interpreter','latex', ...
    'Color',colors{4})
text(10^0.1,-1.5e-8,['$\mathbf{\ r_{RLS \ Lag}= }$',num2str(r_RLS_Lag)],'Interpreter','latex', ...
    'Color',colors{4})

set(gca,'fontsize',12,'fontweight','bold')

subplot(2,3,2)
patch(x_coords, y_coords, [0.9 0.9 0.9],'FaceAlpha', 0.3, ...
    'EdgeColor', 'k'); 
set(gca,'xscale','log')
hold on

b1=semilogx(filtscale./1e3,Thm_Eulerian,'LineWidth',1.5,'Color',colors{1});
hold on
r_cg_eul = findXatYZero(filtscale./1e3, Thm_Eulerian);

b2=semilogx(1./kf_E(dstr:end)./1e3.*2.*pi,SpecFlux_E(dstr:end), ...
    'LineWidth',1.5, ...
'Color','k');
r_RLS_eul = findXatYZero(1./kf_E(dstr:end)./1e3.*2.*pi,SpecFlux_E(dstr:end));

B1=semilogx(filtscale./1e3,Th_Lag625,'LineWidth',1.5,'Color',colors{5}, ...
    'LineStyle','--')
r_cg_Lag = findXatYZero(filtscale./1e3,Th_Lag625);

B2=semilogx(1./kf_Lag(1:dend).*2.*pi./1e3,SpecFlux_Lag625(1:dend),'LineWidth',1.5,'Color',colors{5})
r_RLS_Lag = findXatYZero(1./kf_Lag(1:dend).*2.*pi./1e3,SpecFlux_Lag625(1:dend));

fill(x_fill, [CI_SpecFlux625(1,1:dend), fliplr(CI_SpecFlux625(2,1:dend))], ...
    colors_rgb{5}, 'FaceAlpha', 0.1, 'EdgeColor', 'none');

grid on
ylim([-2e-8,4e-8]);
xlim([1e3,1e6]./1e3);
xtick([1,4,1e1,1e2,200,1e3])

xlabel('$$\mathbf{r \ [km]}$$','Interpreter','latex')
ylabel('$$\mathbf{F(r) \ [m^{2}/s^{3}]}$$','Interpreter','latex')
% % title('Hf: CG vs sf3-RLS flux')
% legend([b1,b2],{['CG (',win,')'],'SF3-fitting (RLS)'})
legend([b1,b2,B1,B2],{['Eul CG (',win,')'],'SF3-fitting (RLS)',['Lag CG (',win,')'],'Lag SF3-fitting (RLS)'})
text(0.03, 0.95, 'b) Hf P625', 'Units', 'normalized', ...
     'FontSize', 12, 'FontWeight', 'bold', ...
     'BackgroundColor', [1, 1, 0.8, 0.6], ... % 半透明背景
     'EdgeColor', [.7,.7,.7], ... % 边框颜色
     'Margin', 3, ... % 边距
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')
text(10^0.1,-0.6e-8,['$\mathbf{\ r_{CG \ Eul}= }$',num2str(r_cg_eul)],'Interpreter','latex', ...
    'Color',colors{1})
text(10^0.1,-0.9e-8,['$\mathbf{\ r_{RLS \ Eul}= }$',num2str(r_RLS_eul)],'Interpreter','latex')
text(10^0.1,-1.2e-8,['$\mathbf{\ r_{CG \ Lag}= }$',num2str(r_cg_Lag)],'Interpreter','latex', ...
    'Color',colors{5})
text(10^0.1,-1.5e-8,['$\mathbf{\ r_{RLS \ Lag}= }$',num2str(r_RLS_Lag)],'Interpreter','latex', ...
    'Color',colors{5})
set(gca,'fontsize',12,'fontweight','bold')

subplot(2,3,3)
patch(x_coords, y_coords, [0.9 0.9 0.9],'FaceAlpha', 0.3, ...
    'EdgeColor', 'k'); 
set(gca,'xscale','log')
hold on

b1=semilogx(filtscale./1e3,Thm_Eulerian,'LineWidth',1.5,'Color',colors{1});
hold on
r_cg_eul = findXatYZero(filtscale./1e3, Thm_Eulerian);
b2=semilogx(1./kf_E(dstr:end)./1e3.*2.*pi,SpecFlux_E(dstr:end), ...
    'LineWidth',1.5, ...
'Color','k');
r_RLS_eul = findXatYZero(1./kf_E(dstr:end)./1e3.*2.*pi,SpecFlux_E(dstr:end));

B1=semilogx(filtscale./1e3,Th_Lag2500,'LineWidth',1.5,'Color',colors{6}, ...
    'LineStyle','--')
r_cg_Lag = findXatYZero(filtscale./1e3,Th_Lag2500);
B2=semilogx(1./kf_Lag(1:dend).*2.*pi./1e3,SpecFlux_Lag2500(1:dend),'LineWidth',1.5,'Color',colors{6})
r_RLS_Lag = findXatYZero(1./kf_Lag(1:dend).*2.*pi./1e3,SpecFlux_Lag2500(1:dend));

fill(x_fill, [CI_SpecFlux2500(1,1:dend), fliplr(CI_SpecFlux2500(2,1:dend))], ...
    colors_rgb{6}, 'FaceAlpha', 0.1, 'EdgeColor', 'none');

grid on
ylim([-2e-8,4e-8]);
xlim([1e3,1e6]./1e3);
xtick([1,4,1e1,1e2,200,1e3])
xlabel('$$\mathbf{r \ [km]}$$','Interpreter','latex')
ylabel('$$\mathbf{F(r) \ [m^{2}/s^{3}]}$$','Interpreter','latex')
% % title('Hf: CG vs sf3-RLS flux')
% legend([b1,b2],{['CG (',win,')'],'SF3-fitting (RLS)'})
legend([b1,b2,B1,B2],{['Eul CG (',win,')'],'SF3-fitting (RLS)',['Lag CG (',win,')'],'Lag SF3-fitting (RLS)'})
text(0.03, 0.95, 'c) Hf P2500', 'Units', 'normalized', ...
     'FontSize', 12, 'FontWeight', 'bold', ...
     'BackgroundColor', [1, 1, 0.8, 0.6], ... % 半透明背景
     'EdgeColor', [.7,.7,.7], ... % 边框颜色
     'Margin', 3, ... % 边距
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')
text(10^0.1,-0.6e-8,['$\mathbf{\ r_{CG \ Eul}= }$',num2str(r_cg_eul)],'Interpreter','latex', ...
    'Color',colors{1})
text(10^0.1,-0.9e-8,['$\mathbf{\ r_{RLS \ Eul}= }$',num2str(r_RLS_eul)],'Interpreter','latex')
text(10^0.1,-1.2e-8,['$\mathbf{\ r_{CG \ Lag}= }$',num2str(r_cg_Lag)],'Interpreter','latex', ...
    'Color',colors{6})
text(10^0.1,-1.5e-8,['$\mathbf{\ r_{RLS \ Lag}= }$',num2str(r_RLS_Lag)],'Interpreter','latex', ...
    'Color',colors{6})

set(gca,'fontsize',12,'fontweight','bold')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% nowave 
Case='nowave'
timerange=1:2148;
param_Eul_cg_sf3fk
lambda=[10,1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6, ...
    1e-7,1e-8,1e-9,1e-10,1e-11,1e-12];
timerange=1:720;
[SpecFlux_Lag289,Vt_Lag289,ebs_Lag289,~,...
    ~,CI_ebs289,CI_Vt289,...
    CI_SpecFlux289,Th_Lag289]=get_param_Lag_sf3fk_bootstrap(Case,289,lambda,timerange);
timerange=1:720;
[SpecFlux_Lag625,Vt_Lag625,ebs_Lag625,~,...
    ~,CI_ebs625,CI_Vt625,...
    CI_SpecFlux625,Th_Lag625]=get_param_Lag_sf3fk_bootstrap(Case,625,lambda,timerange);
timerange=1:720;
[SpecFlux_Lag2500,Vt_Lag2500,ebs_Lag2500,kf_Lag,...
    lf_Lag,CI_ebs2500,CI_Vt2500,...
    CI_SpecFlux2500,Th_Lag2500]=get_param_Lag_sf3fk_bootstrap(Case,2500,lambda,timerange);
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


subplot(2,3,4)
patch(x_coords, y_coords, [0.9 0.9 0.9],'FaceAlpha', 0.3, ...
    'EdgeColor', 'k'); 
set(gca,'xscale','log')
hold on
b1=semilogx(filtscale./1e3,Thm_Eulerian,'LineWidth',1.5,'Color',colors{2});
hold on
r_cg_eul = findXatYZero(filtscale./1e3, Thm_Eulerian);

b2=semilogx(1./kf_E(dstr:end)./1e3.*2.*pi,SpecFlux_E(dstr:end), ...
    'LineWidth',1.5, ...
'Color','k');
r_RLS_eul = findXatYZero(1./kf_E(dstr:end)./1e3.*2.*pi,SpecFlux_E(dstr:end));

B1=semilogx(filtscale./1e3,Th_Lag289,'LineWidth',1.5,'Color',colors{4}, ...
    'LineStyle','--')
r_cg_Lag = findXatYZero(filtscale./1e3,Th_Lag289);

B2=semilogx(1./kf_Lag(1:dend).*2.*pi./1e3,SpecFlux_Lag289(1:dend),'LineWidth',1.5,'Color',colors{4})
r_RLS_Lag = findXatYZero(1./kf_Lag(1:dend).*2.*pi./1e3,SpecFlux_Lag289(1:dend));

fill(x_fill, [CI_SpecFlux289(1,1:dend), fliplr(CI_SpecFlux289(2,1:dend))], ...
    colors_rgb{4}, 'FaceAlpha', 0.1, 'EdgeColor', 'none');

grid on
ylim([-2e-8,4e-8]);
xlim([1e3,1e6]./1e3);
xtick([1,4,1e1,1e2,200,1e3])
xlabel('$$\mathbf{r \ [km]}$$','Interpreter','latex')
ylabel('$$\mathbf{F(r) \ [m^{2}/s^{3}]}$$','Interpreter','latex')
% % title('Hf: CG vs sf3-RLS flux')
% legend([b1,b2],{['CG (',win,')'],'SF3-fitting (RLS)'})
legend([b1,b2,B1,B2],{['Eul CG (',win,')'],'SF3-fitting (RLS)',['Lag CG (',win,')'],'Lag SF3-fitting (RLS)'})
text(0.03, 0.95, 'd) Sm P289', 'Units', 'normalized', ...
     'FontSize', 12, 'FontWeight', 'bold', ...
     'BackgroundColor', [1, 1, 0.8, 0.6], ... % 半透明背景
     'EdgeColor', [.7,.7,.7], ... % 边框颜色
     'Margin', 3, ... % 边距
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')
text(10^0.1,-0.6e-8,['$\mathbf{\ r_{CG \ Eul}= }$',num2str(r_cg_eul)],'Interpreter','latex', ...
    'Color',colors{2})
text(10^0.1,-0.9e-8,['$\mathbf{\ r_{RLS \ Eul}= }$',num2str(r_RLS_eul)],'Interpreter','latex')
text(10^0.1,-1.2e-8,['$\mathbf{\ r_{CG \ Lag}= }$',num2str(r_cg_Lag)],'Interpreter','latex', ...
    'Color',colors{4})
text(10^0.1,-1.5e-8,['$\mathbf{\ r_{RLS \ Lag}= }$',num2str(r_RLS_Lag)],'Interpreter','latex', ...
    'Color',colors{4})
set(gca,'fontsize',12,'fontweight','bold')

subplot(2,3,5)
patch(x_coords, y_coords, [0.9 0.9 0.9],'FaceAlpha', 0.3, ...
    'EdgeColor', 'k'); 
set(gca,'xscale','log')
hold on
b1=semilogx(filtscale./1e3,Thm_Eulerian,'LineWidth',1.5,'Color',colors{2});
r_cg_eul = findXatYZero(filtscale./1e3, Thm_Eulerian);

b2=semilogx(1./kf_E(dstr:end)./1e3.*2.*pi,SpecFlux_E(dstr:end), ...
    'LineWidth',1.5, 'Color','k');
r_RLS_eul = findXatYZero(1./kf_E(dstr:end)./1e3.*2.*pi,SpecFlux_E(dstr:end));

B1=semilogx(filtscale./1e3,Th_Lag625,'LineWidth',1.5,'Color',colors{5}, ...
    'LineStyle','--')
r_cg_Lag = findXatYZero(filtscale./1e3,Th_Lag625);
B2=semilogx(1./kf_Lag(1:dend).*2.*pi./1e3,SpecFlux_Lag625(1:dend),'LineWidth',1.5,'Color',colors{5})
r_RLS_Lag = findXatYZero(1./kf_Lag(1:dend).*2.*pi./1e3,SpecFlux_Lag625(1:dend));

fill(x_fill, [CI_SpecFlux625(1,1:dend), fliplr(CI_SpecFlux625(2,1:dend))], ...
    colors_rgb{5}, 'FaceAlpha', 0.1, 'EdgeColor', 'none');

grid on
ylim([-2e-8,4e-8]);
xlim([1e3,1e6]./1e3);
xtick([1,4,1e1,1e2,200,1e3])
xlabel('$$\mathbf{r \ [km]}$$','Interpreter','latex')
ylabel('$$\mathbf{F(r) \ [m^{2}/s^{3}]}$$','Interpreter','latex')

legend([b1,b2,B1,B2],{['Eul CG (',win,')'],'SF3-fitting (RLS)',['Lag CG (',win,')'],'Lag SF3-fitting (RLS)'})
text(0.03, 0.95, 'e) Sm P625', 'Units', 'normalized', ...
     'FontSize', 12, 'FontWeight', 'bold', ...
     'BackgroundColor', [1, 1, 0.8, 0.6], ... % 半透明背景
     'EdgeColor', [.7,.7,.7], ... % 边框颜色
     'Margin', 3, ... % 边距
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')
text(10^0.1,-0.6e-8,['$\mathbf{\ r_{CG \ Eul}= }$',num2str(r_cg_eul)],'Interpreter','latex', ...
    'Color',colors{2})
text(10^0.1,-0.9e-8,['$\mathbf{\ r_{RLS \ Eul}= }$',num2str(r_RLS_eul)],'Interpreter','latex')
text(10^0.1,-1.2e-8,['$\mathbf{\ r_{CG \ Lag}= }$',num2str(r_cg_Lag)],'Interpreter','latex', ...
    'Color',colors{5})
text(10^0.1,-1.5e-8,['$\mathbf{\ r_{RLS \ Lag}= }$',num2str(r_RLS_Lag)],'Interpreter','latex', ...
    'Color',colors{5})
set(gca,'fontsize',12,'fontweight','bold')

subplot(2,3,6)
patch(x_coords, y_coords, [0.9 0.9 0.9],'FaceAlpha', 0.3, ...
    'EdgeColor', 'k'); 
set(gca,'xscale','log')
hold on
b1=semilogx(filtscale./1e3,Thm_Eulerian,'LineWidth',1.5,'Color',colors{2});
r_cg_eul = findXatYZero(filtscale./1e3, Thm_Eulerian);

b2=semilogx(1./kf_E(dstr:end)./1e3.*2.*pi,SpecFlux_E(dstr:end), ...
    'LineWidth',1.5, 'Color','k');
r_RLS_eul = findXatYZero(1./kf_E(dstr:end)./1e3.*2.*pi,SpecFlux_E(dstr:end));

B1=semilogx(filtscale./1e3,Th_Lag2500,'LineWidth',1.5,'Color',colors{6}, ...
    'LineStyle','--')
r_cg_Lag = findXatYZero(filtscale./1e3,Th_Lag2500);
B2=semilogx(1./kf_Lag(1:dend).*2.*pi./1e3,SpecFlux_Lag2500(1:dend),'LineWidth',1.5,'Color',colors{6})
r_RLS_Lag = findXatYZero(1./kf_Lag(1:dend).*2.*pi./1e3,SpecFlux_Lag2500(1:dend));

fill(x_fill, [CI_SpecFlux2500(1,1:dend), fliplr(CI_SpecFlux2500(2,1:dend))], ...
    colors_rgb{6}, 'FaceAlpha', 0.1, 'EdgeColor', 'none');

grid on
ylim([-2e-8,4e-8]);
xlim([1e3,1e6]./1e3);
xtick([1,4,1e1,1e2,200,1e3])
xlabel('$$\mathbf{r \ [km]}$$','Interpreter','latex')
ylabel('$$\mathbf{F(r) \ [m^{2}/s^{3}]}$$','Interpreter','latex')

legend([b1,b2,B1,B2],{['Eul CG (',win,')'],'SF3-fitting (RLS)',['Lag CG (',win,')'],'Lag SF3-fitting (RLS)'})
text(0.03, 0.95, 'f) Sm P2500', 'Units', 'normalized', ...
     'FontSize', 12, 'FontWeight', 'bold', ...
     'BackgroundColor', [1, 1, 0.8, 0.6], ... % 半透明背景
     'EdgeColor', [.7,.7,.7], ... % 边框颜色
     'Margin', 3, ... % 边距
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')
text(10^0.1,-0.6e-8,['$\mathbf{\ r_{CG \ Eul}= }$',num2str(r_cg_eul)],'Interpreter','latex', ...
    'Color',colors{2})
text(10^0.1,-0.9e-8,['$\mathbf{\ r_{RLS \ Eul}= }$',num2str(r_RLS_eul)],'Interpreter','latex')
text(10^0.1,-1.2e-8,['$\mathbf{\ r_{CG \ Lag}= }$',num2str(r_cg_Lag)],'Interpreter','latex', ...
    'Color',colors{6})
text(10^0.1,-1.5e-8,['$\mathbf{\ r_{RLS \ Lag}= }$',num2str(r_RLS_Lag)],'Interpreter','latex', ...
    'Color',colors{6})
set(gca,'fontsize',12,'fontweight','bold')
saveas(gcf,['Eul_CG_vs_RLS_',win],'png')
