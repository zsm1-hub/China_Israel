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
timerange=1:1920;
[SpecFlux_Lag289,Vt_Lag289,ebs_Lag289,~,...
    ~,CI_ebs289,CI_Vt289,...
    CI_SpecFlux289,Th_Lag289]=get_param_Lag_sf3fk_bootstrap(Case,289,lambda,timerange);
timerange=1:1920;
[SpecFlux_Lag625,Vt_Lag625,ebs_Lag625,~,...
    ~,CI_ebs625,CI_Vt625,...
    CI_SpecFlux625,Th_Lag625]=get_param_Lag_sf3fk_bootstrap(Case,625,lambda,timerange);
timerange=1:720;

[SpecFlux_Lag2500,Vt_Lag2500,ebs_Lag2500,kf_Lag,...
    lf_Lag,CI_ebs2500,CI_Vt2500,...
    CI_SpecFlux2500,Th_Lag2500]=get_param_Lag_sf3fk_bootstrap(Case,2500,lambda,timerange);
x_fill = [1./kf_Lag./1e3.*2.*pi, fliplr(1./kf_Lag./1e3.*2.*pi)];

if timerange(end)<500
    if strcmpi(Case, 'wave')
        load wave_CG_Eul_T480_rough.mat
        
    end
    if strcmpi(Case, 'nowave')
        load nowave_CG_Eul_T480_rough.mat
    end
    Thm_Eulerian=nanmean(Th_all480,2);
end


% b1=semilogx(filtscale./1e3,Thm_Eulerian,'LineWidth',1.5,'Color',colors{1});
% r_cg = findXatYZero(filtscale./1e3, Thm_Eulerian);

% b2=semilogx(1./kf_E(dstr:end)./1e3.*2.*pi,SpecFlux_E(dstr:end), ...
%     'LineWidth',1.5, 'Color','k');
% r_fit = findXatYZero(1./kf_E(dstr:end)./1e3.*2.*pi, SpecFlux_E(dstr:end));

subplot(2,3,1)
b1=semilogx(filtscale./1e3,Thm_Eulerian,'LineWidth',1.5,'Color',colors{1});
hold on
r_cg = findXatYZero(filtscale./1e3, Thm_Eulerian);

b2=semilogx(1./kf_E(dstr:end)./1e3.*2.*pi,SpecFlux_E(dstr:end), ...
    'LineWidth',1.5, ...
'Color','k');

B1=semilogx(filtscale./1e3,Th_Lag289,'LineWidth',1.5,'Color',colors{4}, ...
    'LineStyle','--')
% hold on

B2=semilogx(1./kf_Lag.*2.*pi./1e3,SpecFlux_Lag289,'LineWidth',1.5,'Color',colors{4})

fill(x_fill, [CI_SpecFlux289(1,:), fliplr(CI_SpecFlux289(2,:))], ...
    colors_rgb{4}, 'FaceAlpha', 0.1, 'EdgeColor', 'none');

grid on
ylim([-2e-8,4e-8]);
xlim([1e3,1e6]./1e3);
xlabel('$$\mathbf{r \ [km]}$$','Interpreter','latex')
ylabel('$$\mathbf{F(r) \ [m^{2}/s^{3}]}$$','Interpreter','latex')
% % title('Hf: CG vs sf3-RLS flux')
% legend([b1,b2],{['CG (',win,')'],'SF3-fitting (RLS)'})
legend([b1,b2,B1,B2],{['Eul CG (',win,')'],'SF3-fitting (RLS)',['Lag CG (',win,')'],'Lag SF3-fitting (RLS)'})
text(0.05, 0.95, 'a) Hf P289', 'Units', 'normalized', ...
     'FontSize', 12, 'FontWeight', 'bold', ...
     'BackgroundColor', [1, 1, 0.8, 0.6], ... % 半透明背景
     'EdgeColor', [.7,.7,.7], ... % 边框颜色
     'Margin', 3, ... % 边距
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')
% text(10^0.1,-1.2e-8,['$\mathbf{CG \ r_{F=0}= }$',num2str(r_cg)],'Interpreter','latex', ...
%     'Color',colors{1})
% text(10^0.1,-1.5e-8,['$\mathbf{Fitting \ r_{F=0}= }$',num2str(r_fit)],'Interpreter','latex')

set(gca,'fontsize',12,'fontweight','bold')

subplot(2,3,2)
b1=semilogx(filtscale./1e3,Thm_Eulerian,'LineWidth',1.5,'Color',colors{1});
hold on
r_cg = findXatYZero(filtscale./1e3, Thm_Eulerian);

b2=semilogx(1./kf_E(dstr:end)./1e3.*2.*pi,SpecFlux_E(dstr:end), ...
    'LineWidth',1.5, ...
'Color','k');
B1=semilogx(filtscale./1e3,Th_Lag625,'LineWidth',1.5,'Color',colors{5}, ...
    'LineStyle','--')

B2=semilogx(1./kf_Lag.*2.*pi./1e3,SpecFlux_Lag625,'LineWidth',1.5,'Color',colors{5})

fill(x_fill, [CI_SpecFlux625(1,:), fliplr(CI_SpecFlux625(2,:))], ...
    colors_rgb{5}, 'FaceAlpha', 0.1, 'EdgeColor', 'none');

grid on
ylim([-2e-8,4e-8]);
xlim([1e3,1e6]./1e3);
xlabel('$$\mathbf{r \ [km]}$$','Interpreter','latex')
ylabel('$$\mathbf{F(r) \ [m^{2}/s^{3}]}$$','Interpreter','latex')
% % title('Hf: CG vs sf3-RLS flux')
% legend([b1,b2],{['CG (',win,')'],'SF3-fitting (RLS)'})
legend([b1,b2,B1,B2],{['Eul CG (',win,')'],'SF3-fitting (RLS)',['Lag CG (',win,')'],'Lag SF3-fitting (RLS)'})
text(0.05, 0.95, 'a) Hf P625', 'Units', 'normalized', ...
     'FontSize', 12, 'FontWeight', 'bold', ...
     'BackgroundColor', [1, 1, 0.8, 0.6], ... % 半透明背景
     'EdgeColor', [.7,.7,.7], ... % 边框颜色
     'Margin', 3, ... % 边距
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')
% text(10^0.1,-1.2e-8,['$\mathbf{CG \ r_{F=0}= }$',num2str(r_cg)],'Interpreter','latex', ...
%     'Color',colors{1})
% text(10^0.1,-1.5e-8,['$\mathbf{Fitting \ r_{F=0}= }$',num2str(r_fit)],'Interpreter','latex')

set(gca,'fontsize',12,'fontweight','bold')

subplot(2,3,3)
b1=semilogx(filtscale./1e3,Thm_Eulerian,'LineWidth',1.5,'Color',colors{1});
hold on
r_cg = findXatYZero(filtscale./1e3, Thm_Eulerian);

b2=semilogx(1./kf_E(dstr:end)./1e3.*2.*pi,SpecFlux_E(dstr:end), ...
    'LineWidth',1.5, ...
'Color','k');

B1=semilogx(filtscale./1e3,Th_Lag2500,'LineWidth',1.5,'Color',colors{6}, ...
    'LineStyle','--')

B2=semilogx(1./kf_Lag.*2.*pi./1e3,SpecFlux_Lag2500,'LineWidth',1.5,'Color',colors{6})

fill(x_fill, [CI_SpecFlux2500(1,:), fliplr(CI_SpecFlux2500(2,:))], ...
    colors_rgb{6}, 'FaceAlpha', 0.1, 'EdgeColor', 'none');

grid on
ylim([-2e-8,4e-8]);
xlim([1e3,1e6]./1e3);
xlabel('$$\mathbf{r \ [km]}$$','Interpreter','latex')
ylabel('$$\mathbf{F(r) \ [m^{2}/s^{3}]}$$','Interpreter','latex')
% % title('Hf: CG vs sf3-RLS flux')
% legend([b1,b2],{['CG (',win,')'],'SF3-fitting (RLS)'})
legend([b1,b2,B1,B2],{['Eul CG (',win,')'],'SF3-fitting (RLS)',['Lag CG (',win,')'],'Lag SF3-fitting (RLS)'})
text(0.05, 0.95, 'a) Hf P2500', 'Units', 'normalized', ...
     'FontSize', 12, 'FontWeight', 'bold', ...
     'BackgroundColor', [1, 1, 0.8, 0.6], ... % 半透明背景
     'EdgeColor', [.7,.7,.7], ... % 边框颜色
     'Margin', 3, ... % 边距
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')
% text(10^0.1,-1.2e-8,['$\mathbf{CG \ r_{F=0}= }$',num2str(r_cg)],'Interpreter','latex', ...
%     'Color',colors{1})
% text(10^0.1,-1.5e-8,['$\mathbf{Fitting \ r_{F=0}= }$',num2str(r_fit)],'Interpreter','latex')

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
timerange=1:480;
[SpecFlux_Lag2500,Vt_Lag2500,ebs_Lag2500,kf_Lag,...
    lf_Lag,CI_ebs2500,CI_Vt2500,...
    CI_SpecFlux2500,Th_Lag2500]=get_param_Lag_sf3fk_bootstrap(Case,2500,lambda,timerange);
x_fill = [1./kf_Lag./1e3.*2.*pi, fliplr(1./kf_Lag./1e3.*2.*pi)];

if timerange(end)<500
    if strcmpi(Case, 'wave')
        load wave_CG_Eul_T480_rough.mat
        
    end
    if strcmpi(Case, 'nowave')
        load nowave_CG_Eul_T480_rough.mat
    end
    Thm_Eulerian=nanmean(Th_all480,2);
end


subplot(2,3,4)
b1=semilogx(filtscale./1e3,Thm_Eulerian,'LineWidth',1.5,'Color',colors{2});
hold on
r_cg = findXatYZero(filtscale./1e3, Thm_Eulerian);

b2=semilogx(1./kf_E(dstr:end)./1e3.*2.*pi,SpecFlux_E(dstr:end), ...
    'LineWidth',1.5, ...
'Color','k');

B1=semilogx(filtscale./1e3,Th_Lag289,'LineWidth',1.5,'Color',colors{4}, ...
    'LineStyle','--')
% hold on

B2=semilogx(1./kf_Lag.*2.*pi./1e3,SpecFlux_Lag289,'LineWidth',1.5,'Color',colors{4})

fill(x_fill, [CI_SpecFlux289(1,:), fliplr(CI_SpecFlux289(2,:))], ...
    colors_rgb{4}, 'FaceAlpha', 0.1, 'EdgeColor', 'none');

grid on
ylim([-2e-8,4e-8]);
xlim([1e3,1e6]./1e3);
xlabel('$$\mathbf{r \ [km]}$$','Interpreter','latex')
ylabel('$$\mathbf{F(r) \ [m^{2}/s^{3}]}$$','Interpreter','latex')
% % title('Hf: CG vs sf3-RLS flux')
% legend([b1,b2],{['CG (',win,')'],'SF3-fitting (RLS)'})
legend([b1,b2,B1,B2],{['Eul CG (',win,')'],'SF3-fitting (RLS)',['Lag CG (',win,')'],'Lag SF3-fitting (RLS)'})
text(0.05, 0.95, 'a) Smooth P289', 'Units', 'normalized', ...
     'FontSize', 12, 'FontWeight', 'bold', ...
     'BackgroundColor', [1, 1, 0.8, 0.6], ... % 半透明背景
     'EdgeColor', [.7,.7,.7], ... % 边框颜色
     'Margin', 3, ... % 边距
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')
% text(10^0.1,-1.2e-8,['$\mathbf{CG \ r_{F=0}= }$',num2str(r_cg)],'Interpreter','latex', ...
%     'Color',colors{1})
% text(10^0.1,-1.5e-8,['$\mathbf{Fitting \ r_{F=0}= }$',num2str(r_fit)],'Interpreter','latex')

set(gca,'fontsize',12,'fontweight','bold')

subplot(2,3,5)
b1=semilogx(filtscale./1e3,Thm_Eulerian,'LineWidth',1.5,'Color',colors{2});
hold on
r_cg = findXatYZero(filtscale./1e3, Thm_Eulerian);

b2=semilogx(1./kf_E(dstr:end)./1e3.*2.*pi,SpecFlux_E(dstr:end), ...
    'LineWidth',1.5, ...
'Color','k');
B1=semilogx(filtscale./1e3,Th_Lag625,'LineWidth',1.5,'Color',colors{5}, ...
    'LineStyle','--')

B2=semilogx(1./kf_Lag.*2.*pi./1e3,SpecFlux_Lag625,'LineWidth',1.5,'Color',colors{5})

fill(x_fill, [CI_SpecFlux625(1,:), fliplr(CI_SpecFlux625(2,:))], ...
    colors_rgb{5}, 'FaceAlpha', 0.1, 'EdgeColor', 'none');

grid on
ylim([-2e-8,4e-8]);
xlim([1e3,1e6]./1e3);
xlabel('$$\mathbf{r \ [km]}$$','Interpreter','latex')
ylabel('$$\mathbf{F(r) \ [m^{2}/s^{3}]}$$','Interpreter','latex')
% % title('Hf: CG vs sf3-RLS flux')
% legend([b1,b2],{['CG (',win,')'],'SF3-fitting (RLS)'})
legend([b1,b2,B1,B2],{['Eul CG (',win,')'],'SF3-fitting (RLS)',['Lag CG (',win,')'],'Lag SF3-fitting (RLS)'})
text(0.05, 0.95, 'a) Smooth P625', 'Units', 'normalized', ...
     'FontSize', 12, 'FontWeight', 'bold', ...
     'BackgroundColor', [1, 1, 0.8, 0.6], ... % 半透明背景
     'EdgeColor', [.7,.7,.7], ... % 边框颜色
     'Margin', 3, ... % 边距
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')
% text(10^0.1,-1.2e-8,['$\mathbf{CG \ r_{F=0}= }$',num2str(r_cg)],'Interpreter','latex', ...
%     'Color',colors{1})
% text(10^0.1,-1.5e-8,['$\mathbf{Fitting \ r_{F=0}= }$',num2str(r_fit)],'Interpreter','latex')

set(gca,'fontsize',12,'fontweight','bold')

subplot(2,3,6)
b1=semilogx(filtscale./1e3,Thm_Eulerian,'LineWidth',1.5,'Color',colors{2});
hold on
r_cg = findXatYZero(filtscale./1e3, Thm_Eulerian);

b2=semilogx(1./kf_E(dstr:end)./1e3.*2.*pi,SpecFlux_E(dstr:end), ...
    'LineWidth',1.5, ...
'Color','k');

B1=semilogx(filtscale./1e3,Th_Lag2500,'LineWidth',1.5,'Color',colors{6}, ...
    'LineStyle','--')

B2=semilogx(1./kf_Lag.*2.*pi./1e3,SpecFlux_Lag2500,'LineWidth',1.5,'Color',colors{6})

fill(x_fill, [CI_SpecFlux2500(1,:), fliplr(CI_SpecFlux2500(2,:))], ...
    colors_rgb{6}, 'FaceAlpha', 0.1, 'EdgeColor', 'none');

grid on
ylim([-2e-8,4e-8]);
xlim([1e3,1e6]./1e3);
xlabel('$$\mathbf{r \ [km]}$$','Interpreter','latex')
ylabel('$$\mathbf{F(r) \ [m^{2}/s^{3}]}$$','Interpreter','latex')
% % title('Hf: CG vs sf3-RLS flux')
% legend([b1,b2],{['CG (',win,')'],'SF3-fitting (RLS)'})
legend([b1,b2,B1,B2],{['Eul CG (',win,')'],'SF3-fitting (RLS)',['Lag CG (',win,')'],'Lag SF3-fitting (RLS)'})
text(0.05, 0.95, 'a) Smooth P2500', 'Units', 'normalized', ...
     'FontSize', 12, 'FontWeight', 'bold', ...
     'BackgroundColor', [1, 1, 0.8, 0.6], ... % 半透明背景
     'EdgeColor', [.7,.7,.7], ... % 边框颜色
     'Margin', 3, ... % 边距
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')
% text(10^0.1,-1.2e-8,['$\mathbf{CG \ r_{F=0}= }$',num2str(r_cg)],'Interpreter','latex', ...
%     'Color',colors{1})
% text(10^0.1,-1.5e-8,['$\mathbf{Fitting \ r_{F=0}= }$',num2str(r_fit)],'Interpreter','latex')

set(gca,'fontsize',12,'fontweight','bold')
