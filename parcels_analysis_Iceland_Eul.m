
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
win='kaiser'
param_Eul_cg_sf3fk

% figure(1)
screenSize = get(0, 'ScreenSize');
% 创建figure并设置位置
figure('Position', [0, 0, screenSize(3), screenSize(4)]);


subplot(2,3,1)
b1=semilogx(filtscale./1e3,Thm_Eulerian,'LineWidth',1.5,'Color',colors{1});
hold on
r_cg = findXatYZero(filtscale./1e3, Thm_Eulerian);

b2=semilogx(1./kf_E(dstr:end)./1e3.*2.*pi,SpecFlux_E(dstr:end), ...
    'LineWidth',1.5, ...
'Color','k');
r_fit = findXatYZero(1./kf_E(dstr:end)./1e3.*2.*pi, SpecFlux_E(dstr:end));


grid on
ylim([-2e-8,2e-8]);
xlim([1e3,1e6]./1e3);
xlabel('$$\mathbf{[km]}$$','Interpreter','latex')
ylabel('$$\mathbf{F(r) \ [m^{2}/s^{3}]}$$','Interpreter','latex')
% title('Hf: CG vs sf3-RLS flux')
legend([b1,b2],{['CG (',win,')'],'SF3-fitting (RLS)'})
text(0.05, 0.95, 'a) Hf', 'Units', 'normalized', ...
     'FontSize', 12, 'FontWeight', 'bold', ...
     'BackgroundColor', [1, 1, 0.8, 0.6], ... % 半透明背景
     'EdgeColor', [.7,.7,.7], ... % 边框颜色
     'Margin', 3, ... % 边距
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')
text(10^0.1,-1.2e-8,['$\mathbf{CG \ r_{F=0}= }$',num2str(r_cg)],'Interpreter','latex', ...
    'Color',colors{1})
text(10^0.1,-1.5e-8,['$\mathbf{Fitting \ r_{F=0}= }$',num2str(r_fit)],'Interpreter','latex')

set(gca,'fontsize',12,'fontweight','bold')


subplot(2,3,2)
a1=semilogx(r(2:end)./1e3,SF3(2:end)./r(2:end),'LineWidth',1.5,'Color',colors{1});
hold on
a2=semilogx(lf_E./1e3,Vt_E./lf_E','LineWidth',1.5,'LineStyle','none', ...
    'Marker','+','Color','k');

grid on
xlim([1e3,1e6]./1e3);
ylim([-3,4].*1e-8)
xlabel('$$\mathbf{[km]}$$','Interpreter','latex')
ylabel('$$\mathbf{D3(r)/r \ [m^{2}/s^{3}]}$$','Interpreter','latex')
% title('Hf: SF3/r and fitting (RLS)')
% legend([a1,a2],{'SF3/r','fitting (RLS)'})
legend([a1,a2],{['CG (',win,')'],'SF3-fitting (RLS)'})
text(0.05, 0.95, 'b) Hf', 'Units', 'normalized', ...
     'FontSize', 12, 'FontWeight', 'bold', ...
     'BackgroundColor', [1, 1, 0.8, 0.6], ... % 半透明背景
     'EdgeColor', [.7,.7,.7], ... % 边框颜色
     'Margin', 3, ... % 边距
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')

set(gca,'fontsize',12,'fontweight','bold')



subplot(2,3,3)
c1=semilogx(kc_mid_s.*1e3,dPidk_s,'LineWidth',1.5,'Color',colors{1});
hold on
c2=semilogx(kf_E.*1e3./2./pi,ebs_E(1:end-1).*dk_E','LineWidth',1.5, ...
    'Color','k')
grid on
ylim([-1,1].*1e-8)
xlim([1e-6,1e-3].*1e3)
xlabel('$$\mathbf{[1/km]}$$','Interpreter','latex')
ylabel('$$\mathbf{\epsilon_j*dk_j \ [m^{2}/s^{3}]}$$','Interpreter','latex')
legend([c1,c2],{['CG (',win,')'],'SF3-fitting (RLS)'})
text(0.05, 0.95, 'c) Hf', 'Units', 'normalized', ...
     'FontSize', 12, 'FontWeight', 'bold', ...
     'BackgroundColor', [1, 1, 0.8, 0.6], ... % 半透明背景
     'EdgeColor', [.7,.7,.7], ... % 边框颜色
     'Margin', 3, ... % 边距
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')

% title('Hf: $$dk_j.\epsilon_j$$','Interpreter','latex')
set(gca,'fontsize',12,'fontweight','bold')


%%%%%%%%%%%%%%%%%%nowave%%%%%%%%%%%%%%%%%
Case='nowave'; % wave
param_Eul_cg_sf3fk

subplot(2,3,4)
b1=semilogx(filtscale./1e3,Thm_Eulerian,'LineWidth',1.5,'Color',colors{2});
hold on
r_cg = findXatYZero(filtscale./1e3, Thm_Eulerian);

b2=semilogx(1./kf_E(dstr:end)./1e3.*2.*pi,SpecFlux_E(dstr:end), ...
    'LineWidth',1.5, ...
'Color','k');
r_fit = findXatYZero(1./kf_E(dstr:end)./1e3.*2.*pi,SpecFlux_E(dstr:end));
if strcmpi(win, 'hann')
    r_cg=r_fit;
end
grid on
ylim([-2e-8,2e-8]);
xlim([1e3,1e6]./1e3);
xlabel('$$\mathbf{[km]}$$','Interpreter','latex')
ylabel('$$\mathbf{F(k) \ [m^{2}/s^{3}]}$$','Interpreter','latex')
% title('Hf: CG vs sf3-RLS flux')
legend([b1,b2],{['CG (',win,')'],'SF3-fitting (RLS)'})
text(0.05, 0.95, 'd) Smooth', 'Units', 'normalized', ...
     'FontSize', 12, 'FontWeight', 'bold', ...
     'BackgroundColor', [1, 1, 0.8, 0.6], ... % 半透明背景
     'EdgeColor', [.7,.7,.7], ... % 边框颜色
     'Margin', 3, ... % 边距
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')
text(10^0.1,-1.2e-8,['$\mathbf{CG \ r_{F=0}= }$',num2str(r_cg)],'Interpreter','latex', ...
    'Color',colors{2})
text(10^0.1,-1.5e-8,['$\mathbf{Fitting \ r_{F=0}= }$',num2str(r_fit)],'Interpreter','latex')

set(gca,'fontsize',12,'fontweight','bold')


subplot(2,3,5)
a1=semilogx(r(2:end)./1e3,SF3(2:end)./r(2:end),'LineWidth',1.5,'Color',colors{2});
hold on
a2=semilogx(lf_E./1e3,Vt_E./lf_E','LineWidth',1.5,'LineStyle','none', ...
    'Marker','+','Color','k');

grid on
xlim([1e3,1e6]./1e3);
ylim([-3,4].*1e-8)
xlabel('$$\mathbf{[km]}$$','Interpreter','latex')
ylabel('$$\mathbf{D3(r)/r \ [m^{2}/s^{3}]}$$','Interpreter','latex')
% title('Hf: SF3/r and fitting (RLS)')
% legend([a1,a2],{'SF3/r','fitting (RLS)'})
legend([a1,a2],{['CG (',win,')'],'SF3-fitting (RLS)'})
text(0.05, 0.95, 'e) Smooth', 'Units', 'normalized', ...
     'FontSize', 12, 'FontWeight', 'bold', ...
     'BackgroundColor', [1, 1, 0.8, 0.6], ... % 半透明背景
     'EdgeColor', [.7,.7,.7], ... % 边框颜色
     'Margin', 3, ... % 边距
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')

set(gca,'fontsize',12,'fontweight','bold')



subplot(2,3,6)
c1=semilogx(kc_mid_s.*1e3,dPidk_s,'LineWidth',1.5,'Color',colors{2});
hold on
c2=semilogx(kf_E.*1e3./2./pi,ebs_E(1:end-1).*dk_E','LineWidth',1.5, ...
    'Color','k')
grid on
ylim([-1,1].*1e-8)
xlim([1e-6,1e-3].*1e3)
xlabel('$$\mathbf{[1/km]}$$','Interpreter','latex')
ylabel('$$\mathbf{\epsilon_j*dk_j \ [m^{2}/s^{3}]}$$','Interpreter','latex')
legend([c1,c2],{['CG (',win,')'],'SF3-fitting (RLS)'})
text(0.05, 0.95, 'f) Smooth', 'Units', 'normalized', ...
     'FontSize', 12, 'FontWeight', 'bold', ...
     'BackgroundColor', [1, 1, 0.8, 0.6], ... % 半透明背景
     'EdgeColor', [.7,.7,.7], ... % 边框颜色
     'Margin', 3, ... % 边距
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')
% title('Hf: $$dk_j.\epsilon_j$$','Interpreter','latex')
set(gca,'fontsize',12,'fontweight','bold')

% close(figure(1))
saveas(gcf,['Eul_CG_vs_RLS_',win],'png')
