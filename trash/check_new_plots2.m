% [x_fill,y_fill]=get_shadow(x,y)
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
ini='_rough'
% npoint=9;

% ini='_rough'
npoint=7;
npoint2=12;


inv='RLS'
yy1=-0.6e-8;
yy2=-1.0e-8;
yy3=-1.4e-8;
yy4=-1.8e-8;

timerange=1:1940;
% timerange=1:720;

% %500 m
% param_Eul_cg_sf3fk_500m
% SF3_E=SF3_mean';
% r_E=r;
% ymin_fk=-3e-8;ymax_fk=10e-8;
% ymin_sf=-3e-7;ymax_sf=1e-7;
% % ymin_ebsdk=-3e-8;ymax_ebsdk=3e-8;
% ymin_ebsdk=-2e-4;ymax_ebsdk=3e-4;
% ymin_right = -4e-8;    % 
% ymax_right = 6e-8;  % 
% xmin=0.5;
% xmax=1e3;

% % 2km param
param_Eul_cg_sf3fk
% param_Eul_cg_sf3fk_no2pi

SF3_E=SF3(2:end)';
r_E=r(2:end)';
ymin_fk=-4e-8;ymax_fk=4e-8;
% ymin_fk=-3e-8;ymax_fk=6e-8;
ymin_sf=-1e-7;ymax_sf=1e-7;
% ymin_ebsdk=-3e-8;ymax_ebsdk=3e-8;
ymin_ebsdk=-2e-4;ymax_ebsdk=3e-4;
xmin=1;
xmax=1e3;

% 2km param rot
% param_Eul_cg_sf3fk_rot
% SF3_E=nanmean(SF3(2:end,:),2);
% r_E=r(2:end)';
% ymin_fk=-3e-8;ymax_fk=3e-8;
% ymin_sf=-1e-7;ymax_sf=1e-7;
% % ymin_ebsdk=-3e-8;ymax_ebsdk=3e-8;
% ymin_ebsdk=-2e-4;ymax_ebsdk=3e-4;
% ymin_right = -4e-8;    % 
% ymax_right = 6e-8;  % 
% xmin=1;
% xmax=1e3;



if strcmpi(Case, 'wave')
    ymin_right = -4e-8;    % 
    ymax_right = 6e-8;  % 
    Casemean='Hf';
    if strcmpi(ini,'_rough')
        lambda=[1e-7,1e-8,1e-9];
        % str1=15;
        % end1=28;
        str1=1;
        end1=19;
        range1=1e3;
        range2=300e3;
    end
    if strcmpi(ini,'_roughsmall') || strcmpi(ini,'_roughsmall_rot') || strcmpi(ini,'_roughLASER') || strcmpi(ini,'_roughsmall_2mon') 
        % lambda=[10,1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6, ...
        % 1e-7,1e-8,1e-9,1e-10,1e-11,2e-11,1e-12,1e-13,1e-14];
        lambda=[1e-7,1e-8,1e-9];
        str1=1;
        end1=19;
        range1=1e3;
        range2=500e3;
    end
    if  strcmpi(ini,'_roughLASER') || strcmpi(ini,'_roughsmall_2mon') 
        % lambda=[10,1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6, ...
        % 1e-7,1e-8,1e-9,1e-10,1e-11,2e-11,1e-12,1e-13,1e-14];
        lambda=[1e-7,1e-8,1e-9];
        str1=1;
        end1=19;
        range1=2.5e3;
        range2=300e3;
    end
    if  strcmpi(ini,'_roughsmall_1mon') 
        % lambda=[10,1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6, ...
        % 1e-7,1e-8,1e-9,1e-10,1e-11,2e-11,1e-12,1e-13,1e-14];
        lambda=[1e-7,1e-8,1e-9];
        str1=1;
        end1=19;
        range1=2.5e3;
        range2=300e3;
        timerange=1:720;
    end
    if strcmpi(ini,'_roughsmall_500m')
        lambda=[1e-8];
        str1=1;
        end1=24;
        range1=2e3;
        range2=300e3;
    end
    if strcmpi(ini,'_rough_500m') || strcmpi(ini,'_roughbox100g_500m')
        lambda=[1e-7];
        str1=1;
        end1=24;
        % end1=17;
        range1=1e3;
        range2=300e3;
    end
    load HF_SF3_check.mat
    % range2=300e3;
end
if strcmpi(Case, 'nowave')
    ymin_right = -1e-8;    % 
    ymax_right = 1.5e-8;  % 
    Casemean='Sm';
    % lambda=[10,1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6, ...
    % 1e-7,1e-8,1e-9,1e-10,1e-11,2e-11,1e-12,1e-13,1e-14,1e-15,1e-16];
    if strcmpi(ini,'_rough')
        % lambda=[10,1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6, ...
        % 1e-7,1e-8,1e-9,1e-10,1e-11,2e-11,1e-12];
        lambda=[1e-7,1e-8,1e-9];
        str1=1;
        end1=19;
        range1=2.5e3;
        range2=300e3;
    end
    if strcmpi(ini,'_roughsmall') || strcmpi(ini,'_roughsmall_rot') || strcmpi(ini,'_roughLASER') || strcmpi(ini,'_roughsmall_2mon') 
     %     lambda=[10,1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6, ...
     % 1e-7,1e-8,1e-9,1e-10,1e-11,2e-11,1e-12,1e-13,7e-14];
         lambda=[1e-7,1e-8,1e-9];
        str1=1;
        end1=18;
        range1=2.5e3;
        range2=300e3;
    end
    if strcmpi(ini,'_roughsmall_1mon') 
     %     lambda=[10,1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6, ...
     % 1e-7,1e-8,1e-9,1e-10,1e-11,2e-11,1e-12,1e-13,7e-14];
         lambda=[1e-7,1e-8,1e-9,1e-10,1e-11];
        str1=1;
        end1=18;
        range1=2.5e3;
        range2=300e3;
        timerange=1:720;
    end

    if strcmpi(ini,'_roughsmall_500m')
        lambda=[1e-8];
        str1=1;
        end1=24;
        range1=1e3;
        range2=300e3;
    end
    if strcmpi(ini,'_rough_500m') || strcmpi(ini,'_roughbox100g_500m')
        lambda=[5e-9];
        str1=1;
        % end1=24;
        end1=15;
        range1=1e3;
        range2=300e3;
    end
   
    load SM_SF3_check.mat

end


range1=2.5e3;
range2=400e3;

[dist_axis,SF3_Lag289,SF3_Lag289_std,optimal_fac289,optimal_po289,...
    Th_Lag289,Th_Lag289_std,bootstrap_results]=Fk_fitting_SF3_Bayesian_RLS_Lcurve_uncertainty4(Case,289,...
    timerange,ini,range1,range2,npoint,500,900);
% x_fill289_kf,y_fill289
kf_Lag=bootstrap_results.kf_final;
SpecFlux_Lag289=nanmean(bootstrap_results.Fk_all,2);
Fk_error289=std(bootstrap_results.Fk_all,0,2);
Vt_Lag289=nanmean(bootstrap_results.Vt_all,2);
eps_Lag289=nanmean(bootstrap_results.eps_all,2);
eps_error289=std(bootstrap_results.eps_all,0,2);

[dist_axis,SF3_Lag625,SF3_Lag625_std,optimal_fac625,optimal_po625,...
    Th_Lag625,Th_Lag625_std,bootstrap_results]=Fk_fitting_SF3_Bayesian_RLS_Lcurve_uncertainty4(Case,625,...
    timerange,ini,range1,range2,npoint,500,900);

SpecFlux_Lag625=nanmean(bootstrap_results.Fk_all,2);
Fk_error625=std(bootstrap_results.Fk_all,0,2);
Vt_Lag625=nanmean(bootstrap_results.Vt_all,2);
eps_Lag625=nanmean(bootstrap_results.eps_all,2);
eps_error625=std(bootstrap_results.eps_all,0,2);

[dist_axis,SF3_Lag1089,SF3_Lag1089_std,optimal_fac1089,optimal_po1089,...
    Th_Lag1089,Th_Lag1089_std,bootstrap_results]=Fk_fitting_SF3_Bayesian_RLS_Lcurve_uncertainty4(Case,1089,...
    timerange,ini,range1,range2,npoint,500,900);

SpecFlux_Lag1089=nanmean(bootstrap_results.Fk_all,2);
Fk_error1089=std(bootstrap_results.Fk_all,0,2);
Vt_Lag1089=nanmean(bootstrap_results.Vt_all,2);
eps_Lag1089=nanmean(bootstrap_results.eps_all,2);
eps_error1089=std(bootstrap_results.eps_all,0,2);
dist_axis=dist_axis';

if strcmpi(Case, 'nowave')
    colors_rgb{1}=colors_rgb{2};
    colors{1}=colors{2};
end
rescale=2*pi;
ds=1;
de=5;
% rescale=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%plot%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% patch(x_coords, y_coords, [0.9 0.9 0.9], 'EdgeColor', 'none'); % 中灰色，无边线
screenSize = get(0, 'ScreenSize');
% 创建figure并设置位置
figure('Position', [0, 0, screenSize(3), screenSize(4)]);

%%%%%%%%%%%%%%%%%%%%%%%%%debug%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,3,1)
b1=semilogx(dist_axis./1e3,...
    SF3_Lag289./dist_axis,'LineWidth',1.5);
hold on
[x_fill,y_fill]=get_shadow(dist_axis',(SF3_Lag289+SF3_Lag289_std)',...
    (SF3_Lag289-SF3_Lag289_std)');
fill(x_fill./1e3,y_fill./x_fill, ...
     'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
b2=semilogx(dist_axis./1e3,Vt_Lag289./dist_axis,'LineStyle','-','LineWidth',1.5);

% b3=semilogx(r_E(1:108)./1e3,...
%     SF3_E(1:108)./r_E(1:108),'LineWidth',1.5,'Color','k');

b3=semilogx(r_E./1e3,...
    SF3_E./r_E,'LineWidth',1.5,'Color','k');
% b4=semilogx(R_check./1e3,...
%     SF3_check./R_check,'LineWidth',1.5,'Color','m');
grid on
xlim([xmin,xmax]);
xticks([1,4,1e1,1e2,200,1e3])
% xticklabels({'10^{0}','4','10^{1}','10^{2}','200','10^{3}'})
ylim([ymin_sf,ymax_sf]);
xlim([xmin,xmax]);
xticks([1,4,1e1,1e2,200,1e3])
xlabel('$$\mathbf{r \ [km]}$$','Interpreter','latex')
ylabel('$$\mathbf{D3(r)/r \ [m^{2}/s^{3}]}$$','Interpreter','latex')

% legend([b1,b2,b3],{['Lag D3(r)/r'],'fitting',['Eul D3(r)/r']},...
%     'Location','southeast')

pos = [0.25, 0.88, 0.1, 0.04];  % [x, y, width, height]
% pos = [0.9, 0.88, 0.1, 0.04];  % [x, y, width, height]
legend([b1,b2,b3], {'Lag D3(r)/r', 'fitting', 'Eul D3(r)/r'}, ...
    'Position', pos);

text(0.03, 0.95, ['d) ',Casemean,' P289'], 'Units', 'normalized', ...
     'FontSize', 12, 'FontWeight', 'bold', ...
     'BackgroundColor', [1, 1, 0.8, 0.6], ... % 半透明背景
     'EdgeColor', [.7,.7,.7], ... % 边框颜色
     'Margin', 3, ... % 边距
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')
set(gca,'fontsize',12,'fontweight','bold')




% %%%%%%%%%%%%%%%%%%%%%%%%%debug%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,3,2)
semilogx(dist_axis./1e3,...
    SF3_Lag625./dist_axis,'LineWidth',1.5);
hold on
[x_fill,y_fill]=get_shadow(dist_axis',(SF3_Lag625+SF3_Lag625_std)',...
    (SF3_Lag625-SF3_Lag625_std)');

fill(x_fill./1e3,y_fill./x_fill, ...
     'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
semilogx(dist_axis./1e3,Vt_Lag625./dist_axis,'LineStyle','-','LineWidth',1.5);
% semilogx(r_E(1:108)./1e3,...
%     SF3_E(1:108)./r_E(1:108),'LineWidth',1.5,'Color','k');
semilogx(r_E./1e3,...
    SF3_E./r_E,'LineWidth',1.5,'Color','k');
% semilogx(R_check./1e3,...
%     SF3_check./R_check,'LineWidth',1.5,'Color','m');
grid on
xlim([xmin,xmax]);
xticks([1,4,1e1,1e2,200,1e3])
% xticklabels({'10^{0}','4','10^{1}','10^{2}','200','10^{3}'})
ylim([ymin_sf,ymax_sf]);
xlim([xmin,xmax]);
xticks([1,4,1e1,1e2,200,1e3])
xlabel('$$\mathbf{r \ [km]}$$','Interpreter','latex')
ylabel('$$\mathbf{D3(r)/r \ [m^{2}/s^{3}]}$$','Interpreter','latex')
text(0.03, 0.95, ['e) ',Casemean,' P625'], 'Units', 'normalized', ...
     'FontSize', 12, 'FontWeight', 'bold', ...
     'BackgroundColor', [1, 1, 0.8, 0.6], ... % 半透明背景
     'EdgeColor', [.7,.7,.7], ... % 边框颜色
     'Margin', 3, ... % 边距
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')
set(gca,'fontsize',12,'fontweight','bold')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,3,3)
semilogx(dist_axis./1e3,...
    SF3_Lag1089./dist_axis,'LineWidth',1.5);
hold on
[x_fill,y_fill]=get_shadow(dist_axis',(SF3_Lag1089+SF3_Lag1089_std)',...
    (SF3_Lag1089-SF3_Lag1089_std)');
fill(x_fill,y_fill, ...
     'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
semilogx(dist_axis./1e3,Vt_Lag1089./dist_axis,'LineStyle','-','LineWidth',1.5);
% semilogx(r_E(1:108)./1e3,...
%     SF3_E(1:108)./r_E(1:108),'LineWidth',1.5,'Color','k');
semilogx(r_E./1e3,...
    SF3_E./r_E,'LineWidth',1.5,'Color','k');
% semilogx(R_check./1e3,...
%     SF3_check./R_check,'LineWidth',1.5,'Color','m');
grid on
xlim([xmin,xmax]);
xticks([1,4,1e1,1e2,200,1e3])
% xticklabels({'10^{0}','4','10^{1}','10^{2}','200','10^{3}'})
ylim([ymin_sf,ymax_sf]);
xlim([xmin,xmax]);
xticks([1,4,1e1,1e2,200,1e3])
xlabel('$$\mathbf{r \ [km]}$$','Interpreter','latex')
ylabel('$$\mathbf{D3(r)/r \ [m^{2}/s^{3}]}$$','Interpreter','latex')

% legend([b1,b2],{['Lag D3(r)/r'],'fitting'})
text(0.03, 0.95, ['f) ',Casemean,' P1089'], 'Units', 'normalized', ...
     'FontSize', 12, 'FontWeight', 'bold', ...
     'BackgroundColor', [1, 1, 0.8, 0.6], ... % 半透明背景
     'EdgeColor', [.7,.7,.7], ... % 边框颜色
     'Margin', 3, ... % 边距
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')
set(gca,'fontsize',12,'fontweight','bold')

%%
ini='_roughsmall'
timerange=1:1940;
[dist_axis,SF3_Lag289,SF3_Lag289_std,optimal_fac289,optimal_po289,...
    Th_Lag289,Th_Lag289_std,bootstrap_results]=Fk_fitting_SF3_Bayesian_RLS_Lcurve_uncertainty4(Case,289,...
    timerange,ini,range1,range2,npoint,500,900);
% x_fill289_kf,y_fill289
kf_Lag=bootstrap_results.kf_final;
SpecFlux_Lag289=nanmean(bootstrap_results.Fk_all,2);
Fk_error289=std(bootstrap_results.Fk_all,0,2);
Vt_Lag289=nanmean(bootstrap_results.Vt_all,2);
eps_Lag289=nanmean(bootstrap_results.eps_all,2);
eps_error289=std(bootstrap_results.eps_all,0,2);

[dist_axis,SF3_Lag625,SF3_Lag625_std,optimal_fac625,optimal_po625,...
    Th_Lag625,Th_Lag625_std,bootstrap_results]=Fk_fitting_SF3_Bayesian_RLS_Lcurve_uncertainty4(Case,625,...
    timerange,ini,range1,range2,npoint,500,900);

SpecFlux_Lag625=nanmean(bootstrap_results.Fk_all,2);
Fk_error625=std(bootstrap_results.Fk_all,0,2);
Vt_Lag625=nanmean(bootstrap_results.Vt_all,2);
eps_Lag625=nanmean(bootstrap_results.eps_all,2);
eps_error625=std(bootstrap_results.eps_all,0,2);

[dist_axis,SF3_Lag1089,SF3_Lag1089_std,optimal_fac1089,optimal_po1089,...
    Th_Lag1089,Th_Lag1089_std,bootstrap_results]=Fk_fitting_SF3_Bayesian_RLS_Lcurve_uncertainty4(Case,1089,...
    timerange,ini,range1,range2,npoint,500,900);

SpecFlux_Lag1089=nanmean(bootstrap_results.Fk_all,2);
Fk_error1089=std(bootstrap_results.Fk_all,0,2);
Vt_Lag1089=nanmean(bootstrap_results.Vt_all,2);
eps_Lag1089=nanmean(bootstrap_results.eps_all,2);
eps_error1089=std(bootstrap_results.eps_all,0,2);
dist_axis=dist_axis';

subplot(2,3,4)
b1=semilogx(dist_axis./1e3,...
    SF3_Lag289./dist_axis,'LineWidth',1.5);
hold on
[x_fill,y_fill]=get_shadow(dist_axis',(SF3_Lag289+SF3_Lag289_std)',...
    (SF3_Lag289-SF3_Lag289_std)');
fill(x_fill./1e3,y_fill./x_fill, ...
     'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
b2=semilogx(dist_axis./1e3,Vt_Lag289./dist_axis,'LineStyle','-','LineWidth',1.5);

% b3=semilogx(r_E(1:108)./1e3,...
%     SF3_E(1:108)./r_E(1:108),'LineWidth',1.5,'Color','k');

b3=semilogx(r_E./1e3,...
    SF3_E./r_E,'LineWidth',1.5,'Color','k');
% b4=semilogx(R_check./1e3,...
%     SF3_check./R_check,'LineWidth',1.5,'Color','m');
grid on
xlim([xmin,xmax]);
xticks([1,4,1e1,1e2,200,1e3])
% xticklabels({'10^{0}','4','10^{1}','10^{2}','200','10^{3}'})
ylim([ymin_sf,ymax_sf]);
xlim([xmin,xmax]);
xticks([1,4,1e1,1e2,200,1e3])
xlabel('$$\mathbf{r \ [km]}$$','Interpreter','latex')
ylabel('$$\mathbf{D3(r)/r \ [m^{2}/s^{3}]}$$','Interpreter','latex')

% legend([b1,b2,b3],{['Lag D3(r)/r'],'fitting',['Eul D3(r)/r']},...
%     'Location','southeast')
text(0.03, 0.95, ['d) ',Casemean,' P289'], 'Units', 'normalized', ...
     'FontSize', 12, 'FontWeight', 'bold', ...
     'BackgroundColor', [1, 1, 0.8, 0.6], ... % 半透明背景
     'EdgeColor', [.7,.7,.7], ... % 边框颜色
     'Margin', 3, ... % 边距
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')
set(gca,'fontsize',12,'fontweight','bold')




% %%%%%%%%%%%%%%%%%%%%%%%%%debug%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,3,5)
semilogx(dist_axis./1e3,...
    SF3_Lag625./dist_axis,'LineWidth',1.5);
hold on
[x_fill,y_fill]=get_shadow(dist_axis',(SF3_Lag625+SF3_Lag625_std)',...
    (SF3_Lag625-SF3_Lag625_std)');

fill(x_fill./1e3,y_fill./x_fill, ...
     'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
semilogx(dist_axis./1e3,Vt_Lag625./dist_axis,'LineStyle','-','LineWidth',1.5);
% semilogx(r_E(1:108)./1e3,...
%     SF3_E(1:108)./r_E(1:108),'LineWidth',1.5,'Color','k');
semilogx(r_E./1e3,...
    SF3_E./r_E,'LineWidth',1.5,'Color','k');
% semilogx(R_check./1e3,...
%     SF3_check./R_check,'LineWidth',1.5,'Color','m');
grid on
xlim([xmin,xmax]);
xticks([1,4,1e1,1e2,200,1e3])
% xticklabels({'10^{0}','4','10^{1}','10^{2}','200','10^{3}'})
ylim([ymin_sf,ymax_sf]);
xlim([xmin,xmax]);
xticks([1,4,1e1,1e2,200,1e3])
xlabel('$$\mathbf{r \ [km]}$$','Interpreter','latex')
ylabel('$$\mathbf{D3(r)/r \ [m^{2}/s^{3}]}$$','Interpreter','latex')
text(0.03, 0.95, ['e) ',Casemean,' P625'], 'Units', 'normalized', ...
     'FontSize', 12, 'FontWeight', 'bold', ...
     'BackgroundColor', [1, 1, 0.8, 0.6], ... % 半透明背景
     'EdgeColor', [.7,.7,.7], ... % 边框颜色
     'Margin', 3, ... % 边距
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')
set(gca,'fontsize',12,'fontweight','bold')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,3,6)
semilogx(dist_axis./1e3,...
    SF3_Lag1089./dist_axis,'LineWidth',1.5);
hold on
[x_fill,y_fill]=get_shadow(dist_axis',(SF3_Lag1089+SF3_Lag1089_std)',...
    (SF3_Lag1089-SF3_Lag1089_std)');
fill(x_fill,y_fill, ...
     'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
semilogx(dist_axis./1e3,Vt_Lag1089./dist_axis,'LineStyle','-','LineWidth',1.5);
% semilogx(r_E(1:108)./1e3,...
%     SF3_E(1:108)./r_E(1:108),'LineWidth',1.5,'Color','k');
semilogx(r_E./1e3,...
    SF3_E./r_E,'LineWidth',1.5,'Color','k');
% semilogx(R_check./1e3,...
%     SF3_check./R_check,'LineWidth',1.5,'Color','m');
grid on
xlim([xmin,xmax]);
xticks([1,4,1e1,1e2,200,1e3])
% xticklabels({'10^{0}','4','10^{1}','10^{2}','200','10^{3}'})
ylim([ymin_sf,ymax_sf]);
xlim([xmin,xmax]);
xticks([1,4,1e1,1e2,200,1e3])
xlabel('$$\mathbf{r \ [km]}$$','Interpreter','latex')
ylabel('$$\mathbf{D3(r)/r \ [m^{2}/s^{3}]}$$','Interpreter','latex')

% legend([b1,b2],{['Lag D3(r)/r'],'fitting'})
text(0.03, 0.95, ['f) ',Casemean,' P1089'], 'Units', 'normalized', ...
     'FontSize', 12, 'FontWeight', 'bold', ...
     'BackgroundColor', [1, 1, 0.8, 0.6], ... % 半透明背景
     'EdgeColor', [.7,.7,.7], ... % 边框颜色
     'Margin', 3, ... % 边距
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')
set(gca,'fontsize',12,'fontweight','bold')

