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
% ini='_rough_500m'
ini='_rough'

% str1=15;
% end1=28;
% str1=1;
% end1=19;
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
SF3_E=SF3(2:end)';
r_E=r(2:end)';
ymin_fk=-3e-8;ymax_fk=6e-8;
ymin_sf=-1e-7;ymax_sf=1e-7;
% ymin_ebsdk=-3e-8;ymax_ebsdk=3e-8;
ymin_ebsdk=-2e-4;ymax_ebsdk=3e-4;
ymin_right = -4e-8;    % 
ymax_right = 6e-8;  % 
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
    Casemean='Hf';
    if strcmpi(ini,'_rough')
        lambda=[1e-7,1e-8,1e-9];
        % str1=15;
        % end1=28;
        str1=1;
        end1=19;
        range1=2.5e3;
        range2=300e3;
    end
    if strcmpi(ini,'_roughsmall') || strcmpi(ini,'_roughsmall_rot') || strcmpi(ini,'_roughLASER') || strcmpi(ini,'_roughsmall_2mon') 
        % lambda=[10,1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6, ...
        % 1e-7,1e-8,1e-9,1e-10,1e-11,2e-11,1e-12,1e-13,1e-14];
        lambda=[1e-7,1e-8,1e-9];
        str1=1;
        end1=19;
        range1=2.5e3;
        range2=300e3;
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
   
    % range2=300e3;
end
if strcmpi(Case, 'nowave')
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
   

end


% 480 h bootstrap
% timerange=1:1940;
% timerange=1:720;
% timerange=1:1920;
% [SpecFlux_Lag289,Vt_Lag289,ebs_Lag289,kf_Lag,...
%     lf_Lag,CI_ebs289,CI_Vt289,...
%     CI_SpecFlux289,Th_Lag289,SF3_289,dist_axis]=get_param_Lag_sf3fk_bootstrap2(Case,ini, ...
%     289,lambda,timerange,inv,range1,range2);
range1=1e3;
range2=300e3;
[dist_axis,SF3_Lag289,SF3_Lag289_std,kf_Lag,dkf,...
    Vt_Lag289,eps_Lag289,SpecFlux_Lag289,Fk_error289,eps_error289,...
    optimal_fac289,optimal_po289,filtscale,...
    Th_Lag289,Th_Lag289_std]=Fk_fitting_SF3_Bayesian_RLS_Lcurve(Case,289,...
    timerange,ini,range1,range2,2);
% x_fill289_kf,y_fill289

[dist_axis,SF3_Lag625,SF3_Lag625_std,kf_Lag,dkf,...
    Vt_Lag625,eps_Lag625,SpecFlux_Lag625,Fk_error625,eps_error625,...
    optimal_fac625,optimal_po625,filtscale,...
    Th_Lag625,Th_Lag625_std]=Fk_fitting_SF3_Bayesian_RLS_Lcurve(Case,625,...
    timerange,ini,range1,range2,2);

[dist_axis,SF3_Lag1089,SF3_Lag1089_std,kf_Lag,dkf,...
    Vt_Lag1089,eps_Lag1089,SpecFlux_Lag1089,Fk_error1089,eps_error1089,...
    optimal_fac1089,optimal_po1089,filtscale,...
    Th_Lag1089,Th_Lag1089_std]=Fk_fitting_SF3_Bayesian_RLS_Lcurve(Case,1089,...
    timerange,ini,range1,range2,2);
% 

if strcmpi(Case, 'nowave')
    colors_rgb{1}=colors_rgb{2};
    colors{1}=colors{2};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%plot%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% patch(x_coords, y_coords, [0.9 0.9 0.9], 'EdgeColor', 'none'); % 中灰色，无边线
screenSize = get(0, 'ScreenSize');
% 创建figure并设置位置
figure('Position', [0, 0, screenSize(3), screenSize(4)]);
subplot(3,3,1)

% patch(x_coords, y_coords, [0.9 0.9 0.9],'FaceAlpha', 0.3, ...
%     'EdgeColor', 'k'); 
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

B2=semilogx(1./kf_Lag(1:end).*2.*pi./1e3,SpecFlux_Lag289(1:end),'LineWidth',1.5,'Color',colors{4})
r_RLS_Lag = findXatYZero(1./kf_Lag(1:end).*2.*pi./1e3,SpecFlux_Lag289(1:end));
[x_fill,y_fill]=get_shadow((1./kf_Lag.*2.*pi./1e3)',...
    (SpecFlux_Lag289+Fk_error289)',(SpecFlux_Lag289-Fk_error289)')
fill(x_fill, y_fill, ...
    colors_rgb{4}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');


grid on
ylim([ymin_fk,ymax_fk]);
xlim([xmin,xmax]);
% xticks([1,4,1e1,1e2,200,1e3])
xticks([1,4,1e1,1e2,200,1e3])

% xticklabels({'10^{0}','4','10^{1}','10^{2}','200','10^{3}'})

xlabel('$$\mathbf{r \ [km]}$$','Interpreter','latex')
ylabel('$$\mathbf{F(r) \ [m^{2}/s^{3}]}$$','Interpreter','latex')
% % title('Hf: CG vs sf3-RLS flux')
% legend([b1,b2],{['CG (',win,')'],'SF3-fitting (RLS)'})
legend([b1,b2,B1,B2],{['Eul CG (',win,')'],['Eul ',inv,'-fitting'], ...
    ['Lag CG (',win,')'],['Lag ',inv,'-fitting']})
text(0.03, 0.95, ['a) ',Casemean,' P289'], 'Units', 'normalized', ...
     'FontSize', 12, 'FontWeight', 'bold', ...
     'BackgroundColor', [1, 1, 0.8, 0.6], ... % 半透明背景
     'EdgeColor', [.7,.7,.7], ... % 边框颜色
     'Margin', 3, ... % 边距
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')
text(10^0.1,yy1,['$\mathbf{\ r_{CG \ Eul}= }$',num2str(r_cg_eul)],'Interpreter','latex', ...
    'Color',colors{1})
text(10^0.1,yy2,['$\mathbf{\ r_{fit \ Eul}= }$',num2str(r_RLS_eul)],'Interpreter','latex')
text(10^0.1,yy3,['$\mathbf{\ r_{CG \ Lag}= }$',num2str(r_cg_Lag)],'Interpreter','latex', ...
    'Color',colors{4})
text(10^0.1,yy4,['$\mathbf{\ r_{fit \ Lag}= }$',num2str(r_RLS_Lag)],'Interpreter','latex', ...
    'Color',colors{4})

set(gca,'fontsize',12,'fontweight','bold')

%%%%%%%%%%%%%%%%%%%%%%%%%debug%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,3,4)
b1=semilogx(dist_axis./1e3,...
    SF3_Lag289./dist_axis,'LineWidth',1.5);
hold on
[x_fill,y_fill]=get_shadow(dist_axis',(SF3_Lag289+SF3_Lag289_std)',...
    (SF3_Lag289-SF3_Lag289_std)');
fill(x_fill./1e3,y_fill./x_fill, ...
     'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
b2=semilogx(dist_axis./1e3,Vt_Lag289./dist_axis,'LineStyle','-','LineWidth',1.5);

b3=semilogx(r_E./1e3,...
    SF3_E./r_E,'LineWidth',1.5,'Color','k');
b4=semilogx(R_check./1e3,...
    SF3_check./R_check,'LineWidth',1.5,'Color','m');
grid on
xlim([xmin,xmax]);
xticks([1,4,1e1,1e2,200,1e3])
% xticklabels({'10^{0}','4','10^{1}','10^{2}','200','10^{3}'})
ylim([ymin_sf,ymax_sf]);
xlim([xmin,xmax]);
xticks([1,4,1e1,1e2,200,1e3])
xlabel('$$\mathbf{r \ [km]}$$','Interpreter','latex')
ylabel('$$\mathbf{D3(r)/r \ [m^{2}/s^{3}]}$$','Interpreter','latex')

legend([b1,b2,b3,b4],{['Lag D3(r)/r'],'fitting',['Eul D3(r)/r'],['check Eul D3(r)/r']},...
    'Location','southeast')
text(0.03, 0.95, ['d) ',Casemean,' P289'], 'Units', 'normalized', ...
     'FontSize', 12, 'FontWeight', 'bold', ...
     'BackgroundColor', [1, 1, 0.8, 0.6], ... % 半透明背景
     'EdgeColor', [.7,.7,.7], ... % 边框颜色
     'Margin', 3, ... % 边距
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')
set(gca,'fontsize',12,'fontweight','bold')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subplot(3,3,7)
% set(gca,'xscale','log')
% hold on
% 
% c4=semilogx(1./kf_Lag(1:end-2)./1e3.*2.*pi,ebs_Lag289(1:end-3).*dk_L(1:end-2)','LineWidth', ...
%     1.5,'Color',colors{4});
% hold on
% 
% fill(x_fillebs,y_fillebs_289, ...
%      colors_rgb{4}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% b1=semilogx(1./kf_E(1:end-2)./1e3.*2.*pi,ebs_E(1:end-3).*dk_E(1:end-2)',...
%     'LineWidth', 1.5,'Color','k');
% b2=semilogx(1./kc_mid_s./1e3,dPidk_s,'LineWidth',1.5,'Color',colors{1});
% 
% % semilogx(1./(kc_mid_s)./1e3,dPidk_s,'LineWidth', ...
% %     1.5,'Color',colors{4},'LineStyle','--');
% % semilogx(1./dk_E./1e3,dPidk_s,'LineWidth', ...
% %     1.5,'Color',colors{1},'LineStyle','-');
% % semilogx(rebs_c,ebs_Lag_cg289,'LineWidth', ...
% %     1.5,'Color',colors{4},'LineStyle','--');
% 
% 
% grid on
% ylim([ymin_ebsdk,ymax_ebsdk]);
% % xlim([1e-6,1e-3].*1e3)
% % xticks([1e-3,1/200,1e-2,1e-1,1/4,1])
% % xticklabels({'1/1000','1/200','1/100','1/10','1/4','1'})
% xlim([xmin,xmax]);
% xticks([1,4,1e1,1e2,200,1e3])
% xlabel('$$\mathbf{r \ [km]}$$','Interpreter','latex')
% % xlabel('$$\mathbf{k \ [1/km]}$$','Interpreter','latex')
% ylabel('$$\mathbf{\epsilon_j*dk_j  \ [m^{2}/s^{3}]}$$','Interpreter','latex')
% % legend([c1,c2,c3,c4],{['Eul CG (',win,')'],'Eul SF3-fitting (RLS)', ...
% %     ['Lag CG (',win,')'],'Lag SF3-fitting (RLS)'})
% text(0.03, 0.95, ['g) ',Casemean,' P289'], 'Units', 'normalized', ...
%      'FontSize', 12, 'FontWeight', 'bold', ...
%      'BackgroundColor', [1, 1, 0.8, 0.6], ... % 半透明背景
%      'EdgeColor', [.7,.7,.7], ... % 边框颜色
%      'Margin', 3, ... % 边距
%      'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')
% set(gca,'fontsize',12,'fontweight','bold')
% 
subplot(3,3,7)
set(gca,'xscale','log')
hold on

% 保存左侧y轴的颜色和线型设置
c4 = semilogx(1./kf_Lag./1e3.*2.*pi, eps_Lag289(2:end), ...
    'LineWidth', 1.5, 'Color', colors{4});

hold on
[x_fill2,y_fill2]=get_shadow(fliplr((1./kf_Lag./1e3.*2.*pi)'),...
    fliplr((eps_Lag289(2:end)+eps_error289(2:end))'),...
    fliplr((eps_Lag289(2:end)-eps_error289(2:end))'));
fill(x_fill2, y_fill2, colors_rgb{4}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% fill(x_fillebs, y_fillebs_289, colors_rgb{4}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
b1 = semilogx(1./kf_E./1e3.*2.*pi, ebs_E(1:end-1), ...
    'LineWidth', 1.5, 'Color', 'k');

grid on
ylim([ymin_ebsdk, ymax_ebsdk]);
xlim([0.5e3, 1e6]./1e3);
xticks([1, 4, 1e1, 1e2, 200, 1e3])
xlabel('$$\mathbf{r \ [km]}$$', 'Interpreter', 'latex')
ylabel('$$\mathbf{\epsilon_j  \ [m^{3}/s^{3}]}$$', 'Interpreter', 'latex')

% 设置左侧y轴的属性
left_ax = gca;
left_ax.YColor = 'k';  % 左侧y轴颜色
left_ax.YLim = [ymin_ebsdk, ymax_ebsdk];

% 创建右侧y轴
yyaxis right

% 
w1 = semilogx(1./K1D./1e3,ww_mean , 'LineWidth', 1.5, 'Color', ...
    [0.7, 0.5, 0.3], 'LineStyle', '-');

% 
right_ax = gca;
right_ax.YColor = [0.7, 0.5, 0.3];  % 
right_ax.YLim = [ymin_right, ymax_right];

% right_ax.YLim = [ymin_right, ymax_right];
ylabel('$\overline{\hat{\tau_{x}}^{*}\hat{u}+\hat{\tau_{y}}^{*}\hat{v}} \ [m^{3}/s^{3}]$', ...
    'Interpreter', 'latex', 'Color', [0.7, 0.5, 0.3])

% 
yyaxis left

%
text(0.03, 0.95, ['g) ', Casemean, ' P289'], 'Units', 'normalized', ...
     'FontSize', 12, 'FontWeight', 'bold', ...
     'BackgroundColor', [1, 1, 0.8, 0.6], ...
     'EdgeColor', [.7,.7,.7], ...
     'Margin', 3, ...
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')
set(gca, 'fontsize', 12, 'fontweight', 'bold')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%625
subplot(3,3,2)

% patch(x_coords, y_coords, [0.9 0.9 0.9],'FaceAlpha', 0.3, ...
%     'EdgeColor', 'k'); 
set(gca,'xscale','log')
hold on
b1=semilogx(filtscale./1e3,Thm_Eulerian,'LineWidth',1.5,'Color',colors{1});
hold on
r_cg_eul = findXatYZero(filtscale./1e3, Thm_Eulerian);

b2=semilogx(1./kf_E./1e3.*2.*pi,SpecFlux_E, ...
    'LineWidth',1.5, ...
'Color','k');
r_RLS_eul = findXatYZero(1./kf_E./1e3.*2.*pi,SpecFlux_E);

B1=semilogx(filtscale./1e3,Th_Lag625,'LineWidth',1.5,'Color',colors{5}, ...
    'LineStyle','--')
% hold on
r_cg_Lag = findXatYZero(filtscale./1e3,Th_Lag625);

B2=semilogx(1./kf_Lag(1:end).*2.*pi./1e3,SpecFlux_Lag625(1:end), ...
    'LineWidth',1.5,'Color',colors{5})
r_RLS_Lag = findXatYZero(1./kf_Lag.*2.*pi./1e3,SpecFlux_Lag625);
[x_fill,y_fill]=get_shadow((1./kf_Lag.*2.*pi./1e3)',...
    (SpecFlux_Lag625+Fk_error625)',(SpecFlux_Lag625-Fk_error625)')
fill(x_fill, y_fill, ...
    colors_rgb{5}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% fill(x_fill, [CI_SpecFlux625(1,1:end), fliplr(CI_SpecFlux625(2,1:end))], ...
%     colors_rgb{5}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
grid on
ylim([ymin_fk,ymax_fk]);
xlim([xmin,xmax]);
xticks([1,4,1e1,1e2,200,1e3])
% xticklabels({'10^{0}','4','10^{1}','10^{2}','200','10^{3}'})

xlabel('$$\mathbf{r \ [km]}$$','Interpreter','latex')
ylabel('$$\mathbf{F(r) \ [m^{2}/s^{3}]}$$','Interpreter','latex')
legend([B1,B2],{['Lag CG (',win,')'],['Lag ',inv,'-fitting']})
text(0.03, 0.95, ['b) ',Casemean,' P625'], 'Units', 'normalized', ...
     'FontSize', 12, 'FontWeight', 'bold', ...
     'BackgroundColor', [1, 1, 0.8, 0.6], ... % 半透明背景
     'EdgeColor', [.7,.7,.7], ... % 边框颜色
     'Margin', 3, ... % 边距
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')
text(10^0.1,yy1,['$\mathbf{\ r_{CG \ Eul}= }$',num2str(r_cg_eul)],'Interpreter','latex', ...
    'Color',colors{1})
text(10^0.1,yy2,['$\mathbf{\ r_{fit \ Eul}= }$',num2str(r_RLS_eul)],'Interpreter','latex')
text(10^0.1,yy3,['$\mathbf{\ r_{CG \ Lag}= }$',num2str(r_cg_Lag)],'Interpreter','latex', ...
    'Color',colors{5})
text(10^0.1,yy4,['$\mathbf{\ r_{fit \ Lag}= }$',num2str(r_RLS_Lag)],'Interpreter','latex', ...
    'Color',colors{5})

set(gca,'fontsize',12,'fontweight','bold')
% %%%%%%%%%%%%%%%%%%%%%%%%%debug%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,3,5)
semilogx(dist_axis./1e3,...
    SF3_Lag625./dist_axis,'LineWidth',1.5);
hold on
[x_fill,y_fill]=get_shadow(dist_axis',(SF3_Lag625+SF3_Lag625_std)',...
    (SF3_Lag625-SF3_Lag625_std)');

fill(x_fill./1e3,y_fill./x_fill, ...
     'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
semilogx(dist_axis./1e3,Vt_Lag625./dist_axis,'LineStyle','-','LineWidth',1.5);
semilogx(r_E./1e3,...
    SF3_E./r_E,'LineWidth',1.5,'Color','k');
semilogx(R_check./1e3,...
    SF3_check./R_check,'LineWidth',1.5,'Color','m');
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
% subplot(3,3,8)
% set(gca,'xscale','log')
% hold on
% 
% c4=semilogx(1./kf_Lag(1:end-2)./1e3.*2.*pi,ebs_Lag625(1:end-3).*dk_L(1:end-2)','LineWidth', ...
%     1.5,'Color',colors{5});
% hold on
% b1=semilogx(1./kf_E(1:end-2)./1e3.*2.*pi,ebs_E(1:end-3).*dk_E(1:end-2)',...
%     'LineWidth', 1.5,'Color','k');
% b2=semilogx(1./kc_mid_s./1e3,dPidk_s,'LineWidth',1.5,'Color',colors{1});
% fill(x_fillebs,y_fillebs_625, ...
%      colors_rgb{5}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% grid on
% ylim([ymin_ebsdk,ymax_ebsdk]);
% % xlim([1e-6,1e-3].*1e3)
% % xticks([1e-3,1/200,1e-2,1e-1,1/4,1])
% % xticklabels({'1/1000','1/200','1/100','1/10','1/4','1'})
% xlim([xmin,xmax]);
% xticks([1,4,1e1,1e2,200,1e3])
% xlabel('$$\mathbf{r \ [km]}$$','Interpreter','latex')
% % xlabel('$$\mathbf{k \ [1/km]}$$','Interpreter','latex')
% ylabel('$$\mathbf{\epsilon_j*dk_j  \ [m^{2}/s^{3}]}$$','Interpreter','latex')
% % legend([c1,c2,c3,c4],{['Eul CG (',win,')'],'Eul SF3-fitting (RLS)', ...
% %     ['Lag CG (',win,')'],'Lag SF3-fitting (RLS)'})
% text(0.03, 0.95, ['h) ',Casemean,' P625'], 'Units', 'normalized', ...
%      'FontSize', 12, 'FontWeight', 'bold', ...
%      'BackgroundColor', [1, 1, 0.8, 0.6], ... % 半透明背景
%      'EdgeColor', [.7,.7,.7], ... % 边框颜色
%      'Margin', 3, ... % 边距
%      'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')
% set(gca,'fontsize',12,'fontweight','bold')
subplot(3,3,8)
set(gca,'xscale','log')
hold on

% 保存左侧y轴的颜色和线型设置
c4 = semilogx(1./kf_Lag./1e3.*2.*pi, eps_Lag625(2:end), ...
    'LineWidth', 1.5, 'Color', colors{5});

hold on
% fill(x_fill, y_fillEBS_625, colors_rgb{5}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% fill(x_fillebs, y_fillebs_289, colors_rgb{4}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
b1 = semilogx(1./kf_E./1e3.*2.*pi, ebs_E(1:end-1), ...
    'LineWidth', 1.5, 'Color', 'k');

grid on
ylim([ymin_ebsdk, ymax_ebsdk]);
xlim([0.5e3, 1e6]./1e3);
xticks([1, 4, 1e1, 1e2, 200, 1e3])
xlabel('$$\mathbf{r \ [km]}$$', 'Interpreter', 'latex')
ylabel('$$\mathbf{\epsilon_j  \ [m^{3}/s^{3}]}$$', 'Interpreter', 'latex')

% 设置左侧y轴的属性
left_ax = gca;
left_ax.YColor = 'k';  % 左侧y轴颜色
left_ax.YLim = [ymin_ebsdk, ymax_ebsdk];

% 创建右侧y轴
yyaxis right

% 
w1 = semilogx(1./K1D./1e3,ww_mean , 'LineWidth', 1.5, 'Color', ...
    [0.7, 0.5, 0.3], 'LineStyle', '-');

% 
right_ax = gca;
right_ax.YColor = [0.7, 0.5, 0.3];  % 
right_ax.YLim = [ymin_right, ymax_right];

% right_ax.YLim = [ymin_right, ymax_right];
ylabel('$\overline{\hat{\tau_{x}}^{*}\hat{u}+\hat{\tau_{y}}^{*}\hat{v}} \ [m^{3}/s^{3}]$', ...
    'Interpreter', 'latex', 'Color', [0.7, 0.5, 0.3])

% 
yyaxis left

%
text(0.03, 0.95, ['h) ', Casemean, ' P625'], 'Units', 'normalized', ...
     'FontSize', 12, 'FontWeight', 'bold', ...
     'BackgroundColor', [1, 1, 0.8, 0.6], ...
     'EdgeColor', [.7,.7,.7], ...
     'Margin', 3, ...
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')
set(gca, 'fontsize', 12, 'fontweight', 'bold')


% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 1089
subplot(3,3,3)
% 
% patch(x_coords, y_coords, [0.9 0.9 0.9],'FaceAlpha', 0.3, ...
%     'EdgeColor', 'k'); 
set(gca,'xscale','log')
hold on
b1=semilogx(filtscale./1e3,Thm_Eulerian,'LineWidth',1.5,'Color',colors{1});
hold on
r_cg_eul = findXatYZero(filtscale./1e3, Thm_Eulerian);

b2=semilogx(1./kf_E./1e3.*2.*pi,SpecFlux_E, ...
    'LineWidth',1.5, ...
'Color','k');
r_RLS_eul = findXatYZero(1./kf_E./1e3.*2.*pi,SpecFlux_E);

B1=semilogx(filtscale./1e3,Th_Lag1089,'LineWidth',1.5,'Color',colors{6}, ...
    'LineStyle','--')
% hold on
r_cg_Lag = findXatYZero(filtscale./1e3,Th_Lag1089);

B2=semilogx(1./kf_Lag.*2.*pi./1e3,SpecFlux_Lag1089,'LineWidth', ...
    1.5,'Color',colors{6})
r_RLS_Lag = findXatYZero(1./kf_Lag.*2.*pi./1e3,SpecFlux_Lag1089);
[x_fill,y_fill]=get_shadow((1./kf_Lag.*2.*pi./1e3)',...
    (SpecFlux_Lag1089+Fk_error1089)',(SpecFlux_Lag1089-Fk_error1089)')
fill(x_fill, y_fill, ...
    colors_rgb{6}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

grid on
ylim([ymin_fk,ymax_fk]);
xlim([xmin,xmax]);
xticks([1,4,1e1,1e2,200,1e3])
legend([B1,B2],{['Lag CG (',win,')'],['Lag ',inv,'-fitting']})

xlabel('$$\mathbf{r \ [km]}$$','Interpreter','latex')
ylabel('$$\mathbf{F(r) \ [m^{2}/s^{3}]}$$','Interpreter','latex')

text(0.03, 0.95, ['c) ',Casemean,' P1089'], 'Units', 'normalized', ...
     'FontSize', 12, 'FontWeight', 'bold', ...
     'BackgroundColor', [1, 1, 0.8, 0.6], ... % 半透明背景
     'EdgeColor', [.7,.7,.7], ... % 边框颜色
     'Margin', 3, ... % 边距
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')
text(10^0.1,yy1,['$\mathbf{\ r_{CG \ Eul}= }$',num2str(r_cg_eul)],'Interpreter','latex', ...
    'Color',colors{1})
text(10^0.1,yy2,['$\mathbf{\ r_{fit \ Eul}= }$',num2str(r_RLS_eul)],'Interpreter','latex')

text(10^0.1,yy3,['$\mathbf{\ r_{CG \ Lag}= }$',num2str(r_cg_Lag)],'Interpreter','latex', ...
    'Color',colors{6})
text(10^0.1,yy4,['$\mathbf{\ r_{fit \ Lag}= }$',num2str(r_RLS_Lag)],'Interpreter','latex', ...
    'Color',colors{6})

set(gca,'fontsize',12,'fontweight','bold')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,3,6)
semilogx(dist_axis./1e3,...
    SF3_Lag1089./dist_axis,'LineWidth',1.5);
hold on
[x_fill,y_fill]=get_shadow(dist_axis',(SF3_Lag1089+SF3_Lag1089_std)',...
    (SF3_Lag1089-SF3_Lag1089_std)');
fill(x_fill2,y_fill2, ...
     'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
semilogx(dist_axis./1e3,Vt_Lag1089./dist_axis,'LineStyle','-','LineWidth',1.5);
semilogx(r_E./1e3,...
    SF3_E./r_E,'LineWidth',1.5,'Color','k');
semilogx(R_check./1e3,...
    SF3_check./R_check,'LineWidth',1.5,'Color','m');
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subplot(3,3,9)
% set(gca,'xscale','log')
% hold on
% 
% c4=semilogx(1./kf_Lag(1:end-2)./1e3.*2.*pi,ebs_Lag1089(1:end-3).*dk_L(1:end-2)','LineWidth', ...
%     1.5,'Color',colors{6});
% hold on
% b1=semilogx(1./kf_E(1:end-2)./1e3.*2.*pi,ebs_E(1:end-3).*dk_E(1:end-2)',...
%     'LineWidth', 1.5,'Color','k');
% b2=semilogx(1./kc_mid_s./1e3,dPidk_s,'LineWidth',1.5,'Color',colors{1});
% fill(x_fillebs,y_fillebs_1089, ...
%      colors_rgb{6}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% grid on
% ylim([ymin_ebsdk,ymax_ebsdk]);
% % xlim([1e-6,1e-3].*1e3)
% % xticks([1e-3,1/200,1e-2,1e-1,1/4,1])
% % xticklabels({'1/1000','1/200','1/100','1/10','1/4','1'})
% xlim([xmin,xmax]);
% xticks([1,4,1e1,1e2,200,1e3])
% xlabel('$$\mathbf{r \ [km]}$$','Interpreter','latex')
% % xlabel('$$\mathbf{k \ [1/km]}$$','Interpreter','latex')
% ylabel('$$\mathbf{\epsilon_j*dk_j  \ [m^{2}/s^{3}]}$$','Interpreter','latex')
% % legend([c1,c2,c3,c4],{['Eul CG (',win,')'],'Eul SF3-fitting (RLS)', ...
% %     ['Lag CG (',win,')'],'Lag SF3-fitting (RLS)'})
% text(0.03, 0.95, ['i) ',Casemean,' P1089'], 'Units', 'normalized', ...
%      'FontSize', 12, 'FontWeight', 'bold', ...
%      'BackgroundColor', [1, 1, 0.8, 0.6], ... % 半透明背景
%      'EdgeColor', [.7,.7,.7], ... % 边框颜色
%      'Margin', 3, ... % 边距
%      'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')
% set(gca,'fontsize',12,'fontweight','bold')

subplot(3,3,9)
set(gca,'xscale','log')
hold on

% 保存左侧y轴的颜色和线型设置
c4 = semilogx(1./kf_Lag./1e3.*2.*pi, eps_Lag1089(2:end), ...
    'LineWidth', 1.5, 'Color', colors{6});

hold on
% fill(x_fill, y_fillEBS_1089, colors_rgb{6}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% fill(x_fillebs, y_fillebs_289, colors_rgb{4}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
b1 = semilogx(1./kf_E./1e3.*2.*pi, ebs_E(1:end-1), ...
    'LineWidth', 1.5, 'Color', 'k');

grid on
ylim([ymin_ebsdk, ymax_ebsdk]);
xlim([0.5e3, 1e6]./1e3);
xticks([1, 4, 1e1, 1e2, 200, 1e3])
xlabel('$$\mathbf{r \ [km]}$$', 'Interpreter', 'latex')
ylabel('$$\mathbf{\epsilon_j  \ [m^{3}/s^{3}]}$$', 'Interpreter', 'latex')

% 设置左侧y轴的属性
left_ax = gca;
left_ax.YColor = 'k';  % 左侧y轴颜色
left_ax.YLim = [ymin_ebsdk, ymax_ebsdk];

% 创建右侧y轴
yyaxis right

% 
w1 = semilogx(1./K1D./1e3,ww_mean , 'LineWidth', 1.5, 'Color', ...
    [0.7, 0.5, 0.3], 'LineStyle', '-');

% 
right_ax = gca;
right_ax.YColor = [0.7, 0.5, 0.3];  % 
right_ax.YLim = [ymin_right, ymax_right];

% right_ax.YLim = [ymin_right, ymax_right];
ylabel('$\overline{\hat{\tau_{x}}^{*}\hat{u}+\hat{\tau_{y}}^{*}\hat{v}} \ [m^{3}/s^{3}]$', ...
    'Interpreter', 'latex', 'Color', [0.7, 0.5, 0.3])

% 
yyaxis left

%
text(0.03, 0.95, ['i) ', Casemean, ' P1089'], 'Units', 'normalized', ...
     'FontSize', 12, 'FontWeight', 'bold', ...
     'BackgroundColor', [1, 1, 0.8, 0.6], ...
     'EdgeColor', [.7,.7,.7], ...
     'Margin', 3, ...
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')
set(gca, 'fontsize', 12, 'fontweight', 'bold')

% 
% 
% saveas(gcf,['Version3_Bayesian_Eul_Lag_CG', ...
%     win,'_vs_RLS_Fr',ini,'_',inv,'_',Casemean],'png')
