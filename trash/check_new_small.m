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
% ini='_roughsmall_1mon'
ini='_roughsmall_div'
% ini='_roughsmall_2mon'
% str1=15;
% end1=28;
% str1=1;
% end1=19;
inv='RLS'
yy1=-0.6e-8;
yy2=-1.0e-8;
yy3=-1.4e-8;
yy4=-1.8e-8;
npoint2=12

timerange=1:1940;
% timerange=1:720;

% %500 m
% param_Eul_cg_sf3fk_500m
% SF3_E=SF3_mean';
% r_E=r;
% ymin_fk=-10e-8;ymax_fk=10e-8;
% ymin_sf=-2e-7;ymax_sf=3e-7;
% % ymin_ebsdk=-3e-8;ymax_ebsdk=3e-8;
% ymin_ebsdk=-12e-4;ymax_ebsdk=6e-4;
% ymin_right = -6e-8;    % 
% ymax_right = 3e-8;  % 
% xmin=0.5;
% xmax=1e3;
% % 2km param
% param_Eul_cg_sf3fk
% SF3_E=SF3(2:end)';
% r_E=r(2:end)';
% ymin_fk=-3e-8;ymax_fk=6e-8;
% ymin_sf=-1e-7;ymax_sf=1e-7;
% % ymin_ebsdk=-3e-8;ymax_ebsdk=3e-8;
% ymin_ebsdk=-2e-4;ymax_ebsdk=3e-4;
% ymin_right = -4e-8;    % 
% ymax_right = 6e-8;  % 
% xmin=1;
% xmax=1e3;

% 2km param rot
param_Eul_cg_sf3fk_div
SF3_E=nanmean(SF3(2:end,:),2);
r_E=r(2:end)';
ymin_fk=-3e-8;ymax_fk=3e-8;
ymin_sf=-1e-7;ymax_sf=1e-7;
% ymin_ebsdk=-3e-8;ymax_ebsdk=3e-8;
ymin_ebsdk=-2e-4;ymax_ebsdk=3e-4;
ymin_right = -4e-8;    % 
ymax_right = 6e-8;  % 
xmin=1;
xmax=1e3;



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
    if strcmpi(ini,'_roughsmall') || strcmpi(ini,'_roughsmall_rot') || strcmpi(ini,'_roughsmall_div') || strcmpi(ini,'_roughLASER') || strcmpi(ini,'_roughsmall_2mon') 
        % lambda=[10,1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6, ...
        % 1e-7,1e-8,1e-9,1e-10,1e-11,2e-11,1e-12,1e-13,1e-14];
        % lambda=[1e-7,1e-8,1e-9];
        lambda=[1e-5];
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
        range1=2.2e3;
        range2=300e3;
        timerange=1:720;
    end
    if strcmpi(ini,'_roughsmall_500m')
        lambda=[1e-10];
        str1=1;
        end1=24;
        range1=1e3;
        range2=200e3;
    end
    if strcmpi(ini,'_rough_500m') || strcmpi(ini,'_roughbox100g_500m')
        lambda=[1e-7];
        str1=1;
        end1=24;
        % end1=17;
        range1=1e3;
        range2=200e3;
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
    if strcmpi(ini,'_roughsmall') || strcmpi(ini,'_roughsmall_rot') || strcmpi(ini,'_roughsmall_div') || strcmpi(ini,'_roughLASER') || strcmpi(ini,'_roughsmall_2mon') 
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
        lambda=[1e-10];
        str1=1;
        end1=24;
        range1=1e3;
        range2=200e3;
    end
    if strcmpi(ini,'_rough_500m') || strcmpi(ini,'_roughbox100g_500m')
        lambda=[1e-8];
        str1=1;
        % end1=24;
        end1=24;
        range1=1e3;
        range2=100e3;
    end
   

end


% 480 h bootstrap
timerange=1:1940;
% timerange=1:720;
% timerange=1:1920;
[SpecFlux_Lag289,Vt_Lag289,ebs_Lag289,kf_Lag,...
    lf_Lag,CI_ebs289,CI_Vt289,...
    CI_SpecFlux289,Th_Lag289,SF3_289,dist_axis]=get_param_Lag_sf3fk_bootstrap2(Case,ini, ...
    289,lambda,timerange,inv,range1,range2);
dk_L=u2rho_2d(abs(diff(kf_Lag)));


SF3_Lag289=nanmean(SF3_289,2);
for i = 1:size(SF3_Lag289,1)
    CI_SF3_Lag289(:,i) = prctile(SF3_289(i,:), [95, 5]); 
end


% [SpecFlux_Lag625,Vt_Lag625,ebs_Lag625,~,...
%     ~,CI_ebs625,CI_Vt625,...
%     CI_SpecFlux625,Th_Lag625,SF3_625,dist_axis]=get_param_Lag_sf3fk_bootstrap2(Case,ini, ...
%     625,lambda,timerange,inv,range1,range2);
% dk_L=u2rho_2d(abs(diff(kf_Lag)));
% 
% SF3_Lag625=nanmean(SF3_625,2);
% for i = 1:size(SF3_Lag625,1)
%     CI_SF3_Lag625(:,i) = prctile(SF3_625(i,:), [95, 5]); 
% end
% % 
% % % % 
% [SpecFlux_Lag1089,Vt_Lag1089,ebs_Lag1089,~,...
%     ~,CI_ebs1089,CI_Vt1089,...
%     CI_SpecFlux1089,Th_Lag1089,SF3_1089,dist_axis]=get_param_Lag_sf3fk_bootstrap2(Case,ini, ...
%     1089,lambda,timerange,inv,range1,range2);
% % 
% dk_L=u2rho_2d(abs(diff(kf_Lag)));
% 
% 
% SF3_Lag1089=nanmean(SF3_1089,2);
% for i = 1:size(SF3_Lag1089,1)
%     CI_SF3_Lag1089(:,i) = prctile(SF3_1089(i,:), [95, 5]); 
% end

dd1=dist_axis(1:end-1);
dd1=dd1(str1:end1);
x_fill2 = [dd1./1e3, fliplr(dd1)./1e3];

cisf3_289=CI_SF3_Lag289(:,1:end-1);
cisf3_289=cisf3_289(:,str1:end1);
y_fillSF3_289=[cisf3_289(1,:),...
    fliplr(cisf3_289(2,:))]./x_fill2./1e3;
y_fillebs_289=[CI_ebs289(1,1:end-3).*dk_L(1:end-2)...
    fliplr(CI_ebs289(2,1:end-3).*dk_L(1:end-2))];
y_fillEBS_289=[CI_ebs289(1,1:end-1)...
    fliplr(CI_ebs289(2,1:end-1))];

r1 = generate_r1_from_r(1./kf_Lag(1:end-2)./1e3.*2.*pi);
rebs_c=0.5*(r1(1:end-1)+r1(2:end));
var=fliplr(interp1(filtscale./1e3,Th_Lag289,r1));
ebs_Lag_cg289=var(2:end)-var(1:end-1);

% 
% cisf3_625=CI_SF3_Lag625(:,1:end-1);
% cisf3_625=cisf3_625(:,str1:end1);
% y_fillSF3_625=[cisf3_625(1,:),...
%     fliplr(cisf3_625(2,:))]./x_fill2./1e3;
% y_fillebs_625=[CI_ebs625(1,1:end-3).*dk_L(1:end-2)...
%     fliplr(CI_ebs625(2,1:end-3).*dk_L(1:end-2))];
% y_fillEBS_625=[CI_ebs625(1,1:end-1)...
%     fliplr(CI_ebs625(2,1:end-1))];
% % 
% 
% cisf3_1089=CI_SF3_Lag1089(:,1:end-1);
% cisf3_1089=cisf3_1089(:,str1:end1);
% y_fillSF3_1089=[cisf3_1089(1,:),...
%     fliplr(cisf3_1089(2,:))]./x_fill2./1e3;
% y_fillebs_1089=[CI_ebs1089(1,1:end-3).*dk_L(1:end-2)...
%     fliplr(CI_ebs1089(2,1:end-3).*dk_L(1:end-2))];
% y_fillEBS_1089=[CI_ebs1089(1,1:end-1)...
%     fliplr(CI_ebs1089(2,1:end-1))];

dk_L=u2rho_2d(abs(diff(kf_Lag)));
x_fill = [1./kf_Lag(1:end)./1e3.*2.*pi, fliplr(1./kf_Lag(1:end)./1e3.*2.*pi)];
x_fillebs = [1./kf_Lag(1:end-2)./1e3.*2.*pi, fliplr(1./kf_Lag(1:end-2)./1e3.*2.*pi)];


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

fill(x_fill, [CI_SpecFlux289(1,1:end), fliplr(CI_SpecFlux289(2,1:end))], ...
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
b1=semilogx(dist_axis(str1:end1)./1e3,...
    SF3_Lag289(str1:end1)./dist_axis(str1:end1)','LineWidth',1.5);
hold on
fill(x_fill2,y_fillSF3_289, ...
     'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
b2=semilogx(lf_Lag./1e3,Vt_Lag289./lf_Lag','LineStyle','-','LineWidth',1.5);
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

legend([b1,b2,b3],{['Lag D3(r)/r'],'fitting',['Eul D3(r)/r']},...
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
c4 = semilogx(1./kf_Lag./1e3.*2.*pi, ebs_Lag289(1:end-1), ...
    'LineWidth', 1.5, 'Color', colors{4});

hold on
fill(x_fill, y_fillEBS_289, colors_rgb{4}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
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



% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %625
% subplot(3,3,2)
% 
% % patch(x_coords, y_coords, [0.9 0.9 0.9],'FaceAlpha', 0.3, ...
% %     'EdgeColor', 'k'); 
% set(gca,'xscale','log')
% hold on
% b1=semilogx(filtscale./1e3,Thm_Eulerian,'LineWidth',1.5,'Color',colors{1});
% hold on
% r_cg_eul = findXatYZero(filtscale./1e3, Thm_Eulerian);
% 
% b2=semilogx(1./kf_E(dstr:end)./1e3.*2.*pi,SpecFlux_E(dstr:end), ...
%     'LineWidth',1.5, ...
% 'Color','k');
% r_RLS_eul = findXatYZero(1./kf_E(dstr:end)./1e3.*2.*pi,SpecFlux_E(dstr:end));
% 
% B1=semilogx(filtscale./1e3,Th_Lag625,'LineWidth',1.5,'Color',colors{5}, ...
%     'LineStyle','--')
% % hold on
% r_cg_Lag = findXatYZero(filtscale./1e3,Th_Lag625);
% 
% B2=semilogx(1./kf_Lag(1:end).*2.*pi./1e3,SpecFlux_Lag625(1:end), ...
%     'LineWidth',1.5,'Color',colors{5})
% r_RLS_Lag = findXatYZero(1./kf_Lag(1:end).*2.*pi./1e3,SpecFlux_Lag625(1:end));
% 
% fill(x_fill, [CI_SpecFlux625(1,1:end), fliplr(CI_SpecFlux625(2,1:end))], ...
%     colors_rgb{5}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% 
% 
% grid on
% ylim([ymin_fk,ymax_fk]);
% xlim([xmin,xmax]);
% xticks([1,4,1e1,1e2,200,1e3])
% % xticklabels({'10^{0}','4','10^{1}','10^{2}','200','10^{3}'})
% 
% xlabel('$$\mathbf{r \ [km]}$$','Interpreter','latex')
% ylabel('$$\mathbf{F(r) \ [m^{2}/s^{3}]}$$','Interpreter','latex')
% legend([B1,B2],{['Lag CG (',win,')'],['Lag ',inv,'-fitting']})
% text(0.03, 0.95, ['b) ',Casemean,' P625'], 'Units', 'normalized', ...
%      'FontSize', 12, 'FontWeight', 'bold', ...
%      'BackgroundColor', [1, 1, 0.8, 0.6], ... % 半透明背景
%      'EdgeColor', [.7,.7,.7], ... % 边框颜色
%      'Margin', 3, ... % 边距
%      'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')
% text(10^0.1,yy1,['$\mathbf{\ r_{CG \ Eul}= }$',num2str(r_cg_eul)],'Interpreter','latex', ...
%     'Color',colors{1})
% text(10^0.1,yy2,['$\mathbf{\ r_{fit \ Eul}= }$',num2str(r_RLS_eul)],'Interpreter','latex')
% text(10^0.1,yy3,['$\mathbf{\ r_{CG \ Lag}= }$',num2str(r_cg_Lag)],'Interpreter','latex', ...
%     'Color',colors{5})
% text(10^0.1,yy4,['$\mathbf{\ r_{fit \ Lag}= }$',num2str(r_RLS_Lag)],'Interpreter','latex', ...
%     'Color',colors{5})
% 
% set(gca,'fontsize',12,'fontweight','bold')
% % %%%%%%%%%%%%%%%%%%%%%%%%%debug%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subplot(3,3,5)
% semilogx(dist_axis(str1:end1)./1e3,...
%     SF3_Lag625(str1:end1)./dist_axis(str1:end1)','LineWidth',1.5);
% hold on
% fill(x_fill2,y_fillSF3_625, ...
%      'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% semilogx(lf_Lag./1e3,Vt_Lag625./lf_Lag','LineStyle','-','LineWidth',1.5);
% semilogx(r_E./1e3,...
%     SF3_E./r_E,'LineWidth',1.5,'Color','k');
% % semilogx(R_check./1e3,...
% %     SF3_check./R_check,'LineWidth',1.5,'Color','m');
% grid on
% xlim([xmin,xmax]);
% xticks([1,4,1e1,1e2,200,1e3])
% % xticklabels({'10^{0}','4','10^{1}','10^{2}','200','10^{3}'})
% ylim([ymin_sf,ymax_sf]);
% xlim([xmin,xmax]);
% xticks([1,4,1e1,1e2,200,1e3])
% xlabel('$$\mathbf{r \ [km]}$$','Interpreter','latex')
% ylabel('$$\mathbf{D3(r)/r \ [m^{2}/s^{3}]}$$','Interpreter','latex')
% text(0.03, 0.95, ['e) ',Casemean,' P625'], 'Units', 'normalized', ...
%      'FontSize', 12, 'FontWeight', 'bold', ...
%      'BackgroundColor', [1, 1, 0.8, 0.6], ... % 半透明背景
%      'EdgeColor', [.7,.7,.7], ... % 边框颜色
%      'Margin', 3, ... % 边距
%      'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')
% set(gca,'fontsize',12,'fontweight','bold')
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % subplot(3,3,8)
% % set(gca,'xscale','log')
% % hold on
% % 
% % c4=semilogx(1./kf_Lag(1:end-2)./1e3.*2.*pi,ebs_Lag625(1:end-3).*dk_L(1:end-2)','LineWidth', ...
% %     1.5,'Color',colors{5});
% % hold on
% % b1=semilogx(1./kf_E(1:end-2)./1e3.*2.*pi,ebs_E(1:end-3).*dk_E(1:end-2)',...
% %     'LineWidth', 1.5,'Color','k');
% % b2=semilogx(1./kc_mid_s./1e3,dPidk_s,'LineWidth',1.5,'Color',colors{1});
% % fill(x_fillebs,y_fillebs_625, ...
% %      colors_rgb{5}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% % grid on
% % ylim([ymin_ebsdk,ymax_ebsdk]);
% % % xlim([1e-6,1e-3].*1e3)
% % % xticks([1e-3,1/200,1e-2,1e-1,1/4,1])
% % % xticklabels({'1/1000','1/200','1/100','1/10','1/4','1'})
% % xlim([xmin,xmax]);
% % xticks([1,4,1e1,1e2,200,1e3])
% % xlabel('$$\mathbf{r \ [km]}$$','Interpreter','latex')
% % % xlabel('$$\mathbf{k \ [1/km]}$$','Interpreter','latex')
% % ylabel('$$\mathbf{\epsilon_j*dk_j  \ [m^{2}/s^{3}]}$$','Interpreter','latex')
% % % legend([c1,c2,c3,c4],{['Eul CG (',win,')'],'Eul SF3-fitting (RLS)', ...
% % %     ['Lag CG (',win,')'],'Lag SF3-fitting (RLS)'})
% % text(0.03, 0.95, ['h) ',Casemean,' P625'], 'Units', 'normalized', ...
% %      'FontSize', 12, 'FontWeight', 'bold', ...
% %      'BackgroundColor', [1, 1, 0.8, 0.6], ... % 半透明背景
% %      'EdgeColor', [.7,.7,.7], ... % 边框颜色
% %      'Margin', 3, ... % 边距
% %      'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')
% % set(gca,'fontsize',12,'fontweight','bold')
% subplot(3,3,8)
% set(gca,'xscale','log')
% hold on
% 
% % 保存左侧y轴的颜色和线型设置
% c4 = semilogx(1./kf_Lag./1e3.*2.*pi, ebs_Lag625(1:end-1), ...
%     'LineWidth', 1.5, 'Color', colors{5});
% 
% hold on
% fill(x_fill, y_fillEBS_625, colors_rgb{5}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% % fill(x_fillebs, y_fillebs_289, colors_rgb{4}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% b1 = semilogx(1./kf_E./1e3.*2.*pi, ebs_E(1:end-1), ...
%     'LineWidth', 1.5, 'Color', 'k');
% 
% grid on
% ylim([ymin_ebsdk, ymax_ebsdk]);
% xlim([0.5e3, 1e6]./1e3);
% xticks([1, 4, 1e1, 1e2, 200, 1e3])
% xlabel('$$\mathbf{r \ [km]}$$', 'Interpreter', 'latex')
% ylabel('$$\mathbf{\epsilon_j  \ [m^{3}/s^{3}]}$$', 'Interpreter', 'latex')
% 
% % 设置左侧y轴的属性
% left_ax = gca;
% left_ax.YColor = 'k';  % 左侧y轴颜色
% left_ax.YLim = [ymin_ebsdk, ymax_ebsdk];
% 
% % 创建右侧y轴
% yyaxis right
% 
% % 
% w1 = semilogx(1./K1D./1e3,ww_mean , 'LineWidth', 1.5, 'Color', ...
%     [0.7, 0.5, 0.3], 'LineStyle', '-');
% 
% % 
% right_ax = gca;
% right_ax.YColor = [0.7, 0.5, 0.3];  % 
% right_ax.YLim = [ymin_right, ymax_right];
% 
% % right_ax.YLim = [ymin_right, ymax_right];
% ylabel('$\overline{\hat{\tau_{x}}^{*}\hat{u}+\hat{\tau_{y}}^{*}\hat{v}} \ [m^{3}/s^{3}]$', ...
%     'Interpreter', 'latex', 'Color', [0.7, 0.5, 0.3])
% 
% % 
% yyaxis left
% 
% %
% text(0.03, 0.95, ['h) ', Casemean, ' P625'], 'Units', 'normalized', ...
%      'FontSize', 12, 'FontWeight', 'bold', ...
%      'BackgroundColor', [1, 1, 0.8, 0.6], ...
%      'EdgeColor', [.7,.7,.7], ...
%      'Margin', 3, ...
%      'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')
% set(gca, 'fontsize', 12, 'fontweight', 'bold')
% 
% 
% % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % 1089
% subplot(3,3,3)
% % 
% % patch(x_coords, y_coords, [0.9 0.9 0.9],'FaceAlpha', 0.3, ...
% %     'EdgeColor', 'k'); 
% set(gca,'xscale','log')
% hold on
% b1=semilogx(filtscale./1e3,Thm_Eulerian,'LineWidth',1.5,'Color',colors{1});
% hold on
% r_cg_eul = findXatYZero(filtscale./1e3, Thm_Eulerian);
% 
% b2=semilogx(1./kf_E(dstr:end)./1e3.*2.*pi,SpecFlux_E(dstr:end), ...
%     'LineWidth',1.5, ...
% 'Color','k');
% r_RLS_eul = findXatYZero(1./kf_E(dstr:end)./1e3.*2.*pi,SpecFlux_E(dstr:end));
% 
% B1=semilogx(filtscale./1e3,Th_Lag1089,'LineWidth',1.5,'Color',colors{6}, ...
%     'LineStyle','--')
% % hold on
% r_cg_Lag = findXatYZero(filtscale./1e3,Th_Lag1089);
% 
% B2=semilogx(1./kf_Lag(1:end).*2.*pi./1e3,SpecFlux_Lag1089(1:end),'LineWidth', ...
%     1.5,'Color',colors{6})
% r_RLS_Lag = findXatYZero(1./kf_Lag(1:end).*2.*pi./1e3,SpecFlux_Lag1089(1:end));
% 
% fill(x_fill, [CI_SpecFlux1089(1,1:end), fliplr(CI_SpecFlux1089(2,1:end))], ...
%     colors_rgb{6}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% 
% grid on
% ylim([ymin_fk,ymax_fk]);
% xlim([xmin,xmax]);
% xticks([1,4,1e1,1e2,200,1e3])
% legend([B1,B2],{['Lag CG (',win,')'],['Lag ',inv,'-fitting']})
% 
% xlabel('$$\mathbf{r \ [km]}$$','Interpreter','latex')
% ylabel('$$\mathbf{F(r) \ [m^{2}/s^{3}]}$$','Interpreter','latex')
% 
% text(0.03, 0.95, ['c) ',Casemean,' P1089'], 'Units', 'normalized', ...
%      'FontSize', 12, 'FontWeight', 'bold', ...
%      'BackgroundColor', [1, 1, 0.8, 0.6], ... % 半透明背景
%      'EdgeColor', [.7,.7,.7], ... % 边框颜色
%      'Margin', 3, ... % 边距
%      'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')
% text(10^0.1,yy1,['$\mathbf{\ r_{CG \ Eul}= }$',num2str(r_cg_eul)],'Interpreter','latex', ...
%     'Color',colors{1})
% text(10^0.1,yy2,['$\mathbf{\ r_{fit \ Eul}= }$',num2str(r_RLS_eul)],'Interpreter','latex')
% 
% text(10^0.1,yy3,['$\mathbf{\ r_{CG \ Lag}= }$',num2str(r_cg_Lag)],'Interpreter','latex', ...
%     'Color',colors{6})
% text(10^0.1,yy4,['$\mathbf{\ r_{fit \ Lag}= }$',num2str(r_RLS_Lag)],'Interpreter','latex', ...
%     'Color',colors{6})
% 
% set(gca,'fontsize',12,'fontweight','bold')
% 
% 
% subplot(3,3,6)
% semilogx(dist_axis(str1:end1)./1e3,...
%     SF3_Lag1089(str1:end1)./dist_axis(str1:end1)','LineWidth',1.5);
% hold on
% fill(x_fill2,y_fillSF3_1089, ...
%      'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% semilogx(lf_Lag./1e3,Vt_Lag1089./lf_Lag','LineStyle','-','LineWidth',1.5);
% semilogx(r_E./1e3,...
%     SF3_E./r_E,'LineWidth',1.5,'Color','k');
% % semilogx(R_check./1e3,...
% %     SF3_check./R_check,'LineWidth',1.5,'Color','m');
% grid on
% xlim([xmin,xmax]);
% xticks([1,4,1e1,1e2,200,1e3])
% % xticklabels({'10^{0}','4','10^{1}','10^{2}','200','10^{3}'})
% ylim([ymin_sf,ymax_sf]);
% xlim([xmin,xmax]);
% xticks([1,4,1e1,1e2,200,1e3])
% xlabel('$$\mathbf{r \ [km]}$$','Interpreter','latex')
% ylabel('$$\mathbf{D3(r)/r \ [m^{2}/s^{3}]}$$','Interpreter','latex')
% 
% % legend([b1,b2],{['Lag D3(r)/r'],'fitting'})
% text(0.03, 0.95, ['f) ',Casemean,' P1089'], 'Units', 'normalized', ...
%      'FontSize', 12, 'FontWeight', 'bold', ...
%      'BackgroundColor', [1, 1, 0.8, 0.6], ... % 半透明背景
%      'EdgeColor', [.7,.7,.7], ... % 边框颜色
%      'Margin', 3, ... % 边距
%      'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')
% set(gca,'fontsize',12,'fontweight','bold')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % subplot(3,3,9)
% % set(gca,'xscale','log')
% % hold on
% % 
% % c4=semilogx(1./kf_Lag(1:end-2)./1e3.*2.*pi,ebs_Lag1089(1:end-3).*dk_L(1:end-2)','LineWidth', ...
% %     1.5,'Color',colors{6});
% % hold on
% % b1=semilogx(1./kf_E(1:end-2)./1e3.*2.*pi,ebs_E(1:end-3).*dk_E(1:end-2)',...
% %     'LineWidth', 1.5,'Color','k');
% % b2=semilogx(1./kc_mid_s./1e3,dPidk_s,'LineWidth',1.5,'Color',colors{1});
% % fill(x_fillebs,y_fillebs_1089, ...
% %      colors_rgb{6}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% % grid on
% % ylim([ymin_ebsdk,ymax_ebsdk]);
% % % xlim([1e-6,1e-3].*1e3)
% % % xticks([1e-3,1/200,1e-2,1e-1,1/4,1])
% % % xticklabels({'1/1000','1/200','1/100','1/10','1/4','1'})
% % xlim([xmin,xmax]);
% % xticks([1,4,1e1,1e2,200,1e3])
% % xlabel('$$\mathbf{r \ [km]}$$','Interpreter','latex')
% % % xlabel('$$\mathbf{k \ [1/km]}$$','Interpreter','latex')
% % ylabel('$$\mathbf{\epsilon_j*dk_j  \ [m^{2}/s^{3}]}$$','Interpreter','latex')
% % % legend([c1,c2,c3,c4],{['Eul CG (',win,')'],'Eul SF3-fitting (RLS)', ...
% % %     ['Lag CG (',win,')'],'Lag SF3-fitting (RLS)'})
% % text(0.03, 0.95, ['i) ',Casemean,' P1089'], 'Units', 'normalized', ...
% %      'FontSize', 12, 'FontWeight', 'bold', ...
% %      'BackgroundColor', [1, 1, 0.8, 0.6], ... % 半透明背景
% %      'EdgeColor', [.7,.7,.7], ... % 边框颜色
% %      'Margin', 3, ... % 边距
% %      'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')
% % set(gca,'fontsize',12,'fontweight','bold')
% 
% subplot(3,3,9)
% set(gca,'xscale','log')
% hold on
% 
% % 保存左侧y轴的颜色和线型设置
% c4 = semilogx(1./kf_Lag./1e3.*2.*pi, ebs_Lag1089(1:end-1), ...
%     'LineWidth', 1.5, 'Color', colors{6});
% 
% hold on
% fill(x_fill, y_fillEBS_1089, colors_rgb{6}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% % fill(x_fillebs, y_fillebs_289, colors_rgb{4}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% b1 = semilogx(1./kf_E./1e3.*2.*pi, ebs_E(1:end-1), ...
%     'LineWidth', 1.5, 'Color', 'k');
% 
% grid on
% ylim([ymin_ebsdk, ymax_ebsdk]);
% xlim([0.5e3, 1e6]./1e3);
% xticks([1, 4, 1e1, 1e2, 200, 1e3])
% xlabel('$$\mathbf{r \ [km]}$$', 'Interpreter', 'latex')
% ylabel('$$\mathbf{\epsilon_j  \ [m^{3}/s^{3}]}$$', 'Interpreter', 'latex')
% 
% % 设置左侧y轴的属性
% left_ax = gca;
% left_ax.YColor = 'k';  % 左侧y轴颜色
% left_ax.YLim = [ymin_ebsdk, ymax_ebsdk];
% 
% % 创建右侧y轴
% yyaxis right
% 
% % 
% w1 = semilogx(1./K1D./1e3,ww_mean , 'LineWidth', 1.5, 'Color', ...
%     [0.7, 0.5, 0.3], 'LineStyle', '-');
% 
% % 
% right_ax = gca;
% right_ax.YColor = [0.7, 0.5, 0.3];  % 
% right_ax.YLim = [ymin_right, ymax_right];
% 
% % right_ax.YLim = [ymin_right, ymax_right];
% ylabel('$\overline{\hat{\tau_{x}}^{*}\hat{u}+\hat{\tau_{y}}^{*}\hat{v}} \ [m^{3}/s^{3}]$', ...
%     'Interpreter', 'latex', 'Color', [0.7, 0.5, 0.3])
% 
% % 
% yyaxis left
% 
% %
% text(0.03, 0.95, ['i) ', Casemean, ' P1089'], 'Units', 'normalized', ...
%      'FontSize', 12, 'FontWeight', 'bold', ...
%      'BackgroundColor', [1, 1, 0.8, 0.6], ...
%      'EdgeColor', [.7,.7,.7], ...
%      'Margin', 3, ...
%      'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')
% set(gca, 'fontsize', 12, 'fontweight', 'bold')
% 
% % 
% % 
% saveas(gcf,['Version2_Eul_Lag_CG', ...
%     win,'_vs_RLS_Fr',ini,'_',inv,'_',Casemean],'png')
% 
% 
% figure(2)
% B2=semilogx(1./kf_Lag(1:end).*2.*pi./1e3,SpecFlux_Lag289(1:end),'LineWidth',1.5,'Color',colors{4});hold on
% B3=semilogx(1./kf_Lag(1:end).*2.*pi./1e3,SpecFlux_Lag625(1:end),'LineWidth',1.5,'Color',colors{5})
% B4=semilogx(1./kf_Lag(1:end).*2.*pi./1e3,SpecFlux_Lag1089(1:end),'LineWidth',1.5,'Color',colors{6})
