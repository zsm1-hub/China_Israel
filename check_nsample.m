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
ini='_rough_500m'
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
param_Eul_cg_sf3fk_500m
SF3_E=SF3_mean';
r_E=r;

if strcmpi(Case, 'wave')
    Casemean='Hf';
    if strcmpi(ini,'rough')
        lambda=[10,1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6, ...
        1e-7,1e-8,1e-9,1e-10,1e-11,2e-11,1e-12];
        str1=15;
        end1=28;
        range1=2.5e3;
        range2=500e3;
    end
    if strcmpi(ini,'roughsmall')
        lambda=[10,1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6, ...
        1e-7,1e-8,1e-9,1e-10,1e-11,2e-11,1e-12,1e-13,1e-14];
        str1=1;
        end1=19;
        range1=2.5e3;
        range2=500e3;
    end
    if strcmpi(ini,'_roughsmall_500m')
         lambda=[1e-8];
        str1=1;
        end1=24;
        range1=1.5e3;
        range2=300e3;
    end
    if strcmpi(ini,'_rough_500m')
        lambda=[1e-7];
        str1=1;
        end1=24;
        range1=1e3;
        range2=300e3;
    end
   
    % range2=300e3;
end
if strcmpi(Case, 'nowave')
    Casemean='Sm';
    % lambda=[10,1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6, ...
    % 1e-7,1e-8,1e-9,1e-10,1e-11,2e-11,1e-12,1e-13,1e-14,1e-15,1e-16];
    if strcmpi(ini,'rough')
        lambda=[10,1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6, ...
        1e-7,1e-8,1e-9,1e-10,1e-11,2e-11,1e-12];
        str1=15;
        end1=28;
        range1=2.5e3;
        range2=300e3;
    end
    if strcmpi(ini,'roughsmall')
         lambda=[10,1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6, ...
     1e-7,1e-8,1e-9,1e-10,1e-11,2e-11,1e-12,1e-13,7e-14];
        str1=1;
        end1=19;
        range1=2.5e3;
        range2=300e3;
    end

    if strcmpi(ini,'_roughsmall_500m')
        lambda=[1e-9,1e-10];
        str1=1;
        end1=24;
        range1=1e3;
        range2=300e3;
    end
    if strcmpi(ini,'_rough_500m')
         lambda=[1e-8];
        str1=1;
        end1=24;
        range1=1.5e3;
        range2=300e3;
    end
   

end


timerange=1:1940;

fname=[Case,'_pars_P289T',num2str(timerange(end)),ini,'bootstrap.mat'];
load(fname);
nsample289w=nsample;
nsample289_distriw=nsample289w./nansum(nsample289w);

fname=[Case,'_pars_P625T',num2str(timerange(end)),ini,'bootstrap.mat'];
load(fname);
nsample625w=nsample;
nsample625_distriw=nsample625w./nansum(nsample625w);

fname=[Case,'_pars_P1089T',num2str(timerange(end)),ini,'bootstrap.mat'];
load(fname);
nsample1089w=nsample;
nsample1089_distriw=nsample1089w./nansum(nsample1089w);

%%%%%%%%%%%
Case='nowave'; 
fname=[Case,'_pars_P289T',num2str(timerange(end)),ini,'bootstrap.mat'];
load(fname);
nsample289nw=nsample;
nsample289_distrinw=nsample289nw./nansum(nsample289nw);

fname=[Case,'_pars_P625T',num2str(timerange(end)),ini,'bootstrap.mat'];
load(fname);
nsample625nw=nsample;
nsample625_distrinw=nsample625nw./nansum(nsample625nw);

fname=[Case,'_pars_P1089T',num2str(timerange(end)),ini,'bootstrap.mat'];
load(fname);
nsample1089nw=nsample;
nsample1089_distrinw=nsample1089nw./nansum(nsample1089nw);


%%%%%%%%%%%%%%% plot
figure(2)
semilogx(dist_axis./1e3,...
    nsample289_distriw,'LineWidth',1.5);
hold on
semilogx(dist_axis./1e3,...
    nsample289_distrinw,'LineWidth',1.5);



grid on
xlim([0.5e3,1e6]./1e3);
xticks([1,4,1e1,1e2,200,1e3])
% xticklabels({'10^{0}','4','10^{1}','10^{2}','200','10^{3}'})
% ylim([-3e-7,2e-7]);
xlim([0.5e3,1e6]./1e3);
xticks([1,4,1e1,1e2,200,1e3])
xlabel('$$\mathbf{r \ [km]}$$','Interpreter','latex')
ylabel('$$\mathbf{nsample/tot \ [m^{2}/s^{3}]}$$','Interpreter','latex')
text(0.03, 0.95, ['e) ',Casemean,' P625'], 'Units', 'normalized', ...
     'FontSize', 12, 'FontWeight', 'bold', ...
     'BackgroundColor', [1, 1, 0.8, 0.6], ... % 半透明背景
     'EdgeColor', [.7,.7,.7], ... % 边框颜色
     'Margin', 3, ... % 边距
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')
set(gca,'fontsize',12,'fontweight','bold')

%%%%%%%%%%%%%%%%%%%%%%
figure(3)
semilogx(dist_axis./1e3,...
    nsample289w,'LineWidth',1.5);
hold on
semilogx(dist_axis./1e3,...
    nsample289nw,'LineWidth',1.5);



grid on
xlim([0.5e3,1e6]./1e3);
xticks([1,4,1e1,1e2,200,1e3])
% xticklabels({'10^{0}','4','10^{1}','10^{2}','200','10^{3}'})
% ylim([-3e-7,2e-7]);
xlim([0.5e3,1e6]./1e3);
xticks([1,4,1e1,1e2,200,1e3])
xlabel('$$\mathbf{r \ [km]}$$','Interpreter','latex')
ylabel('$$\mathbf{nsample/tot \ [m^{2}/s^{3}]}$$','Interpreter','latex')
text(0.03, 0.95, ['e) ',Casemean,' P625'], 'Units', 'normalized', ...
     'FontSize', 12, 'FontWeight', 'bold', ...
     'BackgroundColor', [1, 1, 0.8, 0.6], ... % 半透明背景
     'EdgeColor', [.7,.7,.7], ... % 边框颜色
     'Margin', 3, ... % 边距
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')
set(gca,'fontsize',12,'fontweight','bold')

