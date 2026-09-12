clear all;close all;clc
addpath('/meddy/simingzhang/Analysis/matlab/Parcels_SF/')
addpath('/meddy/simingzhang/Data/Parcels_data')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          1. Basic setup and read data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Case='wave'; % wave
np=[289,625,2500]; % numbers of particles
days=89.5;  % days
dt=3600; % s  Advection_RK4 delta_t drift时间间隔
% input_dir='D:\LIN2023\model\RoyBarkan\LLC4320/'; % drift所在文件夹
ini='_rough'
% lambda=1e-9;

if strcmpi(Case, 'wave')
    cgname='s2sflux_spec_hf_hann_corr.0002.nc';
    cgstd='wavecase_modified_cg_tukey1_RodivofstrainofCG_Eul_std_.mat'
end

if strcmpi(Case, 'nowave')
    cgname='s2sflux_spec_smooth.0002.nc';
    cgstd='nowavecase_modified_cg_tukey1_RodivofstrainofCG_Eul_std_.mat'
end

% cgname='s2sflux_spec_smooth.0002.nc';
ncdisp(cgname)
Thm_Eulerian=ncread(cgname,'Thm');
filtscale=ncread(cgname,'filtscale');
filtscale=filtscale(2:end);
Thm_Eulerian=Thm_Eulerian(2:end);

colors={'#0072BD','#D95319','#EDB120','#7E2F8E'};
colors_rgb = {...
    [0, 114/255, 189/255], ...   % #0072BD
    [217/255, 83/255, 25/255], ... % #D95319
    [237/255, 177/255, 32/255], ... % #EDB120
    [126/255, 47/255, 142/255] ... % #7E2F8E
};
kf1=fliplr(1./filtscale'.*2.*pi);
A={'a1','a2','a3','a4'}
B={'b1','b2','b3','b4'}
lambda=[10,1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6, ...
    1e-7,1e-8,1e-9,1e-10,1e-11,1e-12];
for ii=1:length(np)
    nparticles=np(ii);
    
    outputname=[Case,'_pars_P',num2str(nparticles),'T',num2str(days),...
        'bootstrap.mat']
    load(outputname)
    dstr=1;

    for tt=1:size(SF3,2)
         [SpecFlux(:,tt),Vt(:,tt),ebs(:,tt),kf,lf]=Fk_fitting_SF3_Lcurve(SF3(dstr:end,tt),dist_axis(dstr:end), ...
             3e3,500e3,'log','RLS',lambda,kf1,0);
    end
    
    % [residual_norms, solution_norms,...
    % lambda_opt_idx]=Fk_fitting_SF3_Lcurve(SF3_mean(dstr:end),dist_axis(dstr:end), ...
    % 1,300e3, 'fuc','RLS',[10,1,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6, ...
    %     1e-7,1e-8,1e-9,1e-10,1e-11,1e-12],kf1);
    clear CI_ebs CI_Vt CI_SpecFlux CI_SF3
    
    % Convert parameters to variance preserving form for plotting
    ebs_varp = ebs;
    ebs_varp(1:end-1,:) = ebs(1:end-1,:).*kf';
    
    for i = 1:size(ebs,1)
        CI_ebs(:,i) = prctile(ebs_varp(i,:), [95, 5]); 
    end
    
    for i = 1:size(Vt,1)
        CI_Vt(:,i) = prctile(Vt(i,:), [95,5]);
    end
    
    for i = 1:size(SF3,1)
        CI_SF3(:,i) = prctile(SF3(i,:), [95,5]);
    end
    
    for i = 1:size(SpecFlux,1)    
        CI_SpecFlux(:,i) = prctile(SpecFlux(i,:), [95,5]);
    end
    
    median_ebs = median(ebs_varp,2); 
    mean_ebs = nanmean(ebs_varp,2); 
    median_Vt = median(Vt,2);
    mean_Vt = nanmean(Vt,2);
    
    %median_S = median(S,2);
    mean_SF3 = nanmean(SF3,2);
    
    median_SpecFlux = median(SpecFlux,2);
    mean_SpecFlux = nanmean(SpecFlux,2);
    
    
    figure(1)
    semilogx(dist_axis./1e3,nanmean(SF3,2)./dist_axis','Color',colors{ii}, ...
    'LineWidth',1.5);
    hold on
    semilogx(lf./1e3,mean_Vt./lf','Color',colors{ii}, ...
    'LineWidth',1.5,'linestyle','--');
    grid on
    xlim([1e3,1e6]./1e3);
    xlabel('km')
    ylabel('m^{3}/s^{3}')
    title([Case,': sf3/r'])
    set(gca,'fontsize',16,'FontWeight','b')

    figure(2)
    A{ii}=semilogx(1./kf./1e3.*2.*pi,mean_SpecFlux,'Color',colors{ii}, ...
    'LineWidth',1.5);
    hold on
    x_fill = [1./kf./1e3.*2.*pi, fliplr(1./kf./1e3.*2.*pi)];
    y_fill = [CI_SpecFlux(1,:), fliplr(CI_SpecFlux(2,:))];
    % y_fill = [(nanmean(SpecFlux,2)+std1)', fliplr((nanmean(SpecFlux,2)-std1)')];
    fill(x_fill, y_fill, colors_rgb{ii}, 'FaceAlpha', 0.1, 'EdgeColor', 'none');

    eval(['load ',Case,'_pars_P',num2str(nparticles),'T89.5daysCG_Lag_spectukey', ...
    ini,'.mat']);
    
    Th_all=Th_all(2:end,:);
    Th=Th(2:end);
    std1=std(Th_all, 0, 2, 'omitnan');
    x_fill = [filtscale'./1e3, fliplr(filtscale'./1e3)];
    y_fill = [(Th'+std1)', fliplr((Th'-std1)')];
    B{ii}=semilogx(filtscale./1e3,Th','LineWidth',1.5,'Color',colors{ii},'LineStyle','--');
    fill(x_fill, y_fill, colors_rgb{ii}, 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    grid on
    xlim([1e3,1e6]./1e3);
    ylim([-0.4e-7,0.4e-7]);
    xlabel('km')
    ylabel('m^{2}/s^{3}')
    if ii==3
    eval(['load ',Case,'_Eulerian_SF3.mat'])
    N=287;dstr=1;dot=202;
    
    [r,SF3,S3L1,S3T1]=calc_radial(S3L1,S3T1,N,xscale);
    clear SpecFlux_E;clear kf_E;
    % [SpecFlux_E,Vt_E,ebs_E,kf_E,lf_E]=Fk_fitting_SF3(SF3(3:dot)', ...
    %     r(3:dot),1,300e3,'log','RLS',1e-10,kf1);
    [SpecFlux_E,Vt_E,ebs_E,kf_E,lf_E]=Fk_fitting_SF3(SF3(1:dot)', ...
        r(1:dot),1,300e3,'fuc','RLS',1e-10,kf1);
    c1=semilogx(1./kf_E(dstr:end)./1e3.*2.*pi,(SpecFlux_E(dstr:end)),'LineWidth',1.5,'Color','k','LineStyle','-');
    % legend([A{1}, B{1}, A{2}, B{2}, A{3}, B{3},c1], ...
    % {'RLS-P289','CG-P289','RLS-P625','CG-P625', ...
    % 'RLS-P2500','CG-P2500','Eul fk'})
    end
    title([Case,': Moving Bootstrap RLS-fk'])
    set(gca,'fontsize',16,'FontWeight','b')

    % figure(3)
    % A{ii}=semilogx(1./kf./1e3,mean_ebs(1:end-1),'Color',colors{ii}, ...
    % 'LineWidth',1.5);
    % hold on
    % grid on

    figure(4)
    subplot(2,3,ii)
    A{ii}=semilogx(1./kf./1e3.*2.*pi,mean_SpecFlux,'Color',colors{ii}, ...
    'LineWidth',1.5);
    hold on
    
    B{ii}=semilogx(filtscale./1e3,Th', ...
        'LineWidth',1.5,'Color',colors{ii},'LineStyle','--');
    grid on
    ylim([-5e-8,5e-8]);
    

    figure(5)
    ebs1=nanmean(ebs,2);
    semilogx(kf.*1e3./(2.*pi),ebs1(1:end-1).*abs(mean(diff(kf))),'LineWidth',1.5, ...
    'Color',colors{ii});
    % abs(mean(diff(kf)))
    hold on
    grid on
end
