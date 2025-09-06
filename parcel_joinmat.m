% join mat
clear
Case='wave'; % wave
nparticles=15376; % numbers of particles
np=15376;
days=89.5;  % days
dt=3600; % s  Advection_RK4 delta_t drift时间间隔
% input_dir='D:\LIN2023\model\RoyBarkan\LLC4320/'; % drift所在文件夹
input_dir='/meddy/simingzhang/Data/Parcels_data/';
% timerange=24*10:24*11-6; % 计算结构函数用的时间范围
timerange=1:2140;

if strcmpi(Case, 'wave')
    fname=['wave_pars_P',num2str(nparticles),'T',num2str(days),'days.nc'];
end

if strcmpi(Case, 'nowave')
    fname=['nowave_pars_P',num2str(nparticles),'T',num2str(days),'days.nc'];
end


colors={'#0072BD'};
for jj=1:length(np)
    npart=np(jj);
    SF1l_time=0;
    SF1t_time=0;
    SF2ll_time=0;
    SF2tt_time=0;
    SF3ltt_time=0;
    SF3lll_time=0;
    nvaild_time=0;

    % 集合之前分开计算的chunk
    ii=1;
    eval(['load ',Case,'_pars_P',num2str(npart),'T89.5dayschunk',num2str(ii,'%02d'),'SF123.mat']);
    num_chunk=floor(total_steps./time_batch);
    remainder=mod(total_steps,time_batch);

    for ii = 1:num_chunk
        eval(['load ',Case,'_pars_P',num2str(npart),'T89.5dayschunk',num2str(ii,'%02d'),'SF123.mat'])

        SF1t_time=SF1t_time+SF1t;
        SF1l_time=SF1l_time+SF1l;
        SF2tt_time=SF2tt_time+SF2tt;
        SF2ll_time=SF2ll_time+SF2ll;
        SF3ltt_time=SF3ltt_time+SF3ltt;
        SF3lll_time=SF3lll_time+SF3lll;
        nvaild_time=nvaild_time+nvaild;
        disp(ii);
    end

    if remainder>0
        eval(['load ',Case,'_pars_P',num2str(npart),'T89.5dayschunk',num2str(ii+1,'%02d'),'SF123.mat'])
        % SF3ltt(isnan(SF3ltt))=0;
        % SF3lll(isnan(SF3lll))=0;
        SF1t_time=SF1t_time+SF1t;
        SF1l_time=SF1l_time+SF1l;
        SF2tt_time=SF2tt_time+SF2tt;
        SF2ll_time=SF2ll_time+SF2ll;
        SF3ltt_time=SF3ltt_time+SF3ltt;
        SF3lll_time=SF3lll_time+SF3lll;

        nvaild_time=nvaild_time+nvaild;
    end

    % 计算平均值
    SF1t_time(SF1t_time==0)=nan;
    SF1l_time(SF1l_time==0)=nan;


    SF2tt_time(SF2tt_time==0)=nan;
    SF2ll_time(SF2ll_time==0)=nan;

    SF3ltt_time(SF3ltt_time==0)=nan;
    SF3lll_time(SF3lll_time==0)=nan;

    SF1l_time=SF1l_time./nvaild_time;
    SF1t_time=SF1t_time./nvaild_time;

    SF2ll_time=SF2ll_time./nvaild_time;
    SF2tt_time=SF2tt_time./nvaild_time;

    SF3lll_time=SF3lll_time./nvaild_time;
    SF3ltt_time=SF3ltt_time./nvaild_time;


    % semilogx(dist_axis, (SF3lll_time+SF3ltt_time), 'linewidth',2,...
    % 'Color',colors{jj})
    % hold on
    % loglog(dist_axis, abs(SF3lll_time+SF3ltt_time), 'linewidth',2, ...
    %     'Color',colors{jj})
    % hold on
    % loglog(dist_axis, (SF3lll_time+SF3ltt_time), '+', 'linewidth',2, ...
    %     'Color',colors{jj})
    % loglog(dist_axis, (-SF3lll_time-SF3ltt_time), 'o', 'linewidth',2, ...
    %     'Color',colors{jj})
    % text(dist_axis(1), abs(SF3lll_time(1)+SF3ltt_time(1)),num2str(np(jj)))
end
save([fname(1:end-3),'SF123.mat'],'dist_axis', ...
    'SF1l_time','SF1t_time','SF2ll_time','SF2tt_time','SF3ltt_time','SF3lll_time','-v7.3');

disp('done!!')

%%
% colors={'#0072BD','#D95319','#EDB120','#7E2F8E'};
% jj=1;
% load wave_pars_P289T89.5daysSF123.mat
% loglog(dist_axis, abs(SF3lll+SF3ltt), 'linewidth',2,'Color',colors{jj})
% hold on
% loglog(dist_axis, SF3lll+SF3ltt, '+', 'linewidth',2,'Color',colors{jj})
% loglog(dist_axis, -SF3lll-SF3ltt, 'o', 'linewidth',2,'Color',colors{jj})
% 
% jj=jj+1;
% load wave_pars_P625T89.5daysSF123.mat
% loglog(dist_axis, abs(SF3lll+SF3ltt), 'linewidth',2,'Color',colors{jj})
% loglog(dist_axis, SF3lll+SF3ltt, '+', 'linewidth',2,'Color',colors{jj})
% loglog(dist_axis, -SF3lll-SF3ltt, 'o', 'linewidth',2,'Color',colors{jj})
% 
% jj=jj+1;
% load wave_pars_P2500T89.5daysSF123.mat
% loglog(dist_axis, abs(SF3lll_time+SF3ltt_time), 'linewidth',2, ...
%         'Color',colors{jj})
%     hold on
% loglog(dist_axis, (SF3lll_time+SF3ltt_time), '+', 'linewidth',2, ...
%     'Color',colors{jj})
% loglog(dist_axis, (-SF3lll_time-SF3ltt_time), 'o', 'linewidth',2, ...
%     'Color',colors{jj})
% 
% jj=jj+1;
% load wave_pars_P15376T89.5daysSF123.mat
% loglog(dist_axis, abs(SF3lll_time+SF3ltt_time), 'linewidth',2, ...
%         'Color',colors{jj})
%     hold on
% loglog(dist_axis, (SF3lll_time+SF3ltt_time), '+', 'linewidth',2, ...
%     'Color',colors{jj})
% loglog(dist_axis, (-SF3lll_time-SF3ltt_time), 'o', 'linewidth',2, ...
%     'Color',colors{jj})

%%%%%%%%%%%%%%%%%%%      wave          %%%%%%%%%%%%%%%%%%%%
clear
clf
colors={'#0072BD','#D95319','#EDB120','#7E2F8E'};
jj=1;
load wave_pars_P289T89.5daysSF123.mat
a1=semilogx(dist_axis./1e3, (SF3lll+SF3ltt)./dist_axis, 'linewidth',2,'Color',colors{jj})
hold on


jj=jj+1;
load wave_pars_P625T89.5daysSF123.mat
a2=semilogx(dist_axis./1e3, (SF3lll+SF3ltt)./dist_axis, 'linewidth',2,'Color',colors{jj})

jj=jj+1;
load wave_pars_P2500T89.5daysSF123.mat
a3=semilogx(dist_axis./1e3, (SF3lll_time+SF3ltt_time)./dist_axis, 'linewidth',2, ...
        'Color',colors{jj})

jj=jj+1;
load wave_pars_P15376T89.5daysSF123.mat
a4=semilogx(dist_axis./1e3, (SF3lll_time+SF3ltt_time)./dist_axis, 'linewidth',2, ...
        'Color',colors{jj})

legend(['a1','a2','a3','a4'],{'287','625','2500','15376'},'Location', 'northwest')
grid on
xlim([1e3,1e6]./1e3);
ylim([-2.5e-7,1.5e-7]);
xlabel('km')
ylabel('m^{2}/s^{3}')
title('HF case: Lagrangian SF3/r')
set(gca,'fontsize',16,'FontWeight','b')

%%%%%%%%%%%%%%%%%%%      nowave          %%%%%%%%%%%%%%%%%%%%
clear
colors={'#0072BD','#D95319','#EDB120','#7E2F8E'};
jj=1;
load nowave_pars_P289T89.5daysSF123.mat
a1=semilogx(dist_axis./1e3, (SF3lll+SF3ltt), 'linewidth',2,'Color',colors{jj})
hold on


jj=jj+1;
load nowave_pars_P625T89.5daysSF123.mat
a2=semilogx(dist_axis./1e3, (SF3lll+SF3ltt), 'linewidth',2,'Color',colors{jj})

jj=jj+1;
load nowave_pars_P2500T89.5daysSF123.mat
a3=semilogx(dist_axis./1e3, (SF3lll_time+SF3ltt_time), 'linewidth',2, ...
        'Color',colors{jj})

jj=jj+1;
load nowave_pars_P15376T89.5daysSF123.mat
a4=semilogx(dist_axis./1e3, (SF3lll_time+SF3ltt_time), 'linewidth',2, ...
        'Color',colors{jj})

legend(['a1','a2','a3','a4'],{'287','625','2500','15376'},'Location', 'northwest')
grid on
xlim([1e3,1e6]./1e3);
ylim([-5e-3,15e-3]);
xlabel('km')
ylabel('m^{3}/s^{3}')
title('Smooth case: Lagrangian SF3')
set(gca,'fontsize',14,'FontWeight','b')
