%   data source: Iceland (wave and no wave case)
%   utility: import parcels data to calc SF2,SF3
%   doesn't use bootstrap to resample,insteadly, calc time-mean SF2 and SF3 directly
%   code writer: zsm, modified from Balwada 2022 sciadv supplyment
clear all;close all;clc
% addpath('D:\LIN2023\model\RoyBarkan\LLC4320/')
% addpath('D:\LIN2023\crocotools\Preprocessingtools') % add function "spheric_dist.m"
% 
addpath('/meddy/simingzhang/Analysis/matlab/Parcels_SF/')
addpath('/meddy/simingzhang/Data/Parcels_data')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          1. Basic setup and read data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Case='nowave'; % wave
nparticles=15376; % numbers of particles
days=89.5;  % days
dt=3600; % s  Advection_RK4 delta_t drift时间间隔
% input_dir='D:\LIN2023\model\RoyBarkan\LLC4320/'; % drift所在文件夹
input_dir='/meddy/simingzhang/Data/Parcels_data/onetime_spectukey/';
% timerange=24*10:24*11-6; % 计算结构函数用的时间范围
time_batch=2;

if strcmpi(Case, 'wave')
    fname=[input_dir,'wave_pars_P',num2str(nparticles),'T',num2str(days),'days.nc'];
end

if strcmpi(Case, 'nowave')
    fname=[input_dir,'nowave_pars_P',num2str(nparticles),'T',num2str(days),'days.nc'];
end

lon=ncread(fname,'lon');
lat=ncread(fname,'lat');

ue=ncread(fname,'ue');
ve=ncread(fname,'ve');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              2. Calc Lagrangian Velocity and save *traj.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for t=1:size(lon,1)-1
    U(t,:)=(spheric_dist(lat(t,:),lat(t,:),lon(t,:),lon(t+1,:)))./dt.*...
        sign(lon(t+1,:)-lon(t,:));
    V(t,:)=(spheric_dist(lat(t+1,:),lat(t,:),lon(t,:),lon(t,:)))./dt.*...
        sign(lat(t+1,:)-lat(t,:));
end
lon(end,:)=[];lat(end,:)=[];

lon=lon;
lat=lat;
u=U;
v=V;
%%%%%%%%%%%%%%%%%%%%%there is a 2D experiment, So I assuming H=-501
%%%%%%%%%%%%%%%%%%%%%H 没有意义, 只是在Balwada的code里面只采样了500米以上的粒子

traj=struct();
traj.trajmat_X=lon;traj.trajmat_Y=lat;
traj.trajmat_U=u;traj.trajmat_V=v;
traj.H=-520.*ones(size(v,1),size(v,2));
traj.T_axis=linspace(dt, (size(v,1))*dt, size(v,1))./86400;

% save([fname(1:end-3),'traj.mat'],'traj');

%%%%%%%%%%%%%%%%%%% 时间维度按time_batchsize分块，防止内存溢出
TRAJ=traj;
total_steps = size(traj.trajmat_X, 1);  % 511
num_chunks = floor(total_steps / time_batch);  % 10
remainder = mod(total_steps, time_batch);      % 11
clear traj;


for i = 1:num_chunks
    % 当前分块的索引范围
    start_idx = (i-1)*time_batch + 1;
    end_idx = i*time_batch;
    
    % 提取当前分块的数据
    traj = struct();
    traj.trajmat_X = TRAJ.trajmat_X(start_idx:end_idx, :);
    traj.trajmat_Y = TRAJ.trajmat_Y(start_idx:end_idx, :);
    traj.trajmat_U = TRAJ.trajmat_U(start_idx:end_idx, :);
    traj.trajmat_V = TRAJ.trajmat_V(start_idx:end_idx, :);
    traj.H = TRAJ.H(start_idx:end_idx, :);
    traj.T_axis = TRAJ.T_axis(start_idx:end_idx);
    
    % 保存为MAT文件
    % save([fname(1:end-3),'traj.mat'],'traj');
    save([fname(1:end-3),'chunk',num2str(i,'%02d'),'traj.mat'], 'traj');
    fprintf('Saved chunk %d (steps %d-%d)\n', i, start_idx, end_idx);
end

% 处理余数部分
if remainder > 0
    start_idx = num_chunks*time_batch + 1;
    end_idx = total_steps;
    
    traj = struct();
    traj.trajmat_X = TRAJ.trajmat_X(start_idx:end_idx, :);
    traj.trajmat_Y = TRAJ.trajmat_Y(start_idx:end_idx, :);
    traj.trajmat_U = TRAJ.trajmat_U(start_idx:end_idx, :);
    traj.trajmat_V = TRAJ.trajmat_V(start_idx:end_idx, :);
    traj.H = TRAJ.H(start_idx:end_idx, :);
    traj.T_axis = TRAJ.T_axis(start_idx:end_idx);
    
    save([fname(1:end-3),'chunk',num2str(i+1,'%02d'),'traj.mat'], 'traj');
    fprintf('Saved remainder (steps %d-%d)\n', start_idx, end_idx);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 3. Calc pairs of particles' variables 
%                    pair_time 储存不同粒子对的变量
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% follow Balwada's routine : "trajectories2binnedpairs_vectorized.m"
% need code "dist_rx.m","dist_ry.m","dist_geo.m","dist_du.m"


%% Break up by time and separation bins
% took 58.9s for GLAD
% took 146s for LASER 
% parfor ii=1:num_chunks+1
%     eval(['load ',[fname(1:end-3),'chunk',num2str(ii,'%02d'),'traj.mat']])
%     Htraj=traj.H;
%     tic 
%     for i=1:length(traj.T_axis)
%         disp(i)
%         %id = find(~isnan(traj.trajmat_X(i,:))); % find non-NaN
% 
%         %id = find(~isnan(traj.trajmat_X(i,:)) & Htraj(i,:)<-500); % find non-Nan and deep
%         % id = find(~isnan(traj.trajmat_X(i,:)) & Htraj(i,:)<-500 & ...
%         %     traj.trajmat_X(i,:)>=-91 & traj.trajmat_X(i,:)<=-84 & ...
%         %     traj.trajmat_Y(i,:)>=24); % find non-Nan and deep and in similar region to GLAD
%         id = find(Htraj(i,:)<-500);
%         % if mod(i,300)==0
%         %     disp(i)
%         % end
% 
%         if length(id)>1
%             X = traj.trajmat_X(i,id)';
%             Y = traj.trajmat_Y(i,id)';
%             U = traj.trajmat_U(i,id)';
%             V = traj.trajmat_V(i,id)';
% 
%             Xvec = [X, Y];
% 
%             pairs_time(i).dist = pdist(Xvec, @dist_geo);
% 
%             rx = pdist(Xvec, @dist_rx);
%             ry = pdist(Xvec, @dist_ry);
% 
%             magr = sqrt(rx.^2 + ry.^2);
% 
%             rx = rx./magr; ry = ry./magr;
% 
%             dux = pdist(U, @dist_du);
%             duy = pdist(V, @dist_du);
% 
%             pairs_time(i).dul = dux.*rx + duy.*ry;
%             pairs_time(i).dut = duy.*rx - dux.*ry;
%         else
%             pairs_time(i).dul = NaN;
%             pairs_time(i).dut = NaN;
%             pairs_time(i).dist = NaN;
%         end
%     end
%     toc
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %              4. Calc 2-order and 3-order structure function (time-mean)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % follow Balwada's routine : "pairtimes2SF.m"
% 
%     tpts = length(pairs_time);
%     %
%     npairs = zeros(tpts,1);
%     for i = 1:tpts
%         % npairs(i) = length(find(~isnan(pairs_time(i).dul)));
%         npairs(i) = length(pairs_time(i).dist);;
%     end
% 
%     %% % Align pairs in a single vector
%     % 
%     dul = zeros(sum(npairs),1);
%     dut = zeros(sum(npairs),1);
%     dist = zeros(sum(npairs),1);
% 
%     %
%     % estimate num of pairs
% 
%     empty1 = 1;
%     for i = 1:tpts % time loop
%         if npairs(i) == 0
%             continue
%         end
% 
%         if npairs(i) == 1
%             dist(empty1) = pairs_time(i).dist;
%             dul(empty1) = pairs_time(i).dul;
%             dut(empty1) = pairs_time(i).dut;
%             empty1 = empty1+1;
%         end
% 
%         if npairs(i) >1
%             dist(empty1: empty1+npairs(i)-1) = pairs_time(i).dist;
%             dul(empty1: empty1+npairs(i)-1) = pairs_time(i).dul;
%             dut(empty1: empty1+npairs(i)-1) = pairs_time(i).dut;
%             empty1 = empty1+npairs(i);
%         end
%         disp(i)
%     end
% 
%     %%
%     clear pairs_time
% 
%     gamma = 1.5;
% 
%     dist_bin(1) = 10; % in m
%     dist_bin = gamma.^[0:100]*dist_bin(1);
%     id = find(dist_bin>1000*10^3,1);
%     dist_bin = dist_bin(1:id);
%     dist_bin(2:end+1) = dist_bin(1:end);
%     dist_bin(1) = 0;
%     dist_axis = 0.5*(dist_bin(1:end-1) + dist_bin(2:end));
%     % dist_axis=dist_axis(dist_axis>2e3);
%     % Generate vel axis
%     % vel_bins = linspace(-2, 2, 50);
%     % vel_axis = 0.5*(vel_bins(1:end-1) + vel_bins(2:end));
% 
% 
%     tic
%     for i = 1:length(dist_axis)
%         disp(i)
%         id = find(dist>= dist_bin(i) & dist<dist_bin(i+1));
% 
%         pairs_sep(i).dul = dul(id);
%         pairs_sep(i).dut = dut(id);
% 
%         pairs_per_bin(i) = length(id);
%         SF1l(i) = nansum(pairs_sep(i).dul.^1);
%         SF1t(i) = nansum(pairs_sep(i).dut.^1);
% 
%         SF2ll(i) = nansum(pairs_sep(i).dul.^2);
%         SF2tt(i) = nansum(pairs_sep(i).dut.^2);
%         SF2lt(i) = nansum(pairs_sep(i).dut.*pairs_sep(i).dul);
% 
%         SF3lll(i) = nansum(pairs_sep(i).dul.^3);
%         SF3ltt(i) = nansum(pairs_sep(i).dul.*pairs_sep(i).dut.^2);
%         nvaild(i)=sum(~isnan(pairs_sep(i).dul));
%     end
%     toc
%     % save([fname(1:end-3),'chunk',num2str(ii,'%02d'),'SForigin.mat'], ...
%     %     'dist_axis','pairs_sep', ...
%     % '-v7.3');
%     save([fname(1:end-3),'chunk',num2str(ii,'%02d'),'SF123.mat'],'dist_axis','nvaild' ,...
%          'SF1l','SF1t','SF2ll','SF2tt','SF3ltt','SF3lll','time_batch','total_steps','-v7.3');
%     clear pairs_per_bin;
%     clear pairs_sep
%     clear pairs_time
%     clear dist
%     clear dul
%     clear dut
%     clear npairs
% end
% 初始化并行池
% 初始化并行池

gamma = 1.5;
dist_bin = gamma.^(0:100) * 10; % 初始距离箱
id_stop = find(dist_bin > 1000e3, 1); % 找到第一个超过1000km的箱子
if ~isempty(id_stop)
    dist_bin = dist_bin(1:id_stop);
end
dist_bin = [0, dist_bin]; % 在开头插入0
dist_axis = 0.5 * (dist_bin(1:end-1) + dist_bin(2:end));

num_workers = 20; % 例如使用4个工作节点

if isempty(gcp('nocreate'))
    parpool('local', num_workers);
end

% 将外层循环改为 parfor
parfor ii = 1:num_chunks+1
    % 使用函数加载数据以避免 eval 在并行环境中的问题
    traj = loadTrajData(fname, ii);
    
    Htraj = traj.H;
    tic 
    
    % 预分配 pairs_time 数组
    pairs_time = struct('dist', cell(1, length(traj.T_axis)), ...
                       'dul', cell(1, length(traj.T_axis)), ...
                       'dut', cell(1, length(traj.T_axis)));
    
    for i = 1:length(traj.T_axis)
        % 使用逻辑索引替代 find 提高效率
        id = Htraj(i,:) < -500;
        
        if sum(id) > 1
            X = traj.trajmat_X(i,id)';
            Y = traj.trajmat_Y(i,id)';
            U = traj.trajmat_U(i,id)';
            V = traj.trajmat_V(i,id)';
            
            Xvec = [X, Y];
            
            % 计算距离
            pairs_time(i).dist = pdist(Xvec, @dist_geo);
            
            rx = pdist(Xvec, @dist_rx);
            ry = pdist(Xvec, @dist_ry);
            
            magr = sqrt(rx.^2 + ry.^2);
            
            rx = rx./magr; 
            ry = ry./magr;
            
            dux = pdist(U, @dist_du);
            duy = pdist(V, @dist_du);
            
            pairs_time(i).dul = dux.*rx + duy.*ry;
            pairs_time(i).dut = duy.*rx - dux.*ry;
        else
            pairs_time(i).dul = NaN;
            pairs_time(i).dut = NaN;
            pairs_time(i).dist = NaN;
        end
    end
    toc

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %              4. Calc 2-order and 3-order structure function (time-mean)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tpts = length(pairs_time);
    
    % 预分配 npairs 数组
    npairs = zeros(tpts, 1);
    for i = 1:tpts
        if ~isempty(pairs_time(i).dist)
            npairs(i) = length(pairs_time(i).dist);
        end
    end
    
    %% % Align pairs in a single vector
    total_pairs = sum(npairs);
    dul = zeros(total_pairs, 1);
    dut = zeros(total_pairs, 1);
    dist = zeros(total_pairs, 1);
    
    idx = 1;
    for i = 1:tpts
        if npairs(i) > 0
            end_idx = idx + npairs(i) - 1;
            dist(idx:end_idx) = pairs_time(i).dist;
            dul(idx:end_idx) = pairs_time(i).dul;
            dut(idx:end_idx) = pairs_time(i).dut;
            idx = end_idx + 1;
        end
    end
    

    
    %% 计算结构函数
    % 预分配结果数组
    pairs_per_bin = zeros(length(dist_axis), 1);
    SF1l = zeros(length(dist_axis), 1);
    SF1t = zeros(length(dist_axis), 1);
    SF2ll = zeros(length(dist_axis), 1);
    SF2tt = zeros(length(dist_axis), 1);
    SF2lt = zeros(length(dist_axis), 1);
    SF3lll = zeros(length(dist_axis), 1);
    SF3ltt = zeros(length(dist_axis), 1);
    nvaild = zeros(length(dist_axis), 1);
    
    tic
    for i = 1:length(dist_axis)
        % 找出距离在[dist_bin(i), dist_bin(i+1)]之间的点对
        idx_bin = dist >= dist_bin(i) & dist < dist_bin(i+1);
        nvaild(i) = sum(idx_bin);
        
        if nvaild(i) > 0
            dul_bin = dul(idx_bin);
            dut_bin = dut(idx_bin);
            
            SF1l(i) = nansum(dul_bin);
            SF1t(i) = nansum(dut_bin);
            SF2ll(i) = nansum(dul_bin.^2);
            SF2tt(i) = nansum(dut_bin.^2);
            SF2lt(i) = nansum(dut_bin .* dul_bin);
            SF3lll(i) = nansum(dul_bin.^3);
            SF3ltt(i) = nansum(dul_bin .* dut_bin.^2);
        end
    end
    toc
    
    % 使用辅助函数保存结果
    saveResults(fname, ii, dist_axis, nvaild, ...
                SF1l, SF1t, SF2ll, SF2tt, ...
                SF3ltt, SF3lll, time_batch, total_steps);
end
delete(gcp('nocreate'));

%% join
np=nparticles;
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
    eval(['load ',input_dir,Case,'_pars_P',num2str(npart),'T89.5dayschunk',num2str(ii,'%02d'),'SF123.mat']);
    num_chunk=floor(total_steps./time_batch);
    remainder=mod(total_steps,time_batch);

    for ii = 1:num_chunk
        eval(['load ',input_dir,Case,'_pars_P',num2str(npart),'T89.5dayschunk',num2str(ii,'%02d'),'SF123.mat'])

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
        eval(['load ',input_dir,Case,'_pars_P',num2str(npart),'T89.5dayschunk',num2str(ii+1,'%02d'),'SF123.mat'])
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
system('rm *chunk*');

% clear
% load nowave_pars_P289T89.5daysSF123.mat
% SF3_=(SF3lll+SF3ltt)';
% 
% load nowave_pars_P289T89.5daysSF123test.mat
% SF3_test=(SF3lll_time+SF3ltt_time)';
% 
% semilogx(dist_axis,SF3_);hold on
% semilogx(dist_axis,SF3_test);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %              5. plot 2-order and 3-order structure function
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 前面计算的SF123都是同一时刻或者同一dist_bin的累加，还没有取平均
% %这里需要平均下
% 
% % 
% % clear
% np=[289,625,2500,15376];
% colors={'#0072BD','#D95319','#EDB120','#7E2F8E'};
% for jj=1:length(np)
%     npart=np(jj);
%     SF1l_time=0;
%     SF1t_time=0;
%     SF2ll_time=0;
%     SF2tt_time=0;
%     SF3ltt_time=0;
%     SF3lll_time=0;
%     nvaild_time=0;
% 
%     % 集合之前分开计算的chunk
%     ii=1;
%     eval(['load wave_pars_P',num2str(npart),'T89.5dayschunk',num2str(ii,'%02d'),'SF123.mat']);
%     num_chunk=floor(total_steps./time_batch);
%     remainder=mod(total_steps,time_batch);
% 
%     for ii = 1:num_chunk
%         eval(['load wave_pars_P',num2str(npart),'T89.5dayschunk',num2str(ii,'%02d'),'SF123.mat'])
% 
%         SF1t_time=SF1t_time+SF1t;
%         SF1l_time=SF1l_time+SF1l;
%         SF2tt_time=SF2tt_time+SF2tt;
%         SF2ll_time=SF2ll_time+SF2ll;
%         SF3ltt_time=SF3ltt_time+SF3ltt;
%         SF3lll_time=SF3lll_time+SF3lll;
%         nvaild_time=nvaild_time+nvaild;
%         disp(ii);
%     end
% 
%     if remainder>0
%         eval(['load wave_pars_P',num2str(npart),'T89.5dayschunk',num2str(ii+1,'%02d'),'SF123.mat'])
%         % SF3ltt(isnan(SF3ltt))=0;
%         % SF3lll(isnan(SF3lll))=0;
%         SF1t_time=SF1t_time+SF1t;
%         SF1l_time=SF1l_time+SF1l;
%         SF2tt_time=SF2tt_time+SF2tt;
%         SF2ll_time=SF2ll_time+SF2ll;
%         SF3ltt_time=SF3ltt_time+SF3ltt;
%         SF3lll_time=SF3lll_time+SF3lll;
% 
%         nvaild_time=nvaild_time+nvaild;
%     end
% 
%     % 计算平均值
%     SF1t_time(SF1t_time==0)=nan;
%     SF1l_time(SF1l_time==0)=nan;
% 
% 
%     SF2tt_time(SF2tt_time==0)=nan;
%     SF2ll_time(SF2ll_time==0)=nan;
% 
%     SF3ltt_time(SF3ltt_time==0)=nan;
%     SF3lll_time(SF3lll_time==0)=nan;
% 
%     SF1l_time=SF1l_time./nvaild_time;
%     SF1t_time=SF1t_time./nvaild_time;
% 
%     SF2ll_time=SF2ll_time./nvaild_time;
%     SF2tt_time=SF2tt_time./nvaild_time;
% 
%     SF3lll_time=SF3lll_time./nvaild_time;
%     SF3ltt_time=SF3ltt_time./nvaild_time;
% 
% 
%     % semilogx(dist_axis, (SF3lll_time+SF3ltt_time), 'linewidth',2,...
%     % 'Color',colors{jj})
%     % hold on
%     loglog(dist_axis, abs(SF3lll_time+SF3ltt_time), 'linewidth',2, ...
%         'Color',colors{jj})
%     hold on
%     loglog(dist_axis, (SF3lll_time+SF3ltt_time), '+', 'linewidth',2, ...
%         'Color',colors{jj})
%     loglog(dist_axis, (-SF3lll_time-SF3ltt_time), 'o', 'linewidth',2, ...
%         'Color',colors{jj})
%     % text(dist_axis(1), abs(SF3lll_time(1)+SF3ltt_time(1)),num2str(np(jj)))
% end
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %              6. 对5这种分块计算和上个版本不分块计算的验证
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% clear
% clf
% %%%%%%%%%%%%%%%%%%%%%  分块计算的case
% np=[2500];
% colors={'#D95319'};
% for jj=1:length(np)
%     npart=np(jj);
%     SF1l_time=0;
%     SF1t_time=0;
%     SF2ll_time=0;
%     SF2tt_time=0;
%     SF3ltt_time=0;
%     SF3lll_time=0;
%     nvaild_time=0;
% 
%     % 集合之前分开计算的chunk，并计算SF123的平均
%     ii=1;
%     eval(['load wave_pars_P',num2str(npart),'T89.5dayschunk',num2str(ii,'%02d'),'SF123.mat']);
%     num_chunk=floor(total_steps./time_batch);
%     remainder=mod(total_steps,time_batch);
% 
%     for ii = 1:num_chunk
%         eval(['load wave_pars_P',num2str(npart),'T89.5dayschunk',num2str(ii,'%02d'),'SF123.mat'])
% 
%         SF1t_time=SF1t_time+SF1t;
%         SF1l_time=SF1l_time+SF1l;
%         SF2tt_time=SF2tt_time+SF2tt;
%         SF2ll_time=SF2ll_time+SF2ll;
%         SF3ltt_time=SF3ltt_time+SF3ltt;
%         SF3lll_time=SF3lll_time+SF3lll;
%         nvaild_time=nvaild_time+nvaild;
%         disp(ii);
%     end
% 
%     if remainder>0
%         eval(['load wave_pars_P',num2str(npart),'T89.5dayschunk',num2str(ii+1,'%02d'),'SF123.mat'])
%         % SF3ltt(isnan(SF3ltt))=0;
%         % SF3lll(isnan(SF3lll))=0;
%         SF1t_time=SF1t_time+SF1t;
%         SF1l_time=SF1l_time+SF1l;
%         SF2tt_time=SF2tt_time+SF2tt;
%         SF2ll_time=SF2ll_time+SF2ll;
%         SF3ltt_time=SF3ltt_time+SF3ltt;
%         SF3lll_time=SF3lll_time+SF3lll;
% 
%         nvaild_time=nvaild_time+nvaild;
%     end
% 
%     SF1t_time(SF1t_time==0)=nan;
%     SF1l_time(SF1l_time==0)=nan;
% 
% 
%     SF2tt_time(SF2tt_time==0)=nan;
%     SF2ll_time(SF2ll_time==0)=nan;
% 
%     SF3ltt_time(SF3ltt_time==0)=nan;
%     SF3lll_time(SF3lll_time==0)=nan;
% 
%     SF1l_time=SF1l_time./nvaild_time;
%     SF1t_time=SF1t_time./nvaild_time;
% 
%     SF2ll_time=SF2ll_time./nvaild_time;
%     SF2tt_time=SF2tt_time./nvaild_time;
% 
%     SF3lll_time=SF3lll_time./nvaild_time;
%     SF3ltt_time=SF3ltt_time./nvaild_time;
% 
% 
%     a1=semilogx(dist_axis, (SF3lll_time+SF3ltt_time), 'linewidth',0.5,...
%     'Color',colors{jj})
%     hold on
%     % loglog(dist_axis, abs(SF3lll_time+SF3ltt_time)./nvaild_time, 'linewidth',2, ...
%     %     'Color',colors{jj})
%     % hold on
%     % loglog(dist_axis, (SF3lll_time+SF3ltt_time)./nvaild_time, '+', 'linewidth',2, ...
%     %     'Color',colors{jj})
%     % loglog(dist_axis, (-SF3lll_time-SF3ltt_time)./nvaild_time, 'o', 'linewidth',2, ...
%     %     'Color',colors{jj})
% end
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%% 直接计算的case
% load wave_pars_P2500T89.5daysSF123.mat
% a2=semilogx(dist_axis, (SF3lll+SF3ltt), 'linewidth',2,'LineStyle',':');hold on
% legend(['a1','a2'],{'chunk','direct'})
% 
% 
% disp(['error chunk-direct :      ',num2str(SF3lll_time+SF3ltt_time-SF3lll-SF3ltt)]);