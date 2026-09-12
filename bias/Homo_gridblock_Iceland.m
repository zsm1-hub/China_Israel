%% ============================================================
% Iceland Lagrangian structure functions
% + grid-index spatial homogeneity
% + optional original block bootstrap
%
% Spatial homogeneity:
%
% H_L^n(r) =
% std_space(local_SF_n) /
% (abs(mean_space(local_SF_n)) + rms_space(local_SF_n))
%
% local_SF_n is first calculated in each spatial block.
%
% Pair distance and direction definitions retain the original:
%   dist_geo.m
%   dist_rx.m
%   dist_ry.m
%   dist_du.m
%% ============================================================

clear;
close all;
clc;

%% ============================================================
% 1. Paths and settings
%% ============================================================

addpath( ...
    '/meddy/simingzhang/Analysis/matlab/Parcels_SF/');

addpath( ...
    '/meddy/simingzhang/Data/Parcels_data');

addpath(genpath( ...
    '/meddy/simingzhang/Data/RB_iceland_data'));

Case = 'wave';             % 'wave' or 'nowave'
nparticles = 289;
days = 89.5;
dt = 3600;                 % seconds

ini = '_roughsmall';
timerange = 1:1940;

% Spatial-block settings
nblock_I = 4;
nblock_J = 4;
min_pairs = 500;

% At least this many valid spatial blocks are required for H
min_valid_blocks = 8;

% true: also perform the original bootstrap calculation
% false: calculate only overall SF and spatial homogeneity
do_bootstrap = true;

num_boot = 1000;

%% ============================================================
% 2. Select input directory
%% ============================================================

if strcmpi(ini, '_grid')

    input_dir = ...
        '/meddy/simingzhang/Data/Parcels_data/' + ...
        "tranV_onetime_spectukey/";

    xscale = [2:18,21:3:48,54:6:114];

elseif strcmpi(ini, '_rough') || ...
       strcmpi(ini, '_rough_1mon') || ...
       strcmpi(ini, '_rough_2mon')

    input_dir = ...
        '/meddy/simingzhang/Data/Parcels_data/' + ...
        "tranV_onetime_roughdistr_tukey/";

    xscale = [2:18,21:3:48,54:6:114];

elseif strcmpi(ini, '_roughsmall') || ...
       strcmpi(ini, '_roughsmall_1mon') || ...
       strcmpi(ini, '_roughsmall_2mon') || ...
       strcmpi(ini, '_roughsmall_3mon')

    input_dir = ...
        '/meddy/simingzhang/Data/Parcels_data/' + ...
        "tranV_onetime_roughsmallregion/";

    xscale = [2:18,21:3:48,54:6:114];

elseif strcmpi(ini, '_roughLASER')

    input_dir = ...
        '/meddy/simingzhang/Data/Parcels_data/' + ...
        "tranV_onetime_roughLASER/";

    xscale = [2:18,21:3:48,54:6:114];

elseif strcmpi(ini, '_roughsmall_rot')

    input_dir = ...
        '/meddy/simingzhang/Data/Parcels_data/' + ...
        "tranV_onetime_roughsmallregion_rot/";

    xscale = [2:18,21:3:48,54:6:114];

elseif strcmpi(ini, '_roughsmall_div')

    input_dir = ...
        '/meddy/simingzhang/Data/Parcels_data/' + ...
        "tranV_onetime_roughsmallregion_div/";

    xscale = [2:18,21:3:48,54:6:114];

elseif strcmpi(ini, '_cruise')

    input_dir = ...
        '/meddy/simingzhang/Data/Parcels_data/' + ...
        "tranV_cruise_roughsmallregion/";

    xscale = [2:18,21:3:48,54:6:114];

elseif strcmpi(ini, '_roughsmall_500m')

    input_dir = ...
        '/meddy/simingzhang/Data/Parcels_data/' + ...
        "tranV_onetime_roughsmallregion_500m/";

    xscale = ...
        [2:18,21:3:48,54:6:114,120:12:228,240:24:336];

elseif strcmpi(ini, '_rough_500m')

    input_dir = ...
        '/meddy/simingzhang/Data/Parcels_data/' + ...
        "tranV_onetime_rough_500m/";

    xscale = ...
        [2:18,21:3:48,54:6:114,120:12:228,240:24:336];

elseif strcmpi(ini, '_roughbox200g_500m')

    input_dir = ...
        '/meddy/simingzhang/Data/Parcels_data/' + ...
        "tranV_onetime_roughbox200g_500m/";

    xscale = ...
        [2:18,21:3:48,54:6:114,120:12:228,240:24:336];

elseif strcmpi(ini, '_roughbox100g_500m')

    input_dir = ...
        '/meddy/simingzhang/Data/Parcels_data/' + ...
        "tranV_onetime_roughbox100g_500m/";

    xscale = ...
        [2:18,21:3:48,54:6:114,120:12:228,240:24:336];

else

    error('Unknown ini setting: %s',ini);
end

input_dir = char(input_dir);

%% ============================================================
% 3. Build trajectory filename
%% ============================================================

if strcmpi(Case,'wave')

    fname = [ ...
        input_dir, ...
        'wave_pars_P',num2str(nparticles), ...
        'T',num2str(days),'days.nc'];

elseif strcmpi(Case,'nowave')

    fname = [ ...
        input_dir, ...
        'nowave_pars_P',num2str(nparticles), ...
        'T',num2str(days),'days.nc'];

else

    error('Case must be wave or nowave.');
end

fprintf('Input file:\n%s\n',fname);

%% ============================================================
% 4. Read trajectories
%% ============================================================

lon_all = ncread(fname,'lon');
lat_all = ncread(fname,'lat');

ue = ncread(fname,'ue');
ve = ncread(fname,'ve');

lon_all = lon_all(timerange,:);
lat_all = lat_all(timerange,:);

ue = ue(timerange,:);
ve = ve(timerange,:);

nPositionTime = size(lon_all,1);
nParticle = size(lon_all,2);

if nParticle ~= nparticles
    warning( ...
        'Expected %d particles but file contains %d.', ...
        nparticles,nParticle);
end

%% ============================================================
% 5. Read th variables without creating th2, th3, ... arrays
%
% This produces the same Th_all but avoids keeping all 38
% large th* arrays in memory.
%% ============================================================

Th_all = nan(length(xscale),length(timerange));

for iii = 1:length(xscale)

    pistr = ['th',num2str(xscale(iii))];

    th_tmp = ncread(fname,pistr);

    Th_all(iii,:) = mean( ...
        th_tmp(timerange,:), ...
        2, ...
        'omitnan');

    clear th_tmp;
end

%% ============================================================
% 6. Calculate Lagrangian velocity
%
% This is retained from the original calculation.
%% ============================================================

nTime = nPositionTime-1;

Utraj = nan(nTime,nParticle);
Vtraj = nan(nTime,nParticle);

for t = 1:nTime

    Utraj(t,:) = ...
        spheric_dist( ...
            lat_all(t,:), ...
            lat_all(t,:), ...
            lon_all(t,:), ...
            lon_all(t+1,:)) ...
        ./dt ...
        .*sign(lon_all(t+1,:)-lon_all(t,:));

    Vtraj(t,:) = ...
        spheric_dist( ...
            lat_all(t+1,:), ...
            lat_all(t,:), ...
            lon_all(t,:), ...
            lon_all(t,:)) ...
        ./dt ...
        .*sign(lat_all(t+1,:)-lat_all(t,:));
end

% Align position time with velocity time
lon = lon_all(1:nTime,:);
lat = lat_all(1:nTime,:);

clear lon_all lat_all;

%% ============================================================
% 7. Read rotated model grid
%% ============================================================

gname = ...
    '/meddy/simingzhang/Data/RB_iceland_data/' + ...
    "niskin2km_500m_grd.nc";

gname = char(gname);

lon_rho = ncread(gname,'lon_rho');
lat_rho = ncread(gname,'lat_rho');

if ~isequal(size(lon_rho),size(lat_rho))
    error('lon_rho and lat_rho have different sizes.');
end

[nI,nJ] = size(lon_rho);

fprintf('Grid size: %d x %d\n',nI,nJ);

%% ============================================================
% 8. Geographic coordinates -> continuous grid indices
%
% 'none' prevents points outside the model grid from being
% forced into an edge block.
%% ============================================================

[Igrid,Jgrid] = ndgrid(1:nI,1:nJ);

valid_grid = ...
    isfinite(lon_rho) & ...
    isfinite(lat_rho);

FI = scatteredInterpolant( ...
    lon_rho(valid_grid), ...
    lat_rho(valid_grid), ...
    Igrid(valid_grid), ...
    'linear', ...
    'none');

FJ = scatteredInterpolant( ...
    lon_rho(valid_grid), ...
    lat_rho(valid_grid), ...
    Jgrid(valid_grid), ...
    'linear', ...
    'none');

%% ============================================================
% 9. Determine particle-covered grid-index region
%
% This avoids dividing the entire 287x287 model domain when
% particles occupy only a smaller release/analysis region.
%
% The range is determined from all valid particle positions
% during the selected time interval.
%% ============================================================

lon_vector = lon(:);
lat_vector = lat(:);

valid_position = ...
    isfinite(lon_vector) & ...
    isfinite(lat_vector);

I_particle_all = FI( ...
    lon_vector(valid_position), ...
    lat_vector(valid_position));

J_particle_all = FJ( ...
    lon_vector(valid_position), ...
    lat_vector(valid_position));

valid_index = ...
    isfinite(I_particle_all) & ...
    isfinite(J_particle_all);

if ~any(valid_index)
    error('Particle positions cannot be mapped to the model grid.');
end

I_min = max(1,floor(min(I_particle_all(valid_index))));
I_max = min(nI,ceil(max(I_particle_all(valid_index))));

J_min = max(1,floor(min(J_particle_all(valid_index))));
J_max = min(nJ,ceil(max(J_particle_all(valid_index))));

clear lon_vector lat_vector;
clear I_particle_all J_particle_all;
clear valid_position valid_index;

I_edges = linspace( ...
    I_min,I_max+eps(I_max),nblock_I+1);

J_edges = linspace( ...
    J_min,J_max+eps(J_max),nblock_J+1);

nBlock = nblock_I*nblock_J;

fprintf( ...
    'Particle-covered I range: %.1f to %.1f\n', ...
    I_min,I_max);

fprintf( ...
    'Particle-covered J range: %.1f to %.1f\n', ...
    J_min,J_max);

fprintf( ...
    'Blocks: %d x %d = %d\n', ...
    nblock_I,nblock_J,nBlock);

%% ============================================================
% 10. Use the original distance-bin construction
%% ============================================================

gamma = 1.3;

is_500m = ...
    strcmpi(ini,'_roughsmall_500m') || ...
    strcmpi(ini,'_rough_500m') || ...
    strcmpi(ini,'_roughbox200g_500m') || ...
    strcmpi(ini,'_roughbox100g_500m');

if is_500m
    first_edge = 1000;
else
    first_edge = 4000;
end

dist_bin = first_edge*gamma.^(0:100);

first_large = find(dist_bin>600e3,1);

if isempty(first_large)
    error('Distance-bin vector does not reach 600 km.');
end

dist_bin = dist_bin(1:first_large-1);

dist_bin(2:end+1) = dist_bin(1:end);
dist_bin(1) = 0;

dist_axis = ...
    0.5*(dist_bin(1:end-1)+dist_bin(2:end));

nScale = length(dist_axis);

%% ============================================================
% 11. Streaming accumulators
%
% H calculation retains only nScale x nBlock arrays.
%% ============================================================

sum_SF1_block = zeros(nScale,nBlock);
sum_SF2_block = zeros(nScale,nBlock);
sum_SF3_block = zeros(nScale,nBlock);

pair_count = zeros(nScale,nBlock);

% Overall structure functions
sum_SF1_total = zeros(nScale,1);
sum_SF2ll_total = zeros(nScale,1);
sum_SF2tt_total = zeros(nScale,1);
sum_SF3lll_total = zeros(nScale,1);
sum_SF3ltt_total = zeros(nScale,1);
count_total = zeros(nScale,1);

%% ============================================================
% 12. Optional pair storage for original bootstrap
%
% Only enabled when do_bootstrap=true.
%% ============================================================

if do_bootstrap

    pairs_sep = repmat( ...
        struct('dul',[],'dut',[]), ...
        nScale,1);
end

%% ============================================================
% 13. Main pair calculation
%
% The original pdist definitions are retained exactly.
%% ============================================================

tic;

for it = 1:nTime

    if mod(it,100)==0
        fprintf('Time %d / %d\n',it,nTime);
    end

    id = find( ...
        isfinite(lon(it,:)) & ...
        isfinite(lat(it,:)) & ...
        isfinite(Utraj(it,:)) & ...
        isfinite(Vtraj(it,:)));

    nCurrent = length(id);

    if nCurrent<2
        continue;
    end

    X = lon(it,id)';
    Y = lat(it,id)';

    U_now = Utraj(it,id)';
    V_now = Vtraj(it,id)';

    Xvec = [X,Y];

    % Original pair calculations
    dist_now = pdist(Xvec,@dist_geo);
    rx = pdist(Xvec,@dist_rx);
    ry = pdist(Xvec,@dist_ry);

    magr = sqrt(rx.^2+ry.^2);

    dux = pdist(U_now,@dist_du);
    duy = pdist(V_now,@dist_du);

    % Convert all pdist outputs to columns
    dist_now = dist_now(:);
    rx = rx(:);
    ry = ry(:);
    magr = magr(:);
    dux = dux(:);
    duy = duy(:);

    % nchoosek order is:
    % (1,2), (1,3), ..., (1,N), (2,3), ...
    % which matches pdist pair order.
    pair_local = nchoosek(1:nCurrent,2);

    pair_a = pair_local(:,1);
    pair_b = pair_local(:,2);

    if size(pair_local,1)~=length(dist_now)
        error('Pair-index order/length does not match pdist.');
    end

    good_pair = ...
        isfinite(dist_now) & ...
        isfinite(rx) & ...
        isfinite(ry) & ...
        isfinite(magr) & ...
        isfinite(dux) & ...
        isfinite(duy) & ...
        magr>0;

    if ~any(good_pair)
        continue;
    end

    dist_now = dist_now(good_pair);

    rx = rx(good_pair)./magr(good_pair);
    ry = ry(good_pair)./magr(good_pair);

    dux = dux(good_pair);
    duy = duy(good_pair);

    pair_a = pair_a(good_pair);
    pair_b = pair_b(good_pair);

    % Same longitudinal/transverse definitions as original
    dul = dux.*rx+duy.*ry;
    dut = duy.*rx-dux.*ry;

    % Pair midpoint in longitude/latitude
    lon_mid = ...
        0.5*(X(pair_a)+X(pair_b));

    lat_mid = ...
        0.5*(Y(pair_a)+Y(pair_b));

    % Map midpoint into continuous rotated-grid indices
    I_mid = FI(lon_mid,lat_mid);
    J_mid = FJ(lon_mid,lat_mid);

    % Separation bin
    scale_id = discretize(dist_now,dist_bin);

    % Spatial block in grid-index space
    block_I = discretize(I_mid,I_edges);
    block_J = discretize(J_mid,J_edges);

    block_id = ...
        (block_I-1)*nblock_J+block_J;

    good = ...
        isfinite(scale_id) & ...
        isfinite(block_id) & ...
        isfinite(dul) & ...
        isfinite(dut) & ...
        scale_id>=1 & ...
        scale_id<=nScale & ...
        block_id>=1 & ...
        block_id<=nBlock;

    if ~any(good)
        continue;
    end

    scale_id = scale_id(good);
    block_id = block_id(good);

    dul = dul(good);
    dut = dut(good);

    %% Samples used by the structure functions

    sf1_sample = dul;
    sf2ll_sample = dul.^2;
    sf2tt_sample = dut.^2;
    sf3lll_sample = dul.^3;
    sf3ltt_sample = dul.*dut.^2;

    % Appendix-C H_L^3 currently uses the longitudinal moment:
    sf3_homogeneity_sample = sf3lll_sample;

    % If you want H for the full 2-D energy-flux combination,
    % replace the previous line with:
    %
    % sf3_homogeneity_sample = ...
    %     sf3lll_sample+sf3ltt_sample;

    %% Overall SF accumulators

    sum_SF1_total = sum_SF1_total+accumarray( ...
        scale_id,sf1_sample,[nScale,1],@sum,0);

    sum_SF2ll_total = sum_SF2ll_total+accumarray( ...
        scale_id,sf2ll_sample,[nScale,1],@sum,0);

    sum_SF2tt_total = sum_SF2tt_total+accumarray( ...
        scale_id,sf2tt_sample,[nScale,1],@sum,0);

    sum_SF3lll_total = sum_SF3lll_total+accumarray( ...
        scale_id,sf3lll_sample,[nScale,1],@sum,0);

    sum_SF3ltt_total = sum_SF3ltt_total+accumarray( ...
        scale_id,sf3ltt_sample,[nScale,1],@sum,0);

    count_total = count_total+accumarray( ...
        scale_id,ones(size(scale_id)), ...
        [nScale,1],@sum,0);

    %% Spatial block accumulators

    linear_id = sub2ind( ...
        [nScale,nBlock], ...
        scale_id,block_id);

    temp = accumarray( ...
        linear_id,sf1_sample, ...
        [nScale*nBlock,1],@sum,0);

    sum_SF1_block = sum_SF1_block+ ...
        reshape(temp,nScale,nBlock);

    temp = accumarray( ...
        linear_id,sf2ll_sample, ...
        [nScale*nBlock,1],@sum,0);

    sum_SF2_block = sum_SF2_block+ ...
        reshape(temp,nScale,nBlock);

    temp = accumarray( ...
        linear_id,sf3_homogeneity_sample, ...
        [nScale*nBlock,1],@sum,0);

    sum_SF3_block = sum_SF3_block+ ...
        reshape(temp,nScale,nBlock);

    temp = accumarray( ...
        linear_id,ones(size(linear_id)), ...
        [nScale*nBlock,1],@sum,0);

    pair_count = pair_count+ ...
        reshape(temp,nScale,nBlock);

    %% Optional storage for original bootstrap

    if do_bootstrap

        for ir = unique(scale_id(:))'

            in_scale = scale_id==ir;

            pairs_sep(ir).dul = [ ...
                pairs_sep(ir).dul; ...
                dul(in_scale)];

            pairs_sep(ir).dut = [ ...
                pairs_sep(ir).dut; ...
                dut(in_scale)];
        end
    end
end

pair_runtime = toc;

fprintf( ...
    'Pair calculation completed in %.1f seconds.\n', ...
    pair_runtime);

%% ============================================================
% 14. Overall SF calculations
%% ============================================================

SF1 = nan(nScale,1);
SF2ll = nan(nScale,1);
SF2tt = nan(nScale,1);
SF3lll = nan(nScale,1);
SF3ltt = nan(nScale,1);

valid_total = count_total>0;

SF1(valid_total) = ...
    sum_SF1_total(valid_total) ./ ...
    count_total(valid_total);

SF2ll(valid_total) = ...
    sum_SF2ll_total(valid_total) ./ ...
    count_total(valid_total);

SF2tt(valid_total) = ...
    sum_SF2tt_total(valid_total) ./ ...
    count_total(valid_total);

SF3lll(valid_total) = ...
    sum_SF3lll_total(valid_total) ./ ...
    count_total(valid_total);

SF3ltt(valid_total) = ...
    sum_SF3ltt_total(valid_total) ./ ...
    count_total(valid_total);

SF2 = SF2ll+SF2tt;

% Full 2-D third-order structure function
SF3 = SF3lll+SF3ltt;

%% ============================================================
% 15. Local structure functions in spatial blocks
%% ============================================================

local_SF1 = nan(nScale,nBlock);
local_SF2 = nan(nScale,nBlock);
local_SF3 = nan(nScale,nBlock);

valid_local = pair_count>=min_pairs;

local_SF1(valid_local) = ...
    sum_SF1_block(valid_local) ./ ...
    pair_count(valid_local);

local_SF2(valid_local) = ...
    sum_SF2_block(valid_local) ./ ...
    pair_count(valid_local);

local_SF3(valid_local) = ...
    sum_SF3_block(valid_local) ./ ...
    pair_count(valid_local);

%% ============================================================
% 16. Spatial homogeneity
%% ============================================================

[H1,mean_SF1,std_SF1,rms_SF1,nvalid1] = ...
    calculate_H( ...
        local_SF1, ...
        pair_count, ...
        min_pairs, ...
        min_valid_blocks);

[H2,mean_SF2,std_SF2,rms_SF2,nvalid2] = ...
    calculate_H( ...
        local_SF2, ...
        pair_count, ...
        min_pairs, ...
        min_valid_blocks);

[H3,mean_SF3,std_SF3,rms_SF3,nvalid3] = ...
    calculate_H( ...
        local_SF3, ...
        pair_count, ...
        min_pairs, ...
        min_valid_blocks);

%% ============================================================
% 17. Optional original bootstrap
%% ============================================================

if do_bootstrap

    Ttot = nTime*dt;

    Tscale_tot = ...
        dist_axis(:) ./ sqrt(SF2ll+SF2tt);

    dof = ones(nScale,1);

    valid_dof = ...
        isfinite(Tscale_tot) & ...
        Tscale_tot>0;

    dof(valid_dof) = ceil( ...
        Ttot ./ Tscale_tot(valid_dof));

    dof = max(dof,1);

    SF1l = nan(nScale,num_boot);
    SF2l = nan(nScale,num_boot);
    SF3l = nan(nScale,num_boot);

    SF3full_boot = nan(nScale,num_boot);

    nsample = zeros(nScale,1);

    tic;

    for ir = 1:nScale

        fprintf( ...
            'Bootstrap scale %d / %d\n', ...
            ir,nScale);

        dul_bin = pairs_sep(ir).dul;
        dut_bin = pairs_sep(ir).dut;

        nsample(ir) = numel(dul_bin);

        if nsample(ir)<=10
            continue;
        end

        n_blocks = min( ...
            dof(ir), ...
            nsample(ir));

        blocksize = floor( ...
            nsample(ir)/n_blocks);

        if blocksize<1
            continue;
        end

        n_use = n_blocks*blocksize;

        % Randomize once before reshaping.
        % This avoids arbitrary dependence on append order.
        random_order = randperm(nsample(ir),n_use);

        dul_use = dul_bin(random_order);
        dut_use = dut_bin(random_order);

        blocks_dul = reshape( ...
            dul_use, ...
            blocksize,n_blocks)';

        blocks_dut = reshape( ...
            dut_use, ...
            blocksize,n_blocks)';

        SF1l_sample = blocks_dul;
        SF2l_sample = blocks_dul.^2;
        SF3l_sample = blocks_dul.^3;

        SF3full_sample = ...
            blocks_dul.^3 + ...
            blocks_dul.*blocks_dut.^2;

        SF1l(ir,:) = bootstrp( ...
            num_boot, ...
            @(x) mean(x(:),'omitnan'), ...
            SF1l_sample);

        SF2l(ir,:) = bootstrp( ...
            num_boot, ...
            @(x) mean(x(:),'omitnan'), ...
            SF2l_sample);

        SF3l(ir,:) = bootstrp( ...
            num_boot, ...
            @(x) mean(x(:),'omitnan'), ...
            SF3l_sample);

        SF3full_boot(ir,:) = bootstrp( ...
            num_boot, ...
            @(x) mean(x(:),'omitnan'), ...
            SF3full_sample);
    end

    bootstrap_runtime = toc;

    fprintf( ...
        'Bootstrap completed in %.1f seconds.\n', ...
        bootstrap_runtime);

    SF3_mean = mean( ...
        SF3full_boot,2,'omitnan');

    SF3_stderr = std( ...
        SF3full_boot,0,2,'omitnan');

else

    dof = [];
    SF1l = [];
    SF2l = [];
    SF3l = [];
    SF3full_boot = [];
    SF3_mean = [];
    SF3_stderr = [];
    nsample = [];
    bootstrap_runtime = NaN;
end

%% ============================================================
% 18. Save result
%% ============================================================

outputname = [ ...
    input_dir, ...
    Case, ...
    '_pars_P', ...
    num2str(nparticles), ...
    'T', ...
    num2str(timerange(end)), ...
    ini, ...
    'gridBlockHL.mat'];

save( ...
    outputname, ...
    'Case', ...
    'nparticles', ...
    'days', ...
    'dt', ...
    'timerange', ...
    'xscale', ...
    'Th_all', ...
    'dist_axis', ...
    'dist_bin', ...
    'SF1', ...
    'SF2', ...
    'SF2ll', ...
    'SF2tt', ...
    'SF3', ...
    'SF3lll', ...
    'SF3ltt', ...
    'count_total', ...
    'local_SF1', ...
    'local_SF2', ...
    'local_SF3', ...
    'H1', ...
    'H2', ...
    'H3', ...
    'mean_SF1', ...
    'mean_SF2', ...
    'mean_SF3', ...
    'std_SF1', ...
    'std_SF2', ...
    'std_SF3', ...
    'rms_SF1', ...
    'rms_SF2', ...
    'rms_SF3', ...
    'nvalid1', ...
    'nvalid2', ...
    'nvalid3', ...
    'pair_count', ...
    'I_edges', ...
    'J_edges', ...
    'I_min', ...
    'I_max', ...
    'J_min', ...
    'J_max', ...
    'nblock_I', ...
    'nblock_J', ...
    'nBlock', ...
    'min_pairs', ...
    'min_valid_blocks', ...
    'do_bootstrap', ...
    'num_boot', ...
    'dof', ...
    'SF1l', ...
    'SF2l', ...
    'SF3l', ...
    'SF3full_boot', ...
    'SF3_mean', ...
    'SF3_stderr', ...
    'nsample', ...
    'pair_runtime', ...
    'bootstrap_runtime');
