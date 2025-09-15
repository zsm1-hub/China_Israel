function [kr, Er] = calc_ispec2(k, l, E, ndim)
    % 计算方位平均谱 (二维波数谱的径向积分)
    % 输入参数:
    %   k: 波数 x 轴 (1D 数组)
    %   l: 波数 y 轴 (1D 数组)
    %   E: 能谱数据 (2D 或 3D 数组, 维度 [length(l), length(k)] 或 [length(l), length(k), nomg])
    %   ndim: E 的维度 (可选, 默认 2)
    % 输出参数:
    %   kr: 径向波数轴 (1D 数组)
    %   Er: 方位平均后的能谱 (1D 或 2D 数组)
    
    % 处理可选参数 ndim
    if nargin < 4
        ndim = 2;
    end
    
    % 计算波数步长 (假设 k 和 l 是等间距且对称的)
    dk = abs(k(2) - k(1));
    dl = abs(l(2) - l(1));
    
    % 生成波数网格 [K, L]
    [K, L] = meshgrid(k, l);
    wv = sqrt(K.^2 + L.^2); % 计算波数模长
    
    % 确定最大波数 (正半轴范围)
    kmax = max([max(k), max(l)]);
    
    % 分箱设置 (按波数模长分箱)
    % 注意: 此处与原 Python 代码逻辑一致，但可能存在潜在问题
    kr = sqrt(k.^2 + l.^2); % 生成径向波数轴
    dkr = sqrt(dk^2 + dl^2); % 分箱间隔
    
    % 初始化输出数组 (考虑复数)
    if ndim == 3
        nomg = size(E, 3);
    else
        nomg = 1;
    end
    Er = zeros(length(kr), nomg, 'like', complex(1+1i)); % 保持复数类型
    
    % 遍历每个径向波数分箱
    for i = 1:length(kr)
        % 生成掩膜: 选择当前波数区间 [kr(i)-dkr/2, kr(i)+dkr/2]
        mask = (wv >= kr(i) - dkr/2) & (wv <= kr(i) + dkr/2);
        
        % 检查有效数据点
        if sum(mask(:)) > 0
            % 计算角度步长 (近似积分)
            dtheta = pi / (sum(mask(:)) - 1); % 避免除以零
            
            % 计算积分 (考虑不同维度)
            if ndim == 2
                Er(i) = sum(E(mask) .* wv(mask) * dtheta, 'all');
            else
                % 3D 情况: 沿前两维积分，保持第三维
                E_slice = reshape(E(mask, :), [], nomg);
                wv_slice = repmat(wv(mask), 1, nomg);
                Er(i, :) = sum(E_slice .* wv_slice * dtheta, 1);
            end
        end
    end
    
    % 去除冗余维度并取实数部分
    Er = real(squeeze(Er));
end