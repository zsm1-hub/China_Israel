function [po_idx, fac0_idx] = find_L_curve_corner(residual_norms, solution_norms, n_po, n_fac0)
    % 修正后的拐点检测函数
    % 输入:
    %   residual_norms: 残差范数矩阵 (n_po × n_fac0)
    %   solution_norms: 解范数矩阵 (n_po × n_fac0)
    %   n_po: po参数的数量
    %   n_fac0: fac0参数的数量
    % 输出:
    %   po_idx: 最优po参数的索引
    %   fac0_idx: 最优fac0参数的索引
    
    % 转换为对数尺度
    log_res = log10(residual_norms(:));
    log_sol = log10(solution_norms(:));
    
    % 移除无穷大值
    valid_idx = isfinite(log_res) & isfinite(log_sol);
    log_res = log_res(valid_idx);
    log_sol = log_sol(valid_idx);
    
    if length(log_res) < 3
        % 如果没有足够的数据点，使用默认值
        po_idx = round(n_po/2);
        fac0_idx = round(n_fac0/2);
        fprintf('警告：有效数据点不足，使用默认索引\n');
        return;
    end
    
    % 计算曲率
    curvature = zeros(size(log_res));
    for i = 2:length(log_res)-1
        % 计算曲率的近似值
        dx1 = log_res(i) - log_res(i-1);
        dy1 = log_sol(i) - log_sol(i-1);
        dx2 = log_res(i+1) - log_res(i);
        dy2 = log_sol(i+1) - log_sol(i);
        
        denominator = (dx1^2 + dy1^2)^1.5;
        if denominator > 0
            curvature(i) = abs(dx1*dy2 - dx2*dy1) / denominator;
        else
            curvature(i) = 0;
        end
    end
    
    % 找到最大曲率点
    [~, max_idx] = max(curvature);
    
    % 将一维索引转换为二维索引
    % 注意：MATLAB使用列优先存储
    po_idx = mod(max_idx-1, n_po) + 1;
    fac0_idx = floor((max_idx-1) / n_po) + 1;
    
    % 确保索引在有效范围内
    po_idx = max(1, min(n_po, po_idx));
    fac0_idx = max(1, min(n_fac0, fac0_idx));
    
    fprintf('曲率最大点: 一维索引=%d, 二维索引=(%d, %d)\n', max_idx, po_idx, fac0_idx);
end