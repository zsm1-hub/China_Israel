% 辅助函数：计算L曲线的曲率并找到拐点
function [k, idx_max] = max_curvature(log_res, log_sol)
    % 确保输入是行向量
    if size(log_res, 1) > size(log_res, 2)
        log_res = log_res';
    end
    if size(log_sol, 1) > size(log_sol, 2)
        log_sol = log_sol';
    end
    
    % 计算一阶导数
    dlog_res = gradient(log_res);
    dlog_sol = gradient(log_sol);
    
    % 计算二阶导数
    d2log_res = gradient(dlog_res);
    d2log_sol = gradient(dlog_sol);
    
    % 计算曲率
    numerator = abs(d2log_sol .* dlog_res - d2log_res .* dlog_sol);
    denominator = (dlog_res.^2 + dlog_sol.^2).^(3/2);
    k = numerator ./ denominator;
    
    % 找到最大曲率点
    [~, idx_max] = max(k);
    
    % 处理特殊情况
    if isempty(idx_max) || isnan(idx_max)
        idx_max = round(length(log_res)/2);
        warning('Could not find curvature maximum, using midpoint.');
    end
end