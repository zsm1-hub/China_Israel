function [optimal_fac, optimal_po, lcurve_results] = L_curve_analysis(Y, W, A, kbins, fac0_range, po_range)
    % 修正后的L曲线分析函数
    % 输入:
    %   Y: 观测数据
    %   W: 权重矩阵
    %   A: 设计矩阵
    %   kbins: 波数分箱
    %   fac0_range: ε_u先验方差范围
    %   po_range: ε_j先验方差范围
    % 输出:
    %   optimal_fac: 最优的ε_u先验方差
    %   optimal_po: 最优的ε_j先验方差
    %   lcurve_results: 分析结果

    nk = length(kbins);
    n_po = length(po_range);
    n_fac0 = length(fac0_range);

    % 初始化存储矩阵
    norm2_eps = zeros(n_po, n_fac0);
    norm2_epsu = zeros(n_po, n_fac0);
    norm2_res = zeros(n_po, n_fac0);
    norm2_err = zeros(n_po, n_fac0);
    norm2_erru = zeros(n_po, n_fac0);
    norm2_P0 = zeros(n_po, n_fac0);

    % 提取权重
    wgts = sqrt(diag(W));

    % 显示进度
    fprintf('执行L曲线分析 (%d x %d = %d 次反演)...\n', n_po, n_fac0, n_po*n_fac0);

    % 双重循环测试所有参数组合
    for n = 1:n_po
        for m = 1:n_fac0
            % 构建P矩阵
            P_diag = [fac0_range(m), ones(1, nk) * po_range(n)];
            P_test = diag(P_diag);

            try
                % 执行RLS反演
                [eps_test, sf3_test, res_test, cxx_test] = RLS(Y, W, P_test, A);
                eps_err_test = sqrt(diag(cxx_test));

                % 计算各种范数
                norm2_eps(n, m) = norm(eps_test(2:end), 2);      % ε_j的L2范数
                norm2_epsu(n, m) = eps_test(1)^2;                % ε_u的平方
                norm2_res(n, m) = norm(res_test ./ wgts, 2);     % 加权残差范数
                norm2_err(n, m) = norm(eps_err_test(2:end), 2);  % ε_j误差的L2范数
                norm2_erru(n, m) = eps_err_test(1)^2;            % ε_u误差的平方
                norm2_P0(n, m) = norm(sqrt(P_diag(2:end)), 2);  % 先验的L2范数
            catch ME
                % 如果反演失败，设置默认值
                fprintf('反演失败: po=%.2e, fac0=%.2e\n', po_range(n), fac0_range(m));
                norm2_eps(n, m) = Inf;
                norm2_epsu(n, m) = Inf;
                norm2_res(n, m) = Inf;
                norm2_err(n, m) = Inf;
                norm2_erru(n, m) = Inf;
                norm2_P0(n, m) = Inf;
            end
        end
        fprintf('完成 %d/%d\n', n, n_po);
    end

    % 找到最优参数
    [optimal_po_idx, optimal_fac0_idx] = find_L_curve_corner(norm2_res, norm2_eps, n_po, n_fac0);

    % 检查索引是否在有效范围内
    if optimal_po_idx < 1 || optimal_po_idx > n_po || optimal_fac0_idx < 1 || optimal_fac0_idx > n_fac0
        fprintf('警告：自动选择的索引超出范围，使用默认值\n');
        optimal_po_idx = round(n_po/2);
        optimal_fac0_idx = round(n_fac0/2);
    end

    optimal_fac = fac0_range(optimal_fac0_idx);
    optimal_po = po_range(optimal_po_idx);

    % 打包结果
    lcurve_results.norm2_eps = norm2_eps;
    lcurve_results.norm2_epsu = norm2_epsu;
    lcurve_results.norm2_res = norm2_res;
    lcurve_results.norm2_err = norm2_err;
    lcurve_results.norm2_erru = norm2_erru;
    lcurve_results.norm2_P0 = norm2_P0;
    lcurve_results.fac0_range = fac0_range;
    lcurve_results.po_range = po_range;
    lcurve_results.optimal_po_idx = optimal_po_idx;
    lcurve_results.optimal_fac0_idx = optimal_fac0_idx;

    fprintf('最优参数: po_idx=%d, fac0_idx=%d\n', optimal_po_idx, optimal_fac0_idx);
    fprintf('最优值: po=%.2e, fac0=%.2e\n', optimal_po, optimal_fac);
end
