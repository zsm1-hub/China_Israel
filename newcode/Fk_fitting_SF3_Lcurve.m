function [SpecFlux, Vt, ebs, kf, lf, lambda_opt,kfo] = Fk_fitting_SF3_Lcurve(SF3, ...
    dist_axis, mindist, maxdist, kftype, inv_style, lambda_vec, kf1,plt)
    % Fk_fitting_SF3_Lcurve - 计算谱通量和能量注入密度，仅绘制L曲线
    %
    % 输入参数:
    %   SF3: 结构函数数据 (r × time)
    %   dist_axis: 距离向量 (1 × r)
    %   mindist: 最小距离 (m)
    %   maxdist: 最大距离 (m)
    %   kftype: 波数类型 ('log' 或 'linear')
    %   inv_style: 反演方法 ('NNLS', 'LS', 'RLS')
    %   lambda_vec: 正则化参数向量 (用于RLS方法)
    %   kf1: 自定义波数向量 (当kftype不为'log'时使用)
    %
    % 输出参数:
    %   SpecFlux: 谱通量 (m²/s³)
    %   Vt: 重建的结构函数
    %   ebs: 能量注入密度 (m³/s³)
    %   kf: 波数向量 (m⁻¹)
    %   lf: 距离向量 (m)
    %   lambda_opt: 最优正则化参数
    
    % 显示输入参数信息
    nsamps = size(SF3, 2);
    r = dist_axis;
    Nr = length(r);
    
    % 选择合理的距离范围
    ns = find(dist_axis >= mindist, 1);
    ne = find(dist_axis <= maxdist, 1, 'last');
    
    if isempty(ns) || isempty(ne)
        error('No valid points found in the specified distance range.');
    end
    
    R = r(ns:ne); 
    NR = length(R);
    lf = R;
    
    % 创建波数向量
    if strcmp(kftype, 'log')
        % kf = logspace(log10(1/max(R)), log10(1/min(R)), 8);
        kf= logspace(log10(1/max(R)), log10(1/min(R)), 11).*2.*pi;

    else
        kf = kf1;
    end
    % kf=kf(1:28:end);
    kfo=kf;
    dk = diff(kf);
    kf = 0.5*(kf(1:end-1) + kf(2:end));
    kf = fliplr(kf); 
    dk = fliplr(dk);
    
    Nk = length(kf);
    
    % 初始化输出矩阵
    ebs = zeros(Nk+1, nsamps);
    S = zeros(Nr, nsamps);
    Vt = zeros(NR, nsamps);
    SpecFlux = zeros(Nk, nsamps); 
    
    % 初始化最优λ
    lambda_opt = [];
    
    % 主循环处理每个样本
    for n = 1:nsamps
        % 获取当前样本的SF3数据
        S(:, n) = SF3(:, n);
        V = S(ns:ne, n)';
        
        % 构建矩阵A
        A = zeros(NR, Nk+1);
        for j = 1:Nk
            A(:, j) = -4/kf(j) * besselj(1, kf(j)*R)' * dk(j);
        end
        A(:, end) = 2*R';
        
        % 归一化处理
        norm_flag = 1;
        if norm_flag == 1
            for j = 1:NR
                A(j, :) = A(j, :) ./ abs(R(j));
            end
            V = V ./ abs(R);
        end
        
        % 根据反演方法求解
        if strcmp(inv_style, 'NNLS')
            ebs(:, n) = lsqnonneg(A, V');
        elseif strcmp(inv_style, 'LS')
            ebs(:, n) = A \ V';
        elseif strcmp(inv_style, 'RLS')
            % ===== L曲线分析 =====
            if n == 1 % 只对第一个样本进行L曲线分析
                % 初始化存储
                residual_norms = zeros(size(lambda_vec));
                solution_norms = zeros(size(lambda_vec));
                
                % 对每个λ进行计算
                for idx = 1:length(lambda_vec)
                    lambda_current = lambda_vec(idx);
                    n_cols = size(A, 2);
                    A_aug = [A; sqrt(lambda_current) * eye(n_cols)];
                    V_aug = [V'; zeros(n_cols, 1)];
                    x = A_aug \ V_aug;
                    
                    % 计算残差范数
                    residual = A * x - V';
                    residual_norms(idx) = norm(residual);
                    
                    % 计算解范数
                    solution_norms(idx) = norm(x);
                end
                
                % 找到L曲线拐点（曲率最大点）
                [~, lambda_opt_idx] = max_curvature(log10(residual_norms), log10(solution_norms));
                lambda_opt = lambda_vec(lambda_opt_idx);

                if plt==1
                    % 绘制L曲线
                    figure(20);
                    loglog(residual_norms, solution_norms, '-o', 'LineWidth', 1.5);
                    hold on;
                    plot(residual_norms(lambda_opt_idx), solution_norms(lambda_opt_idx), ...
                        'ro', 'MarkerSize', 10, 'LineWidth', 2);
                    xlabel('Residual Norm $\|A\mathbf{x} - \mathbf{b}\|_2$', ...
                        'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
                    ylabel('Solution Norm $\|\mathbf{x}\|_2$', ...
                        'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
                    title(sprintf('L-curve for Regularization Parameter Selection'), ...
                        'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
                    grid on;
                    legend('L-curve', 'Optimal \lambda', 'Location', 'best');
                    % 只显示首、尾、拐点附近的λ值
                    show_idx = [1, max(1, lambda_opt_idx-2), lambda_opt_idx, min(length(lambda_vec), lambda_opt_idx+2), length(lambda_vec)];
                    show_idx = unique(show_idx);
                    
                    for i = show_idx
                        text(residual_norms(i) * 1.1, solution_norms(i) * 0.9, ...
                            sprintf('λ=%.1e', lambda_vec(i)), ...
                            'FontSize', 10, 'HorizontalAlignment', 'left');
                    end
                    % 其他设置...
                end

            %     if plt==1
            %         % 绘制L曲线
            %         figure(20);
            %         loglog(residual_norms, solution_norms, '-o', 'LineWidth', 1.5);
            %         hold on;
            %         plot(residual_norms(lambda_opt_idx), solution_norms(lambda_opt_idx), ...
            %             'ro', 'MarkerSize', 10, 'LineWidth', 2);
            % 
            %         % 添加标签和标题
            %         xlabel('Residual Norm $\|A\mathbf{x} - \mathbf{b}\|_2$', ...
            %             'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
            %         ylabel('Solution Norm $\|\mathbf{x}\|_2$', ...
            %             'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
            %         title(sprintf('L-curve for Regularization Parameter Selection'), ...
            %             'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
            %         grid on;
            %         legend('L-curve', 'Optimal \lambda', 'Location', 'best');
            % 
            % 
            %         % 添加λ值标注
            %         for i = 1:length(lambda_vec)
            %             text(residual_norms(i), solution_norms(i), ...
            %                 sprintf('λ=%.1e', lambda_vec(i)), ...
            %                 'FontSize', 10, 'HorizontalAlignment', 'left');
            %         end
            %     end
            %     % set(gca, 'FontSize', 14, 'FontWeight', 'bold');
            % end
            % 使用最优λ求解
            n_cols = size(A, 2);
            A_aug = [A; sqrt(lambda_opt) * eye(n_cols)];
            V_aug = [V'; zeros(n_cols, 1)];
            ebs(:, n) = A_aug \ V_aug;
        end
        
        % 重建V
        Vt(:, n) = 2 * ebs(end, n) * R;
        for j = 1:Nk
            Vt(:, n) = Vt(:, n) - 4 * ebs(j, n) ./ kf(j) .* besselj(1, kf(j)*R)' * dk(j);
        end
        
        % 计算谱通量
        SpecFlux(1, n) = -ebs(end, n); 
        for j = 2:Nk
            SpecFlux(j, n) = SpecFlux(j-1, n) + ebs(end-j+1, n) * dk(end-j+2);
        end
    end
    
    % 翻转谱通量
    SpecFlux = flipud(SpecFlux);
end