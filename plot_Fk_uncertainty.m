function plot_Fk_uncertainty(kf, Fk, uncertainty)
    % 绘制Fk谱及其不确定性分解
    
    figure;
    
    % 主图：Fk谱与总不确定性
    subplot(2,1,1);
    loglog(kf, Fk, 'b-', 'LineWidth', 2); hold on;
    fill([kf; flipud(kf)], ...
         [Fk + 2*uncertainty.total; flipud(Fk - 2*uncertainty.total)], ...
         'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    ylabel('F(k)');
    title('Fk能谱与总不确定性（±2σ）');
    grid on;
    
    % 子图：不确定性分解
    subplot(2,1,2);
    loglog(kf, uncertainty.total, 'k-', 'LineWidth', 2); hold on;
    loglog(kf, uncertainty.epistemic, 'r--', 'LineWidth', 1.5);
    loglog(kf, uncertainty.aleatoric, 'g--', 'LineWidth', 1.5);
    legend('总不确定性', '认知不确定性', '数据不确定性');
    xlabel('波数 k');
    ylabel('不确定性标准差');
    title('Fk不确定性分解');
    grid on;
end