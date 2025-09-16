function [r,SF3,SFL1,SFT1]=calc_radial(SFL1,SFT1,N,xscale)
[X, Y] = meshgrid(-N/2:N/2-1, -N/2:N/2-1);
% x1=X(end/2+1:end,end/2+1:end);
% y1=Y(end/2+1:end,end/2+1:end);

R = sqrt(X.^2 + Y.^2);  % 注意：已修复 R=0 问题

% 定义径向分箱
r_max = max(R(:));
dr = 1;  % 分箱步长（根据需求调整）
r_edges = 0:dr:r_max;
r_centers = (r_edges(1:end-1) + r_edges(2:end))/2;

% 初始化存储
S3L_radial = zeros(size(r_centers));
S3T_radial = zeros(size(r_centers));
counts = zeros(size(r_centers));

% 径向平均
for i = 1:numel(r_centers)
    r_min = r_edges(i);
    r_max = r_edges(i+1);
    
    % 创建当前分箱的掩膜
    mask = (R >= r_min) & (R < r_max);
    
    % 计算平均值
    S3L_radial(i) = mean(S3L1(mask));
    S3T_radial(i) = mean(S3T1(mask));
    counts(i) = sum(mask(:));
end

% 移除空分箱
valid_mask = counts > 0;
r = r_centers(valid_mask);
S3L1 = S3L_radial(valid_mask);
S3T1 = S3T_radial(valid_mask);
SF3=S3L1+S3T1;
r=r.*xscale;
return