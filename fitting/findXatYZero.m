function x_zero = findXatYZero(x, y)
% FINDXATYZERO 通过线性插值计算 Y=0 时的 X 值
%   x_zero = findXatYZero(x, y) 返回当 y=0 时的 x 值
%   使用线性插值方法，结果保留一位小数
%
%   输入：
%       x - 一维数据点的 x 坐标
%       y - 一维数据点的 y 坐标
%   输出：
%       x_zero - 当 y=0 时的 x 值（保留一位小数）
%               如果找不到零点，返回 NaN

% 检查输入数据有效性
if numel(x) ~= numel(y)
    error('输入向量 x 和 y 必须具有相同的长度');
end

if numel(x) < 2
    error('至少需要两个数据点进行插值');
end

% 确保数据是列向量
x = x(:);
y = y(:);

% 查找符号变化的区间
sign_changes = find(diff(sign(y)) ~= 0);

% 如果没有找到符号变化
if isempty(sign_changes)
    % 检查是否有点正好在零点
    zero_idx = find(y == 0, 1);
    if ~isempty(zero_idx)
        x_zero = round(x(zero_idx(1)), 1); % 保留一位小数
    else
        warning('未找到零点：数据没有跨越 Y=0');
        x_zero = NaN;
    end
    return;
end

% 对于每个符号变化区间进行插值
x_zeros = zeros(size(sign_changes));
for i = 1:length(sign_changes)
    idx = sign_changes(i);
    
    % 获取相邻两点
    x1 = x(idx);
    x2 = x(idx+1);
    y1 = y(idx);
    y2 = y(idx+1);
    
    % 线性插值公式：x = x1 + (0 - y1)*(x2 - x1)/(y2 - y1)
    x_zeros(i) = x1 - y1 * (x2 - x1) / (y2 - y1);
end

% 选择第一个零点（或根据需求选择其他）
x_zero = x_zeros(1);

% 四舍五入保留一位小数
x_zero = round(x_zero, 1);
end