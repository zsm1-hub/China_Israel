function triplets=particle_setup(grid_min_x,grid_max_x,grid_min_y,grid_max_y,...
    Nums,a1,a2)
if mod(Nums,9)~=0
    error('you should release nums of 9 beishu')
else
    Nums=Nums/9;
end
disp(['circle radius = ',num2str(a1)]);
disp(['triplets length = ',num2str(a2)]);

% 生成 t 的值，从 0 到 2*pi
t = linspace(0, 2*pi, 1000);

% 定义四叶草的笛卡尔坐标方程
x = cos(2*t) .* cos(t);
y = cos(2*t) .* sin(t);

% 旋转 45 度：应用旋转矩阵
theta = pi / 4;  % 45度的弧度
x_rot = (cos(theta) * x - sin(theta) * y);  % 新的 x 坐标
y_rot = (sin(theta) * x + cos(theta) * y);  % 新的 y 坐标

% 确定四叶草图形的边界
x_min = min(x_rot);
x_max = max(x_rot);
y_min = min(y_rot);
y_max = max(y_rot);

% 网格范围 输入
% grid_min_x = 125;
% grid_max_x = 175;
% grid_min_y = 125;
% grid_max_y = 175;

% 计算缩放因子
scale_x = (grid_max_x - grid_min_x) / (x_max - x_min);  % x 方向的缩放因子
scale_y = (grid_max_y - grid_min_y) / (y_max - y_min);  % y 方向的缩放因子

% 使用不同的缩放因子
scale_factor_x = scale_x;
scale_factor_y = scale_y;

% 计算偏移量，使图形居中
offset_x = (grid_max_x + grid_min_x) / 2 - (x_max + x_min) * scale_factor_x / 2;
offset_y = (grid_max_y + grid_min_y) / 2 - (y_max + y_min) * scale_factor_y / 2;

% 应用缩放和偏移
x_final = scale_factor_x * x_rot + offset_x;
y_final = scale_factor_y * y_rot + offset_y;

t_random = 2 * pi * rand(1, Nums);  % 生成 500 个 0 到 2*pi 之间的随机数

% 计算对应的 x 和 y 坐标
x_random = cos(2*t_random) .* cos(t_random);
y_random = cos(2*t_random) .* sin(t_random);

% 旋转这些随机点
x_random_rot = cos(theta) * x_random - sin(theta) * y_random;
y_random_rot = sin(theta) * x_random + cos(theta) * y_random;

% 缩放和偏移
x_random_final = scale_factor_x * x_random_rot + offset_x;
y_random_final = scale_factor_y * y_random_rot + offset_y;


% 创建图形
figure;
plot(x_final, y_final, 'LineWidth', 2);  % 绘制四叶草曲线
hold on;
plot(x_random_final, y_random_final, 'ro', 'MarkerSize', 5);  % 绘制随机选取的点
% 
% a1=1/4;%%500m dx=2000m
% a1=1;
for ii=1:size(x_random_final,2)
    vertices((ii-1)*3+1:(ii-1)*3+3,:)=equilateral_triangle_on_circle(x_random_final(ii),...
        y_random_final(ii), a1);    
end

% a2=1/10;
for ii=1:size(vertices,1)
    triplets((ii-1)*3+1:(ii-1)*3+3,:) = get_equilateral_triangle_vertices_from_point(vertices(ii,1),...
        vertices(ii,2), a2);
end


plot(vertices(:,1), vertices(:,2), 'bx', 'MarkerSize', 5);  % 绘制随机选取的点
plot(triplets(:,1), triplets(:,2), 'gx', 'MarkerSize', 5);  % 绘制随机选取的点

%%%%%%%%%%%%%%%%%%%%%%%%这么丢的
% scatter(lon(py,px),lat(py,px))
% I=1:size(lon,2);J=1:size(lon,1);