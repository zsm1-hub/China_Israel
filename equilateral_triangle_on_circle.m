function vertices=equilateral_triangle_on_circle(x0, y0, a1)
    % x0, y0: 圆心坐标
    % a1: 圆的半径
    % theta0: 等边三角形第一个顶点的角度（弧度制）
    theta0 = -pi/2 + pi * rand;

    % 计算三角形的三个顶点角度
    angles = theta0 + (0:2) * 2*pi / 3;  % 每个顶点相差120度
    
    % 计算三角形的三个顶点坐标
    vertices = zeros(3, 2);
    for i = 1:3
        vertices(i, 1) = x0 + a1 * cos(angles(i));  % x坐标
        vertices(i, 2) = y0 + a1 * sin(angles(i));  % y坐标
    end
    disp(['radius=',num2str(a1)])
    return