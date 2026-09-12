function [x_fill,y_fill]=get_shadow(x,y1,y2)
[N,M]=size(x);
[N1,M]=size(y1);

if N>1
    error('x should be a column vector!');
elseif N1>1
    error('y should be a column vector!');
end
x_fill = [x, fliplr(x)];
y_fill=[y1,fliplr(y2)];
return