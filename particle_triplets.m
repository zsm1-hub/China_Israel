clear
% a1=1/4;%0.5km
% a2=1/20;%100m;

a1=1/2;%1km
a2=1/4;%500m;

Nums=900;
grid_min_x = 135;
grid_max_x = 145;
grid_min_y = 135;
grid_max_y = 145;
triplets=particle_setup(grid_min_x,grid_max_x,grid_min_y,grid_max_y,...
    Nums,a1,a2);
px_temp=[triplets(:,1)]';
py_temp=[triplets(:,2)]';
save('triplets.mat','px_temp','py_temp');

data = scipy.io.loadmat('triplets.mat')

# 获取 myData 数据
my_data = data['myData']