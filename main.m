% 读取 Excel 文件
data1 = readtable('data1.xlsx','Range',[1 1 28 3]); % 替换为你的文件名
data2 = readtable('data2.xlsx');
save('data.mat','data2','data1');

% 定义地理信息
geo_info.latitude = 367/12; % 纬度，单位：度
geo_info.longitude = 6859/60; % 经度，单位：度
geo_info.year = 2023;
geo_info.month = 5;
geo_info.day = 23;

% 定义太阳能板信息
solar_panel.azimuth = 0; % 太阳能板方位角，单位：度（正南为0度）
solar_panel.tilt = 30; % 太阳能板水平仰角，单位：度
solar_panel.area = 1; %太阳能板面积，单位：平方米

% 定义太阳辐射值
% I1,I2,I3;

% [Sun_altitude_angle, Sun_azimuth_deg]=solar_angles 计算给定时间和地点的太阳高度角和方位角
% 输入参数：
%   year: 计算日期的年份
%   month: 计算日期的月份
%   day: 计算日期的天
%   hour: 计算时间的小时
%   minute: 计算时间的分钟
%   lon: 观测地点的经度（东经为正）
%   lat: 观测地点的纬度（北纬为正）
% 输出参数：
%   Sun_altitude_angle: 太阳高度角（单位：度）
%   Sun_azimuth_deg: 太阳方位角（单位：度）

I0 = 1334;
depth=1000;

decline = zeros(27,1); %衰减系数
Sun_altitude_angle = zeros(27,1);
l = zeros(27,1); %距离
time = zeros(27,1);
I1 = data1.I1;
for i = 1 : size(data1,1)
    [Sun_altitude_angle(i), ~] = solar_angles(geo_info.year,geo_info.month, geo_info.day,data1.hour(i),data1.minute(i), geo_info.longitude,geo_info.latitude);
    l(i) = depth/sin(Sun_altitude_angle(i)/180*pi);
    % decline(i) = (I0 -  data1.I1(i));
    % time(i) = data1.hour(i) + data1.minute(i)/60;
end
l = [0;l];
I1 = [I0;I1];

%对l和I1进行指数拟合，I(l)=I0*exp(-5.7714*e-4*l)

%1334 * exp(b * l)
a = @(x) -loss(x,l,I1);
options = optimoptions('particleswarm', ...
    'Display', 'iter'); % 显示迭代信息
b = particleswarm(a,1,-0.1,0,options);



%%
%2025年每月15日，在晴天条件下，该城区一块面积为1m2的光伏板朝向正南方且水平倾角分别为20、40、60时受到的最大太阳直射强度和太阳直射辐射总能量
%对于求最大太阳辐射，决策变量：小时与分钟
tilt = [20 40 60];
max_I = zeros(12,3);
sum_I = zeros(12,3);
for i = 1:3
    for j = 1:12
        a = @(x) -cal_I(j,15,x(1),x(2),0,tilt(i));
        % 设置初始粒子位置
        initial_positions = [
            10, 30; % 粒子1的初始位置
            15, 45; % 粒子2的初始位置
            5, 20;  % 粒子3的初始位置
            % 添加更多粒子...
        ];

        % 获取粒子数量和维度
        n_particles = size(initial_positions, 1);
        n_dimensions = size(initial_positions, 2);

        % 设置粒子群算法选项
        options = optimoptions('particleswarm', ...
            'InitialSwarmMatrix', initial_positions, ... % 设置初始粒子位置
            'SwarmSize', n_particles, ... % 粒子数量
            'Display', 'iter'); % 显示迭代信息
        x = particleswarm(a,2,[0,0],[23,59],options);
        max_I(j,i) = -a(x);
    end
end
k = 1:12;
max_I = [k' max_I];

fileID = fopen('1.txt','w');
fprintf(fileID,'%s\n','每个月15日最大太阳直射强度和太阳直射辐射总能量');
fprintf(fileID,'%18s\n','maxI');
fprintf(fileID,'%6s %6s %6s %6s\n','month','20度','40度','60度');
fprintf(fileID,'%7.0f %7.2f %7.2f %7.2f\n',max_I');

%太阳直射辐射总能量
%以dt为步长
mark = 'title20title40title60';
dt = 30;azimuth = 0;
for i = 1:3
    for j = 1:12
        for time = 0:dt:24*60-dt
            hour = mod(time,60);
            minute = time - 60*hour;
            dI.(mark(7*i-6:7*i))(j,time/dt+1) = cal_I(j,15,hour,minute,azimuth,tilt(i));
        end
        sum_I(j,i) = sum(dI.(mark(7*i-6:7*i))(j,:))*dt/60;
    end
end
k = 1:12;
sum_I = [k' sum_I];

fprintf(fileID,'\n%18s\n','sumI');
fprintf(fileID,'%6s %6s %6s %6s\n','month','20度','40度','60度');
fprintf(fileID,'%7.0f %7.2f %7.2f %7.2f\n',sum_I');
fclose(fileID);



%%
%如果光伏板受到的太阳直射辐射总能量最大时，可使路灯蓄电池储电量最大。
%请设计该城区固定安装太阳能光伏板的朝向，使光伏板在晴天条件下受到的太阳直射辐射日均总能量最大
tic
a = @(x) -calBestBand(x(1),x(2));
% % 蒙特卡洛法找初始位置（最小）
% n = 100;
% dimen = 2;
% lb = [0 0];
% ub = [359 90];
% m = 3;
% [initial_positions, min_values] = monteCarlo(lb, ub, n, m, a);

% % 获取粒子数量和维度
% n_particles = size(initial_positions, 1);
% n_dimensions = size(initial_positions, 2);

% 设置粒子群算法选项
options = optimoptions('particleswarm', ...
    'Display', 'iter'); % 显示迭代信息
x = particleswarm(a,2,[-90,0],[90,90],options);

toc

% 提示用户输入一个数值
num = input('是否覆盖原数值,1是0否');
if num == 1
    fileID = fopen('2.txt','w');
    fprintf(fileID,'%s\n','最佳太阳能光伏板朝向');
    fprintf(fileID,'%7s %9.4f %s\n','方位角：',x(1),'度');
    fprintf(fileID,'%7s %9.4f %s\n','水平仰角：',x(2),'度');
    fprintf(fileID,'%7s %9.4f %s\n','最大日均总能量：',-a(x)/365,'Wh');
    fclose(fileID);
end

%%
%综合考虑路灯蓄电池的储电效率高和储电量大这两个目标，请设计出光伏板固定安装的最优朝向
%并计算晴天条件下光伏板受到的太阳直射辐射日均总能量和太阳直射辐射（上午大于 150 W/m2、下午大于 100W/m2）时长
%日均总能量即第二问
%多目标优化问题，使用MOPSO优化算法

% Multi-objective function
f =  @(x,y) -calBestBand(x,y);
g =  @(x,y) -calEffectiveTime(x,y);
MultiObj.fun = @(x) [f(x(:,1),x(:,2)),g(x(:,1),x(:,2))];
MultiObj.nVar = 2;
MultiObj.var_min = [-90,0];
MultiObj.var_max = [90,90];
% MultiObj.truePF = PF;

% Parameters
params.Np = 100;        % 种群规模
params.Nr = 100;        % 存储库规模
params.maxgen = 100;    % 最大代数
params.W = 0.4;         % 惯性权重
params.C1 = 2;          % 个体学习因子
params.C2 = 2;          % 群体学习因子
params.ngrid = 20;      % 每个维度的网格数
params.maxvel = 5;      % 最大速度（百分比）
params.u_mut = 0.5;     % 均匀变异（百分比）

% MOPSO
REP = MOPSO(params,MultiObj);

% 提示用户输入一个数值
num = input('是否覆盖原数值,1是0否');
if num == 1
    save('3.mat','REP');
end

h_fig = figure(1);
h_rep = plot(-REP.pos_fit(:,1)/365,-REP.pos_fit(:,2)/60/365,'ok'); hold on;
%%
load('3.mat')
pos = REP.pos;
pos_fit = REP.pos_fit;

data = pos;
% 计算每个点的密度
[X, Y] = meshgrid(linspace(min(data(:,1)), max(data(:,1)), 100), ...
                  linspace(min(data(:,2)), max(data(:,2)), 100));
pos = [X(:), Y(:)];
kde = ksdensity(data, pos);
density = reshape(kde, size(X));  % KDE 密度矩阵

% 将密度矩阵插值到原始点的位置
density_values = interp2(X, Y, density, data(:,1), data(:,2), 'linear');

% 绘制点并根据密度上色
figure;
scatter(data(:,1), data(:,2), 36, density_values, 'filled');  % 使用密度值作为颜色
colorbar;  % 显示颜色条
xlabel('X');
ylabel('Y');
title('点的分布及密度');

% 使用 kmeans 进行聚类
num_clusters = 3; % 簇的数量
[idx, C] = kmeans(data, num_clusters); % idx 是簇标签，C 是簇中心

% 绘制聚类结果
figure;
scatter(data(:,1), data(:,2), 36, idx, 'filled'); % 根据簇标签上色
hold on;
scatter(C(:,1), C(:,2), 100, 'kx', 'LineWidth', 2); % 绘制簇中心
xlabel('X');
ylabel('Y');
title('kmeans 聚类结果');
grid on;

% 计算每个簇的点数
cluster_counts = histcounts(idx, [1:num_clusters+1]);

% 找到点数最多的簇
[~, densest_cluster_idx] = max(cluster_counts);

% 获取最密集区域的中心
densest_center = C(densest_cluster_idx, :);

% 在图中标记最密集区域的中心
scatter(densest_center(1), densest_center(2), 100, 'r+', 'LineWidth', 2);
legend('数据点', '簇中心', '最密集区域中心');





%%
function R_squared = loss(b,l,I1)
    I2 = 1334 * exp(b * l);
    SST = sum((I1 - mean(I1)).^2); % 总体平方和
    SSE = sum((I1 - I2).^2); % 残差平方和
    R_squared = 1 - (SSE / SST);

end

function sum_I = calBestBand(azimuth,title)
    dayForMonths = [31,28,31,30,31,30,...
        31,31,30,31,30,31];
    dt = 60;sum_I = 0;
    for month = 1:12
        for day = 1:dayForMonths(month)
            for time = 0:dt:24*60-dt
                hour = mod(time,60);
                minute = time - 60*hour;
                dI = cal_I(month,day,hour,minute,azimuth,title) * dt/60;
                sum_I = sum_I + dI;
            end
        end
    end
end

function effectiveTime = calEffectiveTime(azimuth,title)
    dayForMonths = [31,28,31,30,31,30,...
        31,31,30,31,30,31];
    dt = 60;
    effectiveTime=0;
    for month = 1:12
        for day = 1:dayForMonths(month)
            for time = 0:dt:24*60-dt
                hour = mod(time,60);
                minute = time - 60*hour;
                dI = cal_I(month,day,hour,minute,azimuth,title);
                % sum_I = sum_I + dI * dt/60;
                if (dI>150 && hour<12)||(dI>100 && hour>=12)
                    effectiveTime=effectiveTime+dt;
                end
            end
        end
    end
end
