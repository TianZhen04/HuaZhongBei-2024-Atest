function [min_solutions, min_values] = monteCarlo(lb, ub, n, m, objective_function)
    % 蒙特卡洛法寻找最小的 m 个解
    % 输入:
    %   lb - 下限向量
    %   ub - 上限向量
    %   n - 蒙特卡洛法抛针数
    %   m - 需要找到的最小的 m 个解
    %   objective_function - 目标函数句柄
    % 输出:
    %   min_solutions - 最小的 m 个解
    %   min_values - 对应的目标函数值

    % 检查输入参数
    if nargin < 5
        error('需要提供所有输入参数: lb, ub, n, m, objective_function');
    end

    % 检查维度是否一致
    if length(lb) ~= length(ub)
        error('lb 和 ub 的维度不一致');
    end

    % 获取维度
    dimen = length(lb);

    % 生成随机点
    points = lb + (ub - lb) .* rand(n, dimen);

    % 计算每个点的目标函数值
    for i = 1:n
        values(i) = objective_function(points(i,:));
    end
   
    % 找到目标函数值最小的 m 个点的索引
    [sorted_values, sorted_indices] = sort(values);
    min_indices = sorted_indices(1:m);

    % 提取最小的 m 个解和对应的目标函数值
    min_solutions = points(min_indices, :);
    min_values = sorted_values(1:m);
end