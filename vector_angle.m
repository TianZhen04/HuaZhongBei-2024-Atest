function angle = vector_angle(a, b)
    % 检查向量是否为零向量
    if all(a == 0) || all(b == 0)
        error('输入向量不能为零向量');
    end

    % 计算点积
    dot_product = dot(a, b);

    % 计算向量的模
    norm_a = norm(a);
    norm_b = norm(b);

    % 计算夹角的余弦值
    cos_theta = dot_product / (norm_a * norm_b);

    % 确保 cos_theta 在 [-1, 1] 范围内，避免数值误差
    cos_theta = max(min(cos_theta, 1), -1);

    % 计算夹角（弧度）并转换为度
    angle = acos(cos_theta) * 180 / pi;
end
