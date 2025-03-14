function day_of_year = dayofYear(month, day)
    % 每个月的天数（非闰年）
    days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
    
    % 检查输入是否有效
    if month < 1 || month > 12
        error('月份输入错误');
    end
    if day < 1 || day > days_in_month(month)
        error('日期输入错误');
    end
    
    % 计算第几天
    day_of_year = sum(days_in_month(1:month-1)) + day;
end