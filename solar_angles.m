function [Sun_altitude_angle, Sun_azimuth_deg] = solar_angles(year, month, day, hour, minute, lon, lat)
    % 设置经纬度和当地日期(下面时间参数支持向量运算)
    % lat   = 30;       % 北纬(单位：度)
    % lon   = 120;      % 东经(单位：度)
    % year  = 2024;     % 年份(单位：年)
    % month = 6;        % 月份(单位：月)
    % day   = 1;        % 日期(单位：天)
    % hour  = 12;       % 小时(单位：24h)
    % minute= 0;        % 分钟(单位：分钟)
     
    % 输出结果规定
    % 太阳高度角输出结果(单位：度)
    % 地面0 --> 90天空
     
    % 太阳方位角输出结果(单位：度)
    %        南
    %        ^
    %       180
    % 东<--90 270-->西
    %        0
    %        v
    %        北
     
     
    % 计算今年的累积天数sumdays (单位：天)
    date1_num = datenum(year, 1, 1);       % 第一个日期的datenum
    date2_num = datenum(year, month, day); % 第二个日期的datenum
    % 计算两个日期之间的差值
    sumdays   = date2_num-date1_num;
     
     
    % 计算太阳赤纬角sun_declinatio(单位：度)
    %lat_rad为所在地的维度（弧度制）
    lat_rad = deg2rad(lat);
    %Bourges太阳赤纬角算法
    n0=78.801+(0.2422*(year-1969))-round(0.25*(year-1969));
    %b为日角(单位：幅度)
    b=2*pi*(double(sumdays)-n0-1)/365.2422;
    % sun_declinatio为太阳赤纬角（角度制）
    sun_declinatio =0.3723...
        +23.2567*sin(b)...
        +0.1149* sin(2.0*b)...
        -0.1712* sin(3.0*b)...
        -0.758*  cos(b)...
        +0.3656* cos(2.0*b)...
        +0.02010*cos(3.0*b);
     
    %转为（弧度制）
    sun_declinatio_rad =deg2rad(sun_declinatio);
     
     
     
    % 计算太阳太阳时角ts (单位：幅度)
    %Ts真太阳时与平太阳时的时差
    Ts=0.0028-1.9587.*sin(b)+9.9059.*sin(2.*b)-7.0924.*cos(b)-0.6882.*cos(2.*b);
    %Sd为平太阳时，为当地时间(注意这里的当地时间修正)
    Sd = hour  + (minute- (120.0 - lon) * 4.0 ) / 60.0;
    %st为真太阳时=真太阳时时差+平太阳
    st=Sd+Ts/60;
    %ts为所在地的太阳时角
    ts = (st - 12.0) * pi / 12.0;
     
    % 计算太阳高度角Sun_altitude_ra(单位：度)
    Sun_altitude_rad = asin(sin(sun_declinatio_rad) .* sin(lat_rad) + cos(sun_declinatio_rad) .* cos(lat_rad) .* cos(ts));
    %转为度
    Sun_altitude_angle=rad2deg(Sun_altitude_rad);
     
     
     
    % 计算太阳方位角Sun_azimuth_deg(单位：度)
    Sun_azimuth_rad = acos((sin(Sun_altitude_rad) .* sin(lat_rad) - sin(sun_declinatio_rad)) ./ (cos(Sun_altitude_rad) .* cos(lat_rad)));
    %转为度
    Sun_azimuth_deg=rad2deg(Sun_azimuth_rad);
    % 此时方位角是以南为零方位角，转换为以北为零方位角
    Sun_azimuth_deg=180-Sun_azimuth_deg;
    %进行角度范围转换
    if ts > 0
        Sun_azimuth_deg =360-Sun_azimuth_deg;
    end
     
     
    % % 输出结果
    % disp("太阳高度角:");
    % Sun_altitude_angle'
    % disp("太阳方位角:");
    % Sun_azimuth_deg'
end