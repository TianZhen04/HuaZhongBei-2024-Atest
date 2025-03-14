function I = cal_I(month,day,hour,minute,azimuth,tilt)
    data2 = [1	1405
            2	1394
            3	1378
            4	1353
            5	1334
            6	1316
            7	1308
            8	1315
            9	1330
            10	1350
            11	1372
            12	1392];
    geo_info.latitude = 367/12; % 纬度，单位：度
    geo_info.longitude = 6859/60; % 经度，单位：度
    geo_info.year = 2025;
    depth=1000;
    w = -5.7714*(10^-4);
    [altitude_angle, azimuth_deg] = solar_angles(geo_info.year,month, day,hour,minute, geo_info.longitude,geo_info.latitude);
    if altitude_angle <= 0
        I=0;return;
    end
    I = data2(month,2)*exp(w*depth/sind(altitude_angle));
    if tilt~=0
        a = [-sind(azimuth), -cosd(azimuth), tand(90-tilt)];
    else
        a = [0 0 1];
    end
    b = [sind(azimuth_deg), cosd(azimuth_deg), tand(altitude_angle)];
    I = I * cosd(vector_angle(a,b));
end