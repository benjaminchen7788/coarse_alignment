% coarse alignment 20190422 21:48
function coarse_alignment
clc; clear all; close all;

deg2rad = 0.017453292;
rad2deg = 57.29578121;
std_gravity = 9.797645; % Standard gravity at 45 deg latitude 

% compute Earth rotation and local gravuty at NED frame
% (24 43 29 N, 120 54 35 E, 25)  latitude, longtitude, altitude(m)
lat = ( 24 + 43/60 + 29/3600 ) * deg2rad;
lon = ( 120 + 54/60 + 35/3600 ) * deg2rad;
alt = 25;
omega_e =( 7.2921158*1e-5 );
omega_ie_n = [omega_e*cos(lat),  0,  -omega_e*sin(lat)];
g0 = 9.7803267715; g1 = 0.0052790414; g2 = 0.0000232718;
g3 = -0.000003087691089; 
g4 = 0.000000004397731; 
g5 = 0.000000000000721;
gravity = g0*(1+ g1*sin(lat)^2 + g2*sin(lat)^4 )...
        + ( g3 + g4*(sin(lat)^2))*alt + g5*alt^2;
g_n = [0,  0,  gravity];

% fusion for accelerometer
acc_biasX = 0.00231658 * std_gravity;
acc_biasY = 0.00077976 * std_gravity;
acc_biasZ = 0.00058692 * std_gravity;
acc_sfX = -0.0017783;
acc_sfY = -0.0020704;
acc_sfZ = -0.0018665;
acc_non_ortho = [1,  5.1930e-17,   1.0296e-19;
                 9.4007e-5,   1,  41.0484e-19;
                 6.5781e-4,   2.8619e-4,   1];

gyro_biasX = -6.6621e-05*deg2rad;
gyro_biasY = -1.5644e-04*deg2rad;
gyro_biasZ =  5.1580e-04*deg2rad;
gyro_sfX =  -0.0011777;
gyro_sfY =  -9.4293e-04;
gyro_sfZ =   8.4346e-04;
gyro_non_ortho = [1.0000e+00   1.9924e-17   1.1016e-19;
                 -1.3217e-03   1.0000e+00  -1.9704e-19;
                  2.1271e-03  -5.2807e-04   1.0000e+00];

% read the data of sensor
IMU_data=csvread('mDNmE.csv');
size = length(IMU_data);

time = IMU_data(:,1);
% feedback accel meters (g)
% acc_x = IMU_data(:,2)*std_gravity;
% acc_y = IMU_data(:,3)*std_gravity;
% acc_z = IMU_data(:,4)*std_gravity;
% feedback delta meters
acc_x = IMU_data(:,2);
acc_y = IMU_data(:,3);
acc_z = IMU_data(:,4);
gyro_x = IMU_data(:,5)*deg2rad;
gyro_y = IMU_data(:,6)*deg2rad;
gyro_z = IMU_data(:,7)*deg2rad;

% % 1/5 data freq. : 200Hz data to 40Hz data
% SDT_num = fix(size/5);
% SDT_gyro_x = zeros(1,SDT_num);
% SDT_gyro_y = zeros(1,SDT_num);
% SDT_gyro_z = zeros(1,SDT_num);
% for i=1:SDT_num
%   SDT_gyro_x(i)= (gyro_x(5*i-4)+ gyro_x(5*i-3)+ gyro_x(5*i-2)+ gyro_x(5*i-1)+ gyro_x(5*i))/5;
%   SDT_gyro_y(i)= (gyro_y(5*i-4)+ gyro_y(5*i-3)+ gyro_y(5*i-2)+ gyro_y(5*i-1)+ gyro_y(5*i))/5;
%   SDT_gyro_z(i)= (gyro_z(5*i-4)+ gyro_z(5*i-3)+ gyro_z(5*i-2)+ gyro_z(5*i-1)+ gyro_z(5*i))/5;
% end
% % 1/5 data freq. test
% fprintf("\n"); 
% SDT_gyro_x_m = mean(SDT_gyro_x)
% SDT_gyro_y_m = mean(SDT_gyro_y)
% SDT_gyro_z_m = mean(SDT_gyro_z)
% gyro_x_m = mean(gyro_x)
% gyro_y_m = mean(gyro_y)
% gyro_z_m = mean(gyro_z)
% fprintf("\n");

delta_t = 0.002; % 500Hz, 0.002sec, delta_theta = w * delta_t
acc_x_m = mean(acc_x)/delta_t;
acc_y_m = mean(acc_y)/delta_t;
acc_z_m = mean(acc_z)/delta_t;
acc_m = [acc_x_m; acc_y_m; acc_z_m];
% fprintf("acc_m after inv : ");
acc_m = inv(acc_non_ortho)*acc_m;

% fusion
% gyro_x_m = (mean(gyro_x) - gyro_biasX) / (1+gyro_sfX);
% gyro_y_m = (mean(gyro_y) - gyro_biasY) / (1+gyro_sfY);
% gyro_z_m = (mean(gyro_z) - gyro_biasZ) / (1+gyro_sfZ);
% without fusion
gyro_x_m = mean(gyro_x)/delta_t;
gyro_y_m = mean(gyro_y)/delta_t;
gyro_z_m = mean(gyro_z)/delta_t;

gyro_m = [gyro_x_m; gyro_y_m; gyro_z_m];
fprintf("gyro_m after inv : ");
gyro_m = inv(gyro_non_ortho)*gyro_m
gyro_x_m = gyro_m(1); gyro_y_m = gyro_m(2); gyro_z_m = gyro_m(3);

% print
gyro_m_norm = norm(gyro_m)*3600*rad2deg;

% gyro_norm_error = omega_e - norm(gyro_m);
% gyro_m = (1 + gyro_norm_error/norm(gyro_m))*gyro_m;
% norm_gyro_m = norm(gyro_m)*rad2deg*3600

% f X w for [fx, fy, fz; wx, wy, wz; f X w...] matrix
a_cross_w_m = [acc_y_m*gyro_z_m - acc_z_m*gyro_y_m;
               acc_z_m*gyro_x_m - acc_x_m*gyro_z_m;
               acc_x_m*gyro_y_m - acc_y_m*gyro_x_m ];
% sensor data matrix [fx, fy, fz; wx, wy, wz; f X w...];
IMU_matrix_m = [acc_x_m,        acc_y_m,        acc_z_m;
                gyro_x_m,       gyro_y_m,       gyro_z_m;
                a_cross_w_m(1), a_cross_w_m(2), a_cross_w_m(3)];
     
% plot(gyro_x), hold on
% plot(gyro_y),
% plot(gyro_z), hold off 
% legend({'gyrox','gyroy','gyroZ'},'FontSize',16,'Location','NorthEast');


NED_matrix_inv = [-tan(lat)/gravity,  1/(omega_e*cos(lat)),  0;
                   0,                 0,  -1/(gravity*omega_e*cos(lat));
                  -1/gravity,         0,                     0];

% transformation matrix from calibration
C_bn_ans = inv( [   1.8030e-03, -7.0438e-01,  7.0982e-01;
                   -9.5733e-04, -7.0889e-01, -7.0532e-01;
                    1.0000e+00, -5.3132e-04, -7.0610e-05] );
% [roll pitch yaw] = dcm2euler(C_bn_ans);
% fprintf("ans:   roll= %f,   pitch= %f,   yaw= %f \n\n", roll*rad2deg, pitch*rad2deg, yaw*rad2deg);

% transformation matrix from coarse alignment
NED_matrix_inv;
IMU_matrix_m
C_bn_m = NED_matrix_inv * IMU_matrix_m
[roll pitch yaw] = dcm2euler(C_bn_m);
fprintf("C_bn_m:   roll= %f,   pitch= %f,   yaw= %f \n\n\n", roll*rad2deg, pitch*rad2deg, yaw*rad2deg);
% fprintf("norm of C_bn_m = %f\n\n", norm(C_bn_m));

% C_bn_m = C_bn_m + 0.5*(eye(3)-C_bn_m*C_bn_m') * C_bn_m
% [roll pitch yaw] = dcm2euler(C_bn_m);
% norm(C_bn_m)
% fprintf("C_bn_m:   roll= %f,   pitch= %f,   yaw= %f \n", roll*rad2deg, pitch*rad2deg, yaw*rad2deg);

end

