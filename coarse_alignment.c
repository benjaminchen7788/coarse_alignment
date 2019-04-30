// coarse alignment in c program benjamin_2019/04/25_13:55
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ROW 3
#define COLUMN 3

// sub-program for array calcution_benjamin
void inverse_33(double (*matrix)[3] ){
	double det= (*(*(matrix+0)+0)) * (*(*(matrix+1)+1)) * (*(*(matrix+2)+2)) + 
	(*(*(matrix+1)+0)) * (*(*(matrix+2)+1)) * (*(*(matrix+0)+2)) + 
	(*(*(matrix+2)+0)) * (*(*(matrix+0)+1)) * (*(*(matrix+1)+2)) -
	(*(*(matrix+0)+2)) * (*(*(matrix+1)+1)) * (*(*(matrix+2)+0)) -
	(*(*(matrix+1)+2)) * (*(*(matrix+2)+1)) * (*(*(matrix+0)+0)) -
	(*(*(matrix+2)+2)) * (*(*(matrix+0)+1)) * (*(*(matrix+1)+0));
	double a[3][3] = {0};
	a[0][0] = matrix[1][1]*matrix[2][2]-matrix[1][2]*matrix[2][1];
	a[0][1] = -(matrix[0][1]*matrix[2][2]-matrix[0][2]*matrix[2][1]);
	a[0][2] = matrix[0][1]*matrix[1][2]-matrix[0][2]*matrix[1][1];
	a[1][0] = -(matrix[1][0]*matrix[2][2]-matrix[1][2]*matrix[2][0]);
	a[1][1] = matrix[0][0]*matrix[2][2]-matrix[0][2]*matrix[2][0];
	a[1][2] = -(matrix[0][0]*matrix[1][2]-matrix[0][2]*matrix[1][0]);
	a[2][0] = matrix[1][0]*matrix[2][1]-matrix[1][1]*matrix[2][0];
	a[2][1] = -(matrix[0][0]*matrix[2][1]-matrix[0][1]*matrix[2][0]);
	a[2][2] = matrix[0][0]*matrix[1][1]-matrix[0][1]*matrix[1][0];
	for(int i =0;i < ROW;++i){
		for(int j = 0;j < COLUMN;++j)
			matrix[i][j] = (1/det)*a[i][j];
	}
	return;
}

void product_33_31(double (*matrix)[COLUMN], double *vector, double *result){
	for(int i = 0;i < ROW;++i){
		for(int j = 0;j < COLUMN;++j)
			result[i] = result[i] + matrix[i][j]*vector[j];
	}
	return;
}

void product_33_33(double (*matrix_1)[COLUMN], double (*matrix_2)[COLUMN], double *result){
	for(int i = 0;i < ROW;++i){
		for(int j = 0;j < COLUMN;++j){
			*(result + (i*ROW)+j) = matrix_1[i][0]*matrix_2[0][j]
					      + matrix_1[i][1]*matrix_2[1][j]
					      + matrix_1[i][2]*matrix_2[2][j];
		}
	}
	return;
}

void dcm2euler(double *dc){
	const float rad2deg = 57.29578121;
	const float pi = 3.1415926;
	double roll = 0;
	double pitch = 0;
	double yaw = 0;

	if(*(dc+(2*3)+0) <= -0.999){
		roll = 0; // NAN
		yaw = atan2(*(dc+(1*3)+2)-*(dc+(0*3)+1), *(dc+(0*3)+2)+*(dc+(1*3)+1));
	} else if(*(dc+(2*3)+0) >= 0.999){
		roll = 0; // NAN
		yaw = pi + atan2(*(dc+(1*3)+2)-*(dc+(0*3)+1), *(dc+(0*3)+2)+*(dc+(1*3)+1));
	} else{
		roll = atan2(*(dc+(2*3)+1), *(dc+(2*3)+2));
		yaw = atan2(*(dc+(1*3)+0), *(dc+(0*3)+0));
	}
	pitch = atan( -(*(dc+(2*3)+0)) / pow(pow(*(dc+(2*3)+1),2)+pow(*(dc+(2*3)+2),2),0.5));
	roll = roll * rad2deg;
	pitch = pitch * rad2deg;
	yaw = yaw * rad2deg;

	printf("\nroll = %f, pitch = %f, yaw = %f \n\n", roll, pitch, yaw);

	return;
}

int main(int argc, char *argv[]){
	// printf("1 \n"); // test point 1
	// file reading: dynamic memory for time and 6DOF data----------------------
	int size = 60000; // enter the size of data ***
	char *filename = "mDNmE.csv"; // enter the file name here ***
	// enter the parameters from calibration for correcting the sensor value ***
	const float std_gravity = 9.787645;
	double acc_biasX = 0.00231658*std_gravity;
	double acc_biasY = 0.00077976*std_gravity;
	double acc_biasZ = 0.00058692*std_gravity;
	double acc_sfX = -0.0017783;
	double acc_sfY = -0.0020704;
	double acc_sfZ = -0.0018665;
	double acc_non_ortho[ROW][COLUMN] = {
				{1, 5.1930e-17, 1.0296e-19},
				{9.4007e-5, 1, 41.0484e-19},
				{6.5781e-4, 2.8619e-4, 1}
			};
	double gyro_biasX = -6.6621e-05;
	double gyro_biasY = -1.5644e-04;
	double gyro_biasZ =  5.1580e-04;
	double gyro_sfX = -0.0011777;
	double gyro_sfY = -9.4293e-04;
	double gyro_sfZ =  8.4346e-04;
	double gyro_non_ortho[ROW][COLUMN] = {
				{1.0000e+00, 1.9924e-17, 1.1016e-19},
				{-1.3217e-03, 1.0000e+00, -1.9704e-19},
				{2.1271e-03, -5.2807e-04, 1.0000e+00}
			};

	float time_temp = 0.0;
	float a_x = 0; float a_x_temp = 0.0; 
	float a_y = 0; float a_y_temp = 0.0; 
	float a_z = 0; float a_z_temp = 0.0; 
	float w_x = 0; float w_x_temp = 0.0; 
	float w_y = 0; float w_y_temp = 0.0; 
	float w_z = 0; float w_z_temp = 0.0; 
	// file reading
	FILE* fp; 
	fp = fopen(filename, "r");
	if(fp == NULL) { 
		printf("fail to open the file. \n");
		return 1;
	} else {
		// data format: time, accel_x, accel_y, accel_z, ang_x, ang_y, ang_z
		int sensor_data_size = 0;
		for(int i = 0;i < size; ++i){
			fscanf(fp, "%f, %f, %f, %f, %f, %f, %f", &time_temp, 
			&a_x_temp, &a_y_temp, &a_z_temp, &w_x_temp, &w_y_temp, &w_z_temp);

			a_x += a_x_temp; a_y += a_y_temp; a_z += a_z_temp;
			w_x += w_x_temp; w_y += w_y_temp; w_z += w_z_temp;
		}
	}

	// parameters for data type transformation, fusion--------------------------
	const float deg2rad = 0.017453292;
	const float rad2deg = 57.29578121;

	double lat = (24 + 43.0/60.0 + 29.0/3600.0)* deg2rad;
	double lon = (120 + 54.0/60.0 + 35.0/3600.0)* deg2rad;
	double alt = 25;
	double g0 = 9.7803267715;		double g1 = 0.0052790414;
	double g2 = 0.0000232718;		double g3 = -0.000003087691089;
	double g4 = 0.000000004397731;	double g5 = 0.000000000000721;
	double gravity = g0* (1+ g1* pow(sin(lat),2) + g2* pow(sin(lat),4))
						 + (g3 + g4*pow(sin(lat),2))* alt + g5* pow(alt,2);
	double gn[ROW] = {0, 0, gravity};
	// earth rotation
	double omega_e = 7.2921158*1e-5;
	double omega_ie_n[ROW] = {omega_e*cos(lat), 0, -omega_e*sin(lat)};

	// coarse alignment algorithm-----------------------------------------------
	float delta_t = 0.002; // *enter the step of time, 500Hz=0.002sec
	double a_x_m = 0.0; double a_y_m = 0.0; double a_z_m = 0.0;
	double w_x_m = 0.0; double w_y_m = 0.0; double w_z_m = 0.0;

	// sensor data type: delta theta, delta v
	a_x_m = a_x / size / delta_t;
	a_y_m = a_y / size / delta_t;
	a_z_m = a_z / size / delta_t;
	w_x_m = w_x / size / delta_t;
	w_y_m = w_y / size / delta_t;
	w_z_m = w_z / size / delta_t;
	// sensor data without fusion
	double acc_x_m = a_x_m ;
	double acc_y_m = a_y_m ;
	double acc_z_m = a_z_m ;
	double gyro_x_m = w_x_m *deg2rad;
	double gyro_y_m = w_y_m *deg2rad;
	double gyro_z_m = w_z_m *deg2rad;
	// with fusion
	/*double acc_x_m = (a_x_m - acc_biasX) / (1+acc_sfX);
	double acc_y_m = (a_y_m - acc_biasY) / (1+acc_sfY);
	double acc_z_m = (a_z_m - acc_biasZ) / (1+acc_sfZ);
	double gyro_x_m = (w_x_m - gyro_biasX) / (1+gyro_sfX);
	double gyro_y_m = (w_y_m - gyro_biasY) / (1+gyro_sfY);
	double gyro_z_m = (w_z_m - gyro_biasZ) / (1+gyro_sfZ);*/

	// calculate: a[] and w[] with data after fusion
	// accelerometer
	double acc_m[COLUMN] = {acc_x_m, acc_y_m, acc_z_m};
	inverse_33(acc_non_ortho);
	double acc_m_non_ortho[3] = {0};
	product_33_31(acc_non_ortho, acc_m, acc_m_non_ortho);

	// gyro
	double gyro_m[COLUMN] = {gyro_x_m, gyro_y_m, gyro_z_m};
	inverse_33(gyro_non_ortho);
	double gyro_m_non_ortho[3] = {0};
	product_33_31(gyro_non_ortho, gyro_m, gyro_m_non_ortho);
	// test function: for check the value of gyro
	double a_cross_w_m[3] = {acc_y_m*gyro_z_m - acc_z_m*gyro_y_m,
                   			 acc_z_m*gyro_x_m - acc_x_m*gyro_z_m,
                   			 acc_x_m*gyro_y_m - acc_y_m*gyro_x_m };

	//sensor data matrix [fx, fy, fz; wx, wy, wz; f X w...];
	double IMU_matrix_m[3][3] = {acc_x_m,        acc_y_m,        acc_z_m,
                    			gyro_x_m,       gyro_y_m,       gyro_z_m,
                    			a_cross_w_m[0], a_cross_w_m[1], a_cross_w_m[2]};
   	// inverse matrix of theoritical value: a[] and w[]
    double NED_matrix_inv[3][3] = {-tan(lat)/gravity, 1/(omega_e*cos(lat)), 0,
                                    0, 0, -1/(gravity*omega_e*cos(lat)),
                      		   	   -1/gravity, 0,                    0};
    // transformation matrix Cbn
    double C_bn[9] = {0};

    product_33_33(NED_matrix_inv, IMU_matrix_m, C_bn);
    dcm2euler(C_bn);

    // test function for C_bn
	printf("C_bn = \n");
    for(int i = 0;i < ROW;++i){
    	for(int j = 0;j < COLUMN;++j)
    		printf("%-13.6G ", *(C_bn + (i*ROW) + j));
    	printf("\n");
    } printf("\n");

	fclose(fp);
	return 0;
}

