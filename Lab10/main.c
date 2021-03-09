/****** program lab_10, sprawozdanie nr 10 na jego podstawie ******/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define N 200


double d_rand (const double min , const double max)
{
	double r = ( double ) rand () / RAND_MAX ;  	// Przedzial [0 , 1]
	r = r * ( max - min ) + min ; 		// Przeskalowanie do [min , max]
	
	return r;
}

double func(double x, double y)
{	
	return sin(x) * sin(y) - exp( - pow((x + M_PI/2), 2.0) - pow((y - M_PI/2), 2.0) );
}

////////////////////////////////////////////////////////////////////////////////////////////////////


int main()
{

	FILE* file_w0 = fopen("w0.dat", "wr");
	FILE* file_T = fopen("T.dat", "wr");
	
	
	double x_min = -10.0, y_min = -10.0, x_max = 10.0,y_max = 10.0, T, delta_x, delta_y;
	double x[N], y[N];
	double x_0 = 5.0, y_0 = 5.0;
	
////////////////////////////////////////////////////////////////////////////////////////////////////	


	for ( int i = 0; i < N; i++){
		x[i] = x_0;
		y[i] = y_0;
	}

			
	for (int i_T = 0; i_T <= 20; i_T += 1) {
			
		T = 10.0 / pow(2, i_T);
		
		for (int k = 0; k < 100; k += 1) {
			for (int i = 0; i < N; i += 1) {
				
				delta_x = d_rand( -1.0, 1.0);
				delta_y = d_rand( -1.0, 1.0);
					
				while ( fabs(x[i] + delta_x) > x_max || fabs(x[i] + delta_x) < x_min)
					delta_x = d_rand( -1.0, 1.0);
				while ( fabs(y[i] + delta_y) > y_max || fabs(y[i] + delta_y) < y_min)
					delta_y = d_rand( -1.0, 1.0);
					
				if ( func( x[i] + delta_x, y[i] + delta_y ) < func( x[i], y[i]) ) {	
				
					x[i] += delta_x;
					y[i] += delta_y;
				}
						
				else if( d_rand(0.0, 1.0) < exp( -( func( x[i] + delta_x, y[i] + delta_y ) - func( x[i], y[i]) ) / T) ) {
					
					x[i] += delta_x;
					y[i] += delta_y;
				}
			}
			fprintf(file_w0, "%f\n", func(x[0],y[0]));
		}
			

		if( i_T == 0 || i_T == 7 || i_T == 20){
			for ( int i = 0; i < N; i++){
				///printf("%f\t%f\t%f\n", x[i], y[i],func(x[i],y[i]));
				fprintf(file_T, "%f\t%f\n", x[i], y[i]);
			}
			fprintf(file_T, "\n\n");
		}
		
	}
	
	double x_minimum = x[0];
	double y_minimum = y[0];
	double mininum_val = func(x[0], y[0]);
	
	for ( int i = 0; i < N; i++){
		if ( i > 0){
			if ( func(x[i], y[i]) < mininum_val) {
				mininum_val = func(x[i], y[i]);
				x_minimum = x[i];
				y_minimum = y[i];
			}
		}	
	}

	printf("\nMinimum = %f, położenie: (%f, %f)\n\n", mininum_val, x_minimum, y_minimum);
	
	
////////////////////////////////////////////////////////////////////////////////////////////////////	
	
	fclose(file_w0);
	fclose(file_T);
	
	return 0;
}