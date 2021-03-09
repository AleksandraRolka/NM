/****** program lab_11, sprawozdanie nr 11 na jego podstawie ******/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "gsl_errno.h"
#include "gsl_fft_complex.h"

double f_0_bez_szumu(double omega, double t)
{
	return sin( omega * t) + sin(2 * omega * t) + sin(3 * omega * t);
}

double f_0_z_szumem(double omega, double t)
{
	return sin(omega * t) + sin(2 * omega * t) + sin(3 * omega * t) + ( rand() / (RAND_MAX + 1.0) - 0.5);
}

double f_gaussowska( double sigma, double t)
{
	return ( 1.0 / ( sigma * sqrt(2.0 * M_PI)) ) * exp(- (t * t) / ( 2 * sigma * sigma));
}

////////////////////////////////////////////////////////////////////////////////////////////////////


int main()
{
	
	FILE* file_k8 = fopen("k8.dat", "wr");
	FILE* file_k10 = fopen("k10.dat", "wr");
	FILE* file_k12 = fopen("k12.dat", "wr");
	
	double T = 1.0;
	double t_max = 3 * T;
	double sigma = T / 20.0;
	double omega = (2 * M_PI) / T;

	for (int k = 8; k <= 12; k += 2) {
		
		int N = pow(2, k);
		double delta_t = t_max /(double) N;
		
		
		double f_0[2*N];
		double f[2*N];
		double g_1[2*N];
		double g_2[2*N];
		
		for (int i = 0; i < N; i++)
		{
			f_0[2 * i] = f_0_bez_szumu(omega, i * delta_t);
			f[2 * i] = f_0_z_szumem(omega, i * delta_t);
			g_1[2 * i] = f_gaussowska( sigma, i * delta_t);
			g_2[2 * i] = f_gaussowska( sigma, i * delta_t);
			
			f_0[2 * i + 1] = 0.0;
			f[2 * i +1] = 0.0;
			g_1[2 * i + 1] = 0.0;
			g_2[2 * i + 1] = 0.0;
		}
	   
		for (int i = 0; i < N; i++){
			if ( k == 8)
				fprintf(file_k8, "%f\t%f\n", i * delta_t, f[2 * i] );
			if ( k == 10)
				fprintf(file_k10, "%f\t%f\n", i * delta_t, f[2 * i] );
			if ( k == 12)
				fprintf(file_k12, "%f\t%f\n", i * delta_t, f[2 * i] );
	   }
	   
		if ( k == 8)
			fprintf(file_k8, "\n\n");
		if ( k == 10)
			fprintf(file_k10, "\n\n");
		if ( k == 12)
			fprintf(file_k12, "\n\n");
	   
	   
	    gsl_fft_complex_radix2_forward( f, 1, N);
	    gsl_fft_complex_radix2_forward( g_1, 1, N);
	    gsl_fft_complex_radix2_backward( g_2, 1, N);

	
		for (int i = 0; i < N; i++) 
		{
			double a_1 = f[2 * i];							// Re { f(k_1) }
			double b_1 = f[2 * i + 1];						// Im { f(k_1) }
			double a_2 = g_1[2 * i] + g_2[2 * i];			// Re { g(k_1) }
			double b_2 = g_1[2 * i + 1] + g_2[2 * i + 1];	// Im { g(k_1) }
			
			f[2 * i] = a_1 * a_2 - b_1 * b_2;				// Re { f(k_1) * g(k_1) }
			f[2 * i + 1] = a_1 * b_2 + a_2 * b_1;			// IM { f(k_1) * g(k_1) }
		}
	
		gsl_fft_complex_radix2_backward( f, 1, N);
	
	
		double f_max = fabs(f[0]);
	
		for (int i = 0; i < N; i++) {
			if (fabs( f[2 * i] ) > f_max)
				f_max = fabs( f[2 * i] );
		}

		printf("delta_t = %f\n", delta_t);
		printf("f_max= %f\n\n", f_max);

		for (int i = 0; i < N; i++){
			if ( k == 8)
				fprintf(file_k8, "%f\t%f\n", i * delta_t, (f[2 * i] * 2.5) / f_max);
			if ( k == 10)
				fprintf(file_k10, "%f\t%f\n", i * delta_t, (f[2 * i] * 2.5) / f_max);
			if ( k == 12)
				fprintf(file_k12, "%f\t%f\n", i * delta_t, (f[2 * i] * 2.5) / f_max);
	   }
	   
	}	   


////////////////////////////////////////////////////////////////////////////////////////////////////	
	fclose(file_k8);
	fclose(file_k10);
	fclose(file_k12);
	
	return 0;
}