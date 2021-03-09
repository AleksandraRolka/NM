/****** program lab_9, sprawozdanie nr 9 na jego podstawie ******/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define N 201
#define frand() ( (double)rand() ) / ( RAND_MAX + 1.0 )




double func(double x, double x_0, double x_min, double x_max, double sigma)
{
	return sin(( 14 * M_PI * x ) / ( x_max - x_min ) ) * ( exp( -( x - x_0 ) * ( x - x_0 ) / ( 2 * sigma * sigma ) ) + exp( -( x + x_0 ) * ( x + x_0) / ( 2 * sigma * sigma ) ) );
}

////////////////////////////////////////////////////////////////////////////////////////////////////


int main()
{

	FILE* file_Gram = fopen("Gram.dat", "wr");
	FILE* file_pkt = fopen("pkt.dat", "wr");
	FILE* file_approx = fopen("approx.dat", "wr");

	int  m = 50;
	
	double x_min = -4.0, x_max = 4.0;
	double sigma = ( x_max - x_min ) / 16.0;
	double x_0 = 2.0, alfa = 0, beta = 0;
	double x[N], y[N], y_szum[N];
	double phi[m+2][N];
	double c, s;
	
	double delta_x = ( x_max - x_min ) / ( N - 1.0 );
	
////////////////////////////////////////////////////////////////////////////////////////////////////
	
	for (int i = 0; i < N; i++){
			x[i] = x_min + i * delta_x;
			y[i] = func( x[i], x_0, x_min, x_max, sigma);
			y_szum[i] = y[i] + (frand() - 0.5) / 5.0;
			//printf("%d\t%f\t%f\t%f\n",i, x[i],y[i],y_szum[i]);
			fprintf(file_pkt, "%f\t%f\n", x[i],y_szum[i]);
	}

	for (int k = 0; k < N; k++)
				phi[0][k] = 0.0;
	for (int k = 0; k < N; k++)
				phi[1][k] = 1.0;

	for (int j = 1; j < m+1; j++) {
					
		double a_l = 0, a_m = 0, b_l = 0, b_m = 0; 	/// zmienne pomocnicze do obliczania wartosci mianownikow, liczników 
													/// wyrazen na alfe i bete
		for( int i = 0; i < N; i++){
			a_l += x[i] * phi[j][i] * phi[j][i];
			a_m += phi[j][i] * phi[j][i];
			b_l += x[i] * phi[j-1][i] * phi[j][i];
			b_m += phi[j-1][i] * phi[j-1][i];
		}
		alfa = a_l / a_m;
		if ( j == 1)
			beta = 0.0;
		else{
			beta = b_l / b_m ;
		}
		printf("%g\t%g\n",alfa,beta);		// wypisanie do sprawdzenia wartości alfa, beta   
		
		for (int k = 0; k < N; k++) 
				phi[j+1][k] = ( x[k] - alfa ) * phi[j][k] - beta * phi[j-1][k];
	}
	

//	for (int j = 1; j < m+1; j++) {
//		for (int i = 0; i < N; i++)
//			printf("phi[%d][%d] = %g\n",j,i,phi[j][i]);		// wypisanie do sprawdzenia wartości phi[j][i]
//		printf("\n");
//	}


	for (int i = 0; i < N; i++)        
    {        
        fprintf(file_Gram, "%g  ", x[i]);
        for (int j = 1; j < 8; j++)
        {
            fprintf(file_Gram, "%g  ", phi[j][i] / phi[j][0]);
        }
        fprintf(file_Gram, "\n");
    }


	double F, mm[3]={10, 30, 50};
	for (int l = 0; l < 3; l++) {
		m = mm[l];
		for (int k = 0; k < N; k++) {
			F = 0;
			for (int j = 1; j < m+2; j++) {
				c = 0;
				s = 0;

				for(int i = 0; i < N ; i++) {
					c += y_szum[i] * phi[j][i];
					s += phi[j][i] * phi[j][i];
				}
				F += ( c / s )* phi[j][k];
			}

			fprintf(file_approx, "%f\t%f\n", x[k], F);
		}
		fprintf(file_approx, "\n\n");
	}


////////////////////////////////////////////////////////////////////////////////////////////////////	
	
	fclose(file_Gram);
	fclose(file_pkt);
	fclose(file_approx);
	

	return 0;
}