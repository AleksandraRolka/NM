/****** program lab_8, sprawozdanie nr 8 na jego podstawie ******/
#include <stdio.h>
#include <math.h>
#include "nrutil.h"
#include "nrutil.c"
#include "gaussj.c"

/* Dyrektywy dla Taurusa  */
/*
#include "/opt/NR/numerical_recipes.c/nrutil.h"
#include "/opt/NR/numerical_recipes.c/nrutil.c"
#include "/opt/NR/numerical_recipes.c/gaussj.c"
*/



float func_1(float x)
{
	return 1.0 / ( 1.0 + pow(x,2) );
}

float func_2(float x)
{
	return cos(2.0 * x);
}

void wyzM(float *xw,float *yw, float *w, int n, float alfa, float beta)
{
	float **A = matrix(1, n, 1, n);
	float **d = matrix(1, n, 1, n);
	float lambda = 0.5;
	float mi = 1 - lambda;
	float h = ( xw[n] - xw[1] ) / ( n - 1 );

	for (int i = 1; i <= n; i++){
		for (int j = 1; j <= n; j++){
			A[i][j] = 0.0;
			if ( i == j ){
				if ( i ==1 && j == 1 || i == n && j == n){
					A[i][j] = 1.0;
				}
				else{
					A[i][j] = 2.0;
				}
			}
			if (i == (j - 1) )
				A[i][j] = lambda;
			
			if ( i == (j + 1) )
				A[i][j] = mi;
			if ( i == 1 && j ==2 || i == n && j == n - 1)
				A[i][j] = 0.0;
		}
		
		if (i == 1)
			d[i][1] = alfa;
		else if (i == n)
			d[i][1] = beta;
		else{
			d[i][1] = ( 6.0 / (2 * h) ) * ( ( yw[i+1] - yw[i] ) / h - ( yw[i] - yw[i-1] ) / h );
		}
		// printf("%d\t%f\n",i,d[i][1]); 
	}
	
	gaussj(A,n,d,1);
	
	// przypisanie wyniku procedury gaussj z nadpisanego 'wektora' d do wektora w
	for (int i = 1; i <= n; i++)
			w[i] = d[i][1];
	
	free_matrix(A, 1, n, 1, n);
	free_matrix(d, 1, n, 1, n);
}

float wyzSx(float *xw,float *yw, float *m, int n, float x)
{
	float sx = 0.0;
	float A,B;
	float h = ( xw[n] - xw[1] ) / ( n - 1 );
	
	for (int i = 2; i <= n; i++){
		if ( x >= xw[i-1] && x <= xw[i] ) {
			A = ( yw[i] - yw[i-1] ) / h - ( h / 6.0 ) * ( m[i] - m[i-1] );
			B = yw[i-1] - m[i-1] * ( pow( h, 2.0 ) / 6.0 );
		
			sx = m[i-1] * pow(( xw[i] - x ), 3.0) / (6.0 * h) + m[i] * pow( x - xw[i-1], 3.0) / ( 6.0 * h ) + A * ( x - xw[i-1] ) + B;
		}
	}
	
	return sx;
	
}

////////////////////////////////////////////////////////////////////////////////////////////////////


int main()
{

	FILE* file_f1 = fopen("f1.dat", "wr");
	//FILE* file_f2 = fopen("f2.dat", "wr");	// <-- dla drugiej funkcji

	int tab_n[4] = { 5, 8, 21, 10 };
	int n;
	
	for (int k = 0; k < 4; k++){
		n = tab_n[k];
		
		float x_min = -5.0, x_max = 5.0;
		float h = (x_max - x_min) / (n - 1.0);
		float alfa = 0.0, beta = 0.0;
		float delta_x = 0.01;

		float *xw = vector(1, n);
		float *yw = vector(1, n);
		float *m = vector(1, n);
////////////////////////////////////////////////////////////////////////////////////////////////////
	
		for (int i = 1; i <= n; i++){
			
			xw[i] = x_min + (i-1) * h;
			yw[i] = func_1( xw[i] );
			//yw[i] = func_2( xw[i] );		// <-- wersja dla drugiej funkcji
		}
			
		wyzM(xw, yw, m, n, alfa, beta);	

	
		float x = x_min;

		while ( x < x_max ){
			fprintf(file_f1,"%f\t%f\n",x, wyzSx(xw, yw, m, n, x));
		//	fprintf(file_f2,"%f\t%f\n",x, wyzSx(xw, yw, m, n, x));		// <-- wersja dla drugiej funkcji
			x += delta_x;
		}
		fprintf(file_f1,"\n\n");
		//fprintf(file_f2,"\n\n");		// <-- wersja dla drugiej funkcji
		
		// 	<-- odkomentować dla funkcji func_1
		///*
		if (n == 10) {
			printf("%d\n",n);
			FILE* file_poch = fopen("pochodne.dat", "wr");
			
			for (int i = 1; i <= n; i++){
				float iloraz_roznicowy = ( func_1(xw[i] - delta_x) - 2 * func_1(xw[i]) + func_1(xw[i] + delta_x) ) / pow(delta_x, 2.0);
				
				fprintf(file_poch,"%f\t%f\t%f\n", xw[i], m[i], iloraz_roznicowy);
			}
			fclose(file_poch);
		}
		//*/

////////////////////////////////////////////////////////////////////////////////////////////////////
	
		//zwolnienie pamieci
		free_vector(xw, 1, n);
		free_vector(yw, 1, n);
		free_vector(m, 1, n);
	}
	
	fclose(file_f1);
	//fclose(file_f2);	// <-- dla drugiej funkcji

	

	return 0;
}