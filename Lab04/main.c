/****** program lab_4, sprawozdanie na jego podstawie ******/
#include <stdio.h>
#include <math.h>
#include "nrutil.h"
#include "nrutil.c"
#include "tred2.c"
#include "tqli.c"
#include "pythag.c"

/* Dyrektywy dla Taurusa  */
/*
#include "/opt/NR/numerical_recipes.c/nrutil.h"
#include "/opt/NR/numerical_recipes.c/nrutil.c"
#include "/opt/NR/numerical_recipes.c/tred2.c"
#include "/opt/NR/numerical_recipes.c/tqli.c"
#include "/opt/NR/numerical_recipes.c/pythag.c"
*/

int main()
{
	int n_x = 20;
	int n_y = 20;
	int n = n_x * n_y;
	int m = 10;
	float t = -0.021;

	float **H, **Y, **X;
	float *d, *e;
	
	H = matrix(1, n, 1, n);
    Y = matrix(1, n, 1, n);
    X = matrix(1, n, 1, n);
    d = vector(1, n);
    e = vector(1, n);
	
	//wypełnienie macierzy H
	for ( int i = 1; i <= n_x; i++){
		for (int j = 1; j <= n_y; j++){
			int l = j + (i - 1) * n_y;
			for (int k = 1; k <= n; k++)
				H[l][k] = 0.;
			if(i > 1)
				H[l][l - n_y] = t;  //dla i=1 nie ma sasiada z lewej strony
			if(i < n_x)
				H[l][l + n_y] = t; //dla i=n_x nie ma sasiada z prawej strony
			H[l][l] = -4 * t;
			if(j > 1) 
				H[l][l - 1] = t;   //dla j=1 nie ma sasiada ponizej siatki
			if(j < n_y) 
				H[l][l + 1] = t;  //dla j=n_y nie ma sasiada powyzej siatki
		}
	}
	
	//przekształcenie macierzy H do postaci trójdiagonalnej
	tred2(H, n, d, e );
	
	//Y = I, macierz jednostkowa
	for (int i = 1 ; i <= n ; i++){
		for (int j = 1 ; j <= n ; j++){
			if( i == j )
				Y[i][j] = 1.0;
			else
				Y[i][j] = 0.0;
		}
	}
	
	//diagonalizacja macierzy T
	tqli(d, e, n, Y );
	
	//wymnozenie macierzy 
	for (int i = 1 ; i <= n ; i++){
		for (int j = 1 ; j <= n ; j++){
			float sum = 0.0;
			for (int k = 1 ; k <= n ; k++){
				sum += H[i][k] * Y[k][j];
			}
			X[i][j] = sum;
		}
	}
	
	//sortowanie energii oraz indeksów wektorów (tablica indx)
	int indx[n+1];
	float e1, e2, l1, l2;
	for (int l = 1; l <= n; l++) 
		indx[l] = l; 			// inicjalizacja
	for (int l = 1; l <= n - 1; l++){
		for (int k = n; k >= l + 1; k--){
			e1 = d[k-1];
			e2 = d[k];
			l1 = indx[k - 1];
			l2 = indx[k];
			if(e2 < e1){ //wymieniamy energie i indeksy wektorów miejscami
				d[k] = e1;
				d[k - 1] = e2;
				indx[k] = l1;
				indx[k-1] = l2;
			}
		}
	}

	FILE *fp;
	fp = fopen("dane.dat", "w");
	
	for (int i = 1; i <= n_x; i++){
		for (int j = 1; j <= n_y; j++){
			int l = j + (i-1) * n_y;
			fprintf(fp, "%6d %6d ", i, j);
			for (int k = 1;k <= m; k++)
				fprintf(fp, " %12.6f ", X[l][ indx[k] ]);
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	
	//wartości własne
	for (int i = 1; i <= m; i++)
        printf("%.9f\n", d[i]);
	
	//zwolnienie pamieci
	free_matrix(H, 1, n, 1, n);
	free_matrix(Y, 1, n, 1, n);
	free_matrix(X, 1, n, 1, n);
	free_vector(d, 1, n);
	free_vector(e, 1, n);
	
	return 0;
}