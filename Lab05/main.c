/****** program lab_5, sprawozdanie nr 5 na jego podstawie ******/
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

float* matrix_vector_mltpl (float** M, float* v, int size);
float scalar_product (float* v1, float* v2, int size);
float norma( float* v, int size);
float* fun_vec_norm( float* v, float norm, int size);
float** hotelling_reduction(float **W, float lambda, float* v1, float* v2, int size);

int main()
{

	FILE* file = fopen("dane.dat", "wa");
	
	int n = 7;
	int itr = 8;
	float **A = matrix(1, n, 1, n);
	float **W = matrix(1, n, 1, n);
	float *d = vector(1, n);
	float *e = vector(1, n);
	float *x_0 = vector(1, n);
	float *x = vector(1, n);
	float *lambda = vector(1, n);
	

/////////////////////////////////////// metoda I - bezpośrednia ///////////////////////////////////////

	
	//wypełnienie macierzy A
	for (int i = 1 ; i <= n ; i++){
		for (int j = 1 ; j <= n ; j++)
			A[i][j] = sqrt( i + j );
	}
	
	printf("\n");
	for (int i = 1; i <= n; i++)
	{
		for (int j = 1; j <= n; j++)
		{
			printf("%f ", A[j][i]);
		}
		printf("\n");
	}
	printf("\n");
	//
	tred2(A, n, d, e );
	//
	tqli(d, e, n, A );
	
	
	fprintf (file, "Metoda I - bezpośrednia:\n\n");
	
	//wartości własne - zapis do pliku
	fprintf (file, " Wartosci wlasne: \n");
	
	for (int i = 1 ; i <= n ; i++)
		fprintf(file, "  %g\n", d[i]);
	fprintf(file, "\n");
	

////////////////////////////////////// metoda II - iteracyjna //////////////////////////////////////	

	
	//wypełnienie macierzy W
	for (int i = 1 ; i <= n ; i++){
		for (int j = 1 ; j <= n ; j++)
			W[i][j] = sqrt( i + j );
	}
	
	
	for (int k = 1 ; k <= n ; k++){
		//wypełnienie wektora startowego x_0
		for (int i = 1 ; i <= n ; i++)
			x_0[i] = 1.0;
	
	
		for (int i = 1 ; i <= itr ; i++){
			x = matrix_vector_mltpl (W, x_0, n);
			lambda[k] = scalar_product( x, x_0, n) / scalar_product( x_0, x_0, n);
			//x = fun_vec_norm( x, norma( x, n), n );
			
			x_0 = x;
		}
		W = hotelling_reduction(W, lambda[k], x_0, x_0, n);
	}

	fprintf (file, "Metoda II - iteracyjna:\n\n");
	
	//wartości własne - zapis do pliku
	fprintf (file, " Wartosci wlasne: \n");
	
	for (int i = 1 ; i <= n ; i++)
		fprintf(file, "  %g\n", lambda[i]);
	fprintf(file, "\n");

////////////////////////////////////////////////////////////////////////////////////////////////////
	fclose(file);
	
	//zwolnienie pamieci
	free_matrix(A, 1, n, 1, n);
	free_matrix(W, 1, n, 1, n);
	free_vector(d, 1, n);
	free_vector(e, 1, n);
	free_vector(x_0, 1, n);
	free_vector(x, 1, n);
	free_vector(lambda, 1, n);
	
	
	return 0;
}



////////////////////////////////////////////////////////////////////////////////////////////////////
//wykorzystane w programie wlasne funkcje :
////////////////////////////////////////////////////////////////////////////////////////////////////


float* matrix_vector_mltpl (float** M, float* v, int size)
{
    float* res = vector(1, size);
	for (int i = 1 ; i <= size ; i++){
        float sum = 0;
    	for (int j = 1 ; j <= size ; j++)
            sum += M[i][j] * v[j];
        
        res[i] = sum;
    }
    return res;
}

float scalar_product (float* v1, float* v2, int size)
{
	float sum = 0;
	for (int i = 1 ; i <= size ; i++)
		sum += v1[i] * v2[i];
	
	return sum;
}

float norma( float* v, int size)
{
	float res = 0;
	for (int i = 1 ; i <= size ; i++)
		res += v[i] * v[i];
	
	//res = sqrt(scalar_product( v, v, size) );
	res = sqrt(res);
	return res;
}

float* fun_vec_norm( float* v, float norm, int size)
{
	float* res = vector(1, size);
	for (int i = 1 ; i <= size ; i++)
		res[i] = v[i] / norm;
	
	return res;
}

float** hotelling_reduction(float **W, float lambda, float* v1, float* v2, int size)
{
	float **R = matrix(1, size, 1, size);
	float **RH = matrix(1, size, 1, size);

    for (int i = 1 ; i <= size ; i++){
    	for (int j = 1 ; j <= size ; j++)
            R[i][j] = v1[i] * v2[j];
    }
    for (int i = 1 ; i <= size ; i++){
    	for (int j = 1 ; j <= size ; j++)
			RH[i][j] = W[i][j] - lambda * R[i][j];
	}
	return RH;
}