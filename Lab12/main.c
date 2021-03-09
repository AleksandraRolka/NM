/****** program lab_12, sprawozdanie nr 12 na jego podstawie ******/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define n 8


double func(double x)
{
	return log(  pow(x, 3.0) + 3 * pow(x, 2.0) + x + 0.1  ) * sin( 18 * x );
}


void printTab (double **tab)
{
	for (int i = 0; i <= n; i++ )
	{
		for (int j = 0; j < i+1; j++ )
			printf("%.10f\t",tab[i][j]);
		printf("\n");
	}
}


void RichardsonExtrapolation(double** tab)
{
	for (int i = 1; i <= n; i++)
	{
		for (int k = 1; k < i+1; k++)
			tab[i][k] = ( pow(4, k) * tab[i][k-1] - tab[i-1][k-1] ) / ( pow(4, k) - 1 );
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////


int main()
{
	
	double a = 0.0, b = 1.0;
	
	double** D;
	
	D = malloc( (n+1)* sizeof(double*));
	
	for (int i = 0; i <= n; i++)
		D[i] = malloc( (i+1) * sizeof(double));
	
	
	
////////// metoda Simpsona //////////

	for (int w = 0; w <= n; w++)
	{
		double h = (b - a) / pow(2.0, (double)(w + 1));
		int N = pow(2, w + 1);
		double sum = 0.0;
		
		for (int i = 0; i < N/2; i++)
		{
			sum += (h / 3) * ( func(a + (2 * i) * h ) + 4 * func( a + (2 * i + 1) * h )+ func(a + (2 * i + 2) * h) );
		} 
		D[w][0] = sum;
	}
	RichardsonExtrapolation(D);
	
	printf("Całkowanie metodą Simpsona:\n");
	printTab(D);
	
////////// metoda Milne'a //////////


	for (int w = 0; w <= n; w++)
	{
		double h = (b - a) / pow(2.0, (double)(w + 2));
		int N = pow(2, w + 2);
		double sum = 0.0;
		
		for (int i = 0; i <= N / 4 - 1; i++)
		{
			sum += ( (4 * h) / 90 ) * ( 7 * func(a + (4 * i) * h ) + 32 * func( a + (4 * i + 1) * h ) + 12 * func( a + (4 * i + 2) * h ) + 32 * func( a + (4 * i + 3) * h ) + 7 * func( a + (4 * i + 4) * h ));
		} 
		D[w][0] = sum;
	}
	RichardsonExtrapolation(D);
	
	printf("\nCałkowanie metodą Milne'a:\n");
	printTab(D);



////////////////////////////////////////////////////////////////////////////////////////////////////	

// zwolnienie pamięci

	for(int i = 0; i <= n; i++)
		free(D[i]); 
	
    free(D);
	
	return 0;
}