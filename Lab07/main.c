/****** program lab_7, sprawozdanie nr 7 na jego podstawie ******/
#include <stdio.h>
#include <math.h>

#define PI 3.141593


double interpolacja_wielomianowa_Newtona(double x, int n, double* xm, double* fm)
{
	double sum = 0,iloczyn;
	
	for(int j = 0; j <= n; j++){
		iloczyn = 1;
		for (int i = 0; i <= j-1; i++)
			iloczyn *= (x - xm[i]);
		sum += fm[j]  * iloczyn;
	}
	return sum;
}

double func(double x)
{
	return 1.0 / ( 1.0 + pow(x,2) );
}



int main()
{
	FILE* file = fopen("zad_1.dat", "wr");    	// <-- dla wersji z węzłami równoległymi
	//FILE* file = fopen("zad_2.dat", "wr");	// <-- dla wersji z zerami wielomianu Czebyszewa
	
	
	
	double x_min = -5.0, x_max = 5.0;
	
	for( int n = 5; n <= 20; n+=5){
		
		double h = (x_max - x_min) / n;
		double xm[n+1];
		double ym[n+1];
		double fm[n+1][n+1];
		double fm_d[n+1];
		
		for (int i = 0; i < n +1; i++){
			
			xm[i] = x_min + i * h;		// <-- dla wersji z węzłami równoległymi
			
			//xm[i] = 0.5 * ( ( x_min - x_max ) * cos(PI *( (double)( 2 * i + 1 ) / (double)( 2 * n + 2 ) ) ) + x_min + x_max );	// <-- dla wersji z zerami wielomianu Czebyszewa
			
			ym[i] = func(xm[i]);
			fm[i][0] = ym[i];
		}
		
		printf("xm:\n");
		for (int i = 0; i < n +1; i++)
			printf("%f ",xm[i]);
		printf("\n\nym:\n");
		for (int i = 0; i < n +1; i++)
			printf("%.7f ",ym[i]);
		printf("\n\nfm:\n");
		printf("\n\n\n");
		
		for(int j = 1; j <= n; j++){
			for(int i = j; i <= n; i++){
				fm[i][j] = ( fm[i][j-1] - fm[i-1][j-1] ) / ( xm[i] - xm[i-j] );
			}
		}


		for(int i = 0; i <= n; i++){
			fm_d[i] = fm[i][i];
		}
	
		double x = x_min, k = 0.01;
	
		while ( x < x_max ){
			printf("%f\t%f\n",x, interpolacja_wielomianowa_Newtona(x,n,xm,fm_d));
			fprintf(file,"%f\t%f\n",x, interpolacja_wielomianowa_Newtona(x,n,xm,fm_d));
			x += k;
		}
		
		printf("\n\n");
		fprintf(file,"\n\n");
	}
	
	fclose(file);
	
	return 0;
}