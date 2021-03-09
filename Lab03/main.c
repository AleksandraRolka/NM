/****** program lab_3, sprawozdanie na jego podstawie ******/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>


int main()
{
	
//nr w nazwie do zmiany przy zmianie wersji danych beta, F_0

	FILE* file = fopen("output_1.dat", "w");
    
//DANE :
	int N = 2000;
	int itr = 0;
	
// odkomentowac odpowiednia wersje
	double beta = 0.0, F_0 = 0.0, omega = 0.8;
//	double beta = 0.4, F_0 = 0.0, omega = 0.8;	
//	double beta = 0.4, F_0 = 0.1, omega = 0.8;
	
	double V_0 = 0.0;
	double x_o = 1.0;
	double w = 1.0;
	double h = 0.02;
	
	double a_1 = 1.0;
	double a_2 = w*w * h*h - 2.0 - beta*h;
	double a_3 = 1 + beta*h;
	

//alokacja wektorow
	double *b = malloc((N+1)*sizeof(double));
	double *d_0 = malloc((N+1)*sizeof(double));
	double *d_1 = malloc((N+1)*sizeof(double));
	double *d_2 = malloc((N+1)*sizeof(double));
	double *x_s = malloc((N+1)*sizeof(double));
	double *x_n = malloc((N+1)*sizeof(double));
	
//wypelnienie wektorow	

	for(int i = 2; i<=N;i++){
		d_0[i] = a_3;
		d_1[i] = a_2;
		d_2[i] = a_1;
		b[i] = F_0 * sin(omega * h * (i - 1.)) * h * h;
	}
	
	d_0[0] = d_0[1] = 1.0;
	d_1[0] = 0.0;
	d_1[1] = -1.0;
	d_2[0] = d_2[1] = 0.0;
	
	b[0] = 1.0;
	b[1] = 0.0;
	
	for(int i = 0; i<=N; i++)
		x_s[i] = 1.0;
	
	double S_s, S_n;
	
	while( itr < 100000){
		++itr;
		S_n = 0.0;
		S_s = 0.0;
		
		for(int i=0;i<=N;i++){
			
			if ( i == 0 )
				x_n[i] = (1 / d_0[i]) * (b[i]);
			else if ( i == 1 )
				x_n[i] = (1 / d_0[i]) * (b[i] - d_1[i] * x_n[i - 1]);
			else
				x_n[i] = (1 / d_0[i]) * ( b[i] - d_1[i] * x_n[i-1] - d_2[i] * x_n[i-2]);
			S_n += x_n[i]*x_n[i];
			S_s += x_s[i]*x_s[i];
		}
		
		
		if( fabs(S_n-S_s) <  1e-6 ){
			break;
		}
		
		
		for(int i=0;i<=N;i++)
			x_s[i] = x_n[i];
		
	}
	
	for(int i=0;i<=N;i++)
		printf("%4f  %4f\n",i*h,x_n[i]);
	
	printf("Liczba iteracji: %d \n",itr);
	
	for(int i=0;i<=N;i++){
        fprintf(file,"%4f %4f\n", i * h, x_n[i]);
    }
    fclose(file);
	
	free(b);
	free(d_0);
	free(d_1);
	free(d_2);
	free(x_s);
	free(x_n);
	return 0;
}