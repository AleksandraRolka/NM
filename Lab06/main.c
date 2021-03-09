/****** program lab_6, sprawozdanie nr 6 na jego podstawie ******/
#include <stdio.h>
#include <math.h>

#define N 5
#define IT_MAX 30

double licz_r(double* a, double* b, int n, double xj)
{
	double rj;
	b[n] = 0;
	
	for (int i = n - 1; i >=0; i--)
		b[i] = a[i+1] + xj * b[i+1];
	
	rj = a[0] + xj * b[0];
	
	return rj;
}

int main()
{
	FILE* file = fopen("results.dat", "wr");
	
	double a[N+1] = {240, -196, -92, 33, 14, 1};
	double b[N+1];
	double c[N];
	double x0, x1, Rj, Rj_p;
	int n;
	
	printf("\nL\t it\tx_it\t\tR_it\t\tR_it'\n\n");
	
	for (int L = 1; L <= N; L++){
		n = N - L + 1;
		x0 = 0;
		
		for (int it = 1; it <= IT_MAX; it++){
			Rj = licz_r (a, b, n, x0);
			Rj_p = licz_r (b, c, n-1, x0);
			x1 = x0 - ( Rj / Rj_p );
			
			printf("%d%10d\t%10g\t%10g\t%10g\n", L, it, x1, Rj, Rj_p);
			fprintf(file,"%d%10d\t%10g\t%10g\t%10g\n", L, it, x1, Rj, Rj_p);
			
			if(fabs( x1 - x0 ) < 1.0e-7)
				break;
			
			x0 = x1;
		}
		
		for (int i = 0; i <= (n-1); i++)
			a[i] = b[i];

		printf("\n");
	}
	
	fclose(file);
	
	return 0;
}