/****** Dokończony program z lab_2, sprawozdanie na jego podstawie ******/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_linalg.h>

//kompilacja z: -lgsl -lgslcblas

#define N 4


int main(){

    gsl_matrix* A = gsl_matrix_calloc(N,N);
    gsl_matrix* A_orginal = gsl_matrix_calloc(N,N);
    gsl_matrix* A_L = gsl_matrix_calloc(N,N);
    gsl_matrix* A_U = gsl_matrix_calloc(N,N);
	
    gsl_permutation* p = gsl_permutation_calloc(N);

    int sigma = 2;
    int signum;
	
	
	printf("---------------------------------------------------------------------------\n");

    for(int i=0; i<N; i++) {
		for(int j=0; j<N; j++) {
            gsl_matrix_set(A,i,j,1.0/(i+j+sigma));
            gsl_matrix_set(A_orginal,i,j,1.0/(i+j+sigma));
    }
    }

    printf("Macierz A: \n");
    for(int i=0; i<N; i++) {
		for(int j=0; j<N; j++) {
            printf("%10f ",gsl_matrix_get(A,i,j));
        }
        printf("\n");
    }

    gsl_linalg_LU_decomp(A,p,&signum);

    printf("Macierz LU: \n");
    for(int i=0; i<N; i++) {
		for(int j=0; j<N; j++) {
            printf("%10f ",gsl_matrix_get(A,i,j));
        }
        printf("\n");
    }


    for(int i=0; i<N; i++) {
		for(int j=0; j<N; j++) {
            if(i == j){
                gsl_matrix_set(A_L,i,j,1.0);
                 gsl_matrix_set(A_U,i,j, gsl_matrix_get(A,i,j));
            }
            else if(i < j) {
                gsl_matrix_set(A_L,i,j,0.0);
                gsl_matrix_set(A_U,i,j, gsl_matrix_get(A,i,j));
            }
            else{
                gsl_matrix_set(A_L,i,j, gsl_matrix_get(A,i,j));
                gsl_matrix_set(A_U,i,j,0.0);
            }
        }
    }
    
    printf("Macierz L: \n");
    for(int i=0; i<N; i++) {
		for(int j=0; j<N; j++) {
            printf("%10f ",gsl_matrix_get(A_L,i,j));
        }
        printf("\n");
    }
    printf("Macierz U: \n");
    for(int i=0; i<N; i++) {
		for(int j=0; j<N; j++) {
            printf("%10f ",gsl_matrix_get(A_U,i,j));
        }
        printf("\n");
    }

    printf("\nElementy diagonalne macierzy U: \n");
    for(int i=0; i<N; i++) {
		for(int j=0; j<N; j++) {
            if(i == j)
                printf("%10f ",gsl_matrix_get(A_U,i,j));
        }
    }
    printf("\n");

   
    
    double det_A=1.0;
	
    for(int i=0; i<N; i++) {
		for(int j=0; j<N; j++) {
            if(i == j){
                det_A *=gsl_matrix_get(A_U,i,j);}
        }
    }
	
    det_A *=signum;

    printf("\nWyznacznik macierzy A=%g\n", det_A);
    

//int gsl_linalg_LU_solve(const gsl_matrix * LU, const gsl_permutation * p, const gsl_vector * b, gsl_vector * x)
    
    gsl_vector* b_1 = gsl_vector_calloc(N);
    gsl_vector* b_2 = gsl_vector_calloc(N);
    gsl_vector* b_3 = gsl_vector_calloc(N);
    gsl_vector* b_4 = gsl_vector_calloc(N);

    gsl_vector* x_1 = gsl_vector_calloc(N);
    gsl_vector* x_2 = gsl_vector_calloc(N);
    gsl_vector* x_3 = gsl_vector_calloc(N);
    gsl_vector* x_4 = gsl_vector_calloc(N);
	
	gsl_vector_set_zero(b_1);
	gsl_vector_set_zero(b_2);
	gsl_vector_set_zero(b_3);
	gsl_vector_set_zero(b_4);
	
	gsl_vector_set(b_1,0,1);
	gsl_vector_set(b_2,1,1);
	gsl_vector_set(b_3,2,1);
	gsl_vector_set(b_4,3,1);


    gsl_linalg_LU_solve(A,p,b_1,x_1);
    gsl_linalg_LU_solve(A,p,b_2,x_2);
    gsl_linalg_LU_solve(A,p,b_3,x_3);
    gsl_linalg_LU_solve(A,p,b_4,x_4);

    double val;
    printf("\nMacierz odwrotna A^(-1):\n");

    for(int i=0; i<N; i++) {
            val=gsl_vector_get(x_1,i);
            printf("%15f  ",val);
            val=gsl_vector_get(x_2,i);
            printf("%15f  ",val);
            val=gsl_vector_get(x_3,i);
            printf("%15f  ",val);
            val=gsl_vector_get(x_4,i);
            printf("%15f  ",val);
            printf("\n");
    }
	
	
	gsl_matrix *A_inv = gsl_matrix_calloc(N, N);
	
    for(int i=0; i<N; i++) {
		gsl_matrix_set(A_inv, 0, i, gsl_vector_get(x_1, i));
		gsl_matrix_set(A_inv, 1, i, gsl_vector_get(x_2, i));
		gsl_matrix_set(A_inv, 2, i, gsl_vector_get(x_3, i));
		gsl_matrix_set(A_inv, 3, i, gsl_vector_get(x_4, i));
	}

//macierz C=A*A^(-1)

	gsl_matrix *C = gsl_matrix_calloc(N, N);

	for(int i=0; i<N; i++) {
		for(int j=0; j<N; j++) {
			for(int k=0; k<N; k++) {
				gsl_matrix_set(C, i, j, gsl_matrix_get(C, i, j) + gsl_matrix_get(A_orginal, i, k) * gsl_matrix_get(A_inv, k, j));
			}
		}
	}

	printf("Iloczyn macierzy A*A^(-1): \n");
	for(int i=0; i<N; i++) {
		for(int j=0; j<N; j++) {
            printf("%15g ",gsl_matrix_get(C,i,j));
        }
        printf("\n");
    }

// wskaznik uwarunkowania macierzy:
//          k(A)=||A|| * ||A^(-1)||
// norma:	||A||= max(A_ij) ,   1 < i,j < N

	float norma_A = gsl_matrix_get(A_orginal,1,1);
	float norma_A_1 = gsl_matrix_get(A_inv,1,1);
	
	for(int i=0; i<N; i++) {
		for(int j=0; j<N; j++) {
            if(norma_A < gsl_matrix_get(A_orginal,i,j))
			{
                 norma_A = gsl_matrix_get(A_orginal,i,j);
			}
			if(norma_A_1 < gsl_matrix_get(A_inv,i,j))
			{
                 norma_A_1 = gsl_matrix_get(A_inv,i,j);
			}
		}
	}
	
	float wsk_k_A = norma_A * norma_A_1;
	
	printf("\n||A|| = %f    ||A^-1|| = %f    Wskaznik k(A)= %f\n", norma_A, norma_A_1, wsk_k_A);
	printf("---------------------------------------------------------------------------\n");
	
	

//zwalnianie pamieci

    gsl_matrix_free(A);
    gsl_matrix_free(A_orginal);
    gsl_matrix_free(A_U);
    gsl_matrix_free(A_L);
    gsl_matrix_free(C);
    gsl_matrix_free(A_inv);

    gsl_permutation_free(p);
	
    gsl_vector_free(b_1);
    gsl_vector_free(b_2);
    gsl_vector_free(b_3);
    gsl_vector_free(b_4);
    
    gsl_vector_free(x_1);
    gsl_vector_free(x_2);
    gsl_vector_free(x_3);
    gsl_vector_free(x_4);




    return 0;
}
