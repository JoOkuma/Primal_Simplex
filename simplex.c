#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifndef DEBUG
	// #define DEBUG
#endif

typedef struct {
	double *val; /* Matrix values */
	int nrow; /* Number of rows */
	int ncol; /* Number of columns */
	int* tbcol; /* Lookup table for columns index */
	int n; /* Total number of elements */
} matrix;

/* Lapack function for linear system solving */
void dgesv_(int* N, int* NRHS, double* A, int* LDA, int* IPIV, double* B, int* LDB, int* INFO);

/* Functions prototype */
matrix* create_matrix(int m, int n);
void destroy_matrix(matrix** A);
matrix* copy_matrix(matrix* A);
void print_matrix(char* txt, matrix* mat);
matrix* transpose_matrix(matrix* A);
matrix* inverse_matrix(matrix* A);
double* solve_matrix_array(matrix* A, double* V);
void print_double_array(char* txt, int n, double* arr);
void print_int_array(char* txt, int n, int* arr);
void switch_columns(matrix* A, int A_i, matrix* B, int B_i);
double* cost_coef(matrix* N, double* cN, double* lambda);


int main(void){

	double *ank, *b, *x_b, *lambda, *y, *tmp_dbl_array; 
	double *cost, *phaseICost, *cB, *cN, *cost_cf;
	double e, tmp_dbl, sec;

	matrix *A, *B, *BT, *N;
	
	clock_t start, diff;

	int m, n;	
	int *B_i, *N_i, *aux_i; /* Arrays for matrices columns saving */
	int i, j, nk, e_i, tmp_i;
	int stop, zero, ite;

	/* Read rows and coluns lengths */
	scanf("%d", &m);
	scanf("%d", &n);
	
	/* Alloc array sizes */
	cost = calloc(n, sizeof(double));
	b = calloc(m, sizeof(double));
	
	/* All matrices are saved as transpose because Lapack and Cblas (Fortran) works this way */
	A = create_matrix(n, m);
	
	/* Read input matrices */
	for(i = 0; i < n; i++){
		scanf("%lf", &cost[i]);
	}

	for(i = 0; i < m; i++){
		scanf("%lf", &b[i]);
	}
	
	for(i = 0; i < A->ncol; i++){
		for(j = 0; j < A->nrow; j++){
			scanf("%lf", &A->val[i + A->tbcol[j]]);
		}
	}

	print_matrix("Matrix A", A);
	print_double_array("Cost Array", A->nrow, cost);
	print_double_array("Array b", A->ncol, b);
	
	/* Phase I */

	printf("\n.: Starting Phase I :.\n");
	
	start = clock();

	/* Create inpendent matrix for Phase I */
	B = create_matrix(m, m);
	
	for(i = 0, j = 0; i < B->nrow; i++, j++){
		B->val[i + j * B->nrow] = 1;
	}

	N = copy_matrix(A);
	N_i = calloc(N->nrow, sizeof(int));

	for(i = 0; i < N->nrow; i++){
		N_i[i] = i;
	}

	B_i = calloc(B->nrow, sizeof(int));
	
	for(j = 0; j < B->nrow; j++, i++){
		B_i[j] = i;
	}

	/* All costs are zero except the artificial variables */
	phaseICost = calloc(N->nrow + B->nrow, sizeof(double));
	cN = phaseICost;
	cB = &phaseICost[N->nrow];

	for(i = 0; i < B->nrow; i++){
		cB[i] = 1;
	}

	aux_i = calloc(B->nrow, sizeof(int));
	memcpy(aux_i, B_i, B->nrow * sizeof(int));

	#ifdef DEBUG
		print_matrix("Matrix A", A);
		print_matrix("Matrix B", B);
		print_matrix("Matrix N", N);
	#endif

	stop = 0;
	ite = 0;
	while(1){
		ite++;

		x_b = solve_matrix_array(B, b);
		#ifdef DEBUG
			print_double_array("Array x_b", B->ncol, x_b);
		#endif

		BT = transpose_matrix(B);
		#ifdef DEBUG
			print_matrix("Matrix B Transpose", BT);
		#endif

		lambda = solve_matrix_array(BT, cB);
		#ifdef DEBUG
			print_double_array("Array Lambda", B->nrow, lambda);
		#endif

		cost_cf = cost_coef(N, cN, lambda);
		#ifdef DEBUG
			print_double_array("Array Cost Coef", N->nrow, cost_cf);
		#endif

		nk = 0;
		for(i = 0; i < N->nrow; i++){
			if(cost_cf[nk] > cost_cf[i]) nk = i;
		}

		ank = &N->val[ N->tbcol[nk] ];
		#ifdef DEBUG
			print_double_array("Array aNK", N->ncol, ank);
		#endif

		y = solve_matrix_array(B, ank);
		#ifdef DEBUG
			print_double_array("Array y", N->ncol, y);
		#endif

		e_i = -1;
		e = 1e10;
		for(i = 0; i < B->nrow; i++){
			tmp_dbl = x_b[i]/y[i];
			if(tmp_dbl > 0 && e > tmp_dbl){
				e_i = i;
				e = tmp_dbl;
			}
		}

		if(e_i == -1){
			printf("\nProblem doesn't have a finite optimal solution\n");
			exit(-1);
		}
		
		/* Switch columns to remove artificial variables */
		tmp_i = N_i[nk];
		N_i[nk] = B_i[e_i];
		B_i[e_i] = tmp_i;
		
		tmp_dbl = cN[nk];
		cN[nk] = cB[e_i];
		cB[e_i] = tmp_dbl;
		
		tmp_dbl_array = calloc(N->ncol, sizeof(double));
		memcpy(tmp_dbl_array, &N->val[ N->tbcol[nk] ], N->ncol * sizeof(double));
		memcpy(&N->val[ N->tbcol[nk] ], &B->val[ B->tbcol[e_i] ], N->ncol * sizeof(double));
		memcpy(&B->val[ B->tbcol[e_i] ], tmp_dbl_array, N->ncol * sizeof(double));
		
		#ifdef DEBUG
			print_matrix("Matrix B", B);
			print_matrix("Matrix N", N);
		#endif

		#ifdef DEBUG
			printf("\nBegin to Free...\n");
		#endif
		free(x_b);
		destroy_matrix(&BT);
		free(lambda);
		free(cost_cf);
		free(y);
		free(tmp_dbl_array);
		#ifdef DEBUG
			printf("Free Done!\n\n");
		#endif

		/* Check if the basic matrix still contains artificial variables */
		stop = 1;
		for(i = 0; i < B->nrow; i++){
			for(j = 0; j < B->nrow; j++){
				if(aux_i[i] == B_i[j]){
					stop = 0;
				}
			}
		}

		if(stop) {
			diff = clock() - start;
			sec = diff / CLOCKS_PER_SEC;
			printf("\n.: Phase I Ended in %6.6lf Seconds with %d Iterations :.\n", sec, ite);
			break;
		}
	}

	/* Remove artifical variables from the problem and creates a new non-basic matrix */
	for(i = 0; i < B->nrow; i++){
		for(j = 0; j < N->nrow; j++){
			if(N_i[j] == aux_i[i]){
				memcpy(&N->val[ N->tbcol[j] ], &N->val[ N->tbcol[N->nrow-1] ], N->ncol * sizeof(double));
				N_i[j] = N_i[N->nrow - 1];
				N->n -= N->ncol;
				N->nrow -= 1;
			}
		}
	}

	/* Create basic and non-basic costs */
	cB = calloc(B->ncol, sizeof(double));
	for(i = 0; i < B->nrow; i++){
		cB[i] = cost[B_i[i]];
	}

	cN = calloc(N->ncol, sizeof(double));
	for(i = 0; i < N->nrow; i++){
		cN[i] = cost[N_i[i]];
	}
	
	free(aux_i);
	free(phaseICost);

	/* Phase II */

	#ifdef DEBUG
		print_matrix("Matrix B", B);
		print_double_array("Array cB", B->nrow, cB);

		
		print_matrix("Matrix N", N);
		print_double_array("Array cN", N->nrow, cN);
	#endif

	printf("\n.: Starting Phase II :.\n");
	start = clock();
	
	stop = 0;
	ite = 0;
	while(1){
		ite++;

		x_b = solve_matrix_array(B, b);
		#ifdef DEBUG
			print_double_array("Array x_b", B->ncol, x_b);
		#endif

		BT = transpose_matrix(B);
		#ifdef DEBUG
			print_matrix("Matrix B Transpose", BT);
		#endif

		lambda = solve_matrix_array(BT, cB);
		#ifdef DEBUG
			print_double_array("Array Lambda", B->nrow, lambda);
		#endif

		cost_cf = cost_coef(N, cN, lambda);
		#ifdef DEBUG
			print_double_array("Array Cost Coef", N->nrow, cost_cf);
		#endif

		/* Optimality test, check if any objective function coefficient are positive */
		nk = 0;
		stop = 1;
		for(i = 0; i < N->nrow; i++){
			if(cost_cf[i] < 0) stop = 0;
			if(cost_cf[nk] > cost_cf[i]) nk = i;
		}

		if(stop){
			
			diff = clock() - start;
			sec = diff / CLOCKS_PER_SEC;
			printf("\n.: Phase II Ended in %6.6lf Seconds with %d Iterations :.\n", sec, ite);
			
			printf("\n");
			for(i = 0; i < A->nrow; i++){
				
				zero = 0;

				if(i != A->nrow / 2){
					printf("                   ");
				} else {
					printf("The solution x_b = ");
				}

				for(j = 0; j < B->nrow; j++){
					if(i == B_i[j]){
						zero = 1;
						printf("[%6.3lf]", x_b[j]);
						break;
					}
				}
				if(!zero) printf("[%6.3lf]", (double) zero);

				if(i != A->nrow / 2){
					printf("\n");
				} else {
					printf(" is optimal.\n");
				}
			}
			printf("\n");

			tmp_dbl =  0;
			for(i = 0; i < B->nrow; i++){
				tmp_dbl += cB[i] * x_b[i];
			}
			printf("The optimal result is f(x_b) = %6.3lf\n\n", tmp_dbl);

			#ifdef DEBUG
				printf("\nBegin to Free...\n");
			#endif
			free(x_b);
			destroy_matrix(&BT);
			free(lambda);
			free(cost_cf);
			#ifdef DEBUG
				printf("Free Done!\n\n");
			#endif

			break;
		}


		ank = &N->val[ N->tbcol[nk] ];
		#ifdef DEBUG
			print_double_array("Array aNK", N->ncol, ank);
		#endif

		y = solve_matrix_array(B, ank);
		#ifdef DEBUG
			print_double_array("Array y", N->ncol, y);
		#endif

		e_i = -1;
		e = 1e10;
		for(i = 0; i < B->nrow; i++){
			tmp_dbl = x_b[i]/y[i];
			if(tmp_dbl > 0 && e > tmp_dbl){
				e_i = i;
				e = tmp_dbl;
			}
		}

		if(e_i == -1){
			printf("\nProblem doesn't have a finite optimal solution\n");
			exit(-1);
		}
		
		/* Switch columns to improve solution*/
		tmp_i = N_i[nk];
		N_i[nk] = B_i[e_i];
		B_i[e_i] = tmp_i;
		
		tmp_dbl = cN[nk];
		cN[nk] = cB[e_i];
		cB[e_i] = tmp_dbl;
		
		tmp_dbl_array = calloc(N->ncol, sizeof(double));
		memcpy(tmp_dbl_array, &N->val[ N->tbcol[nk] ], N->ncol * sizeof(double));
		memcpy(&N->val[ N->tbcol[nk] ], &B->val[ B->tbcol[e_i] ], N->ncol * sizeof(double));
		memcpy(&B->val[ B->tbcol[e_i] ], tmp_dbl_array, N->ncol * sizeof(double));
			
		#ifdef DEBUG
			print_matrix("Matrix B", B);
			print_matrix("Matrix N", N);
		#endif

		#ifdef DEBUG
			printf("\nBegin to Free...\n");
		#endif
		free(x_b);
		destroy_matrix(&BT);
		free(lambda);
		free(cost_cf);
		free(y);
		free(tmp_dbl_array);
		#ifdef DEBUG
			printf("Free Done!\n\n");
		#endif

	}

	free(cost);
	free(b);
	destroy_matrix(&A);
	destroy_matrix(&B);
	destroy_matrix(&N);
	free(cB);
	free(B_i);
	free(cN);
	free(N_i);
	
	return 0;
}

matrix* create_matrix(int row, int col){
	matrix* A;
	int i;
	
	A = calloc(1, sizeof(matrix));
	
	A->nrow = row;
	A->ncol = col;
	A->n = row * col;
	
	A->val = calloc(A->n, sizeof(double));
	A->tbcol = calloc(A->nrow, sizeof(int));
	
	for(i = 0; i < A->nrow; i++){
		A->tbcol[i] = i * A->ncol;
	}

	return A;
}
	
void destroy_matrix(matrix** A){
	free((*A)->val);
	free((*A)->tbcol);
	free(*A);
	*A = NULL;
}

matrix* copy_matrix(matrix* A){
	matrix* B = calloc(1, sizeof(matrix));

	B->ncol = A->ncol;
	B->nrow = A->nrow;
	B->n = A->n;
	B->val = calloc(A->n, sizeof(double));
	B->tbcol = calloc(A->nrow, sizeof(double));

	memcpy(B->val, A->val, A->n * sizeof(double));
	memcpy(B->tbcol, A->tbcol, A->nrow * sizeof(int));

	return B;
}

void print_matrix(char* txt, matrix* A){
	int i, j;
	
	printf("\n %s\n", txt);
	for(i = 0; i < A->ncol; i++){
		for(j = 0; j < A->nrow; j++){
			printf("%6.3lf ", A->val[i + A->tbcol[j]]);
		}
		printf("\n");
	}
}

matrix* transpose_matrix(matrix* A){
	matrix *B;
	int i, j, k;

	B = create_matrix(A->ncol, A->nrow);

	/* If matrix big enough compute transpose in parallel
	 * otherwise the overhead of openmp is not worth it */
	
	if(A->n > 1000){
		#pragma omp parallel for
		for(k = 0; k < A->n; k++){
			i = k / A->ncol;
			j = k % A->ncol;
			B->val[k] = A->val[A->nrow*j + i];
		}
	} else {
		for(k = 0; k < A->n; k++){
			i = k / A->ncol;
			j = k % A->ncol;
			B->val[k] = A->val[A->nrow*j + i];
		}
	}

	return B;
}

double* solve_matrix_array(matrix* M, double* V){
	double *X, *A;
	int *IPIV, INFO, NRHS;

	#ifdef DEBUG
		printf("Entering Solve...\n");
	#endif

	A = calloc(M->n, sizeof(double));
	memcpy(A, M->val, M->n * sizeof(double));
	X = calloc(M->nrow, sizeof(double));
	IPIV = calloc(M->nrow, sizeof(int));
	NRHS = 1;

	memcpy(X, V, M->nrow * sizeof(double));
	
	/*
	 * dgesv (from lapack), references:
	 * dgesv_(N, NRHS, A, LDA, IPIV, B, LDB, INFO);
	 * 
	 * N    = Number of linear equations
	 * NRHS = Number of columns of matrix B
	 * A    = Is a double precision array
	 *		  On entry N x N column-major matrix
	 *		  On exit the factor L and U from fatorization A = P*L*U
	 * LDA  = Leading dimension of A
	 * IPIV = On exit is the pivot indices that defines the permutation of matrix P
	 * B    = Is a double precision array
	 * 		  On exit if INFO = 0, the solution of A*x = B
	 * LDB  = Leading dimenion of B
	 * INFO = If INFO = 0, succesful exit
	 */

	dgesv_(&M->nrow, &NRHS, A, &M->nrow, IPIV, X, &M->nrow, &INFO);

	free(IPIV);
	free(A);

	#ifdef DEBUG
		printf("Solve Done!\n");
	#endif

	return X;
}

void print_double_array(char* txt, int n, double* arr){
	int i;
	
	printf("\n %s\n", txt);
	for(i = 0; i < n; i++) printf("%6.3lf ", arr[i]);
	printf("\n");
}

void print_int_array(char* txt, int n, int* arr){
	int i;
	
	printf("\n %s\n", txt);
	for(i = 0; i < n; i++) printf("%6d ", arr[i]);
	printf("\n");
}

void switch_columns(matrix* A, int A_i, matrix* B, int B_i){
	int i;
	double d;

	for(i = 0; i < A->nrow; i++){
		d = A->val[A_i + A->tbcol[i]];
		A->val[A_i + A->tbcol[i]] = B->val[B_i + B->tbcol[i]];
		B->val[B_i + B->tbcol[i]] = d;
	}
}

double* cost_coef(matrix* N, double* cN, double* lambda){
	int i, j;
	double tmp, *cost_cf;

	cost_cf = calloc(N->nrow, sizeof(double));

	for(i = 0; i < N->nrow; i++){
		tmp = 0;
		for(j = 0; j < N->ncol; j++){
			tmp += N->val[j + N->tbcol[i]] * lambda[j];
		}
		cost_cf[i] += cN[i] - tmp;
	}

	return cost_cf;
}
