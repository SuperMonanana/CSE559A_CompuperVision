/* svmCMath.c */

#include "svmCMath.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
//#include <vecLib/clapack.h>
#else
#include "lapack.h"
#endif



#ifndef MIN
#define MIN(x,y) (((x) < (y)) ? (x) : (y))
#endif

/* Fill a given matrix with an n x n identity matrix */
void matrix_ident(int n, double *A) {
    int i, j;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i == j)
                A[i * n + j] = 1.0;
            else
                A[i * n + j] = 0.0;
        }
    }
}

/* Transpose the m x n matrix A and put the result in the n x m matrix AT */
void matrix_transpose(int m, int n, double *A, double *AT) {
    int i, j;
    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
            AT[j * m + i] = A[i * n + j];
}

/* Scale a matrix by a scalar */
void matrix_scale(int m, int n, double *A, double s, double *R) {
    int i;
    int entries = m * n;

    for (i = 0; i < entries; i++) {
        R[i] = A[i] * s;
    }
}

/* Compute the matrix difference R = A - B */
void matrix_diff(int Am, int An, int Bm, int Bn, double *A, double *B, double *R) {
    int r = Am;
    int c = An;
    int n = r * c, i;

    if (Am != Bm || An != Bn) {
        printf("[matrix_sum] Error: mismatched dimensions\n");
        return;
    }

    for (i = 0; i < n; i++) {
        R[i] = A[i] - B[i];
    }    
}

/* Compute the matrix product R = AB */
void matrix_product(int Am, int An, int Bm, int Bn, 
                    const double *A, const double *B, double *R) {
    int r = Am;
    int c = Bn;
    int m = An;

    int i, j, k;

    if (An != Bm) {
        printf("[matrix_product] Error: the number of columns of A and the "
               "number of rows of B must be equal\n");
        return;
    }

    for (i = 0; i < r; i++) {
        for (j = 0; j < c; j++) {
            R[i * c + j] = 0.0;
            for (k = 0; k < m; k++) {
                R[i * c + j] += A[i * An + k] * B[k * Bn + j];
            }
        }
    }
}

/* Compute the matrix product R = A^T B */
void matrix_transpose_product(int Am, int An, int Bm, int Bn, double *A, double *B, double *R) {
    int r = An;
    int c = Bn;
    int m = Am;

    int i, j, k;

    if (Am != Bm) {
        printf("Error: the number of rows of A and the "
               "number of rows of B must be equal\n");
        return;
    }

    for (i = 0; i < r; i++) {
        for (j = 0; j < c; j++) {
            R[i * c + j] = 0.0;
            for (k = 0; k < m; k++) {
                R[i * c + j] += A[k * An + i] * B[k * Bn + j];
            }
        }
    }
}

/* Compute the matrix product R = A B^T */
void matrix_transpose_product2(int Am, int An, int Bm, int Bn, double *A, double *B, double *R) {
    int r = Am;
    int c = Bm;
    int m = An;

    int i, j, k;

    if (An != Bn) {
        printf("Error: the number of columns of A and the "
               "number of columns of B must be equal\n");
        return;
    }
    
    for (i = 0; i < r; i++) {
        for (j = 0; j < c; j++) {
            R[i * c + j] = 0.0;
            for (k = 0; k < m; k++) {
                R[i * c + j] += A[i * An + k] * B[j * Bn + k];
            }
        }
    }
}


/* Compute an RQ factorization of an m by n matrix A */
void dgerqf_driver(int m, int n, double *A, double *R, double *Q)
{
#if 0
    double *AT;
    int lda = m;
    double *tau;
    int tau_dim = MIN(m, n);
    double *work;
    int block_size = 64; /* Just a guess... */
    int lwork = n * block_size;
    int info;
    double *H;
    double *v, *vvT;
    double *Qtmp;

    int i, j;

    /* Transpose A */
    AT = (double *) malloc(sizeof(double) * m * n);
    matrix_transpose(m, n, A, AT);

    /* Call the LAPACK routine */
    tau = (double *) malloc(sizeof(double) * tau_dim);
    work = (double *) malloc(sizeof(double) * lwork);
    dgerqf_(&m, &n, AT, &lda, tau, work, &lwork, &info);

    if (info < 0) {
        printf("[dgeqrf_driver] An error occurred.\n");

        free(AT);
        free(work);
        free(tau);

        return;
    }

    /* Extract the R matrix */
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            if (j < i)
                R[i * n + j] = 0.0;
            else
                R[i * n + j] = AT[(n - m + j) * m + i];
        }
    }


    /* Now extract the Q matrix */
    H = (double *) malloc(sizeof(double) * n * n);
    v = (double *) malloc(sizeof(double) * n);
    vvT = (double *) malloc(sizeof(double) * n * n);
    Qtmp = (double *) malloc(sizeof(double) * n * n);

    for (i = 0; i < tau_dim; i++) {
        matrix_ident(m, H);

        for (j = 0; j < n; j++) {
            if (j > n - tau_dim + i)
                v[j] = 0.0;
            else if (j == n - tau_dim + i)
                v[j] = 1.0;
            else
                v[j] = AT[j * m + (m-tau_dim+i)];
        }

        matrix_transpose_product2(n, 1, n, 1, v, v, vvT);
        matrix_scale(n, n, vvT, tau[i], vvT);
        matrix_diff(n, n, n, n, H, vvT, H);

        if (i == 0) {
            memcpy(Q, H, sizeof(double) * n * n);
        } else {
            matrix_product(n, n, n, n, Q, H, Qtmp);
            memcpy(Q, Qtmp, sizeof(double) * n * n);
        }
    }

    matrix_product(m, n, n, n, R, Q, H);

    free(H);
    free(v);
    free(vvT);
    free(Qtmp);

    free(tau);
    free(work);
    free(AT);
#endif
}

void matrix_invert(int n, double *A, double *Ainv) {

#if 0
    double *At = (double *) malloc(sizeof(double) * n * n);
    int m = n;
    int lda = n;
    int info;
    int *ipiv = (int *) malloc(sizeof(int) * n);
    int lwork = n * 512;
    int i, j;
    double *work = (double *) malloc(sizeof(double) * lwork);

    /* Transpose A info At like FORTRAN likes */
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            At[i * n + j] = A[j * n + i];

    /* Make calls to FORTRAN routines */
    dgetrf_(&m, &n, At, &lda, ipiv, &info);
    dgetri_(&n, At, &lda, ipiv, work, &lwork, &info);

    /* Transpose back into Ainv */
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            Ainv[i * n + j] = At[j * n + i];

    free(At);
    free(ipiv);
    free(work);
#endif  
}
