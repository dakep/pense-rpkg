/*
 * BLAS.h
 * pense
 *
 * Created by David Kepplinger on 2016-22-01.
 * Copyright (c) 2016 David Kepplinger. All rights reserved.
 *
 * This file is part of the R package pense.
 *
 * pense is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * pense is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with R. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef pense_BLAS_h
#define pense_BLAS_h

#define BLAS_CHAR const char * const
#define BLAS_INT const int

#ifdef __cplusplus
#   ifdef RcppArmadillo__RcppArmadilloForward__h
#       error "BLAS.h must be included before RcppArmadillo.h!"
#   endif
#   define ARMA_DONT_USE_FORTRAN_HIDDEN_ARGS  // Don't use hidden args in calls to BLAS/LAPACK
#   include <RcppArmadillo.h>
#   include <R_ext/RS.h>
#   define F77_R_NAME(x) F77_CALL(x)
#else
#   include <R_ext/RS.h>
#   include <R_ext/BLAS.h>
#   include <R_ext/Lapack.h>

#   define F77_R_NAME(x) F77_NAME(x)
#endif


#ifndef BLAS_extern
#define BLAS_extern extern
#endif

#ifndef La_extern
#define La_extern extern
#endif

#define BLAS_R_NAME(x) F77_NAME(x)
#define LAPACK_R_NAME(x) F77_NAME(x)

#define ARMA_TOKENPASTE(x) arma::x
#define ARMA_TOKENPASTE2(x) ARMA_TOKENPASTE(x)

#ifdef ARMA_USE_BLAS
#   define BLAS_ARMA_NAME(x) ARMA_TOKENPASTE2(F77_NAME(x))

    extern "C" {

    /* DGEMM - perform one of the matrix-matrix operations    */
    /* C := alpha*op( A )*op( B ) + beta*C */
    BLAS_extern void
    F77_NAME(dgemm)(const char *transa, const char *transb, const int *m,
            const int *n, const int *k, const double *alpha,
            const double *a, const int *lda,
            const double *b, const int *ldb,
            const double *beta, double *c, const int *ldc);


    /* DSYR - perform the symmetric rank 1 operation A := alpha*x*x' + A */
    BLAS_extern void
    F77_NAME(dsyr)(const char *uplo, const int *n, const double *alpha,
               const double *x, const int *incx,
               double *a, const int *lda);

    /* DSYRK - perform the symmetric rank k operation C := alpha*A*A' + beta*C */
    BLAS_extern void
    F77_NAME(dsyrk)(const char *uplo, const char *trans, const int *n, const int *k,
                    const double *alpha, const double *a, const int *lda,
                    const double *beta, double *c, const int *ldc);

    BLAS_extern void
    F77_NAME(dtrsv)(const char *uplo, const char *trans,
            const char *diag, const int *n,
            const double *a, const int *lda,
            double *x, const int *incx);

    BLAS_extern void
    F77_NAME(dtrsm)(const char *side, const char *uplo,
            const char *transa, const char *diag,
            const int *m, const int *n, const double *alpha,
            const double *a, const int *lda,
            double *b, const int *ldb);
    }

#else
#   define BLAS_ARMA_NAME(x) F77_NAME(x)
#endif



#ifdef ARMA_USE_LAPACK
#   define LAPACK_ARMA_NAME(x) ARMA_TOKENPASTE2(F77_NAME(x))
#   define LAPACK_ARMA_WRAP_CHAR(x) (char *) x
#   define LAPACK_ARMA_WRAP_INT_PTR(x) (int *) &x

    extern "C" {

    /* DSYEVR - compute all eigenvalues and, optionally, eigenvectors   */
    /* of a real symmetric matrix A					   */
    La_extern void
    F77_NAME(dsyevr)(const char *jobz, const char *range, const char *uplo,
             const int *n, double *a, const int *lda,
             const double *vl, const double *vu,
             const int *il, const int *iu,
             const double *abstol, int *m, double *w,
             double *z, const int *ldz, int *isuppz,
             double *work, const int *lwork,
             int *iwork, const int *liwork,
             int *info);
    }
#else
#   define LAPACK_ARMA_NAME(x) F77_NAME(x)
#   define LAPACK_ARMA_WRAP_CHAR(x) x
#   define LAPACK_ARMA_WRAP_INT_PTR(x) &x
#endif



#define BLAS_DGEMV(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)							\
	BLAS_ARMA_NAME(dgemv)(trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy)

#define BLAS_DSYMV(uplo, n, alpha, a, lda, x, incx, beta, y, incy)								\
	BLAS_ARMA_NAME(dsymv)(uplo, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy)

#define BLAS_DTRSV(uplo, trans, diag, n, a, lda, x, incx)                                       \
    BLAS_R_NAME(dtrsv)(uplo, trans, diag, &n, a, &lda, x, &incx)

#define BLAS_DGEMM(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)				\
	BLAS_R_NAME(dgemm)(transa, transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc)

#define BLAS_DDOT(n, x, incx, y, incy)															\
	BLAS_ARMA_NAME(ddot)(&n, x, &incx, y, &incy)

#define BLAS_DSCAL(n, alpha, x, incx)															\
	BLAS_ARMA_NAME(dscal)(&n, &alpha, x, &incx)

#define BLAS_DSYR(uplo, n, alpha, x, incx, a, lda)												\
	BLAS_R_NAME(dsyr)(uplo, &n, &alpha, x, &incx, a, &lda)

#define BLAS_DSYRK(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)                        \
  BLAS_R_NAME(dsyrk)(uplo, trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc)

#define BLAS_DTRSM(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)                       \
    BLAS_R_NAME(dtrsm)(side, uplo, transa, diag, &m, &n, &alpha, a, &lda, b, &ldb)


#define LAPACK_DPOTRF(uplo, n, a, lda, info)                                                    \
    LAPACK_ARMA_NAME(dpotrf)(LAPACK_ARMA_WRAP_CHAR(uplo), LAPACK_ARMA_WRAP_INT_PTR(n),          \
                             a, LAPACK_ARMA_WRAP_INT_PTR(lda), &info)

#define LAPACK_DSYEVR_ALL_VECTORS(uplo, n, a, lda, abstol, nevalues, evalues, evectors,         \
                                  ldevectors, isuppevectors, work, lwork, iwork, liwork, info)  \
    LAPACK_R_NAME(dsyevr)("V", "A", uplo, &n, a, &lda, NULL, NULL, NULL, NULL, &abstol,         \
                     &nevalues, evalues, evectors, &ldevectors, isuppevectors, work, &lwork,    \
                     iwork, &liwork, &info)


#define LAPACK_DSYEVR_RANGE(uplo, n, a, lda, rangeLower, rangeUpper, abstol, nevalues,          \
                            evalues, evectors, ldevectors, isuppevectors, work, lwork, iwork,   \
                            liwork, info)                                                       \
    LAPACK_R_NAME(dsyevr)("V", "V", uplo, &n, a, &lda, &rangeLower, &rangeUpper, NULL, NULL,    \
                     &abstol, &nevalues, evalues, evectors, &ldevectors, isuppevectors, work,   \
                     &lwork, iwork, &liwork, &info)


#define LAPACK_DGELS(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)                      \
    LAPACK_R_NAME(dgels)(trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info)


#endif
